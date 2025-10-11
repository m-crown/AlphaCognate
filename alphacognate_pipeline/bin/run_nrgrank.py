import os
from pathlib import Path
import subprocess
import time
from nrgrank import generate_conformers, process_target, process_ligands,nrgrank_main
from clean_cif import filter_cif_by_plddt
import sys
import pandas as pd
import io
import argparse
import gzip
import logging
from gemmi import cif

_log = logging.getLogger(__name__)

def convert_to_mol2(input_path, output_path):
    open_babel_command = f"obabel \"{input_path}\" -O \"{output_path}\" ---errorlevel 1"
    print(f'obabel command: {open_babel_command}')
    subprocess.run(open_babel_command, shell=True, check=True)


def write_centers_pdb(coord_list, name, path, ligand_names, extra_info):
    if not os.path.isdir(os.path.dirname(path)):
        os.mkdir(os.path.dirname(path))
    if not os.path.isdir(path):
        os.mkdir(path)
    textfile = open(os.path.join(path, name + ".pdb"), 'w')
    counter = 0
    if extra_info is not None:
        for line in extra_info:
            textfile.write(line)
    for line in coord_list:
        textfile.write("ATOM   {:>4}  C   SPH X   1{:>12} {:>7} {:>7}  1.00  0.10 \n".format(
            str(counter),
            str(round(line[0], 3)),
            str(round(line[1], 3)),
            str(round(line[2], 3)))
        )
        counter += 1
    textfile.close()
    return os.path.join(path, name + ".pdb")


def prep_bd_site(transplant_table, output_dir, bd_site_offset):
    cluster_mapping = (
        transplant_table.loc[
            (transplant_table["cluster"].notna()) &
            (transplant_table["cluster"] != "") &
            (transplant_table["cluster_center"].notna()) &
            (transplant_table["cluster_center"] != ""),
            ["cluster", "cluster_center"]
        ]
        .drop_duplicates()
    ) #temporary - find out why gemmi is giving nan cluster centers.
    cluster_mapping["cluster"] = (
        cluster_mapping["cluster"].astype(float).astype(int)
    )
    cluster_dict = dict(zip(
        cluster_mapping["cluster"],
        cluster_mapping["cluster_center"]
    ))

    _log.error(cluster_mapping)
    if len(cluster_mapping) == 0:
        _log.error("No cluster centers found in transplant table. Exiting.")
        sys.exit(1)
    output_file_list = []
    for cluster_int in cluster_dict.keys():
        cluster_center_coords = cluster_dict[cluster_int]
        cluster_center = cluster_center_coords.split(',')
        cluster_center = [float(x) for x in cluster_center]
        new_list = [
            [cluster_center[0] + bd_site_offset, cluster_center[1], cluster_center[2]],
            [cluster_center[0] - bd_site_offset, cluster_center[1], cluster_center[2]],
            [cluster_center[0], cluster_center[1] + bd_site_offset, cluster_center[2]],
            [cluster_center[0], cluster_center[1] - bd_site_offset, cluster_center[2]],
            [cluster_center[0], cluster_center[1], cluster_center[2] + bd_site_offset],
            [cluster_center[0], cluster_center[1], cluster_center[2] - bd_site_offset],
        ]
        file_path = write_centers_pdb(new_list, f'bd_site_{cluster_int}', output_dir,None, None)
        output_file_list.append(file_path)
    return output_file_list

from rdkit import Chem

def has_r_group_or_dummy(smiles: str) -> bool:
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False  # invalid SMILES
    
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        atomic_num = atom.GetAtomicNum()
        if symbol == '*' or symbol.startswith('R') or atomic_num == 0:
            return True
    return False


def main():
    args = argparse.ArgumentParser(description="Run NRGRank on a target CIF file.")
    args.add_argument('--target_path_cif', type=str, help='Path to the target CIF file.')
    args.add_argument('--output_prefix', type=str, help='Name of the output file (without extension).')
    args.add_argument('--output_dir', type=str, default='.', help='Directory to save output files.')
    args = args.parse_args()

    start_time = time.time()

    target_path_cif = args.target_path_cif
    if not Path(target_path_cif).exists():
        _log.error(f"Target CIF file {target_path_cif} does not exist.")
        sys.exit(1)

    #make sure the output directory exists
    Path(f"{args.output_dir}").mkdir(parents=True, exist_ok=True)
    
    #load the af structure to which ligands will be transplanted
    af_structure = cif.read(target_path_cif)
    structure_block = af_structure.sole_block()
    structure_table = pd.DataFrame(structure_block.get_mmcif_category('_alphacognate_structure.'))
    if structure_table["num_transplants"][0] == "0":
        _log.error(f"No transplants found in structure {target_path_cif}. Exiting.")
        #now need to produce the required output files to satisfy snakemake.
        #TODO: There must be a nicer way to validate col schemas etc and make empty files - pandera? pydantic?
        with open(f"{args.output_dir}/{args.output_prefix}_cognate_ranking.csv", 'w') as f:
            f.write("Name,CF,Binding site,pose_file,nrgrank_runtime\n")
        _log.info(f"No transplants in structure. Empty results saved to {args.output_dir}/{args.output_prefix}_cognate_ranking.csv")
        return
    transplant_table = pd.DataFrame(structure_block.get_mmcif_category('_alphacognate_transplants.'))
    
    cognate_table = pd.DataFrame(structure_block.get_mmcif_category('_alphacognate_cognate_mapping.'))

    # Filter the cognate table to only include ligands without R groups or dummy atoms.
    cognate_table["valid_smiles"] = cognate_table["cognate_mapping_smiles"].apply(has_r_group_or_dummy)
    cognate_table = cognate_table[cognate_table["valid_smiles"] == False]
    cognate_table = cognate_table.drop(columns=["valid_smiles"])

    if len(cognate_table) == 0:
        _log.error("No valid cognate ligands found in cognate mapping table. Exiting.")
        with open(f"{args.output_dir}/{args.output_prefix}_cognate_ranking.csv", 'w') as f:
            f.write("Name,CF,Binding site,pose_file,nrgrank_runtime\n")
        _log.info(f"No valid cognate ligands found. Empty results saved to {args.output_dir}/{args.output_prefix}_cognate_ranking.csv")
        return
    
    analysis_file_path = Path(f"{args.output_dir}/analysis")
    Path(analysis_file_path).mkdir(parents=True, exist_ok=True)

    bd_site_file_list = prep_bd_site(transplant_table, analysis_file_path, 5)

    # Filter residues in structure based on pLDDT and store as a new CIF file and mol2 file. Store in output directory
    filtered_plddt_filename = Path(target_path_cif)
    filtered_plddt_basename = filtered_plddt_filename.name
    target_path_cif_filtered_plddt = analysis_file_path / f"{filtered_plddt_basename}_filtered_plddt.cif"
    target_path_mol2 = analysis_file_path / f"{filtered_plddt_basename}.mol2"
    filter_cif_by_plddt(target_path_cif, target_path_cif_filtered_plddt)
    if not os.path.exists(target_path_mol2):
        convert_to_mol2(target_path_cif_filtered_plddt, target_path_mol2)
    output_all_bds = []
    for count , bd_site_file in enumerate(bd_site_file_list):
        if count == 0:
            first_binding_site = True
        bd_site_id = os.path.basename(os.path.splitext(bd_site_file)[0]).split('_')[-1]
        ### Target prep ###
        processed_target_path = process_target(str(target_path_mol2), str(bd_site_file), ignore_distance_sphere=True, overwrite=True, USE_CLASH=False)
        ### Ligand prep ###
        output_folder_path = os.path.dirname(processed_target_path)
        smiles_df = cognate_table.drop_duplicates(subset='cognate_mapping_smiles')[
            ['cognate_mapping_smiles', 'cognate_mapping_name']]
        smiles_dict = {'Smiles': smiles_df['cognate_mapping_smiles'], 'Name': smiles_df['cognate_mapping_name']}
        if first_binding_site:
            generate_conformers(smiles_dict, output_folder_path, preprocess=False, convert=True)
            process_ligands(os.path.join(output_folder_path, 'conformer.mol2'), 1, output_dir=output_folder_path)
        processed_ligand_path = os.path.join(output_folder_path, 'preprocessed_ligands_1_conf')
        ### Srceening ###
        _, output_dict = nrgrank_main(
            bd_site_id,
            processed_target_path, 
            processed_ligand_path, 
            args.output_dir,
            write_info=False,
            USE_CLASH=False, 
            write_csv=False,
            file_separator="\t",
            output_dictionary=True,
            unique_run_id=f'bd_site_{bd_site_id}'
        )
        
        #TODO: with a future ver of nrgrank with tsv output, we can skip the below parsing which is made to handle commas in smiles etc
        # Safe parse with pandas
        output_bd = pd.DataFrame(output_dict)
        output_bd["Binding site"] = bd_site_id
        output_all_bds.append(output_bd)

        output_bd["Binding site"] = bd_site_id
        output_all_bds.append(output_bd)
    output_all_bds_df = pd.concat(output_all_bds, ignore_index=True)

    output_all_bds_df["Names"] = output_all_bds_df["Names"].str.strip("'")

    output_all_bds_df["pose_file"] = f"{args.output_dir}/ligand_poses/" + output_all_bds_df["Names"] + "_bd_site_" + output_all_bds_df["Binding site"] + ".cif"

    runtime = time.time() - start_time
    output_all_bds_df["nrgrank_runtime"] = runtime

    output_all_bds_df.to_csv(f"{args.output_dir}/{args.output_prefix}_cognate_ranking.csv", index=False)
    _log.info(f"Results saved to {args.output_dir}/{args.output_prefix}_cognate_ranking.csv")

    _log.info("Converting ligands to cif format")
    for ligand_file in output_all_bds_df["pose_file"].unique():
        ligand_file_in = ligand_file.replace('.cif', '.pdb')
        subprocess.run(f"gemmi convert \"{ligand_file_in}\" \"{ligand_file}\"", shell=True, check=True)

    _log.info(f"NRGRank run completed in {runtime/60:.2f} minutes")

if __name__ == '__main__':
    main()