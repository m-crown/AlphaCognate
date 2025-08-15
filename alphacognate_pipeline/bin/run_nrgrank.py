import os
import subprocess
from nrgrank import generate_conformers,process_target,process_ligands,nrgrank_main
from clean_cif import filter_cif_by_plddt
import sys
import pandas as pd
import io
import argparse
import gzip
import logging

_log = logging.getLogger(__name__)

def convert_to_mol2(input_path, output_path):
    open_babel_command = f"obabel \"{input_path}\" -O \"{output_path}\" ---errorlevel 1"
    print(f'obabel command: {open_babel_command}')
    subprocess.run(open_babel_command, shell=True, check=True)


def get_table_as_df(table_group, table_index):
    data = table_group[table_index]
    data_string = "\n".join(data['Lines'])
    df = pd.read_csv(io.StringIO(data_string), sep=' ', header=None, quotechar="'", names=data['Header'], engine='python')
    return df


def write_pdb(coord_list, name, path, ligand_names, extra_info):
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
    _log.info(transplant_table.head())
    cluster_centers = transplant_table.loc[transplant_table.cluster_center.notna(), 'cluster_center'].unique() #temporary - find out why gemmi is giving nan cluster centers.
    output_file_list = []
    for cluster_center_counter, cluster_center in enumerate(cluster_centers):
        cluster_center = cluster_center.split(',')
        cluster_center = [float(x) for x in cluster_center]
        new_list = [
            [cluster_center[0] + bd_site_offset, cluster_center[1], cluster_center[2]],
            [cluster_center[0] - bd_site_offset, cluster_center[1], cluster_center[2]],
            [cluster_center[0], cluster_center[1] + bd_site_offset, cluster_center[2]],
            [cluster_center[0], cluster_center[1] - bd_site_offset, cluster_center[2]],
            [cluster_center[0], cluster_center[1], cluster_center[2] + bd_site_offset],
            [cluster_center[0], cluster_center[1], cluster_center[2] - bd_site_offset],
        ]
        file_path = write_pdb(new_list, f'bd_site_{cluster_center_counter}', output_dir,None, None)
        output_file_list.append(file_path)
    return output_file_list


def main():
    args = argparse.ArgumentParser(description="Run NRGRank on a target CIF file.")
    args.add_argument('--target_path_cif', type=str, help='Path to the target CIF file.')
    args.add_argument('--output_filename', type=str, help='Name of the output file (without extension).')
    args.add_argument('--output_dir', type=str, default='.', help='Directory to save output files.')
    args = args.parse_args()


    target_path_cif = args.target_path_cif
    if target_path_cif.endswith('.gz'):
        with gzip.open(target_path_cif, 'rt') as f:
            lines = f.readlines()
    else:
        with open(target_path_cif, 'r') as infile:
            lines = infile.readlines()
    group_start = False
    table_group = []
    temp_group = []
    temp_header = []
    table_name = None
    for line in lines:
        if line.startswith('loop_'):
            group_start = True
        if group_start:
            if line != 'loop_\n' and line != '\n':
                temp_line = line.strip()
                if temp_line.startswith('_'):
                    temp_line = temp_line.split('.')
                    table_name = temp_line[0][1:]
                    temp_header.append(temp_line[1])
                else:
                    temp_group.append(temp_line)
        if (line.startswith('\n') and group_start) or (line == lines[-1]):
            group_start = False
            table_group.append({"Name":table_name, 'Header':temp_header, 'Lines':temp_group})
            temp_group = []
            temp_header = []
            table_name=None

    _log.error(f"Found {len(table_group)} tables in the CIF file.")
    transplant_table = get_table_as_df(table_group, -3)
    cognate_table = get_table_as_df(table_group, -2)
    bd_site_file_list = prep_bd_site(transplant_table, os.path.dirname(target_path_cif), 5)
    for bd_site_counter, bd_site_file in enumerate(bd_site_file_list):
        '''
        Target
        '''
        target_path_cif_filtered_plddt = os.path.splitext(target_path_cif)[0] + '_filtered_plddt.cif'
        target_path_mol2 = os.path.splitext(target_path_cif_filtered_plddt)[0] + '.mol2'
        if bd_site_counter == 0:
            filter_cif_by_plddt(target_path_cif, target_path_cif_filtered_plddt)
            if not os.path.exists(target_path_mol2):
                convert_to_mol2(target_path_cif_filtered_plddt, target_path_mol2)
        processed_target_path = process_target(target_path_mol2, bd_site_file, ignore_distance_sphere=True, overwrite=True, USE_CLASH=False)
        '''
        Ligands
        '''
        output_folder_path = os.path.dirname(processed_target_path)
        smiles_df = cognate_table.drop_duplicates(subset='cognate_mapping_smiles')[
            ['cognate_mapping_smiles', 'cognate_mapping_name']]
        smiles_dict = {'Smiles': smiles_df['cognate_mapping_smiles'], 'Name': smiles_df['cognate_mapping_name']}
        if bd_site_counter == 0:
            generate_conformers(smiles_dict, output_folder_path, preprocess=False, convert=True)
            process_ligands(os.path.join(output_folder_path, 'conformer.mol2'), 1, output_dir=output_folder_path)
        processed_ligand_path = os.path.join(output_folder_path, 'preprocessed_ligands_1_conf')
        target = os.path.splitext(os.path.basename(target_path_cif))[0]
        nrgrank_main(
            target, 
            processed_target_path, 
            processed_ligand_path, 
            args.output_dir,
            write_info=False,
            USE_CLASH=False
        )

if __name__ == '__main__':
    main()