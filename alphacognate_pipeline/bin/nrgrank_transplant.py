import logging
from pathlib import Path

_log = logging.getLogger(__name__)

import argparse
import pandas as pd
import gemmi
from gemmi import cif
from alphacognate_transplant import determine_tcs, generate_chain_ids

#TODO: Should type the various output formats for df

def main():
    parser = argparse.ArgumentParser(description="Transplant ligands into AlphaFold structure with configurable parameters.")
    parser.add_argument("--structure_file", type=str, required=True, help="Path to the AlphaFold structure CIF file.")
    parser.add_argument("--ranking_file", type=str, required=True, help="Path to the ligand ranking CSV file (including col with path to ligand pose to be transplanted).")
    parser.add_argument("--output_cif", type=str, required=True, help="Path to output CIF file.")
    parser.add_argument("--output_tsv", type=str, required=True, help="Path to output TSV file (will be gzipped).")
    parser.add_argument("--output_summary", type=str, required=True, help="Path to output structure summary TSV file (will be gzipped).")
    parser.add_argument("--top_ligands_only", action='store_true', help="If set, only the top-ranked ligand for each binding site will be transplanted (all will be included in the loop).")
    parser.add_argument("--transplant_df_info", type=str, required=True, help="Path to the transplant dataframe tsv file (from alphacognate_transplant.py).")
    parser.add_argument('--structure_summary_file', type=str, help='Path to the structure summary TSV file from alphacognate transplants.')
    args = parser.parse_args()
    
    #TODO - Track the time taken to rank cognate ligands, and add this to the alphafold structure summary loop/file

    # Load AlphaFold structure
    af_structure = cif.read(args.structure_file)
    query_block = af_structure.sole_block()
    structure_info = pd.DataFrame(query_block.get_mmcif_category('_alphacognate_structure.'))
    
    if structure_info["num_transplants"][0] == "0":
        _log.info("No transplants found in structure. Exiting.")
        # Create empty output files to satisfy snakemake requirements - we could also evaluate not using a file output earlier in the tool chain and creating this only here. previously ac was the terminal job for the pipeline
        structure_summary_file = pd.read_csv(f"{args.structure_summary_file}", sep="\t")
        structure_summary_file["nrgrank_runtime"] = 0
        structure_summary_file.to_csv(f"{args.output_summary}", sep="\t", index=False)
        
        query_block.set_mmcif_category("_alphacognate_structure.", structure_summary_file.to_dict(orient="list"))
        query_block.write_file(args.output_cif)

        # Create an empty file with no transplants (predominantly because snakemake is expecting it)
        Path(f"{args.output_tsv}").touch()
        return


    query_struct = gemmi.make_structure_from_block(query_block)
    query_struct.merge_chain_parts()
    chem_comp = pd.DataFrame(query_block.get_mmcif_category('_chem_comp.'))
    struct_conf = query_block.get_mmcif_category('_struct_conf.')
    struct_conf_type = query_block.get_mmcif_category('_struct_conf_type.')

    # Load ligand ranking
    ligand_ranking_df = pd.read_csv(args.ranking_file)
    top_ligands_df = ligand_ranking_df.loc[ligand_ranking_df.groupby('Binding site')['Score'].idxmin()]
    ligand_ranking_df["top_ranked"] = ligand_ranking_df.index.isin(top_ligands_df.index)
    existing_chain_ids = [chain.name for chain in query_struct[0]]
    num_transplants = len(ligand_ranking_df) + len(existing_chain_ids) - 1
    transplant_ids = generate_chain_ids(num_transplants)
    transplant_ids = [chain_id for chain_id in transplant_ids if chain_id not in existing_chain_ids]

    for index, row in ligand_ranking_df.iterrows():
        transplant_id = transplant_ids.pop(0)
        ligand_ranking_df.loc[index, "transplanted_chain_id"] = transplant_id
        if args.top_ligands_only and not row["top_ranked"]:
            #skip the actual transplant if top ligands only is set and this is not a top ranked ligand
            ligand_ranking_df.loc[index, "nrgrank_tcs"] = pd.NA
            continue

        cognate_ligand_file = row["pose_file"]
        cognate_ligand = cif.read(cognate_ligand_file)
        ligand_block = cognate_ligand.sole_block()
        ligand_struct = gemmi.make_structure_from_block(ligand_block)
        ligand_struct.merge_chain_parts()
        ligand_chain = ligand_struct[0][0]

        transplant_chain = gemmi.Chain(transplant_id)
        for res in ligand_chain:
            transplant_chain.add_residue(res)

        ligand_struct.ensure_entities()
        query_struct[0].add_chain(transplant_chain)
        
        nrgrank_tcs =  determine_tcs(query_struct, transplant_id)
        ligand_ranking_df.loc[index, "nrgrank_tcs"] = nrgrank_tcs
        transplant_chem_comp = pd.DataFrame([{"id": "LIG", "type": "non-polymer"}])
        chem_comp = pd.concat([chem_comp, transplant_chem_comp]).drop_duplicates(subset="id", keep="first").fillna("?").reset_index(drop=True)

    query_struct.assign_subchains(force=True)
    query_struct.ensure_entities()
    query_struct.update_mmcif_block(query_block)

    struct_asym = pd.DataFrame(query_block.get_mmcif_category('_struct_asym.'))
    struct_asym.loc[struct_asym['entity_id'] == "1", 'id'] = 'A'
    atom_site = pd.DataFrame(query_block.get_mmcif_category('_atom_site.'))
    atom_site.loc[atom_site['label_entity_id'] == '1', 'label_asym_id'] = 'A'
    query_block.set_mmcif_category("_atom_site.", atom_site.to_dict(orient="list"))
    query_block.set_mmcif_category("_struct_asym.", struct_asym.to_dict(orient="list"))
    query_block.set_mmcif_category("_struct_conf.", struct_conf)
    query_block.set_mmcif_category("_struct_conf_type.", struct_conf_type)
    chem_comp = chem_comp.replace(False, "?")
    chem_comp_dict = chem_comp.to_dict(orient="list")
    query_block.set_mmcif_category("_chem_comp.", chem_comp_dict)
    ligand_ranking_df["top_ranked"] = ligand_ranking_df["top_ranked"].astype(int)
    ligand_ranking_df["Names"] = ligand_ranking_df["Names"].str.strip("'") #this may be unnecessary not - investigate
    ligand_ranking_df.rename(columns={"Binding site": "binding_site"}, inplace=True)
    query_block.set_mmcif_category("_nrgrank.", ligand_ranking_df.drop(columns=["pose_file", "nrgrank_runtime"]).to_dict(orient="list"))
    if len(ligand_ranking_df) > 0:
        runtime = ligand_ranking_df["nrgrank_runtime"].values[0] #pass through runtime from the ranking step to the structure summary file
    else:
        runtime = None
    structure_summary_file = pd.read_csv(f"{args.structure_summary_file}", sep="\t")
    structure_summary_file["nrgrank_runtime"] = runtime
    structure_summary_file.to_csv(f"{args.output_summary}", sep="\t", index=False)
    query_block.set_mmcif_category("_alphacognate_structure.", structure_summary_file.to_dict(orient="list"))

    query_block.write_file(args.output_cif)

    #load the output info file from alphacognate transplant step, and combine the information from the ligand ranking step
    transplants_df = pd.read_csv(args.transplant_df_info, sep = "\t")

    combined_transplants_df = transplants_df.merge(ligand_ranking_df, left_on = ["cluster", "cognate_mapping_name"], right_on = ["binding_site", "Names"])
    combined_transplants_df.drop(columns = ["binding_site", "Names", "pose_file"], inplace = True)
    combined_transplants_df.to_csv(args.output_tsv, sep = "\t", index = False, compression = "gzip")

if __name__ == "__main__":
    main()