# /usr/bin/env python3

import pandas as pd
import argparse
import re

def main():
    """example: python3 format_procoggraph.py --cath_domain_ownership /raid/MattC/repos/ProCogGraphData/procoggraph_20240528/flat_files/cath_pdb_residue_interactions.csv.gz --scores_file /raid/MattC/repos/ProCogGraphData/procoggraph_20240528/flat_files/all_parity_calcs.pkl --cognate_ligands /raid/MattC/repos/ProCogGraphData/procoggraph_20240528/flat_files/cognate_ligands_df.pkl"""
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--cath_domain_ownership', metavar='cath_domain_ownership', type=str,
                        help = "path to cath domain ownership file")
    parser.add_argument('--scores_file', type=str, help='The scores file')
    parser.add_argument('--score_cutoff', type=float, default = 0.40, help='The score cutoff for cognate ligands')
    parser.add_argument('--cognate_ligands', type=str, help='The cognate ligands file')
    args = parser.parse_args()
    
    scores = pd.read_pickle(args.scores_file)
    scores_cognate = scores.loc[scores.score >= args.score_cutoff].copy()
    scores_cognate["ec"] = scores_cognate["ec"].str.split(",")
    scores_exploded = scores_cognate.explode("ec")
    scores_mask = scores_exploded.groupby(['ec', 'pdb_ligand'])['score'].transform(max) 
    scores_max = scores_exploded.loc[scores_exploded.score == scores_mask]
    
    cogligs = pd.read_pickle(args.cognate_ligands)
    cogligs.rename(columns = {"uniqueID": "cogliguid"}, inplace = True)

    #load cath domains mapping
    cath_domains = pd.read_csv(f"{args.cath_domain_ownership}", na_values = ["NaN", "None"], keep_default_na = False, dtype = {"bound_ligand_residue_interactions":"str", "bound_entity_pdb_residues": "str", "cath_architecture": "str", "cath_class": "str", "cath_topology": "str", "cath_homologous_superfamily": "str"}, sep = "\t")
    cath_domains["pdb_ec_list"] = cath_domains["pdb_ec_list"].str.split(",")

    #find all bound ligands where domains are all from the same chain, as contigs and af structures are not multichain (yet!) - data to be used for domain interaction tools.
    cath_domains["chainUniqueID"] = cath_domains["chainUniqueID"] + "_" + cath_domains["assembly_chain_id_protein"] #the chainUniqueID is not unique across different chains, so we need to make it unique - potentially update this in the future in the procoggraph pipeline
    cath_assembly_chain_grouped = cath_domains.groupby('uniqueID')['chainUniqueID'].nunique().reset_index()
    cath_assembly_chain_grouped_filtered = cath_assembly_chain_grouped[cath_assembly_chain_grouped['chainUniqueID'] == 1]['uniqueID']
    cath_single_chain = cath_domains[cath_domains['uniqueID'].isin(cath_assembly_chain_grouped_filtered)].copy()
    cath_single_chain_exploded = cath_single_chain.explode("pdb_ec_list")
    cath_single_chain_merged = cath_single_chain_exploded.merge(scores_max, left_on = ["ligand_uniqueID", "ec_list"], right_on = ["pdb_ligand", "ec"], how = "left") #we use a left merge, and transfer pdb ligands. When there is a cognate match we are able to make this annotation, if not we can provide the pdb ligand. 
    cath_merged_nonminor = cath_single_chain_merged.loc[cath_single_chain_merged.domain_ownership != "minor"].copy()
    cath_merged_nonminor["combined_interaction"] = cath_merged_nonminor["cath_code"] + ":" + cath_merged_nonminor["domain_ownership"]
    cath_merged_nonminor["cath_segments_dict"] = cath_merged_nonminor.cath_segments_dict.apply(eval)
    cath_merged_nonminor["cath_min_start"] = cath_merged_nonminor["cath_segments_dict"].apply(lambda x: min([int(re.search(r"START=([-+\d]+)", y.get("SRANGE")).group(1)) for y in x]))
    cath_merged_nonminor_subset = cath_merged_nonminor[["pdb_id", "uniqueID", "hetCode", "pdb_ligand", "bound_entity_pdb_residues", "cognate_ligand", "combined_interaction", "cath_segments_dict", "assembly_chain_id_ligand", "assembly_chain_id_protein", "proteinStructAsymID"]].groupby(["pdb_id", "uniqueID", "hetCode", "pdb_ligand", "bound_entity_pdb_residues", "cognate_ligand", "assembly_chain_id_ligand", "assembly_chain_id_protein", "proteinStructAsymID"]).agg({"combined_interaction":list, "cath_segments_dict": list}).reset_index()#.head(10).T[10:]#[[""]]
    cath_merged_nonminor_subset["combined_interaction"] = cath_merged_nonminor_subset["combined_interaction"].str.join(";")
    cath_merged_nonminor_subset.to_csv("cath_single_chain_domain_interactions.tsv.gz", index = False, sep = "\t")

if __name__ == "__main__":
    main()