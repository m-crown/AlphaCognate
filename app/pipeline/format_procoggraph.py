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
    parser.add_argument('--bound_entity_descriptors', type=str, help='The bound entity descriptors file')
    args = parser.parse_args()
    
    scores = pd.read_pickle(args.scores_file)
    
    #load bound entity and cognate information from individual files so dont keep anything else from scores file
    scores = scores[["ec", "pdb_ligand", "cognate_ligand", "score", "parity_smarts"]]
    
    scores_cognate = scores.loc[scores.score >= args.score_cutoff].copy()
    scores_cognate["ec"] = scores_cognate["ec"].str.split(",")
    scores_exploded = scores_cognate.explode("ec")
    scores_mask = scores_exploded.groupby(['ec', 'pdb_ligand'])['score'].transform(max) 
    scores_max = scores_exploded.loc[scores_exploded.score == scores_mask]
    

    cogligs = pd.read_pickle(args.cognate_ligands)
    cogligs.drop(columns = ["compound_reaction"], inplace = True)
    cogligs.rename(columns = {"uniqueID": "cognate_ligand", "entry": "ec"}, inplace = True)

    scores_max_merged = scores_max.merge(cogligs, on = ["cognate_ligand", "ec"], how = "left", indicator = True)
    assert(len(scores_max_merged.loc[scores_max_merged._merge != "both"]) == 0)
    to_drop = ["_merge"]
    scores_max_merged.drop(columns = to_drop, inplace = True)

    #load cath domains mapping
    cath_domains = pd.read_csv(f"{args.cath_domain_ownership}", na_values = ["NaN", "None"], keep_default_na = False, dtype = {"bound_ligand_residue_interactions":"str", "bound_entity_pdb_residues": "str", "cath_architecture": "str", "cath_class": "str", "cath_topology": "str", "cath_homologous_superfamily": "str"}, sep = "\t")

    #find all bound ligands where domains are all from the same chain, as contigs and af structures are not multichain (yet!) - data to be used for domain interaction tools.
    cath_domains["chainUniqueID"] = cath_domains["chainUniqueID"] + "_" + cath_domains["assembly_chain_id_protein"] #the chainUniqueID is not unique across different chains, so we need to make it unique - potentially update this in the future in the procoggraph pipeline
    cath_assembly_chain_grouped = cath_domains.groupby('uniqueID')['chainUniqueID'].nunique().reset_index()
    cath_assembly_chain_grouped_filtered = cath_assembly_chain_grouped[cath_assembly_chain_grouped['chainUniqueID'] == 1]['uniqueID']
    cath_single_chain = cath_domains[cath_domains['uniqueID'].isin(cath_assembly_chain_grouped_filtered)].copy()
    
    cath_single_chain_nonminor = cath_single_chain.loc[cath_single_chain.domain_ownership != "minor"].copy()
    cath_single_chain_nonminor["combined_interaction"] = cath_single_chain_nonminor["cath_code"] + ":" + cath_single_chain_nonminor["domain_ownership"]
    cath_single_chain_nonminor["combined_positions"] = cath_single_chain_nonminor["uniqueID"] + ":" + cath_single_chain_nonminor["bound_entity_pdb_residues"] + ":" + cath_single_chain_nonminor["assembly_chain_id_ligand"] + ":"  + cath_single_chain_nonminor["assembly_chain_id_protein"] + ":" + cath_single_chain_nonminor["proteinStructAsymID"]
    #commented but retain for future use
    #cath_single_chain_nonminor["cath_segments_dict"] = cath_single_chain_nonminor.cath_segments_dict
    #cath_single_chain_nonminor["cath_min_start"] = cath_single_chain_nonminor["cath_segments_dict"].apply(eval).apply(lambda x: min([int(re.search(r"START=([-+\d]+)", y.get("SRANGE")).group(1)) for y in x]))

    #aggregate into all the interactions for each pdb ligand , and combine the lists into a single string to group into duplicate interactions within a structure - we just take one example of this for each structure to minimise within structure transplants
    cath_single_chain_nonminor_subset = cath_single_chain_nonminor[["pdb_id", "ligand_uniqueID","uniqueID", "hetCode", "combined_interaction", "combined_positions", "pdb_ec_list"]].groupby(
        ["pdb_id", "ligand_uniqueID", "uniqueID", "hetCode", "combined_positions", "pdb_ec_list"]).agg({"combined_interaction":list}).reset_index()
    cath_single_chain_nonminor_subset["combined_interaction"] = cath_single_chain_nonminor_subset["combined_interaction"].apply(lambda x: sorted(x)).str.join(";")

    #aggregate the duplicate pdb ligand interactions within a structure (will select a representative for eahc structure later)
    cath_single_chain_nonminor_subset_grouped = cath_single_chain_nonminor_subset.groupby(["pdb_id", "ligand_uniqueID", "hetCode", "pdb_ec_list", "combined_interaction"]).agg({"combined_positions" : list}).reset_index()
    cath_single_chain_nonminor_subset_grouped["combined_positions"] = cath_single_chain_nonminor_subset_grouped["combined_positions"].str.join(";")
    #expand the ec lists to merge cognate ligands
    cath_single_chain_nonminor_subset_grouped["pdb_ec_list"] = cath_single_chain_nonminor_subset_grouped.pdb_ec_list.str.split(",")
    cath_single_chain_nonminor_subset_grouped_exploded = cath_single_chain_nonminor_subset_grouped.explode("pdb_ec_list")

    print(cath_single_chain_nonminor_subset_grouped_exploded.head(10))

    #merge cogligs and remove duplicate cols
    cath_single_chain_nonminor_subset_grouped_exploded_merged = cath_single_chain_nonminor_subset_grouped_exploded.merge(scores_max_merged, left_on = ["ligand_uniqueID", "pdb_ec_list"], right_on = ["pdb_ligand", "ec"], how = "left") #we use a left merge, and transfer pdb ligands. When there is a cognate match we are able to make this annotation, if not we can provide the pdb ligand.
    cath_single_chain_nonminor_subset_grouped_exploded_merged.drop(columns = ["ec", "pdb_ligand"], inplace = True)
    #aggregate the ec lists into those matching a cognate ligand. i.e. where the cognate ligand is the same in two different ec reactions
    cath_single_chain_nonminor_subset_grouped_exploded_merged_ec_grouped = cath_single_chain_nonminor_subset_grouped_exploded_merged.groupby([col for col in cath_single_chain_nonminor_subset_grouped_exploded_merged if col != "pdb_ec_list"], dropna = False).agg({"pdb_ec_list": list}).reset_index()
    cath_single_chain_nonminor_subset_grouped_exploded_merged_ec_grouped["pdb_ec_list"] = cath_single_chain_nonminor_subset_grouped_exploded_merged_ec_grouped["pdb_ec_list"].str.join(",")
    #take a representative for transplanting
    cath_single_chain_nonminor_subset_grouped_exploded_merged_ec_grouped[["uniqueID", "bound_entity_pdb_residues","assembly_chain_id_ligand", "assembly_chain_id_protein","proteinStructAsymID"]] = cath_single_chain_nonminor_subset_grouped_exploded_merged_ec_grouped["combined_positions"].str.split(";").apply(lambda x: x[0]).str.split(":", expand = True)

    #add bound entity/pdb ligand information
    print(cath_single_chain_nonminor_subset_grouped_exploded_merged_ec_grouped.head(1).T)

    #load bound entity information and merge it with the final dataframe.
    bound_descriptors = pd.read_csv(args.bound_entity_descriptors, sep = "\t", na_values = ["NaN", "None"], keep_default_na = False)

    print(bound_descriptors.columns)

    cath_single_chain_nonminor_subset_grouped_exploded_merged_ec_grouped = cath_single_chain_nonminor_subset_grouped_exploded_merged_ec_grouped.merge(bound_descriptors, left_on = ["ligand_uniqueID"], right_on = ["uniqueID:ID(bd-id)"], how = "left", indicator = True)
    assert(len(cath_single_chain_nonminor_subset_grouped_exploded_merged_ec_grouped.loc[cath_single_chain_nonminor_subset_grouped_exploded_merged_ec_grouped._merge != "both"]) == 0)
    cath_single_chain_nonminor_subset_grouped_exploded_merged_ec_grouped.drop(columns = ["_merge", "uniqueID:ID(bd-id)"], inplace = True)

    #save to a file
    #rename cols for transplanting
    cath_single_chain_nonminor_subset_grouped_exploded_merged_ec_grouped.rename(columns = {"cognate_ligand": "cognate_mapping_id", "compound_name": "cognate_mapping_name", "canonical_smiles": "cognate_mapping_smiles", "ligand_db": "cognate_mapping_xref", "pdb_ec_list": "cognate_mapping_ec_list", "hetCode": "ligand_het_code", "uniqueID": "ligand", "descriptor": "ligand_smiles", "description": "ligand_name", "isCofactor" : "cognate_mapping_cofactor", "score": "cognate_mapping_similarity"}, inplace = True)

    print(cath_single_chain_nonminor_subset_grouped_exploded_merged_ec_grouped.head(10).T)

    cath_single_chain_nonminor_subset_grouped_exploded_merged_ec_grouped.to_csv("cath_single_chain_domain_interactions.tsv.gz", index = False, sep = "\t")

if __name__ == "__main__":
    main()
