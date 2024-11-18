#!/usr/bin/env python3

import pandas as pd
import argparse
import re
from pathlib import Path
import os

def main():
    """example: python3 split_foldseek.py --foldseek_file test_foldseek_file --structure_manifest untitled.txt --output_dir test_split --predicted_structure_domains /raid/MattC/repos/alpha_cognate/cath_alphafold_domains.tsv.gz --procoggraph_data cath_single_chain_domain_interactions.tsv.gz"""
    
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--foldseek_file', metavar='cath_domain_ownership', type=str,
                        help="path to cath domain ownership file")
    parser.add_argument('--structure_manifest', type=str, help='List of structures which were input to foldseek')
    parser.add_argument('--procoggraph_data', type=str, help='ProCogGraph cognate ligand mapping file')
    parser.add_argument('--output_dir', type=str, help='output directory for split foldseek formatted files')
    parser.add_argument('--predicted_structure_domains', type=str, help='file with predicted structure domains')
    parser.add_argument('--foldseek_prob_cutoff', type=float, default = 0.5, help='The score cutoff for foldseek prob')
    parser.add_argument('--foldseek_alntmscore_cutoff', type=float, default = 0.5, help='The score cutoff for foldseek alntmscore')
    parser.add_argument('--foldseek_evalue_cutoff', type=float, default = 0.00001, help='The score cutoff for foldseek evalue')
    parser.add_argument('--max_structures', type=int, default = 250, help='The maximum number of foldseek matches to use for transplants')
    args = parser.parse_args()

    #foldseek results match the af structure to the pre-built pdb database. 
    # We then merege this with procoggraph data to get the structures to transplant

    foldseek_combined = pd.read_csv(args.foldseek_file, sep="\t")
    structure_manifest = pd.read_csv(args.structure_manifest, sep=",", names = ["accession", "file_name",  "structure_dir"])
    predicted_structure_domains = pd.read_csv(args.predicted_structure_domains, sep="\t")

    foldseek_combined["pdb_id"] = foldseek_combined["target"].str.extract("^([A-Za-z0-9]+)_bio-h")
    foldseek_combined["target_chain"] = foldseek_combined["target"].str.extract("^[A-Za-z0-9]+_bio-h_([A-Za-z0-9_]+)")
    foldseek_combined["target_file"] = foldseek_combined["target"].str.extract("^([A-Za-z0-9]+_bio-h)")
    foldseek_combined["query_chain"] = "A" #this is always the same for the alphafold model?
    foldseek_combined = foldseek_combined.merge(structure_manifest, left_on = "query", right_on = "file_name", how = "inner")
    ##adding e value threshold equivalent to 1e-5 here helps to significantly reduce number of hits (using evalue in isolation at this threshold does not reduce enough - need to evaluate this in text
    foldseek_filtered = foldseek_combined.loc[(foldseek_combined.prob >= args.foldseek_prob_cutoff) & (foldseek_combined.alntmscore >= args.foldseek_alntmscore_cutoff) & (foldseek_combined.evalue <= args.foldseek_evalue_cutoff)]
    #foldseek_combined["accession"] = foldseek_combined["query"].apply(lambda x: os.path.basename(x).split(".")[0])
    foldseek_filtered_domains = foldseek_filtered.merge(predicted_structure_domains, on = "accession", how = "left")
    procoggraph = pd.read_csv(args.procoggraph_data, sep="\t")
    procoggraph["bound_entity_chain"] = procoggraph["uniqueID"].str.split("_").apply(lambda x: x[2])

    foldseek_filtered_domains_cognate = foldseek_filtered_domains.merge(procoggraph, left_on = ["pdb_id", "target_chain"], right_on = ["pdb_id", "assembly_chain_id_protein"], how = "inner")

    #check if output dir exists, if not create it
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    for _, group in foldseek_filtered_domains_cognate.groupby("query"):
        structure = group["accession"].values[0]
        #additionally limti at n structures maximum
        group_filtered = group.sort_values("evalue", ascending = True).head(args.max_structures).copy()
        group_filtered.to_csv(f"{args.output_dir}/{structure}_foldseek.tsv.gz", sep="\t", index=False, compression = "gzip")

    for _, structure in structure_manifest.iterrows():
        #can check where it is lost here and create an error output. do this later
        if structure.accession not in foldseek_combined["accession"].values:
            row = {"accession": structure.accession, "file_name": structure.file_name, "structure_dir" : structure.structure_dir, "error":"No foldseek hits found"}
            row_df = pd.DataFrame([row])
            row_df.to_csv(f"{args.output_dir}/{structure.accession}_foldseek.tsv.gz", sep="\t", index=False, compression = "gzip")
        elif structure.accession not in foldseek_filtered["accession"].values:
            row = {"accession": structure.accession, "file_name": structure.file_name, "structure_dir" : structure.structure_dir, "error":"No foldseek hits found above cutoffs"}
            row_df = pd.DataFrame([row])
            row_df.to_csv(f"{args.output_dir}/{structure.accession}_foldseek.tsv.gz", sep="\t", index=False, compression = "gzip")
        elif structure.accession not in foldseek_filtered_domains["accession"].values:
            row = {"accession": structure.accession, "file_name": structure.file_name, "structure_dir" : structure.structure_dir, "error":"No domain mapping found for file"}
            row_df = pd.DataFrame([row])
            row_df.to_csv(f"{args.output_dir}/{structure.accession}_foldseek.tsv.gz", sep="\t", index=False, compression = "gzip")
        elif structure.accession not in foldseek_filtered_domains_cognate["accession"].values:
            row = {"accession": structure.accession, "file_name": structure.file_name, "structure_dir" : structure.structure_dir, "error":"No cognate ligand mapping found for file"}
            row_df = pd.DataFrame([row])
            row_df.to_csv(f"{args.output_dir}/{structure.accession}_foldseek.tsv.gz", sep="\t", index=False, compression = "gzip")

if __name__ == "__main__":
    main()