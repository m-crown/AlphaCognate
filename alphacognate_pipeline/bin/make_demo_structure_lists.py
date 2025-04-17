import pandas as pd
import argparse
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="Get demo structures from a CSV file.")    
    parser.add_argument("combined_results_file", type=str, help="Path to the combined results file from analysis of a subset of structures against the full procoggraph foldseek db")
    parser.add_argument("output_dir", type=str, help="Path to the output directory - where to copy procoggraph structures to")
    
    args = parser.parse_args()

    combined_results_df = pd.read_csv(args.combined_results_file, sep = "\t")
    combined_results_df["uniprot"] = combined_results_df["accession"].str.extract("AF-([A-Z0-9]+)-F1")
    combined_results_df["target_pdb"] = combined_results_df["transplant_structure"].str.extract("(.+)_bio-h")
    
    
    uniprot_ids = combined_results_df["uniprot"].unique()
    target_pdbs = combined_results_df.target_pdb.unique()
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    with open(f"{args.output_dir}/uniprot_ids.txt", "w") as f:
        for uniprot in uniprot_ids:
            f.write(uniprot + "\n")
    with open(f"{args.output_dir}/procoggraph_assemblies.txt", "w") as f:
        for pdb in target_pdbs:
            f.write(pdb + "\n")
    
if __name__ == "__main__":
    main()