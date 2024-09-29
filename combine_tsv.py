import os
import pandas as pd
import argparse

def combine_tsv_files(input_dir, output_file):
    transplant_dfs = []

    for filename in os.listdir(input_dir):
        if filename.endswith(".tsv") or filename.endswith(".tsv.gz"):
            file_path = os.path.join(input_dir, filename)
            # Read the TSV file
            df = pd.read_csv(file_path, sep="\t")
            # Add dataframe to the list
            transplant_dfs.append(df)

    if not transplant_dfs:
        raise ValueError("No transplant TSV files found in the input directory")

    combined_df = pd.concat(transplant_dfs, axis=0, ignore_index=True, sort=False)
    combined_df.to_csv(output_file, sep="\t", index=False, compression="gzip")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine all transplant TSV files in a directory into a single file.")
    parser.add_argument("input_dir", help="Directory containing transplant TSV files to combine")
    parser.add_argument("output_file", help="Path to save the combined transplant TSV file (will be gzipped)")
    args = parser.parse_args()

    combine_tsv_files(args.input_dir, args.output_file)
