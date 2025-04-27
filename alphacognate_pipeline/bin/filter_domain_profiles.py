import gzip
import pandas as pd
import argparse

#add a comment line to top of the file indiciating the source data e.g. ted or other.

def main():
    """
    Filter domain profiles (from alphacognate preprocessing) for uniprot ids based on the accessions provided in the accessions argument.
    Usage:
    python3 filter_cath_af_domains.py --domains_file cath_alphafold_domains.tsv.gz domain--accessions accessions.txt
    """

    parser = argparse.ArgumentParser("Filter domain profiles for uniprot ids based on the accessions provided in the accessions argument.")
    parser.add_argument("--domains_file", type=str, help="Path to the domain profiles file")
    parser.add_argument("--structures_manifest", type=str, help="Path to the structure manifest file")
    parser.add_argument("--output_dir", default=".", type=str, help="Path to the output directory - where to copy procoggraph structures to")

    args = parser.parse_args()


    accessions_file = pd.read_csv(args.structures_manifest, sep=",", names = ["accession", "file_name",  "structure_dir"])

    keep_accessions = set(accessions_file.accession.unique())

    chunk_size = 100_000

    with gzip.open(f"{args.output_dir}/cath_domain_profiles_filtered.tsv.gz", "wt") as f_out: 
        header_written = False
        # write a comment manually
        #f_out.write(f"#Filtered domain profiles generated from {args.domains_file}\n")

        for chunk in pd.read_csv(args.domains_file, sep="\t", chunksize=chunk_size):
            filtered_chunk = chunk[chunk["accession"].isin(keep_accessions)]
            if not header_written:
                filtered_chunk.to_csv(f_out, sep="\t", index=False)
                header_written = True
            else:
                filtered_chunk.to_csv(f_out, sep="\t", index=False, header=False)

if __name__ == "__main__":
    main()
