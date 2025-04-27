
import pandas as pd
import requests
import os
import argparse
from pathlib import Path
from collections import defaultdict

#we should also adapt the ted ooutputs here.
def format_cath_af_domains(output_dir):
    #download the CATH AlphaFold 2022 domain predictions from zenodo
    if not os.path.exists(f"{output_dir}/cath_alphafold_domains.tsv.gz"):
        cath_af_url = "https://zenodo.org/records/7404988/files/cath-v4_3_0.alphafold-v2.2022-11-22.tsv?download=1"
        if not os.path.exists(f"{output_dir}/cath-v4_3_0.alphafold-v2.2022-11-22.tsv"):
            with requests.get(cath_af_url, stream=True) as r:
                r.raise_for_status()  # raise error if download failed
                with open(f"{output_dir}/cath-v4_3_0.alphafold-v2.2022-11-22.tsv", "wb") as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)

        cath_af = pd.read_csv(f"{output_dir}/cath-v4_3_0.alphafold-v2.2022-11-22.tsv", sep = "\t")
        cath_af["af-accession"] = cath_af["domain_ID"].str.extract("^af_([A-Za-z0-9]+)_")
        #specific to the formatting for our alphafold downloads, will add AF- and -F1-model_4 to the accession - normally this owuldnt be done as users would provide the structures.
        cath_af["accession"] = "AF-" + cath_af["af-accession"] + "-F1-model_4"

        cath_af[["domain_start", "domain_end"]] = cath_af["domain_ID"].str.extract("^af_[A-Za-z0-9]+_(\d+)_(\d+)", expand = True)
        cath_af[["domain_start", "domain_end"]] = cath_af[["domain_start", "domain_end"]].astype("int")
        cath_af["domain_profile"] = cath_af["sfam_id"] + ":" + cath_af["domain_start"].astype("str") + "-" + cath_af["domain_end"].astype("str")
        cath_af = cath_af.sort_values(["accession", "domain_start"])
        cath_af_grouped = cath_af.groupby("accession").agg({"domain_profile": list}).reset_index()
        cath_af_grouped["domain_profile"] = cath_af_grouped["domain_profile"].str.join(";")
        cath_af_grouped.to_csv(f"{output_dir}/cath_alphafold_domains.tsv.gz", sep = "\t", compression = "gzip", index = False)#.head(20)

def format_ted_af_domains(output_dir):
    # Download the TED domain predictions
    ted_url = "https://zenodo.org/records/13908086/files/ted_365m.domain_summary.cath.globularity.taxid.tsv.gz?download=1"
    if not os.path.exists(f"{output_dir}/ted_365m.domain_summary.cath.globularity.taxid.tsv.gz"):
        with requests.get(ted_url, stream=True) as r:
            r.raise_for_status()
            with open(f"{output_dir}/ted_365m.domain_summary.cath.globularity.taxid.tsv.gz", "wb") as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)

    columns = [
        "ted_id",
        "md5_domain",
        "consensus_level",
        "chopping",
        "nres_domain",
        "num_segments",
        "plddt",
        "num_helix_strand_turn",
        "num_helix",
        "num_strand",
        "num_helix_strand",
        "num_turn",
        "proteome_id",
        "cath_label",
        "cath_assignment_level",
        "cath_assignment_method",
        "packing_density",
        "norm_rg",
        "tax_common_name",
        "tax_scientific_name",
        "tax_lineage"
    ]

    aggregation = defaultdict(list)
    chunk_size = 100_000
    
    for chunk in pd.read_csv(
        f"{output_dir}/ted_365m.domain_summary.cath.globularity.taxid.tsv.gz",
        sep="\t",
        header=None,
        names=columns,
        chunksize=chunk_size
    ):
        chunk = chunk[chunk["cath_assignment_level"] == "H"]
        chunk["accession"] = chunk["ted_id"].str.extract("^(AF-[A-Za-z0-9]+-F1-model_v4)")
        chunk["accession"] = chunk["accession"].str.replace("-model_v4", "-model_4")
        chunk["domain_start"] = chunk["chopping"].str.extract("(\d+)-").astype("int")
        chunk["domain_profile"] = chunk["cath_label"] + ":" + chunk["chopping"]

        for idx, row in chunk.iterrows():
            if pd.notna(row["accession"]) and pd.notna(row["domain_profile"]):
                aggregation[row["accession"]].append(row["domain_profile"])

    # Now create DataFrame
    cath_af_grouped = pd.DataFrame({
        "accession": list(aggregation.keys()),
        "domain_profile": [";".join(profiles) for profiles in aggregation.values()]
    })

    cath_af_grouped.to_csv(f"{output_dir}/cath_alphafold_domains_ted.tsv.gz", sep="\t", compression="gzip", index=False)

def main():
    """
    Usage:
    
    python3 format_cath_af_domains.py --ted

    """
    
    parser = argparse.ArgumentParser("Download and format CATH domain predictions for alphafold structures.")    
    parser.add_argument("--ted", action="store_true", help="Use domain predictions from ted for predicted structures instead of CATH AlphaFold 2021")
    parser.add_argument("--output_dir", default = ".", type=str, help="Path to the output directory - where to copy procoggraph structures to")

    args = parser.parse_args()
    
    #check the output directory exists and make it if it doesnt
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    if args.ted:
        format_ted_af_domains(args.output_dir)
    else:
        format_cath_af_domains(args.output_dir)

    print("CATH domain predictions downloaded and formatted.")

if __name__ == "__main__":
    main()