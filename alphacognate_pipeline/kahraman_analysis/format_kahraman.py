

import pandas as pd

def format_kahraman():
    """
    Extract and explode uniprot IDs from the Kahraman dataset.
    """
    kahraman_df = pd.read_csv("kahraman_table_1_uniprot_updated.tsv", sep="\t")
    kahraman_df["UniProt"] = kahraman_df["UniProt"].str.split(";")
    kahraman_df = kahraman_df.explode("UniProt").reset_index(drop=True)
    kahraman_df["UniProt"] = kahraman_df["UniProt"].str.strip()
    kahraman_uniprot = kahraman_df["UniProt"].drop_duplicates().reset_index(drop=True)
    kahraman_uniprot.to_csv("kahraman_uniprot_ids.tsv", sep="\t", index=False, header=False)

if __name__ == "__main__":
    format_kahraman()