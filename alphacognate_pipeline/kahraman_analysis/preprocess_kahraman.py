import os
import requests
import tarfile
import gzip
import re
from bs4 import BeautifulSoup

from pathlib import Path
import gzip
import pandas as pd
import shutil

def get_uniprot_from_pdb_entity(pdb_id, entity_id = 1):
    """
    Fetches the UniProt accession and name for a given PDB ID and entity ID using the PDBe API.
    Where multiple mappings exist, the mapping with the highest coverage is returned.
    """
    url = f"https://www.ebi.ac.uk/pdbe/api/v2/pdb/entry/uniprot_mapping/{pdb_id}/{entity_id}"
    response = requests.get(url)
    response.raise_for_status()
    if response.status_code == 200:
        data = response.json()
        if pdb_id in data:
            mapping = [{"accession": item['accession'], "name": item['name'], "coverage": item['residues'][0]['endIndex'] - item['residues'][0]['startIndex']} for item in data[pdb_id]['data']]
            max_mapping = max(mapping, key=lambda x: x['coverage'])
            if max_mapping:
                accession = max_mapping['accession']
                name = max_mapping['name']
                return accession, name
        raise ValueError(f"Unexpected response structure: {data}")
    raise requests.exceptions.HTTPError(f"HTTP error occurred: {response.status_code}")
    # return None

def get_ec_taxon_from_uniprot(uniprot_id):
    if pd.isna(uniprot_id):
        return (pd.NA, pd.NA)
    
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}?fields=ec,organism_id"
    headers = {"Accept": "application/json"}
    response = requests.get(url, headers=headers)
    response.raise_for_status()
    if response.status_code == 200:
        data = response.json()
        protein_description = data.get('proteinDescription', {}).get('recommendedName', {})
        organism = data.get('organism', {})
        if 'ecNumbers' in protein_description:
            ec_values = [ec['value'] for ec in protein_description['ecNumbers']]
        else:
            ec_values = pd.NA
        if 'taxonId' in organism:
            organism_id = organism['taxonId']
        else:
            organism_id = pd.NA
        return (ec_values, organism_id)
    raise ValueError(f"Unexpected response structure: {data}")

def get_alphafold_structure(uniprot_id):
    #Could use the API to get the cifUrl but the URL returned is V6 only. We need v4 for our domain annotations
    url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    try:
        response = requests.get(url)
        response.raise_for_status()
        if response.status_code == 200:
            data = response.json()
            for item in data:
                if item.get("uniprotAccession", None) == uniprot_id:
                    cif_url = item.get("cifUrl", None)
                    if cif_url:
                        return cif_url
        raise ValueError(f"Unexpected response structure: {data}")
    except requests.exceptions.HTTPError as e:
        if response.status_code == 404:
            return None
        else:
            raise e 


def main():
    #First write a function to get the uniprot ID from kahraman dataset (assume that no _1 suffix means entity ID 1)

    kahraman_table_1 = pd.read_csv("kahraman_dataset_table1_updated.tsv", sep = "\t")

    # Getting UniProt accessions and names from PDB
    print("Fetching UniProt accessions and names from PDB...")
    kahraman_table_1[["uniprot_accession", "uniprot_name"]] = kahraman_table_1.apply(lambda row: get_uniprot_from_pdb_entity(row["updated_pdb_id"], row["entity_id"]), axis=1, result_type="expand")

    # Annotating EC numbers and Taxon IDs from UniProt
    print("Fetching EC numbers and Taxon IDs from UniProt...")
    kahraman_table_1[["uniprot_ec", "uniprot_taxon"]] = pd.DataFrame.from_records(
        kahraman_table_1["uniprot_accession"].map(
            lambda x: get_ec_taxon_from_uniprot(x) if pd.notna(x) else (pd.NA, pd.NA)
        ),
        index=kahraman_table_1.index,
        columns=["uniprot_ec", "uniprot_taxon"]
    )

    # print("Fetching AlphaFold structure URLs from UniProt...")
    # kahraman_table_1["af_cif_url"] = kahraman_table_1["uniprot_accession"].apply(lambda x: get_alphafold_structure(x) if pd.notna(x) else None)

    # If using the AF v6 structures from API can use the following code
    # kahraman_table_1.loc[kahraman_table_1["af_cif_url"].notna(), "af_cif_filename"] = kahraman_table_1.loc[kahraman_table_1["af_cif_url"].notna(), "af_cif_url"].apply(lambda x: x.split('/')[-1])

    # Set the cif filename for v4 structures
    kahraman_table_1["af_cif_filename"] = kahraman_table_1["uniprot_accession"].apply(lambda x: f"AF-{x}-F1-model_v4.cif" if pd.notna(x) else pd.NA)

    #make the cif directory if it doesn't exist
    if not os.path.exists("kahraman_af_structures"):
        os.makedirs("kahraman_af_structures")

    # Download the Alphafold CIF archive for all v4 structures.
    af_v4_species_to_download = kahraman_table_1["uniprot_taxon"].unique().tolist()
    print(f"Need to download {len(af_v4_species_to_download)} species tar files from AF v4")
    # search for the species tar directories in the v4 alphafold ftp
    # list the species tars available and match to our species list
    v4_url = "https://ftp.ebi.ac.uk/pub/databases/alphafold/v4/swissprot_cif_v4.tar"
    # resp = requests.get(v4_url)
    # resp.raise_for_status()
    # soup = BeautifulSoup(resp.text, "html.parser")

    #use if doing individual species directories
    # # find all tar files
    # tar_files = [a["href"] for a in soup.find_all("a", href=re.compile(r"\.tar$"))]

    # # match species ids in filenames
    # found_species = {}
    # for sid in af_v4_species_to_download:
    #     print(f"Searching for species {sid}")
    #     pattern = re.compile(rf"UP\d+_{sid}_.*\.tar")
    #     matches = [f for f in tar_files if pattern.search(f)]
    #     if len(matches) > 1:
    #         raise ValueError(f"Error: multiple matches for species {sid}: {matches}")
    #     match = matches[0] if matches else None
    #     if match:
    #         found_species[sid] = match
    #         download_url = v4_url + match
    if not os.path.exists("kahraman_af_structures/swissprot_cif_v4.tar"):
        print(f"Downloading {v4_url}...")
        with requests.get(v4_url, stream=True) as r:
            r.raise_for_status()
            with open("kahraman_af_structures/swissprot_cif_v4.tar", "wb") as out:
                for chunk in r.iter_content(chunk_size=8192):
                    out.write(chunk)
    #     else:
    #         print(f"{sid}: MISSING")

    # missing_species = set(af_v4_species_to_download) - set(found_species)
    # print("Missing species:", missing_species)

    #create list of files to extract from the tar archive
    kahraman_table_1["v4_archive"] = "swissprot_cif_v4.tar"
    
    # v4_files = [file for file in kahraman_table_1.af_cif_filename.dropna().tolist()]

    #extract only the required files from the tar archive
    def gzip_decompress(path, delete=False):
        path = Path(path)
        out_path = path.with_suffix('')  # removes .gz

        with gzip.open(path, 'rb') as f_in, open(out_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

        if delete:
            path.unlink()
        
        return out_path

    output_dir = "kahraman_af_structures/"
    # Open a tar archive (compressed or not)
    with tarfile.open(f"kahraman_af_structures/{tar}", "r") as tar:
        for idx, row in kahraman_table_1.iterrows():
            filename = row["af_cif_filename"]
            tar = row["v4_archive"]
            if not pd.isna(tar) and not pd.isna(filename) and not os.path.exists(os.path.join(output_dir, filename)):
                print(f"Extracting {filename} from {tar}")
                try:
                    tar.extract(filename + ".gz", path=output_dir)
                    #decompress the gz file
                    gzip_decompress("kahraman_af_structures/" + filename + ".gz", delete=False)

                except KeyError:
                    print(f"File {filename} not found in the archive (kahraman_af_structures/{tar}).")

    #now we have a directory of input structures, we need to make the structure manifest for the analysis run - accession, filename and structure directory

    kahraman_table_1["include_in_analysis"] = True
    #check that all files were downloaded
    for id, row in kahraman_table_1.iterrows():
        cif_path = os.path.join("kahraman_af_structures", row["af_cif_filename"])

        if not os.path.exists(cif_path):
            print(f"File {cif_path} not found!")
            kahraman_table_1.loc[id,"include_in_analysis"] = False

    kahraman_table_1.to_csv("kahraman_dataset_table1_preprocessed.tsv", sep="\t", index=False)

    structure_manifest = kahraman_table_1.loc[kahraman_table_1["include_in_analysis"] == True, ["af_cif_filename"]]
    structure_manifest["af_cif_basename"] = structure_manifest["af_cif_filename"].apply(lambda x: x.split(".cif")[0])
    structure_manifest["structure_directory"] = "kahraman_analysis/kahraman_af_structures"
    structure_manifest[["af_cif_basename", "af_cif_filename", "structure_directory"]].to_csv("kahraman_af_structure_manifest.csv", header=None, index=False)

if __name__ == "__main__":
    main()