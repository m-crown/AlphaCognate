import pandas as pd
import argparse
from pathlib import Path
import requests
from requests.adapters import HTTPAdapter, Retry
from gemmi import cif
import math
import gzip
import os 


# Create a session with retries
session = requests.Session()
retries = Retry(
    total=3,                     # Number of total retries
    backoff_factor=1,          # Wait time between retries: 0.5s, 1s, 2s...
    status_forcelist=[502, 503, 504],  # Retry on these HTTP status codes
    allowed_methods=["POST"]     # Make sure POST is included
)
adapter = HTTPAdapter(max_retries=retries)
session.mount('https://', adapter)

def make_modelserver_query(pdb_ids, outdir , chunk_size = 50):
    protonated_id_list = []

    for i in range(0, len(pdb_ids), chunk_size):
        protonated_query_ids = pdb_ids[i:i+chunk_size]
        print(f"Protonated Chunk {i}")
        for id in protonated_query_ids:
            if os.path.exists(f"{outdir}/{id}_bio-h.cif.gz"): 
                print(f"File {id}_bio-h.cif.gz already exists, skipping download.")
                protonated_id_list.append(id)
                protonated_query_ids.remove(id)
        if len(protonated_query_ids) == 0:
            print("All files already downloaded, skipping query")
            continue
        protonated_assembly_query_list = [{"entryId": pdb, "query": "full", "data_source":"pdb-h"} for pdb in protonated_query_ids]
        protonated_assembly_query_json = {"queries": protonated_assembly_query_list}

        # Make the request to the PDBe Model Server
        max_retries = 3
        while max_retries > 0:
            try:
                protonated_assembly_response = session.post('https://www.ebi.ac.uk/pdbe/model-server/v1/query-many',
                                        json=protonated_assembly_query_json,
                                        headers={'accept': 'text/plain', 'Content-Type': 'application/json'})
                break
            except requests.exceptions.RequestException as e:
                print(f"Request failed: {e}. Retrying...")
                max_retries -= 1
                if max_retries == 0:
                    print("Max retries reached. Exiting.")
                    continue

            

        if protonated_assembly_response.status_code == 200:
            file_object = protonated_assembly_response.text
            cif_file = cif.read_string(file_object)
            failed_ids = []
            for block in cif_file:
                if block.find_pair("_model_server_error.message") is not None:
                    print(f"Error in {block.name}:\n{block.find_pair('_model_server_error.message')[1]}")
                    failed_ids.append(block.name) #we need to do something with these ids
                else:
                    pdb_id = block.find_pair("_entry.id")[1].lower()
                    protonated_id_list.append(pdb_id)
                    with gzip.open(f"{outdir}/{pdb_id}_bio-h.cif.gz", "wt") as f:
                        f.write(block.as_string())
            print("Response processed")
        elif protonated_assembly_response.status_code == 503:
            #if the request fails, try again at normal chunk size
            protonated_ids = make_modelserver_query(protonated_query_ids, outdir, chunk_size = chunk_size)
            protonated_id_list.extend(protonated_ids)
        elif protonated_assembly_response.status_code in [502, 504]:
            #process the chunk in smaller groups
            print(f"API Query failed error {protonated_assembly_response.status_code}, trying new chunk size {math.ceil(chunk_size/2)}")
            sub_ids = list(set(protonated_query_ids) - set(protonated_id_list))
            protonated_ids = make_modelserver_query(list(sub_ids), outdir, chunk_size = math.ceil(chunk_size/2))
            protonated_id_list.extend(protonated_ids)
        else:
            # Print an error message if the request failed
            print(f"Request failed with status code {protonated_assembly_response.status_code} on ids: {protonated_query_ids}")

    return protonated_id_list

def main():
    """
    Usage:
    
    python3 bin/get_procoggraph_structures.py data/procoggraph_pdb_ids.txt data/procoggraph_assemblies

    """
    parser = argparse.ArgumentParser(description="Get demo structures from a CSV file.")    
    parser.add_argument("procoggraph_ids", type=str, help="File containing list of uniprot IDs to download AF predictions for")
    parser.add_argument("output_dir", type=str, help="Path to the output directory - where to copy procoggraph structures to")

    args = parser.parse_args()
    
    #check the output directory exists and make it if it doesnt
    with open(args.procoggraph_ids, "r") as f:
        target_pdbs = f.read().splitlines()
    
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    #download the assemblies from the pdbe model server
    protonated_assemblies = make_modelserver_query(target_pdbs, args.output_dir)
    if len(set(target_pdbs) - set(protonated_assemblies)) != 0:
        print(f"Warning: The following PDB IDs were not downloaded: {set(target_pdbs) - set(protonated_assemblies)}")
    #produce empty file to signal download complete to snakemake
    Path(f"{args.output_dir}/download_complete.txt").touch()

if __name__ == "__main__":
    main()
