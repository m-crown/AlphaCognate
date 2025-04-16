import argparse
from pathlib import Path
import requests 

def find_and_download_af_prediction(identifier, outdir = "."):
    #check output dir exists
    outdir_path = Path(outdir)
    outdir_path.mkdir(parents=True, exist_ok=True)
    #check a file containing the identifier itself exists
    files_with_identifier = list(outdir_path.glob(f'*{identifier}*'))
    if files_with_identifier:
        return None, str(files_with_identifier[0])
    
    #specify the API endpoint for the identifier
    endpoint = f"https://alphafold.com/api/prediction/{identifier}"

    response = requests.get(endpoint)

    if response.status_code == 200:
        data = response.json()
        #there should just be one result for an identifier, check this before proceeding
        if len(data) > 1:
            return f"More than one matched result for identifier {identifier}", None
        #download the latest cif file for the prediction
        response = requests.get(data[0]["cifUrl"], stream=True)

        # Check if the request was successful
        cif_path = f'{outdir}/{data[0]["entryId"]}-model_{data[0]["latestVersion"]}.cif'
        if response.status_code == 200:
            with open(cif_path, "wb") as cif_file:
                for chunk in response.iter_content(chunk_size=8192):
                    cif_file.write(chunk)
            return None, cif_path
        else:
            return f'Failed to download file. Status code: {response.status_code}', None
        
    else:
        return f"Failed to retrieve data: {response.status_code}", None

def main():
    """
    Usage:
    
    python3 get_structures.py uniprot_ids.txt target_pdbs.txt data

    """
    parser = argparse.ArgumentParser(description="Get demo structures from a CSV file.")    
    parser.add_argument("uniprot_ids", type=str, help="File containing list of uniprot IDs to download AF predictions for")
    parser.add_argument("output_dir", type=str, help="Path to the output directory - where to copy procoggraph structures to")
    parser.add_argument("complete_file_path", type = str, help = "filepath for output completion file")

    args = parser.parse_args()

    if args.uniprot_ids: 
        with open(args.uniprot_ids, "r") as f:
            uniprot_ids = f.read().splitlines()
        Path(args.output_dir).mkdir(parents=True, exist_ok=True)
        #download the af structures
        for uniprot in uniprot_ids:
            find_and_download_af_prediction(uniprot, args.output_dir)
    #produce empty file to signal download complete to snakemake
    Path(f"{args.complete_file_path}").touch()
if __name__ == "__main__":
    main()