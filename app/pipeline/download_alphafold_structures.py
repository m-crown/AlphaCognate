#!/usr/bin/env python

from pathlib import Path
import requests
import pandas as pd

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

##this should be run before the example cath dataset
#in normal use, the user will provide the directory of structures to be analysed

#loop over input file and download structures according to the af-accession column (required)

def main():
    cath_af_structs = pd.read_csv("/raid/MattC/repos/alpha_cognate/cath_alphafold_domains.tsv.gz", sep = "\t", compression = "gzip")
    cath_af_structs["af_accession"] = cath_af_structs["accession"].str.extract("^AF-([A-Za-z0-9]+)-F1-model_4")
    errors = []
    filepaths = []
    directory = "/raid/MattC/repos/alpha_cognate/alphafold_structures"
    for _, row1 in cath_af_structs.iterrows():
        error, af_structure_file = find_and_download_af_prediction(row1.af_accession, directory)
        if error:
            errors.append({"accession": row1.af_accession, "error": error})
        else:
            filepaths.append({"identifier": row1.af_accession, "directory": directory , "full_fp": af_structure_file})
    
    errors_df = pd.DataFrame(errors)
    errors_df.to_csv("/raid/MattC/repos/alpha_cognate/alphafold_structures/alphafold_download_errors.csv")

    filepaths_df = pd.DataFrame(filepaths)
    #get just the filename from the full filepath
    filepaths_df["file"] = filepaths_df[0].apply(lambda x: Path(x).basename())

    filepaths_df.drop(columns =["full_fp"], inplace = True)
    filepaths_df = filepaths_df[["identifier", "file", "directory"]]
    filepaths_df.to_csv("/raid/MattC/repos/alpha_cognate/alphafold_structures/alphafold_structures.csv", index = False , header = None)

if __name__ == "__main__":
    main()



##snakemake --snakefile pipeline.smk --configfile alphacognate_config.yml --cores 40