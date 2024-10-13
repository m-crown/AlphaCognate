from gemmi import cif
import gzip
import gemmi 
import argparse
import numpy as np
from sklearn.cluster import HDBSCAN
import re 
import pandas as pd

def main():
    """
    example: python3 filter_transplants.py --structure /raid/MattC/repos/alpha_cognate/alphafold_structures_cognate/alphafold_structures_cognate_af_manifest_ho/transplanted_structures/AF-Q9LQ10-F1-model_4_transplants.cif.gz --output_file test_script_transplant.cif
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--structure', metavar='structure with ligand transplants', type=str,
                        help = "path to the transplants structure file")
    parser.add_argument('--output_file', type=str, help='User-specified output file for the filtered transplants, if not specified defaults to [in]_filtered.cif')
    parser.add_argument('--domain_match', type=bool, default = False, help='Filter for domain matches, default is False')
    parser.add_argument('--cognate_match', type=bool, default = False, help='Filter for cognate matches, default is False')
    parser.add_argument('--best_match', type=bool, default = True, help='Filter for best matches, default is True')
    parser.add_argument('--min_cluster_size', type=int, default = 5, help='Minimum cluster size for re-clustering')
    parser.add_argument('--min_samples', type=int, default = 5, help='Minimum samples for a cluster to be defined in re-clustering')
    parser.add_argument('--recluster', type=bool, default = True, help='Whether to recluster the transplants')
    args = parser.parse_args()

    if args.domain_match:
        domain_match = ["True"]
    else:
        domain_match = ["True", "False"]

    if not args.output_file:
        if args.structure.endswith(".cif.gz"):
            args.output_file = args.structure.replace(".cif.gz", "_filtered.cif")
        else:
            args.output_file = args.structure.replace(".cif", "_filtered.cif")
    
    #load in the cif file and get the structure
    af_structure = cif.read(args.structure)
    query_block = af_structure.sole_block()
    query_struct = gemmi.make_structure_from_block(query_block)

    #extract the alphacognate loop and process
    structure_transplants = pd.DataFrame(query_block.get_mmcif_category("_alphacognate"))
    structure_transplants = structure_transplants.replace("", np.nan)
    structure_transplants[["cluster", "tcs"]] = structure_transplants[["cluster", "tcs"]].astype("float")

    #breakdown the summary files into those with errors and those without
    #as standard, we remove rows from the dataframe where there was an error - so nothing was transplanted from that structure
    #then remaining params set by input.
    transplants_to_retain = structure_transplants.loc[(structure_transplants.error.isna()) & (structure_transplants.domain_profile_match.isin(domain_match)) & (structure_transplants.cognateLigand.isna() == args.cognate_match)].copy()

    if args.recluster:
        if not transplants_to_retain["center_of_mass"].isna().all():
            transplants_to_retain["center_of_mass_split"] = transplants_to_retain["center_of_mass"].str.split(",")
            points = transplants_to_retain.loc[(transplants_to_retain.center_of_mass_split.isna() == False), "center_of_mass_split"].apply(lambda x: [float(y) for y in x]).values
            if len(points) < 5:
                transplants_to_retain.loc[(transplants_to_retain.center_of_mass_split.isna() == False), "cluster"] = -1
            else:
                #discuss this in a disccusion section.
                points = np.array([np.array(point) for point in points])
                params = {"min_cluster_size": args.min_cluster_size, "min_samples": args.min_samples} #specifying the default values here for future configuration by user
                hdb = HDBSCAN(**params).fit(points)

                #add cluster labels to the dataframe - clusters with fewer than 5 members will be labelled as noise and given -1 as a cluster value.
                labels = hdb.labels_
                transplants_to_retain.loc[(transplants_to_retain.center_of_mass_split.isna() == False), "cluster"] = labels
        else:
            transplants_to_retain["cluster"] = np.nan

    if args.best_match:
        unclustered_transplants = transplants_to_retain.loc[transplants_to_retain.cluster == -1]
        clustered_transplants = transplants_to_retain.loc[transplants_to_retain.cluster >= 0]
        best_cluster_transplant = clustered_transplants.loc[clustered_transplants.groupby('cluster')['tcs'].idxmin()]
        transplants_to_retain = pd.concat([unclustered_transplants, best_cluster_transplant])

        
    #filter the chains
    transplant_chains = transplants_to_retain.transplanted_ligand_chain.unique().tolist() + ["A"] #keep the protein (a) chain too! 
    for chain in query_struct[0]:
        if chain.name not in transplant_chains:
            query_struct[0].remove_chain(chain.name)
            
    #need to update the num transplants in the output tsv file and set the accession that is not present in the mmcif file the data is loaded from
    transplants_to_retain["num_transplants"] = len(transplants_to_retain)
    print(f"From {len(structure_transplants)} transplants, {len(transplants_to_retain)} were retained")
    #set the accession to be all parts of the filename before _transplants.cif.gz
    transplants_to_retain["accession"] = re.search(r'AF-[\w-]+(?=_transplants)', args.structure).group(0)

    #remove the split center of mass split column and the num transpalnts from the mmcif - it is implicit in the number of rows.
    transplant_dictionary = transplants_to_retain[[col for col in transplants_to_retain.columns if col not in ["center_of_mass_split", "num_transplants", "accession"]]].fillna("").astype("str").to_dict(orient = "list")
    #add the transplant dictionary as an mmcif category to the query block
    query_block.set_mmcif_category("_alphacognate", transplant_dictionary)

    #write the updated file
    query_block.write_file(f"test_transplants.cif")

if __name__ == "__main__":
    main()