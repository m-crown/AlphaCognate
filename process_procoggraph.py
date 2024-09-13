
#this script will process the output files from the procoggraph pipeline and generate a csv file with a subset of information on single chain domain interactions.
#it will also output information on cognate ligands (the best cognate ligand match(es) will be used in this pipeline.)

##maybe we join domains in the format : domain_accession-domain_ownership-cathcode and then aggregate this into lists
#then downstream, we process the transplants through pdbe-arpeggio and use our domain assignment function to determine if the domains have the same interaction mecahsinsms - we can assign various levels of confidence 
#for the transplant based on this.

cath_new = pd.read_csv("/raid/MattC/repos/ProCogGraphData/procoggraph_20240511-16-56-34/contacts/cath_pdb_residue_interactions.csv.gz", sep = "\t")
cath_filtered = cath_new.loc[cath_new.domain_ownership != "minor"]
cath_filtered["num_chains"] = cath_filtered.groupby("uniqueID")["assembly_chain_id_protein"].transform("nunique")
cath_filtered_single_chain = cath_filtered.loc[cath_filtered.num_chains == 1]
cath_filtered_single_chain.groupby(["pdb_id", "ec_list", "uniqueID", "ligand_uniqueID", "bound_entity_pdb_residues", "assembly_chain_id_ligand"]).agg({"domain_accession": list}).reset_index()

