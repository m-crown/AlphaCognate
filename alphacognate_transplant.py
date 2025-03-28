import requests
from pathlib import Path
import subprocess
import pandas as pd
import gemmi
import math
import pandas as pd
from gemmi import cif
from operator import itemgetter
import argparse
import time
import itertools
import string
import numpy as np
from sklearn.cluster import HDBSCAN

def get_alignment_positions(foldseek_result, query_struct, query_chain, target_struct, target_chain):
    query_start = foldseek_result.qstart
    target_start = foldseek_result.tstart 
    
    query_alignment_range = [] 
    target_alignment_range = [] 
    
    query_alignment_range_atom = [] 
    target_alignment_range_atom = []
    
    #start index counters at the 0 indexed start positions of alignment.
    query_idx = query_start - 1
    target_idx = target_start - 1
    #OLD (keeping for ref to error): we use get subchain to get the correct struct asym for the protein sequence (observer issues with using index directly for multiple chains with same name e.g. 4kvm chain H - matches ligand instead of protein)
    query_chain_span = query_struct[0][query_chain].get_polymer() #use get polymer to extract the protein sequence from the auth chain id - instead of using struct asym and get_subchain which doesnt work if struct asym not present
    target_chain_span = target_struct[0][target_chain].get_polymer()
    assert(target_chain_span.check_polymer_type().name in ["PeptideL", "PeptideD"]) #check that the chain is a peptide chain
    qaln = foldseek_result.qaln
    taln = foldseek_result.taln

    #loop over foldseek alignments together
    for q_res, t_res in zip(qaln, taln):
        if q_res != "-" and t_res != "-":
            #we could add more checks in here for example ensuring the one letter code of the residue to be extracted is the same as reported in the alignment.
            if target_idx < len(target_chain_span): #we observe with 1xrr that the diproline bound ligands are included in the alignment at the end - this is incorrect, and results in a tseq that is longer than the target sequence in the mmcif structure. so prevent the indexes going above the length of the target sequence which would result in an index error.
                # if pairs are aligned, extract the residues, using the index counter which considers gaps
                query_res = query_chain_span[query_idx]
                target_res = target_chain_span[target_idx]
                
                query_atoms = []
                target_atoms = []
                
                for atom_name in ["CA"]:  #can expand this down the line if necessary but ca is typical.
                    if atom_name in query_res and atom_name in target_res:
                        query_alignment_range.append(query_res[atom_name][0].pos)
                        target_alignment_range.append(target_res[atom_name][0].pos)
                        query_atoms.append(query_res[atom_name][0])
                        target_atoms.append(target_res[atom_name][0])
                
                query_alignment_range_atom.append({"res": query_res, "atoms": query_atoms})
                target_alignment_range_atom.append({"res": target_res, "atoms": target_atoms})
                
            #increment index counter
            query_idx += 1
            target_idx += 1
            
        #when there is an insertion, we move the opposite sequence index up one.
        elif q_res == "-":
            target_idx += 1
        elif t_res == "-":
            query_idx += 1
    assert(len(query_alignment_range) == len(target_alignment_range)) #check that an equal number of residues are attempting to be aligned on to avoid memory issues in gemmi superposition.
    return query_alignment_range, target_alignment_range, query_alignment_range_atom, target_alignment_range_atom


def contact_search_target(target_struct, target_chain, reslist, contact_search_dist = 6, query_alignment_range = None, target_alignment_range = None, query_alignment_range_atom = None, target_alignment_range_atom = None):

    cs = gemmi.ContactSearch(contact_search_dist)
    ns = gemmi.NeighborSearch(target_struct[0], target_struct.cell, contact_search_dist)
    matched_res = []
    for n_ch, chain in enumerate(target_struct[0]):
        if chain.name == target_chain:
            for n_res, res in enumerate(chain):
                if res.seqid.num in reslist:
                    matched_res.append(res.name)
                    for n_atom, atom in enumerate(res):
                        if not atom.is_hydrogen():
                            ns.add_atom(atom, n_ch, n_res, n_atom)
    if len(matched_res) != len(reslist):
        #observe issue with ac1 and glc in 3poc where two residues have same seqid - need to check this in the future
        print("Failed to match all ligand residues, check chain mappings!")
        return None, None
    
    contacts = cs.find_contacts(ns)
    ca_pos = []

    backbone_atoms = ['CA']
    contact_res_list = []
    for res, res2 in zip(target_alignment_range_atom, query_alignment_range_atom):
        if any(x in [contact.partner1.atom for contact in contacts] for x in res["atoms"]):
            contact_res_list.append(res2["res"].seqid.num) #get a list of contacting residues to check correct domains are in proximity
            res_backbone_atoms = [res["res"][a][0].pos for a in backbone_atoms] #use [a][0] instead of sole_atom(a) to get the first confomer when there are multiple confomers
            #append all of the backbone atoms to the list of ca pos
            ca_pos.extend(res_backbone_atoms)
    #now need to get the index position of the list of ca pos in the mol1 stucture
    mol_range_index_positions = [i for i, e in enumerate(target_alignment_range) if e in ca_pos]
    if len(mol_range_index_positions) < 3:
        return None, None
    
    local_rmsd_atoms_mol1 = itemgetter(*mol_range_index_positions)(query_alignment_range)
    local_rmsd_atoms_mol2 = itemgetter(*mol_range_index_positions)(target_alignment_range)
    
    if isinstance(local_rmsd_atoms_mol1, gemmi.Position) or len(local_rmsd_atoms_mol1) < 3 or isinstance(local_rmsd_atoms_mol2, gemmi.Position) or len(local_rmsd_atoms_mol2) < 3:
        return None, None
    sup = gemmi.superpose_positions(local_rmsd_atoms_mol1, local_rmsd_atoms_mol2) 
    
    return sup, contact_res_list

def generate_chain_ids(n):
    # Start with the alphabet for single-character IDs
    single_alpha = list(string.ascii_uppercase[1:])  # Skip 'A' as it's used for the protein
    single_numerics = [str(x) for x in range(0,10)]
    single_chars = single_alpha + single_numerics
    # If more IDs are needed, generate combinations of two characters
    if n > len(single_chars):
        # Generate combinations of two letters
        double_chars = [''.join(pair) for pair in itertools.product(string.ascii_uppercase + "".join(single_numerics), repeat=2)]
        # Combine single and double character IDs
        all_ids = single_chars + double_chars
    else:
        all_ids = single_chars
    if n  > len(all_ids[:n]):
        raise ValueError("Not enough chain IDs to assign to all transplants")
    else:
        return all_ids[:n]

def create_transplant_chain(target_struct, target_chain, reslist, transplant_chain_name):
    new_chain = gemmi.Chain(transplant_chain_name)
    for chain in target_struct[0]:
        if chain.name == target_chain:
            for res in chain:
                if res.seqid.num in reslist:
                    new_chain.add_residue(res)
    return new_chain

#transplant clash score defined as vdw overlap between all atoms within 4A in the new transplanted structure.
#it seems like partner1 is the contact to query - partner2 is the query.
def determine_tcs(query_struct, transplant_chain_id):
    cs = gemmi.ContactSearch(4)
    ns = gemmi.NeighborSearch(query_struct[0], query_struct.cell,5)
    for n_ch, chain in enumerate(query_struct[0]):
        if chain.name == transplant_chain_id:
            for n_res, res in enumerate(chain):
                for n_atom, atom in enumerate(res):
                    if not atom.is_hydrogen():
                        ns.add_atom(atom, n_ch, n_res, n_atom)

    contacts = cs.find_contacts(ns)
    if len(contacts) == 0:
        tcs = 0
    else:
        sq_dev = [(contact.partner1.atom.element.vdw_r + contact.partner2.atom.element.vdw_r - contact.dist)**2 for contact in contacts]
        mean2_dev = sum(sq_dev) / len(sq_dev)
        tcs = math.sqrt(mean2_dev)
    ##we should also consider providing a method to return a list of clashing atoms like alphafill does in the json
    return tcs

import bisect
def check_domain_profile(af_domain_profiles, residues_in_contact = 0, procoggraph_profile = 0):
    """
    function takes as input a domain profile in list format,
    structure domain:range, and a list of residues that are 
    in contact with a ligand. finally , it takes as input
    a coglig profile. assign a domain to each residue contact
    or none if there are none, and then check that the basic 
    profile of contacts matches pcg mapping.
    """
    domain_range_dict = {}
    domain_starts = []
    for i, domain_profile in enumerate(af_domain_profiles):
        domain,res = domain_profile.split(":")
        res_start, res_end = res.split("-")
        domain_range_dict[i] = {"domain": f"{domain}_{i}" , "start" : int(res_start), "end": int(res_end)}
        domain_starts.append((int(res_start), i))
    domain_starts.sort() #verify the starts are sorted in correct order low to high.
    domain_starts_list = [start for start,idx in domain_starts]
    residue_mappings = {domain_range_dict[idx]['domain']: [] for idx in domain_range_dict}
    residue_mappings.setdefault('no-domain_-1', [])
    for res in residues_in_contact:
        idx = bisect.bisect_right(domain_starts_list, res) - 1
        if idx >= 0 and domain_range_dict[idx]["start"] <= res <= domain_range_dict[idx]["end"]:
            domain = domain_range_dict[idx]['domain']
            residue_mappings[domain].append(res)
        else:
            residue_mappings['no-domain_-1'].append(res)
            
    interacting_domains = [(key.split("_")[0],key.split("_")[1], len(value)) for key,value in residue_mappings.items() if len(value) > 0 and key != "no-domain_-1"]
    procoggraph_map = sorted([domain for domain, idx, count in interacting_domains]) == sorted(procoggraph_profile)
    #format the lists and dictionaries into a flat output to join with dataframe, domains semicolon delimited, information on domains colon delimited
    interacting_domains = ";".join([":".join([str(val) for val in tup]) for tup in interacting_domains])
    residue_mappings = ";".join([":".join([key, ",".join([str(res) for res in value])]) for key, value in residue_mappings.items() if value])
    return residue_mappings, interacting_domains, procoggraph_map


def main():
    """This script takes an input foldseek alignment file and transplants cognate ligands from the aligned structures to an AlphaFold model. 
    The foldseek alignment file is expected to be pre-filtered to include only the best matches that ligands will be transplanted from.
    Example command:
    python3 alphacognate_transplant.py --foldseek_file /raid/MattC/repos/alpha_cognate/combined_foldseek_file.tsv --outdir /raid/MattC/repos/alpha_cognate/results --structure_database_directory /raid/MattC/repos/foldseek/pcg_assemblies/

    The expected input format for predicted structure domains is accession,domain_profile(; delimited list of domain:range)
    """
    parser = argparse.ArgumentParser(description='Transplant cognate ligands from experimentally determined structures to an Predicted Structure Model.')
    parser.add_argument('--foldseek_file', type=str, help='Path to the processed foldseek file of structures to transplant.')
    parser.add_argument('--outdir', type=str, help='Path to the output directory. Created if missing')
    parser.add_argument('--structure_database_directory', type=str, help='Path to the directory containing the structures to transplant ligands from.')
    args = parser.parse_args()

    #check directories exist
    Path(f"{args.outdir}").mkdir(parents=True, exist_ok=True)

    #check if foldseek database directory exists - exit with error if not
    if not Path(args.structure_database_directory).exists():
        print(f"Error: {args.structure_database_directory} does not exist.")
        return

    #from foldseek github c++ code - not using but keeping here for future reference
    threeAA2oneAA = {"ALA":'A', "ARG":'R', "ASN":'N', "ASP":'D', "CYS":'C', "GLN":'Q', "GLU":'E', "GLY":'G', "HIS":'H', "ILE":'I', "LEU":'L', "LYS":'K', "MET":'M', "PHE":'F', "PRO":'P', "SER":'S', "THR":'T', "TRP":'W', "TYR":'Y', "VAL":'V', "MSE":'M', "MLY":'K', "FME":'M', "HYP":'P',"TPO":'T', "CSO":'C', "SEP":'S', "M3L":'K',"HSK":'H', "SAC":'S', "PCA":'E', "DAL":'A',"CME":'C', "CSD":'C', "OCS":'C', "DPR":'P',"B3K":'K', "ALY":'K', "YCM":'C', "MLZ":'K',"4BF":'Y', "KCX":'K', "B3E":'E', "B3D":'D',"HZP":'P', "CSX":'C', "BAL":'A', "HIC":'H',"DBZ":'A', "DCY":'C', "DVA":'V', "NLE":'L',"SMC":'C', "AGM":'R', "B3A":'A', "DAS":'D',"DLY":'K', "DSN":'S', "DTH":'T', "GL3":'G',"HY3":'P', "LLP":'K', "MGN":'Q', "MHS":'H',"TRQ":'W', "B3Y":'Y', "PHI":'F', "PTR":'Y',"TYS":'Y', "IAS":'D', "GPL":'K', "KYN":'W',"SEC":'C'}


    #merge procoggraph data and foldseek data. then we do a foldseek itterrows but we should loop over the structures, where there could be multiple ligands to transplant in the same structure.
    
    foldseek_file = pd.read_csv(args.foldseek_file, sep="\t")
    foldseek_file["fp"] = foldseek_file["structure_dir"] + "/" + foldseek_file["file_name"]
    predicted_structure_id = foldseek_file["accession"].values[0]
    predicted_structure_file = foldseek_file["fp"].values[0]

    if "error" in foldseek_file.columns:
        result = {"num_transplants": 0, "accession": predicted_structure_id, "query_structure": foldseek_file["file_name"].values[0], "transplant_structure": f"{args.outdir}/{predicted_structure_id}_transplants.cif.gz", "foldseek_rmsd": np.nan, "global_rmsd": np.nan, "ligand": np.nan, 'hetCode': np.nan, 'cognateLigand': np.nan, "interaction": np.nan, "local_rmsd": np.nan, "tcs": np.nan, "transplanted_structure_path": np.nan, "transplanted_ligand_chain": np.nan, "transplanted_ligand_residues": np.nan, "domain_residue_contacts": np.nan, "domain_residue_counts": np.nan, "domain_profile_match": np.nan, "center_of_mass": np.nan, "error": foldseek_file.error.values[0]}
        result_df = pd.DataFrame([result])
        result_df.to_csv(f"{args.outdir}/{predicted_structure_id}_transplants.tsv.gz", sep = "\t", index = False, compression = "gzip")

        #despite not doing anything to it, process the predicted structure with gemmi and output so that snakemake is happy.
        af_structure = cif.read(predicted_structure_file)
        query_block = af_structure.sole_block()
        query_block.write_file(f"{args.outdir}/{predicted_structure_id}_transplants.cif")
        return
    
    query = foldseek_file["query"].values[0]
    foldseek_file["domain_profile"] = foldseek_file["domain_profile"].str.split(";")
    foldseek_file["domain"] = foldseek_file["domain_profile"].apply(lambda x: [y.split(":")[0] for y in x])
    #convert bound entity pdb residues to list of ints - force to be a string first in case of places whrre there are no pipe delimiters in files.
    foldseek_file["bound_entity_pdb_residues"] = foldseek_file["bound_entity_pdb_residues"].astype("str").str.split("|").apply(lambda x : [int(y) for y in x])
    foldseek_file["interaction_domains"] = foldseek_file["combined_interaction"].str.split(";")
    foldseek_file["interaction_domains"] = foldseek_file["interaction_domains"].apply(lambda x : [y.split(":")[0] for y in x])
    start_time = time.time()
    transplants = []
    error = None
    
    num_transplants = foldseek_file.uniqueID.nunique() #its the number of ligands, not number of ligands mapped to cogligs that are transplanted
    transplant_chain_ids = generate_chain_ids(num_transplants)

    #load the af structure to which ligands will be transplanted
    #store some of the loops that get lost in processing to be handled later.
    af_structure = cif.read(predicted_structure_file)
    query_block = af_structure.sole_block()
    query_struct = gemmi.make_structure_from_block(query_block)
    query_struct.merge_chain_parts()
    chem_comp = pd.DataFrame(query_block.get_mmcif_category('_chem_comp.'))
    struct_conf = query_block.get_mmcif_category('_struct_conf.')
    struct_conf_type = query_block.get_mmcif_category('_struct_conf_type.')

    #we make a subset of the deduplicated bound entity specific information needed for transplanting and iterate through this.
    #post transplanting, we merge back the cognate information.
    foldseek_transplants = foldseek_file[["target", "assembly_chain_id_protein", "query_chain", "target_file", "uniqueID", "bound_entity_pdb_residues", "assembly_chain_id_ligand", "combined_interaction", "domain_profile", "interaction_domains", "rmsd", "qstart","tstart","qaln","taln"]].drop_duplicates(subset = ["target", "uniqueID", "assembly_chain_id_protein", "query_chain", "assembly_chain_id_ligand"])
    #and the information we merge after transplant, on uniqueID
    foldseek_file["compound_name"] = foldseek_file.compound_name.apply(lambda x: x.split("|")[0] if isinstance(x, str) else np.nan)
    foldseek_other = foldseek_file[["uniqueID", "cognate_ligand", "hetCode", "compound_name", "pdb_ec_list", "score"]].rename(columns = {"uniqueID": "ligand","pdb_ec_list":"ecList", "score": "parityScore", "cognate_ligand": "cognateLigand"})
    if not Path(f"{args.outdir}/{predicted_structure_id}_transplants.tsv.gz", sep = "\t").exists():
        for index, row in foldseek_transplants.iterrows():
            result = {"num_transplants": num_transplants}
            
            target = row["target"]
            target_chain = row["assembly_chain_id_protein"]
            query_chain = row["query_chain"]

            target_doc = cif.read(f"{args.structure_database_directory}/{row.target_file}.cif.gz")
            target_block = target_doc.sole_block()
            target_struct = gemmi.make_structure_from_block(target_block)
            #merge chain parts to help with chain naming issue, use auth ids and get_polymer to identify protein chains
            target_struct.merge_chain_parts()
            query_alignment_range, target_alignment_range, query_alignment_range_atom, target_alignment_range_atom = get_alignment_positions(row, query_struct, query_chain, target_struct, target_chain)

            #superpose the query structure onto the target structure, using the alignment from foldseek result (the global alignment)
            
            
            sup = gemmi.superpose_positions(query_alignment_range, target_alignment_range)
            global_rmsd = sup.rmsd
            target_struct[0].transform_pos_and_adp(sup.transform) #align the target structure onto the query structure

            #get the residues in the assembly that make up the ligand (and the chain they are in)
            res_list = row.bound_entity_pdb_residues
            res_chain = row.assembly_chain_id_ligand

            #identify the contacts to the ligand that should be used to perform the 'local' alignment.
            local_superposition, contact_res_list = contact_search_target(target_struct, res_chain, res_list, 6, query_alignment_range, target_alignment_range, query_alignment_range_atom, target_alignment_range_atom)
            ##we can add a value of matched domain profiles? by checking the residue ranges?
            if local_superposition is not None:
                local_rmsd = local_superposition.rmsd
            else:
                error = "Not enough local contacts"
                result.update({"accession": predicted_structure_id, "query_structure": query, "transplant_structure": target, "foldseek_rmsd": row.rmsd, "global_rmsd": global_rmsd, "ligand": row.uniqueID, "interaction": row.combined_interaction, "error": error, "center_of_mass": np.nan})
                transplants.append(result)
                continue
            
            #check the domain profile of the ligand contacts matches the procoggraph mapping
            #the procoggraph map boolean value describes the match between domain contacts in af structure and expected interacting domains from coglig match. but we do not filter on this, leave it up to end user.
            residue_mappings, interacting_domains, procoggraph_map = check_domain_profile(row.domain_profile, contact_res_list, row.interaction_domains)


            target_struct[0].transform_pos_and_adp(local_superposition.transform)
            
            transplant_chain_id = transplant_chain_ids.pop(0)
            transplant_chain = create_transplant_chain(target_struct, res_chain, res_list, transplant_chain_id)
            transplant_chain_center_of_mass = transplant_chain.calculate_center_of_mass().tolist()
            if all([pd.isna(x) for x in transplant_chain_center_of_mass]):
                transplant_chain_center_of_mass = np.nan
            else:
                transplant_chain_center_of_mass = ",".join([str(x) for x in transplant_chain_center_of_mass])
            
            #get the transplant res names to update the chem_comp_table.
            transplant_res_names = [res.name for res in transplant_chain]
            transplant_chem_comp = pd.DataFrame(target_block.get_mmcif_category("_chem_comp."))
            transplant_chem_comp_filtered = transplant_chem_comp.loc[transplant_chem_comp.id.isin(transplant_res_names)]

            #add any new chem_comp entires from the transplant to the overall chem_comp df
            chem_comp = pd.concat([chem_comp, transplant_chem_comp_filtered]).drop_duplicates(subset="id", keep="first").reset_index(drop=True)
            
            query_struct[0].add_chain(transplant_chain)

            #TODO: explore if determine tcs could run before doing the actual transpalnt, and only pop a new chain id if the tcs is better than any pre-existing transplant of the same bound entity - cognate ligand mapping?
            tcs = determine_tcs(query_struct, transplant_chain_id)

            result.update({"accession": predicted_structure_id, "query_structure": query, "transplant_structure": target, "foldseek_rmsd": row.rmsd, "global_rmsd": global_rmsd, "ligand": row.uniqueID, "interaction": row.combined_interaction, "local_rmsd": local_rmsd, "tcs": tcs, "transplanted_structure_path": f"{args.outdir}/{predicted_structure_id}_transplants.cif.gz", "transplanted_ligand_chain": transplant_chain_id, "transplanted_ligand_residues": ",".join([str(x) for x in res_list]), "domain_residue_contacts": residue_mappings, "domain_residue_counts": interacting_domains, "domain_profile_match": procoggraph_map, "center_of_mass": transplant_chain_center_of_mass})
            transplants.append(result)
        #convert individual transplant dictionaries into a single dataframe
        transplants_df = pd.DataFrame(transplants)
        transplants_df = transplants_df.merge(foldseek_other, on = "ligand", how = "left")
        #calculate cluster membership with default hdbscan parameters for transplants
        if not transplants_df["center_of_mass"].isna().all():
            transplants_df["center_of_mass_split"] = transplants_df["center_of_mass"].str.split(",")
            points = transplants_df.loc[(transplants_df.center_of_mass_split.isna() == False), "center_of_mass_split"].apply(lambda x: [float(y) for y in x]).values
            if len(points) < 5:
                transplants_df.loc[(transplants_df.center_of_mass_split.isna() == False), "cluster"] = -1
            else:
                #discuss this in a disccusion section.
                points = np.array([np.array(point) for point in points])
                params = {"min_cluster_size": 5, "min_samples": 5} #specifying the default values here for future configuration by user
                hdb = HDBSCAN(**params).fit(points)

                #add cluster labels to the dataframe - clusters with fewer than 5 members will be labelled as noise and given -1 as a cluster value.
                labels = hdb.labels_
                transplants_df.loc[(transplants_df.center_of_mass_split.isna() == False), "cluster"] = labels
        else:
            transplants_df["cluster"] = np.nan
        #add runtime to the dataframe and output to a tsv file
        transplants_df["runtime"] = time.time() - start_time

        #remove the split center of mass split column and the num transpalnts from the mmcif - it is implicit in the number of rows.
        transplant_dictionary = transplants_df[[col for col in transplants_df.columns if col not in ["center_of_mass_split", "num_transplants", "accession"]]].fillna("").astype("str").to_dict(orient = "list")
        #add the transplant dictionary as an mmcif category to the query block
        #TODO: Consider adding information on AlphaCognate processing elsewhere in the file, and splitting this table.
        query_block.set_mmcif_category("_alphacognate", transplant_dictionary)
        
        #we have added new entities. Need to make sure theyre present in the _entity. loop.
        query_struct.ensure_entities()

        #update the query block with the new structure.
        query_struct.update_mmcif_block(query_block)
        
        #update loops that get lost or modified in processing.
        query_block.set_mmcif_category("_struct_conf.", struct_conf)
        query_block.set_mmcif_category("_struct_conf_type.", struct_conf_type)
        chem_comp = chem_comp.replace(False, "?")
        chem_comp_dict = chem_comp.to_dict(orient = "list")
        query_block.set_mmcif_category("_chem_comp.", chem_comp_dict)

        #output the AF structure with transplanted ligands
        query_block.write_file(f"{args.outdir}/{predicted_structure_id}_transplants.cif")
        transplants_df.to_csv(f"{args.outdir}/{predicted_structure_id}_transplants.tsv.gz", sep = "\t", index = False, compression = "gzip")

if __name__ == "__main__":
    main()
