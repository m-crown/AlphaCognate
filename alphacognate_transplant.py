import requests
import json
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

    query_chain_span = query_struct[0].get_subchain(query_chain)
    target_chain_span = target_struct[0].get_subchain(target_chain) #we use get subchain to get the correct struct asym for the protein sequence (observer issues with using index directly for multiple chains with same name e.g. 4kvm chain H - matches ligand instead of protein)
    
    qaln = foldseek_result.qaln
    taln = foldseek_result.taln

    #loop over foldseek alignments together
    for q_res, t_res in zip(qaln, taln):
        if q_res != "-" and t_res != "-":
            #we could add more checks in here for example ensuring the one letter code of the residue to be extracted is the same as reported in the alignment.
            if target_idx <= len(target_chain_span.make_one_letter_sequence()) - 1: #we observe with 1xrr that the diproline bound ligands are included in the alignment at the end - this is incorrect, and results in a tseq that is longer than the target sequence in the mmcif structure. so prevent the indexes going above the length of the target sequence which would result in an index error.
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
        print("Failed to match all ligand residues, check chain mappings!")
        return None
    
    contacts = cs.find_contacts(ns)
    ca_pos = []

    backbone_atoms = ['CA']
    contact_res_list = []
    for res in target_alignment_range_atom:
        if any(x in [contact.partner1.atom for contact in contacts] for x in res["atoms"]):
            contact_res_list.append(res["res"].seqid.num) #get a list of contacting residues to check correct domains are in proximity
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

def find_and_download_af_prediction(identifier, outdir = "."):
    #check output dir exists
    outdir_path = Path(outdir)
    outdir_path.mkdir(parents=True, exist_ok=True)
    #check a file containing the identiifer itself exists
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
    python3 alphacognate_transplant.py --af_dir /raid/MattC/repos/alpha_cognate/alphafold_structures --outdir /raid/MattC/repos/alpha_cognate/results --structure_database_directory /raid/MattC/repos/foldseek/pcg_assemblies/ --procoggraph_data /raid/MattC/repos/ProCogGraphData/procoggraph_20240528/flat_files/cath_single_chain_domain_interactions.tsv.gz --cath_af_grouped /raid/MattC/repos/alpha_cognate/cath_alphafold_domains.tsv.gz
    """
    parser = argparse.ArgumentParser(description='Transplant cognate ligands from experimentally determined structures to an AlphaFold model.')
    parser.add_argument('--af_dir', type=str, help='Path to the AlphaFold model structure directory (missing structures downloaded here). Created if missing')
    parser.add_argument('--outdir', type=str, help='Path to the output directory. Created if missing')
    parser.add_argument('--structure_database_directory', type=str, help='Path to the directory containing the structures to transplant ligands from.')
    parser.add_argument('--procoggraph_data', type=str, help='Path to the ProCogGraph data in TSV format.')
    parser.add_argument('--cath_af_grouped', type=str, help='Path to the grouped domain cath af info.')
    args = parser.parse_args()

    #check directories exist
    Path(f"{args.outdir}/foldseek").mkdir(parents=True, exist_ok=True)
    Path(args.af_dir).mkdir(parents=True, exist_ok=True)
    Path(f"{args.outdir}/transplanted_structures").mkdir(parents=True, exist_ok=True)

    #check if foldseek database directory exists - exit with error if not
    if not Path(args.structure_database_directory).exists():
        print(f"Error: {args.structure_database_directory} does not exist.")
        return

    #from foldseek github c++ code
    threeAA2oneAA = {"ALA":'A', "ARG":'R', "ASN":'N', "ASP":'D', "CYS":'C', "GLN":'Q', "GLU":'E', "GLY":'G', "HIS":'H', "ILE":'I', "LEU":'L', "LYS":'K', "MET":'M', "PHE":'F', "PRO":'P', "SER":'S', "THR":'T', "TRP":'W', "TYR":'Y', "VAL":'V', "MSE":'M', "MLY":'K', "FME":'M', "HYP":'P',"TPO":'T', "CSO":'C', "SEP":'S', "M3L":'K',"HSK":'H', "SAC":'S', "PCA":'E', "DAL":'A',"CME":'C', "CSD":'C', "OCS":'C', "DPR":'P',"B3K":'K', "ALY":'K', "YCM":'C', "MLZ":'K',"4BF":'Y', "KCX":'K', "B3E":'E', "B3D":'D',"HZP":'P', "CSX":'C', "BAL":'A', "HIC":'H',"DBZ":'A', "DCY":'C', "DVA":'V', "NLE":'L',"SMC":'C', "AGM":'R', "B3A":'A', "DAS":'D',"DLY":'K', "DSN":'S', "DTH":'T', "GL3":'G',"HY3":'P', "LLP":'K', "MGN":'Q', "MHS":'H',"TRQ":'W', "B3Y":'Y', "PHI":'F', "PTR":'Y',"TYS":'Y', "IAS":'D', "GPL":'K', "KYN":'W',"SEC":'C'}

    #foldseek results match the af structure to the pre-built pdb database. 
    # We then merege this with procoggraph data to get the structures to transplant
    #merge procoggraph data and foldseek data. then we do a foldseek itterrows but we should loop over the structures, where there could be multiple ligands to transplant in the same structure.
    
    cath_af_grouped = pd.read_csv(args.cath_af_grouped, sep="\t")
    cath_af_grouped["domain"] = cath_af_grouped["domain"].str.split(";")
    cath_af_grouped["domain_profile"] = cath_af_grouped["domain_profile"].str.split(";")
    procoggraph = pd.read_csv(args.procoggraph_data, sep="\t")
    procoggraph["bound_entity_pdb_residues"] = procoggraph["bound_entity_pdb_residues"].str.split("|").apply(lambda x : [int(y) for y in x])
    procoggraph["bound_entity_chain"] = procoggraph["uniqueID"].str.split("_").apply(lambda x: x[2])


    all_transplants = []
    for _, row1 in cath_af_grouped[0:1000].iterrows():
        start_time = time.time()
        transplants = []
        error = None
        error, af_structure_file = find_and_download_af_prediction(row1.af_accession, args.af_dir)
        if error:
            result = {"query": row1.af_accession, "error": error}
            transplants.append(result)
        foldseek_file = f"{args.outdir}/foldseek/aln_pcg_{row1.af_accession}"
        command = [
            "foldseek", "easy-search",
            af_structure_file, 
            "/raid/MattC/repos/foldseek/pcgDB", 
            foldseek_file, 
            "/raid/MattC/repos/foldseek/pcgDB_index", '--prefilter-mode' , '1',
            "--format-mode", "4",
            "--format-output", "query,target,u,t,qlen,tlen,alnlen,mismatch,qaln,qstart,qend,alnlen,taln,tstart,tend,qseq,tseq,rmsd,qheader,alntmscore,qtmscore,ttmscore,evalue,prob"
        ]
        
        log_file = f"{args.outdir}/foldseek/{row1.af_accession}_foldseek_output.txt"

        # Run the command using subprocess and capture the output
        if not Path(foldseek_file).exists():
            result = subprocess.run(command, text=True, capture_output=True)

            # Save the stdout to the specified file
            with open(log_file, 'w') as f:
                f.write(result.stdout)

            # Check for errors in stderr
            if result.stderr:
                if not "posix_madvise returned an error (touchMemory)\n" == result.stderr:
                    error = f"Non-tolerated error for {row1['domain'][0]}:\n{result.stderr}"
                    result = {"query": row1.af_accession, "error": error}
                    transplants.append(result)
                    continue  # Stop execution on any other error
                    
        foldseek = pd.read_csv(foldseek_file, sep = "\t")
        foldseek["pdb_id"] = foldseek["target"].str.extract("^([A-Za-z0-9]+)_bio-h")
        foldseek["target_chain"] = foldseek["target"].str.extract("^[A-Za-z0-9]+_bio-h_([A-Za-z0-9]+)")
        foldseek["target_file"] = foldseek["target"].str.extract("^([A-Za-z0-9]+_bio-h)")
        foldseek["query_chain"] = "A" #this is always the same for the alphafold model?
        foldseek["af_accession"] = foldseek["query"].str.extract("AF-([A-Za-z0-9]+)-")
        foldseek_cath = foldseek.merge(cath_af_grouped, on = "af_accession", how = "left")
        foldseek_cath_filtered = foldseek_cath.loc[(foldseek_cath.prob >= 0.5) & (foldseek_cath.alntmscore >= 0.5)].copy()
        if foldseek_cath_filtered is None:
            error = "No structures with high enough alignment scores to transplant from"
            result = {"query": row1.af_accession, "error": error}
            continue

        foldseek_merge = foldseek_cath_filtered.merge(procoggraph, left_on = ["pdb_id", "target_chain"], right_on = ["pdb_id", "assembly_chain_id_protein"], how = "inner")

        if foldseek_merge.empty:
            error = "No structures to transplant ligands from for AF structure"
            result = {"query": row1.af_accession, "error": error}
            continue

        foldseek_merge["interaction_domains"] = foldseek_merge["combined_interaction"].str.split(";")
        foldseek_merge["interaction_domains"] = foldseek_merge["interaction_domains"].apply(lambda x : [y.split(":")[0] for y in x])
        #keep only the matches to structures where all domains in the interaction domains are present in the domain profile
        foldseek_merge_domains = foldseek_merge.loc[foldseek_merge.apply(lambda row: all(domain in row['domain'] for domain in row['interaction_domains']), axis=1)]
        if foldseek_merge_domains.empty:
            error = "No structures with matching domain profiles to transplant from"
            result = {"query": row1.af_accession, "error": error}
            continue
            
        num_transplants = len(foldseek_merge_domains)
        transplant_chain_ids = generate_chain_ids(num_transplants)

        #load the af structure to which ligands will be transplanted
        af_structure = cif.read(af_structure_file)
        query_block = af_structure.sole_block()
        query_struct = gemmi.make_structure_from_block(query_block)
        if not Path(f"{args.outdir}/transplanted_structures/{row1.af_accession}_transplants.tsv", sep = "\t").exists():
            for index, row in foldseek_merge_domains.iterrows():
                result = {"num_transplants": num_transplants}
                query = row["query"]
                target = row["target"]
                target_chain = row["proteinStructAsymID"]
                #target_chain_asym = row[""]
                query_chain = row["query_chain"]

                target_doc = cif.read(f"{args.structure_database_directory}/{row.target_file}.cif.gz")
                target_block = target_doc.sole_block()
                target_struct = gemmi.make_structure_from_block(target_block)
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
                    result.update({"query_structure": query, "transplant_structure": target, "foldseek_rmsd": row.rmsd, "global_rmsd": global_rmsd, "ligand": row.uniqueID, 'hetCode': row.hetCode, 'cognateLigand': row.cognate_ligand, "interaction": row.combined_interaction, "error": error})
                    transplants.append(result)
                    continue
                
                #check the domain profile of the ligand contacts matches the procoggraph mapping
                #the procoggraph map boolean value describes the match between domain contacts in af structure and expected interacting domains from coglig match. but we do not filter on this, leave it up to end user.
                residue_mappings, interacting_domains, procoggraph_map = check_domain_profile(row.domain_profile, contact_res_list, row.interaction_domains)


                target_struct[0].transform_pos_and_adp(local_superposition.transform)
                
                transplant_chain_id = transplant_chain_ids.pop(0)
                transplant_chain = create_transplant_chain(target_struct, res_chain, res_list, transplant_chain_id)
            
                query_struct[0].add_chain(transplant_chain)
                query_struct.update_mmcif_block(query_block)
                
                #explore if determine tcs could run before doing the actual transpalnt, and only pop a new chain id if the tcs is better than any pre-existing transplant of the same bound entity - cognate ligand mapping?
                tcs = determine_tcs(query_struct, transplant_chain_id)


                result.update({"query_structure": query, "transplant_structure": target, "foldseek_rmsd": row.rmsd, "global_rmsd": global_rmsd, "ligand": row.uniqueID, 'hetCode': row.hetCode, 'cognateLigand': row.cognate_ligand, "interaction": row.combined_interaction, "local_rmsd": local_rmsd, "tcs": tcs, "transplanted_structure_path": f"{args.outdir}/transplanted_structures/{row1.af_accession}_cognate_transplants.cif.gz", "transplanted_ligand_chain": transplant_chain_id, "transplanted_ligand_residues": ",".join([str(x) for x in res_list]), "domain_residue_contacts": residue_mappings, "domain_residue_counts": interacting_domains, "domain_profile_match": procoggraph_map})
                transplants.append(result)
                
            #output the AF structure with transplanted ligands
            query_block.write_file(f"{args.outdir}/transplanted_structures/{row1.af_accession}_cognate_transplants.cif")
            #gzip the output file
            subprocess.run(["gzip", f"{args.outdir}/transplanted_structures/{row1.af_accession}_cognate_transplants.cif"])
            #add transplant info to the global df
            transplants_df = pd.DataFrame(transplants)
            transplants_df["runtime"] = time.time() - start_time
            transplants_df.to_csv(f"{args.outdir}/transplanted_structures/{row1.af_accession}_transplants.tsv", sep = "\t", index = False)
        else:
            transplants_df = pd.read_csv(f"{args.outdir}/transplanted_structures/{row1.af_accession}_transplants.tsv", sep = "\t")
            transplants_df["pre-calculated"] = True
        all_transplants.append(transplants_df)
        

    all_transplants_df = pd.concat(all_transplants, ignore_index = True)
    all_transplants_df.to_csv(f"{args.outdir}/all_transplants.tsv", sep = "\t", index = False)
    print(all_transplants_df)

    
    ##we also need to find a way to match the binding pocket residues to domains and confirm that the AF binding pocket matches the domains we expect to be in context.
    ##MAybe we can do a domain rmsd calculation for each of the domains in the transplant structure, to see if they are oriented similarly in the AF sturcture

if __name__ == "__main__":
    main()





