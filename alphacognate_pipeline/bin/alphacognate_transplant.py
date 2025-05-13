from pathlib import Path
import pandas as pd
import gemmi
import math
from gemmi import cif
from operator import itemgetter
import argparse
import time
import itertools
import string
import numpy as np
from typing import Optional, List
from pydantic import BaseModel


#TODO: Consider column choices for transplants.
#TODO: Make errors a defined enum of values.

class CognateMapping(BaseModel):
    cognate_mapping_name: Optional[str] = None
    cognate_mapping_smiles: Optional[str] = None
    cognate_mapping_xref: Optional[str] = None

class TransplantResult(BaseModel):
    accession: Optional[str] = None
    transplant_structure: Optional[str] = None
    foldseek_rmsd: Optional[float] = None
    global_rmsd: Optional[float] = None
    local_rmsd: Optional[float] = None
    ligand: Optional[str] = None
    ligand_het_code: Optional[str] = None
    ligand_name : Optional[str] = None
    ligand_chain: Optional[str] = None
    ligand_residues: Optional[str] = None
    ligand_center_of_mass: Optional[str] = None
    ligand_smiles: Optional[str] = None
    domain_residue_contacts: Optional[str] = None
    domain_residue_counts: Optional[str] = None
    domain_profile_match: Optional[bool] = None
    tcs: Optional[float] = None
    transplant_error: Optional[str] = None

class AlphaCognateStructure(BaseModel):
    accession: str
    runtime: Optional[float] = None
    num_transplants: Optional[int] = None
    num_clusters: Optional[int] = None
    error: Optional[str] = None

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
        #observe issue with ac1 and glc in 3poc where two residues have same seqid - need to check this in the future also 1b5q
        print(f"{target_struct}, {target_chain}, {reslist}, {matched_res}")
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
        #Require a minimum of 3 residues in contact for a valid match - should extend the error message here to be more informative.
        #print(f"Not enough residues in contact with the ligand, only {len(mol_range_index_positions)} residues found.")
        #make this a log thing in future
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

def check_domain_profile(af_domain_profiles, residues_in_contact=0, procoggraph_profile=0):
    domain_range_dict = {}
    if af_domain_profiles == ['']: #sometimes the uniprot id does not have a domain mapping in ted or cath alphafold
        return "", "", False
    
    for i, domain_profile in enumerate(af_domain_profiles):
        domain, res = domain_profile.split(":")
        segments = res.split("_")
        ranges = []
        for segment in segments:
            res_start, res_end = segment.split("-")
            ranges.append((int(res_start), int(res_end)))
        domain_range_dict[f"{domain}_{i}"] = ranges

    residue_mappings = {domain: [] for domain in domain_range_dict}
    residue_mappings.setdefault('no-domain_-1', [])
    
    for res in residues_in_contact:
        found = False
        for domain, ranges in domain_range_dict.items():
            for start, end in ranges:
                if start <= res <= end:
                    residue_mappings[domain].append(res)
                    found = True
                    break
            if found:
                break
        if not found:
            residue_mappings['no-domain_-1'].append(res)
    
    interacting_domains = [(domain.split("_")[0], domain.split("_")[1], len(residues)) for domain, residues in residue_mappings.items() if residues and domain != 'no-domain_-1']
    procoggraph_map = sorted([domain for domain, idx, count in interacting_domains]) == sorted(procoggraph_profile)
    
    interacting_domains = ";".join([":".join([str(val) for val in tup]) for tup in interacting_domains])
    residue_mappings = ";".join([":".join([key, ",".join(map(str, value))]) for key, value in residue_mappings.items() if value])
    
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
    parser.add_argument('--cognate_match', action='store_true', help='Only transplant ligands with cognate ligand mappings (from ProCogGraph database).')
    parser.add_argument('--domain_match', action='store_true', help='Only transplant ligands with matching interactiosn to the ligand mapping (from ProCogGraph database).')
    args = parser.parse_args()

    #check directories exist
    Path(f"{args.outdir}").mkdir(parents=True, exist_ok=True)

    #check if foldseek database directory exists - exit with error if not
    if not Path(args.structure_database_directory).exists():
        print(f"Error: {args.structure_database_directory} does not exist.")
        return

    #merge procoggraph data and foldseek data. then we do a foldseek itterrows but we should loop over the structures, where there could be multiple ligands to transplant in the same structure.
    error = None

    foldseek_file = pd.read_csv(args.foldseek_file, sep="\t", na_values = ["NaN", "None"], keep_default_na = False)
    if "error" in foldseek_file.columns:
        error = foldseek_file.error.values[0]

    foldseek_file["fp"] = foldseek_file["structure_dir"] + "/" + foldseek_file["file_name"]
    
    predicted_structure_id = foldseek_file["accession"].values[0]
    predicted_structure_file = foldseek_file["fp"].values[0]

    if args.cognate_match and not error:
        foldseek_file = foldseek_file[(foldseek_file["cognate_mapping_id"].isna() == False) & (foldseek_file.cognate_mapping_id != "")]
        if len(foldseek_file) == 0:
            error = "No cognate ligands found in foldseek file."

    #load the af structure to which ligands will be transplanted
    af_structure = cif.read(predicted_structure_file)
    query_block = af_structure.sole_block()

    if error: ##eventually want to make this more robust - all foldseek file columns should be the same, the test would become if all the values in foldseek file .error are not na perhaps.

        #set the structure error which is the foldseek file error.
        alphacognate_structure = AlphaCognateStructure(
            accession = predicted_structure_id,
            error = error,
            runtime = 0,
            num_transplants = 0
            )
        
        # Create an empty file with no transplants (predominantly because snakemake is expecting it)
        Path(f"{args.outdir}/{predicted_structure_id}_transplants.tsv.gz").touch()
        af_structure = cif.read(predicted_structure_file)
        query_block = af_structure.sole_block()
        alphacognate_structure_dict = {k:[v] for k,v in alphacognate_structure.model_dump().items()}
        alphacognate_structure_df = pd.DataFrame(alphacognate_structure_dict)
        alphacognate_structure_df.to_csv(f"{args.outdir}/{predicted_structure_id}_structure_summary.tsv.gz", sep = "\t", index = False, compression = "gzip")
        query_block.set_mmcif_category("_alphacognate_structure", alphacognate_structure_dict) #add alphacognate_structure loop with information on the structure overall.
        query_block.write_file(f"{args.outdir}/{predicted_structure_id}_transplants.cif")
        return
    
    #continuing processing foldseek results without errors.
    foldseek_file["domain_profile"] = foldseek_file["domain_profile"].str.split(";")
    foldseek_file["domain"] = foldseek_file["domain_profile"].apply(lambda x: [y.split(":")[0] for y in x])
    #convert bound entity pdb residues to list of ints - force to be a string first in case of places whrre there are no pipe delimiters in files.
    foldseek_file["bound_entity_pdb_residues"] = foldseek_file["bound_entity_pdb_residues"].astype("str").str.split("|").apply(lambda x : [int(y) for y in x])
    foldseek_file["interaction_domains"] = foldseek_file["combined_interaction"].str.split(";")
    foldseek_file["interaction_domains"] = foldseek_file["interaction_domains"].apply(lambda x : [y.split(":")[0] for y in x])
    
    num_potential_transplants = foldseek_file.ligand.nunique() #its the number of ligands, not number of ligands mapped to cogligs that are transplanted
    transplant_chain_ids = generate_chain_ids(num_potential_transplants)

    
    #for now, this is set after handling foldseek errors, but would be good to consolidate this eventually
    transplants: list[TransplantResult] = [] #to store the transplants as we make them
    start_time = time.time() #start the timer for overall transplants.

    #store some of the loops that get lost in processing to be handled later.
    query_struct = gemmi.make_structure_from_block(query_block)
    query_struct.merge_chain_parts()
    chem_comp = pd.DataFrame(query_block.get_mmcif_category('_chem_comp.'))
    struct_conf = query_block.get_mmcif_category('_struct_conf.')
    struct_conf_type = query_block.get_mmcif_category('_struct_conf_type.')

    #we make a subset of the deduplicated bound entity specific information needed for transplanting and iterate through this.
    #post transplanting, we merge back the cognate information.
    foldseek_transplants = foldseek_file[["target", "assembly_chain_id_protein", "query_chain", "target_file", "ligand", "bound_entity_pdb_residues", "assembly_chain_id_ligand", "combined_interaction", "domain_profile", "interaction_domains", "rmsd", "qstart","tstart","qaln","taln", "ligand_het_code", "ligand_smiles", "ligand_name"]].drop_duplicates(
        subset = ["target", "ligand", "assembly_chain_id_protein", "query_chain", "assembly_chain_id_ligand"])
    #and the information we merge after transplant, on ligand
    foldseek_file["cognate_mapping_name"] = foldseek_file.cognate_mapping_name.apply(lambda x: x.split("|")[0] if isinstance(x, str) else np.nan)
    #rename some columns here
    foldseek_other = foldseek_file[["ligand", "cognate_mapping_id", "cognate_mapping_name", "cognate_mapping_ec_list", "cognate_mapping_similarity", "cognate_mapping_smiles", "cognate_mapping_xref"]]

    
    if not Path(f"{args.outdir}/{predicted_structure_id}_transplants.tsv.gz", sep = "\t").exists():
        for index, row in foldseek_transplants.iterrows():
            target = row["target"]
            target_chain = row["assembly_chain_id_protein"]
            query_chain = row["query_chain"]
            #get the residues in the assembly that make up the ligand (and the chain they are in)
            res_list = row.bound_entity_pdb_residues
            res_chain = row.assembly_chain_id_ligand
            #we get a new transplant chain id here, and if it is not used (because there are not enough local contacts)
            transplant_chain_id = transplant_chain_ids.pop(0)

            target_doc = cif.read(f"{args.structure_database_directory}/{row.target_file}.cif.gz")
            target_block = target_doc.sole_block()
            target_struct = gemmi.make_structure_from_block(target_block)
            #merge chain parts to help with chain naming issue, use auth ids and get_polymer to identify protein chains
            target_struct.merge_chain_parts()
            query_alignment_range, target_alignment_range, query_alignment_range_atom, target_alignment_range_atom = get_alignment_positions(row, query_struct, query_chain, target_struct, target_chain)

            #superpose the query structure onto the target structure, using the alignment from foldseek result (the global alignment)
            sup = gemmi.superpose_positions(query_alignment_range, target_alignment_range)
            global_rmsd = sup.rmsd

            #align the target structure onto the query structure
            target_struct[0].transform_pos_and_adp(sup.transform)

            #identify the contacts to the ligand that should be used to perform the 'local' alignment.
            local_superposition, contact_res_list = contact_search_target(target_struct, res_chain, res_list, 6, query_alignment_range, target_alignment_range, query_alignment_range_atom, target_alignment_range_atom)
            ##we can add a value of matched domain profiles? by checking the residue ranges?
            if local_superposition is None:
                # TODO: For now, do not keep failed transplants as it is unnecessary?
                # #initiate a failed transplant result
                # transplantresult = TransplantResult(
                #     accession = predicted_structure_id,
                #     transplant_structure = target,
                #     foldseek_rmsd = row.rmsd,
                #     global_rmsd = global_rmsd,
                #     ligand = row.ligand,
                #     ligand_chain = None
                #     error = "Not enough local contacts"
                #     )
                
                # transplants.append(transplantresult)
                #return the unused transplant chain to the chain list.
                transplant_chain_ids.insert(0, transplant_chain_id)
                continue
            

            #check the domain profile of the ligand contacts matches the procoggraph mapping
            #the procoggraph map boolean value describes the match between domain contacts in af structure and expected interacting domains from coglig match. but we do not filter on this, leave it up to end user.
            residue_mappings, interacting_domains, procoggraph_map = check_domain_profile(row.domain_profile, contact_res_list, row.interaction_domains)
            if args.domain_match and not procoggraph_map:
                #if no domain match is found, we do not want to transplant the ligand.
                transplant_chain_ids.insert(0, transplant_chain_id)
                continue

            target_struct[0].transform_pos_and_adp(local_superposition.transform)
            
            transplant_chain = create_transplant_chain(target_struct, res_chain, res_list, transplant_chain_id)
            transplant_chain_center_of_mass = transplant_chain.calculate_center_of_mass().tolist()

            #TODO: How necessary is this if else statement? Shouldnt 
            if all([pd.isna(x) for x in transplant_chain_center_of_mass]):
                transplant_chain_center_of_mass = ""
            else:
                transplant_chain_center_of_mass = ",".join([str(x) for x in transplant_chain_center_of_mass])
            
            #get the transplant res names to update the chem_comp_table.
            transplant_res_names = [res.name for res in transplant_chain]
            transplant_chem_comp = pd.DataFrame(target_block.get_mmcif_category("_chem_comp."))
            transplant_chem_comp_filtered = transplant_chem_comp.loc[transplant_chem_comp.id.isin(transplant_res_names)]

            #add any new chem_comp entries from the transplant to the overall chem_comp df
            chem_comp = pd.concat([chem_comp, transplant_chem_comp_filtered]).drop_duplicates(subset="id", keep="first").reset_index(drop=True)
            
            query_struct[0].add_chain(transplant_chain)

            #TODO: explore if determine tcs could run before doing the actual transpalnt, and only pop a new chain id if the tcs is better than any pre-existing transplant of the same bound entity - cognate ligand mapping?
            tcs = determine_tcs(query_struct, transplant_chain_id)

            #get the local rmsd from the superposition previously performed
            local_rmsd = local_superposition.rmsd

            #initiate a successful transplant result.
            transplantresult = TransplantResult(
                accession = predicted_structure_id,
                transplant_structure = target,
                foldseek_rmsd = row.rmsd,
                global_rmsd = global_rmsd,
                local_rmsd = local_rmsd,
                ligand = row.ligand,
                ligand_het_code = str(row.ligand_het_code),
                ligand_name = row.ligand_name,
                ligand_smiles = row.ligand_smiles,
                interaction =  row.combined_interaction,
                ligand_chain = transplant_chain_id,
                ligand_residues = ",".join([str(x) for x in res_list]),
                tcs = tcs,
                domain_residue_contacts = residue_mappings,
                domain_residue_counts = interacting_domains,
                domain_profile_match = procoggraph_map,
                ligand_center_of_mass = transplant_chain_center_of_mass
                )
            transplants.append(transplantresult)

        #sometimes no transplants with local contacts will be present i.e. 0 len transplants
        if len(transplants) > 0:
            #convert individual transplant dictionaries into a single dataframe
            #need to consider how to do model dump with list of coglig mappings in future.
            transplants_df = pd.DataFrame([s.model_dump() for s in transplants])
            transplants_df = transplants_df.merge(foldseek_other, on = "ligand", how = "left")
            #calculate cluster membership with default hdbscan parameters for transplants
            if not (transplants_df["ligand_center_of_mass"] == "").all(): #again, is this necessary?
                #this replacement is necessary because we set center of mass as string in transplant model. sometimes gemmi returns np.nan as the center of mass - investigate this.
                #TODO : This has multiple records for same com due to cognate ligand mappings. should remove and merge.
                transplants_df["center_of_mass_split"] = transplants_df["ligand_center_of_mass"].replace("", np.nan).str.split(",")
                points = transplants_df.loc[(transplants_df.center_of_mass_split.isna() == False), "center_of_mass_split"].apply(lambda x: [float(y) for y in x]).values
                if len(points) < 5:
                    transplants_df.loc[(transplants_df.center_of_mass_split.isna() == False), "cluster"] = -1
                else:
                    #discuss this in a disccusion section.
                    points = np.array([np.array(point) for point in points])
                    
                    clustering_method = "agglomerative" #default clustering method - for now not configurable.
                    if clustering_method == "agglomerative":
                        #use agglomerative clustering to cluster the points
                        from sklearn.cluster import AgglomerativeClustering
                        clustering = AgglomerativeClustering(n_clusters=None, distance_threshold=10, linkage = "average").fit(points)
                    elif clustering_method == "hdbscan":
                        #use hdbscan to cluster the points
                        from sklearn.cluster import HDBSCAN
                        if len(points) < 5:
                            labels = [-1] * len(points)
                        else:
                            params = {"min_cluster_size": 5, "min_samples": 5} #specifying the default values here for future configuration by user
                            clustering = HDBSCAN(**params).fit(points)

                    #add cluster labels to the dataframe - clusters with fewer than 5 members will be labelled as noise and given -1 as a cluster value.
                    labels = clustering.labels_
                    transplants_df.loc[(transplants_df.center_of_mass_split.isna() == False), "cluster"] = labels
                    transplants_df.loc[(transplants_df.center_of_mass_split.isna()), "cluster"] = -1
            else:
                transplants_df["cluster"] = 0
            
            #use the assigned cluster labels to generate a cluster center for each cluster.

            cluster_centers = transplants_df[["ligand_chain", "center_of_mass_split", "cluster"]].drop_duplicates(subset = "ligand_chain")
            cluster_centers["center_of_mass_split"] = cluster_centers["center_of_mass_split"].apply(lambda x: [float(y) for y in x])
            cluster_centers = (cluster_centers.groupby('cluster')['center_of_mass_split']
                .apply(lambda coords: np.mean(np.vstack(coords), axis=0))
                .reset_index(name='cluster_center')
            )
            #merge the cluster centers to the df
            transplants_df = transplants_df.merge(cluster_centers, on = "cluster", how = "left")
            transplants_df["cluster_center"] = transplants_df["cluster_center"].apply(lambda x: ",".join([str(y) for y in x]))
            #remove the split center of mass split column and the num transpalnts from the mmcif - it is implicit in the number of rows.
            transplants_df = transplants_df.drop(columns = ["center_of_mass_split"])
            #save to file incl. all cognate mapping
            transplants_df.to_csv(f"{args.outdir}/{predicted_structure_id}_transplants.tsv.gz", sep = "\t", index = False, compression = "gzip")

            #split out the cognate info to separate table.
            cognate_mapping_df = transplants_df[["ligand", "cognate_mapping_name", "cognate_mapping_smiles", "cognate_mapping_xref", "cognate_mapping_ec_list"]].drop_duplicates()
            
            # Keep rows where at least one of the cognate columns is not empty
            mask = (cognate_mapping_df[["cognate_mapping_name", "cognate_mapping_smiles", "cognate_mapping_xref"]]
                    .replace("", np.nan)
                    .notna()
                    .any(axis=1))

            # Apply the mask
            cognate_mapping_df = cognate_mapping_df[mask]

            cognate_mapping_dict = cognate_mapping_df.fillna("").astype("str").to_dict(orient = "list")

            #create a dictionary of the transplant data to add as an mmcif category. can remove duplicates due to multiple coglig mappings here.
            transplants_df = transplants_df[[col for col in transplants_df.columns if col not in ["accession", "cognate_mapping_id", "cognate_mapping_name", "cognate_mapping_smiles", "cognate_mapping_xref", "cognate_mapping_ec_list", "cognate_mapping_similarity"]]].drop_duplicates()
            transplant_dictionary = transplants_df.fillna("").astype("str").to_dict(orient = "list")
            #add the transplant dictionary as an mmcif category to the query block
            query_block.set_mmcif_category("_alphacognate_transplants", transplant_dictionary)
            #add the cognate mapping dictionary as an mmcif category to the query block
            query_block.set_mmcif_category("_alphacognate_cognate_mapping", cognate_mapping_dict)
            
            #we have added new entities. Need to make sure theyre present in the _entity. loop.
            query_struct.assign_subchains(force = True)
            query_struct.ensure_entities()
            #update the query block with the new structure.
            query_struct.update_mmcif_block(query_block)
            
            #manually revert original chain struct asym to A.
            struct_asym = pd.DataFrame(query_block.get_mmcif_category('_struct_asym.'))
            struct_asym.loc[struct_asym['entity_id'] == "1", 'id'] = 'A'
            atom_site = pd.DataFrame(query_block.get_mmcif_category('_atom_site.'))
            atom_site.loc[atom_site['label_entity_id'] == '1', 'label_asym_id'] = 'A'
            #update loops that get lost or modified in processing.
            query_block.set_mmcif_category("_atom_site.", atom_site.to_dict(orient = "list"))
            query_block.set_mmcif_category("_struct_asym.", struct_asym.to_dict(orient = "list"))
            query_block.set_mmcif_category("_struct_conf.", struct_conf)
            query_block.set_mmcif_category("_struct_conf_type.", struct_conf_type)
            chem_comp = chem_comp.replace(False, "?")
            chem_comp_dict = chem_comp.to_dict(orient = "list")
            query_block.set_mmcif_category("_chem_comp.", chem_comp_dict)

        else:
            # Create an empty file with no transplants (predominantly because snakemake is expecting it)
            Path(f"{args.outdir}/{predicted_structure_id}_transplants.tsv.gz").touch()
        
        alphacognate_structure = AlphaCognateStructure(
            accession = predicted_structure_id,
            runtime = time.time() - start_time,
            num_transplants = len(transplants),
            num_clusters = transplants_df["cluster"].nunique()
        )

        #tables need to be keys with lists
        alphacognate_structure_dict = {k:[v] for k,v in alphacognate_structure.model_dump().items()}
        alphacognate_structure_df = pd.DataFrame(alphacognate_structure_dict)
        query_block.set_mmcif_category("_alphacognate_structure", alphacognate_structure_dict) #add alphacognate_structure loop with information on the structure overall.

        #output the AF structure with transplanted ligands
        query_block.write_file(f"{args.outdir}/{predicted_structure_id}_transplants.cif")
        alphacognate_structure_df.to_csv(f"{args.outdir}/{predicted_structure_id}_structure_summary.tsv.gz", sep = "\t", index = False, compression = "gzip")

if __name__ == "__main__":
    main()