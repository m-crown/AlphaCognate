import pandas as pd

#process config options
if "cognate_match" in config:
    if config["cognate_match"] == True:
        config["cognate_match"] = "--cognate_match"
    else:
        config["cognate_match"] = ""
else:
    config["cognate_match"] = ""

if "domain_match" in config:
    if config["domain_match"] == True:
        config["domain_match"] = "--domain_match"
    else:
        config["domain_match"] = ""
else:
    config["domain_match"] = ""

if config.get("demo"):
    config["demo_manifest"] = "demo_manifests/uniprot_ids.txt"
    config["demo_procoggraph_structures"] = "demo_manifests/target_pdbs.txt"
    config["data_dir"] = "demo_data2"
    config["procoggraph_foldseek_db"] = "demo_data2/pcgDB/pcgDB_demo"
    config["procoggraph_foldseek_db_index"] = "demo_data2/pcgDB/pcgDB_demo_index"
else:
    config["procoggraph_foldseek_db"] = "/pcgDB/pcgDB" #use the name of a full pcgdb
    config["procoggraph_foldseek_db_index"] = "/pcgDB/pcgDB_index"

if not config.get("data_dir"):
    config["data_dir"] = "data"

if config.get("predicted_structures_manifest"):
    with open(config["predicted_structures_manifest"], "r") as f:
        predicted_structures = [line.strip() for line in f.readlines()]
else:
    predicted_structures = []

#it is a requirment for running the combined pipeline that a structure manifest is provided.
#and also a file mapping the domains to this . we will offer a preprocessing step to create this from the input manifest and a cath domain dataset.
structures_manifest = pd.read_csv(config["structure_manifest"], sep=",", names = ["accession", "file_name",  "structure_dir"])
structures = structures_manifest.accession.tolist()

#we need to set the sturcutre dir from here
#maybe check that it is a single directory.
structure_dir = structures_manifest.structure_dir.values[0]

#find out if there is a way to check the content of a foldseek database - as we need to know if it is a full or partial db being checked against !!
#if it is demo mode, set the pcg foldseek name to demo, and likewise the index.

rule all:
    input: 
         all_transplant_cifs = expand(config["output_dir"] + "/transplanted_structures/{id}_transplants.cif.gz", id = structures),
         combined_transplant_tsv = config["output_dir"] + "/combined_transplants.tsv.gz",
         combined_structure_summaries = config["output_dir"] + "/combined_structure_summaries.tsv.gz",

rule aggregate_transplants:
    input:
        all_transplant_files = expand(config["output_dir"] + "/transplanted_structures/{id}_transplants.tsv.gz", id = structures), #require the files but dont use in command (for shell max arg limits)
        all_structure_files = expand(config["output_dir"] + "/transplanted_structures/{id}_structure_summary.tsv.gz", id = structures), #require the files but dont use in command (for shell max arg limits)

    output:
        combined_transplant_tsv = config["output_dir"] + "/combined_transplants.tsv.gz",
        combined_structure_summaries = config["output_dir"] + "/combined_structure_summaries.tsv.gz"
    params:
        output_dir = config["output_dir"],
        tsv_dir = config["output_dir"] + "/transplanted_structures/"
    log: config["output_dir"] + "/logs/aggregate_transplant.log"
    shell:
        """
        python3 bin/combine_tsv.py {params.tsv_dir} {params.output_dir}/combined_transplants.tsv.gz {params.output_dir}/combined_structure_summaries.tsv.gz
        """

rule transplant_ligands:
    input:
        foldseek = config["output_dir"] + "/foldseek_split/{id}_foldseek.tsv.gz"
    output:
        config["output_dir"] + "/transplanted_structures/{id}_transplants.cif.gz",
        temp(config["output_dir"] + "/transplanted_structures/{id}_transplants.tsv.gz"),
        temp(config["output_dir"] + "/transplanted_structures/{id}_structure_summary.tsv.gz")
    params:
        output_dir = config["output_dir"] + "/transplanted_structures",
        procoggraph_foldseek_directory = config["procoggraph_foldseek_directory"],
        cognate_match = config["cognate_match"],
        domain_match_only = config["domain_match"]
    log: config["output_dir"] + "/logs/transplants/{id}_transplant.log"
    shell:
        """python3 bin/alphacognate_transplant.py --foldseek_file {input.foldseek} --outdir {params.output_dir} --structure_database_directory {params.procoggraph_foldseek_directory} {params.cognate_match} {params.domain_match_only} && gzip {params.output_dir}/{wildcards.id}_transplants.cif"""

rule split_foldseek:
    input:
        foldseek_file = config["output_dir"] + "/foldseek/combined_foldseek.tsv",
        structure_manifest = config["structure_manifest"],
        predicted_structure_domains = config["predicted_structure_domains"],
        procoggraph_data = config["data_dir"] + "/cath_single_chain_domain_interactions.tsv.gz",
    output:
        expand(config["output_dir"] + "/foldseek_split/{id}_foldseek.tsv.gz", id = structures)
    params:
        output_dir = config["output_dir"] + "/foldseek_split"
    log: config["output_dir"] + "/logs/foldseek_split.log"
    shell:
      """python3 bin/split_foldseek.py --foldseek_file {input.foldseek_file} --structure_manifest {input.structure_manifest} --output_dir {params.output_dir} --predicted_structure_domains {input.predicted_structure_domains} --procoggraph_data {input.procoggraph_data}"""

rule run_foldseek:
    input:
        structures = expand(structure_dir + "/{id}.cif", id = structures),
        procoggraph_foldseek_db = config["procoggraph_foldseek_db"],
    output:
        foldseek_file = config["output_dir"] + "/foldseek/combined_foldseek.tsv"
    threads: workflow.cores
    params: 
        structure_dir = structure_dir,
        procoggraph_foldseek_db_index =  config["procoggraph_foldseek_db_index"]
    log: config["output_dir"] + "/logs/report.log"
    shell:
        """foldseek easy-search {params.structure_dir} {input.procoggraph_foldseek_db} {output.foldseek_file} {params.procoggraph_foldseek_db_index} --threads {threads} --prefilter-mode 1 --format-mode 4 --format-output query,target,u,t,qlen,tlen,alnlen,mismatch,qaln,qstart,qend,taln,tstart,tend,qseq,tseq,rmsd,qheader,alntmscore,qtmscore,ttmscore,evalue,prob > {log} 2> {log}"""

rule build_foldseek_db:
    input:
        config["data_dir"] + "/procoggraph_assemblies/download_complete.txt",
    output:
        config["procoggraph_foldseek_db"],
        directory(config["procoggraph_foldseek_db_index"])
    params:
        output_dir = config["data_dir"],
        input_dir = config["data_dir"] + "/procoggraph_assemblies",
        foldseek_db = config["procoggraph_foldseek_db"],
        foldseek_db_index = config["procoggraph_foldseek_db_index"]
    shell:
        "bin/build_foldseek_database.sh {params.input_dir} {params.output_dir} {params.foldseek_db} {params.foldseek_db_index}"

rule download_procoggraph_structures:
    input:
        procoggraph_manifest = config["demo_procoggraph_structures"] if config.get("demo") else config["data_dir"] + "/procoggraph_pdb_ids.txt",
        all_procoggraph_ids = config["data_dir"] + "/procoggraph_pdb_ids.txt"
    output:
        directory(config["data_dir"] + "/procoggraph_assemblies"),
        config["data_dir"] + "/procoggraph_assemblies/download_complete.txt"
    params:
        output_dir = config["data_dir"] + "/procoggraph_assemblies"
    shell:
        """
        python3 bin/get_procoggraph_structures.py {input.procoggraph_manifest} {params.output_dir} && gzip {params.output_dir}/*.cif
        """

rule format_procoggraph:
    input:
        config["data_dir"] + "/procoggraph_data/bound_descriptors.tsv.gz",
        config["data_dir"] + "/procoggraph_data/cognate_ligands_df.pkl",
        config["data_dir"] + "/procoggraph_data/all_parity_calcs.pkl",
        config["data_dir"] + "/procoggraph_data/cath_pdb_residue_interactions.csv.gz"
    output:
        config["data_dir"] + "/cath_single_chain_domain_interactions.tsv.gz",
        config["data_dir"] + "/procoggraph_pdb_ids.txt"
    params:
        data_dir = config["data_dir"]
    shell:
        """
        python3 bin/format_procoggraph.py --cath_domain_ownership {params.data_dir}/procoggraph_data/cath_pdb_residue_interactions.csv.gz --scores_file {params.data_dir}/procoggraph_data/all_parity_calcs.pkl --cognate_ligands {params.data_dir}/procoggraph_data/cognate_ligands_df.pkl --bound_entity_descriptors {params.data_dir}/procoggraph_data/bound_descriptors.tsv.gz --output_dir {params.data_dir}
        """

rule download_procoggraph_data:
    output:
        config["data_dir"] + "/procoggraph_data/bound_descriptors.tsv.gz",
        config["data_dir"] + "/procoggraph_data/cognate_ligands_df.pkl",
        config["data_dir"] + "/procoggraph_data/all_parity_calcs.pkl",
        config["data_dir"] + "/procoggraph_data/cath_pdb_residue_interactions.csv.gz"
    params:
        data_dir = config["data_dir"]
    shell:
        """
        wget https://zenodo.org/records/15204472/files/alphacognate_data_files.zip?download=1 -O {params.data_dir}/alphacognate_data_files.zip
        unzip {params.data_dir}/alphacognate_data_files.zip -d {params.data_dir}/procoggraph_data
        """

#set default outdir to . for this and for get structures to alphafold

# rule make_manifest:
#     input:
#         config["data_dir"] + "/predicted_structures_download_complete.txt"
#     output:
#         config["data_dir"] + "/alphafold_structures_manifest.csv"
#     params:
#         structure_dir = config["data_dir"] + "/predicted_structures/",
#         output_dir = config["data_dir"]

#     shell:
#         """
#         python3 bin/make_manifest.py {params.structure_dir} {params.output_dir}
#         """

##this one would then live solely as a preprocessing step 

##do this to download demo if demo mode is enabled - add a manifest for the demo structures to demo_manifests.
#think we can get rid of the downlaod complete file here as the check should work on expanded id output
if config.get("demo"):
    rule download_predicted_structures:
        input:
            predicted_structure_manifest = config["demo_structures"]
        output:
            expand(config["data_dir"] + "/predicted_structures/AF-{id}-F1-model_4.cif", id = structures),
            download_complete_file = config["data_dir"] + "/predicted_structures_download_complete.txt"

        params:
            output_dir = config["data_dir"] + "/predicted_structures"
        shell:
            """
            python3 bin/get_alphafold_structures.py {input.demo_structures} {params.output_dir} {output.download_complete_file}
            """