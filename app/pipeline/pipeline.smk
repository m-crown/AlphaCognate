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

structures = pd.read_csv(config["structure_manifest"], sep=",", names = ["accession", "file_name",  "structure_dir"]).accession.tolist()

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
        python3 combine_tsv.py {params.tsv_dir} {params.output_dir}/combined_transplants.tsv.gz {params.output_dir}/combined_structure_summaries.tsv.gz
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
        """python3 alphacognate_transplant.py --foldseek_file {input.foldseek} --outdir {params.output_dir} --structure_database_directory {params.procoggraph_foldseek_directory} {params.cognate_match} {params.domain_match_only} && gzip {params.output_dir}/{wildcards.id}_transplants.cif"""

rule split_foldseek:
    input:
        foldseek_file = config["output_dir"] + "/foldseek/combined_foldseek.tsv",
        structure_manifest = config["structure_manifest"],
        predicted_structure_domains = config["predicted_structure_domains"],
        procoggraph_data = config["procoggraph_data"]
    output:
        expand(config["output_dir"] + "/foldseek_split/{id}_foldseek.tsv.gz", id = structures)
    params:
        output_dir = config["output_dir"] + "/foldseek_split"
    log: config["output_dir"] + "/logs/foldseek_split.log"
    shell:
      """python3 split_foldseek.py --foldseek_file {input.foldseek_file} --structure_manifest {input.structure_manifest} --output_dir {params.output_dir} --predicted_structure_domains {input.predicted_structure_domains} --procoggraph_data {input.procoggraph_data}"""

rule run_foldseek:
    input:
        structures = expand(config["output_dir"] + "/structures/{id}.cif", id = structures),
    output:
        foldseek_file = config["output_dir"] + "/foldseek/combined_foldseek.tsv"
    threads: workflow.cores
    params: 
        structure_dir = config["output_dir"] + "/structures",
        procoggraph_foldseek_db = config["procoggraph_foldseek_db"],
        procoggraph_foldseek_db_index = config["procoggraph_foldseek_db_index"],
    log: config["output_dir"] + "/logs/report.log"
    shell:
        """foldseek easy-search {params.structure_dir} {params.procoggraph_foldseek_db} {output.foldseek_file} {params.procoggraph_foldseek_db_index} --threads {threads} --prefilter-mode 1 --format-mode 4 --format-output query,target,u,t,qlen,tlen,alnlen,mismatch,qaln,qstart,qend,taln,tstart,tend,qseq,tseq,rmsd,qheader,alntmscore,qtmscore,ttmscore,evalue,prob > {log} 2> {log}"""

rule copy_structure:
    input:
        structure_dir = config["structure_dir"]
    output:
        temp(config["output_dir"] + "/structures/{id}.cif")
    shell:
        """
        cp {input.structure_dir}/{wildcards.id}.cif {output}
        """
