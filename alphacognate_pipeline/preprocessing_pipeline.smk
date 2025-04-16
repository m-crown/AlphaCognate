if not config.get("data_dir"):
    config["data_dir"] = "data"

if config.get("demo"):
    config["predicted_structures_manifest"] = "demo_manifests/uniprot_ids.txt"
    config["demo_procoggraph_structures"] = "demo_manifests/target_pdbs.txt"


if config.get("predicted_structures_manifest"):
    with open(config["predicted_structures_manifest"], "r") as f:
        predicted_structures = [line.strip() for line in f.readlines()]
else:
    predicted_structures = []

rule all:
    input:
        config["data_dir"] + "/pcgDB/pcgDB",
        config["data_dir"] + "/predicted_structures_download_complete.txt",
        config["data_dir"] + "/cath_single_chain_domain_interactions.tsv.gz",
        config["data_dir"] + "/alphafold_structures_manifest.csv",
        expand(config["data_dir"] + "/predicted_structures/AF-{id}-F1-model_4.cif", id = predicted_structures)

rule build_foldseek_db:
    input:
        config["data_dir"] + "/procoggraph_assemblies/download_complete.txt",
    output:
        config["data_dir"] + "/pcgDB/pcgDB",
        config["data_dir"] + "/pcgDB/pcgDB.idx"
    params:
        output_dir = config["data_dir"],
        input_dir = config["data_dir"] + "/procoggraph_assemblies"
    shell:
        "preprocessing/build_foldseek_database.sh {params.input_dir} {params.output_dir}"

rule download_procoggraph:
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
        python3 preprocessing/get_procoggraph_structures.py {input.procoggraph_manifest} {params.output_dir} && gzip {params.output_dir}/*.cif
        """

rule make_manifest:
    input:
        config["data_dir"] + "/predicted_structures_download_complete.txt"
    output:
        config["data_dir"] + "/alphafold_structures_manifest.csv"
    params:
        structure_dir = config["data_dir"] + "/predicted_structures/",
        output_dir = config["data_dir"]

    shell:
        """
        python3 preprocessing/make_manifest.py {params.structure_dir} {params.output_dir}
        """

if config.get("predicted_structures_manifest"):
    rule download_predicted_structures:
        #if demo download script is x and if not demo download script is y
        input:
            predicted_structure_manifest = config["predicted_structures_manifest"]
        output:
            expand(config["data_dir"] + "/predicted_structures/AF-{id}-F1-model_4.cif", id = predicted_structures),
            download_complete_file = config["data_dir"] + "/predicted_structures_download_complete.txt"

        params:
            output_dir = config["data_dir"] + "/predicted_structures"
        shell:
            """
            python3 preprocessing/get_alphafold_structures.py {input.predicted_structure_manifest} {params.output_dir} {output.download_complete_file}
            """
else:
    rule check_predicted_structures:
        input:
            expand(config["data_dir"] + "/predicted_structures/{id}.cif", id = predicted_structures)
        output:
            expand(config["data_dir"] + "/predicted_structures/{id}.cif", id = predicted_structures),
        shell:
            """
            echo "Predicted download structure step skipped, all manifests structures present."
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
        python3 preprocessing/format_procoggraph.py --cath_domain_ownership {params.data_dir}/procoggraph_data/cath_pdb_residue_interactions.csv.gz --scores_file {params.data_dir}/procoggraph_data/all_parity_calcs.pkl --cognate_ligands {params.data_dir}/procoggraph_data/cognate_ligands_df.pkl --bound_entity_descriptors {params.data_dir}/procoggraph_data/bound_descriptors.tsv.gz --output_dir {params.data_dir}
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
        unzip {params.data_dir}/alphacognate_data_files.zip -d {params.data_dir}
        """