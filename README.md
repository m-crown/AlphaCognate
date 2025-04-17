# AlphaCognate

AlphaCognate is a tool to identify the cognate ligand binding potential of an AlphaFold structure.

It adapts the AlphaFill methodology first described by Hekkelman et al. in AlphaFill (2021) - updated to use FoldSeek and an expanded set of ligands for transplanting, based upon those with identified cognate ligand matches in the ProCogGraph database.

AlphaCognate processes single chain AlphaFold predictions, and as such utilises a subset of the ProCogGraph database in which domain interactions are detected within a single chain.

AlphaCognate also makes use of the CATH domain annotations of AlphaFold structures from (XYZ) and will be updated to utilise information from The Encyleopedia of Domains (TED) when this becomes available.

## Installation

To install AlphaCognate, clone the repository and install the required dependencies (recommend to use uv for quick installation, but can also use pip):

```bash
git clone https://github.com/m-crown/AlphaCognate.git
cd AlphaCognate
uv venv
uv pip install -r requirements.txt
```

In addition to the Python dependencies, AlphaCognate also requires FoldSeek for structure searching. To install FoldSeek, download the latest appropriate binary file from the [FoldSeek GitHub repo](https://github.com/steineggerlab/foldseek) (example below is shown for Linux with AVX2, binaries are available for Linux ARM64, MacOS and GPU enable - note this is not tested in AlphaCognate).

```bash
#activate the environment
source .venv/bin/activate

#download the latest FoldSeek binary
wget https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz && tar xvzf foldseek-linux-avx2.tar.gz &&  export PATH=$(pwd)/foldseek/bin/:$PATH

#verify FoldSeek is installed
foldseek --version
```

## Pre-requisites

AlphaCognate requires a comma separated manifest file containing the following columns (without a header row):

- structure_basename: The basename of the structure, e.g. AF-A0A023GPK8-F1-model_4
- filename: The filename of the structure, e.g. AF-A0A023GPK8-F1-model_4.cif
- directory: The directory containing the structure, e.g. alphafold_structures/

This manifest fle can be generated using the following command:

```bash
python3 bin/make_manifest.py $STRUCTURE_DIRECTORY $OUTPUT_DIRECTORY
```

You may also need to download the AlphaFold structure files from the AlphaFold database. This can be done using the following command and input file containing the list of UniProt IDs (one per line):

```bash
python3 bin/get_alphafold_structures.py $UNIPROTIDS_FILE $OUTPUTDIRECTORY
```

The manifest file can then be provided in the Snakemake command (see [Running AlphaCognate](#running-the-alphacognate-pipeline)) either as part of a config file (recommended) or directly as a command line argument e.g.:

```snakemake --configfile config.yaml```
or

```snakemake --config manifest=manifest.csv```

## Running The AlphaCognate Pipeline

To run AlphaCognate, the pipeline is executed with Snakemake. The workflow is designed to be run from the command line and requires the manifest file and domain assignments (see [Pre-requisites](#pre-requisites)) as input.

For the TL;DR version, the following command will run the pipeline with default parameters:

```bash 
snakemake -s alphacognate_pipeline.smk --configfile alphacognate_config.yaml --cores 1
```

Read on for more details on execution and options available.

The pipeline will first check the input structure files exist, and then download and process the required ProCogGraph database files (if not already present, total download approx xxMB).

The pipeline will then run the FoldSeek search against the ProCogGraph database, and then run the cognate ligand prediction step.