import gemmi
import pandas as pd
import numpy as np
from pymol import cmd
import matplotlib.pyplot as plt  # For generating distinct colors
import argparse

# Convert the Matplotlib RGBA color format to PyMOL format (hex)
def rgba_to_pymol_color(rgba):
    return f'0x{int(rgba[0]*255):02x}{int(rgba[1]*255):02x}{int(rgba[2]*255):02x}'

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--structure', metavar='structure with ligand transplants', type=str, help="Path to the transplants structure file")
    parser.add_argument('--output_file', type=str, help='User-specified output file for the filtered transplants, if not specified defaults to [in]_coloured.pse')
    args = parser.parse_args()

    if not args.output_file:
        if args.structure.endswith(".cif.gz"):
            args.output_file = args.structure.replace(".cif.gz", "_coloured.pse")
        else:
            args.output_file = args.structure.replace(".cif", "_coloured.pse")

    # Load the structure with gemmi
    af_structure = gemmi.cif.read(args.structure)
    query_block = af_structure.sole_block()
    # Process the data (assuming a loop exists with cluster info)
    structure_transplants = pd.DataFrame(query_block.get_mmcif_category("_alphacognate_transplants."))
    structure_transplants = structure_transplants.replace("", np.nan)
    structure_transplants[["cluster", "tcs"]] = structure_transplants[["cluster", "tcs"]].astype("float")

    # Get the number of unique clusters
    unique_clusters = structure_transplants["cluster"].unique()
    num_clusters = len(unique_clusters)

    # Generate a distinct color for each cluster using Matplotlib's colormap
    colors = plt.cm.get_cmap('tab20', num_clusters)  # 'tab20' colormap can handle up to 20 distinct colors

    # Create a dictionary to map each cluster to a color
    cluster_color_dict = {cluster: colors(i) for i, cluster in enumerate(unique_clusters)}

    

    # Assign colors to clusters and apply to chains
    chain_cluster_dict = {row['ligand_chain']: row['cluster'] for index, row in structure_transplants.iterrows()}

    # Load the structure into PyMOL
    cmd.load(args.structure)
    colour_indices = cmd.get_color_indices()
    # Color chains by cluster in PyMOL
    for chain, cluster in chain_cluster_dict.items():
        pymol_color = rgba_to_pymol_color(cluster_color_dict[cluster])
        cmd.set_color(f'cluster_{int(cluster)}_color', list(cluster_color_dict[cluster][:3]))  # Set custom color
        cmd.color(f'cluster_{int(cluster)}_color', f'chain {chain}')  # Apply color to chain

    # Save the PyMOL session
    print("saving session")
    cmd.save(args.output_file)  # Save the PyMOL session as a .pse file

if __name__ == "__main__":
    main()