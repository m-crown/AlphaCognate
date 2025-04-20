import os
import csv
import argparse

#search directory for cif files to add to a manifest
def find_cif_files_in_directory(directory, output_csv):
    #write outputs to a file
    with open(output_csv, mode='w', newline='') as file:
        writer = csv.writer(file)
        
        #iterate over all the files in the directory and subdirectories
        for root, dirs, files in os.walk(directory):
            for filename in files:
                if filename.endswith('.cif'):
                    print(filename)
                    #extract basename of the file (without path)
                    basename = os.path.splitext(os.path.basename(filename))[0]
                    
                    #write required columns to file - basename, full filename, and directory
                    writer.writerow([basename, filename, root])

def main():
    """
    Script to produce the required manifest for processing alphafold predicted structures in AlphaCognate pipeline.
    """
    parser = argparse.ArgumentParser(description="Get demo structures from a CSV file.")    
    parser.add_argument("structure_dir", type=str, help="Directory to search for structures in")
    parser.add_argument("--output_dir", default = ".", type=str, help="directory to output manifest file to.")

    args = parser.parse_args()
    #specify io
    output_csv_file = f'{args.output_dir}/alphafold_structures_manifest.csv'
    find_cif_files_in_directory(args.structure_dir, output_csv_file)
    print(f".cif file information saved to {output_csv_file}")

if __name__ == "__main__":
    main()