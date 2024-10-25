import os
import csv

#search directory for cif files to add to a manifest
def find_cif_files_in_directory(directory, output_csv):
    #write outputs to a file
    with open(output_csv, mode='w', newline='') as file:
        writer = csv.writer(file)
        
        #iterate over all the files in the directory and subdirectories
        for root, dirs, files in os.walk(directory):
            for filename in files:
                if filename.endswith('.cif'):
                    #extract basename of the file (without path)
                    basename = os.path.splitext(os.path.basename(filename))[0]
                    
                    #write required columns to file - basename, full filename, and directory
                    writer.writerow([basename, filename, root])

#specify io
directory_to_search = '/raid/MattC/repos/alpha_cognate/alphafold_structures'
output_csv_file = 'alphafold_structures_manifest.csv'

find_cif_files_in_directory(directory_to_search, output_csv_file)
print(f".cif file information saved to {output_csv_file}")