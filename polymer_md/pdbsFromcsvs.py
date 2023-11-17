# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 16:51:16 2023

@author: danie
"""

import csv
import subprocess
import sys as sys
from functions.smiles2pdb import *
from functions.molecule_solvater4packmol import *
from functions.packmol_runner import *

def print_error_message():
    print("Error: Incorrect number of arguments.")
    print("This script generates starting structures for individually solvated molecules.")
    print("Usage: python script_name.py <csv_file> <box_size> <solvent_type> <input_directory>")
    print("  <csv_file>: CSV file with names and smilestrings of each compound to solvate.")
    print("  <box_size>: Size of the box for simulations in the form 20.0 (MUST BE A FLOAT).")
    print("  <solvent_type>: Type of the solvent for solvation, e.g., water/ethanol.")
    print("  <input_directory>: Filepath to the solvated monomers input directory.")
    sys.exit(1)

# Check the number of arguments
if len(sys.argv) != 5:
    print_error_message()

def unpack_csv(csv_filename):
    names = []
    smiles = []

    with open(csv_filename, 'r') as csv_file:
        csv_reader = csv.reader(csv_file)
        for row in csv_reader:
            if len(row) >= 2:
                names.append(row[0])
                smiles.append(row[1])

    return names, smiles

# Example usage:
csv_filename = sys.argv[1]
names_list, smiles_list = unpack_csv(csv_filename)

# Now names_list and smiles_list contain the respective lists from the CSV file.
for i in range(len(smiles_list)):
    SmilesToPDB(smiles_list[i], names_list[i], "/scratch/scw1977/dan/polymer_md/pdb_files/molecules")
    
 # Example usage:
box_size = sys.argv[2] # Adjust the box size as needed  
solvent_name = str(sys.argv[3]) # These are searching in the pdb_files/molecules directory for the molecule and solvent

for i in range(len(names_list)):
    calculate_and_construct_input(box_size, names_list[i], solvent_name, sys.argv[4])   # Sys 4 is the final save path of the pdb saved by packmol

packmol_input_directory = "/scratch/scw1977/dan/polymer_md/simulation_systems/packmol_inputs/"
files2packmol = find_matching_files(packmol_input_directory, names_list, solvent_name)

path_to_packmol = "/home/s.983045/bin/packmol/packmol-20.14.0/packmol"
for file in files2packmol:
    print("Running packmol on: ", os.path.basename(file))
    print("Packmol called from: ", path_to_packmol)
    submit_packmol_input(file, path_to_packmol)