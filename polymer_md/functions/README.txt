These functions are called by the main script: csv_to_solvated_pdbs.py

The input to that script is a csv file in the following format:
  name, smiles
  name, smiles
  name, smiles
  name, smiles
  etc..

These are the names and smiles of any molecules to solvate or molecules that will be a solvent (i.e. water/ethanol)

1) smiles2pdb.py
  - This script takes the csv and generates pdb files of each - assigning each molecule (whether it is a solute or solvent) a unique residue so it can be tracked.

2) molecule_solvater4packmol.py
  - This script generates inputs files for packmol of the solvated molecules

3) packmol_runner.py
  - This script runs the packmol program on any constructed input files
