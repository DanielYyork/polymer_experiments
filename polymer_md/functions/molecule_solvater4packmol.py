# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 11:46:49 2023

@author: danie
"""
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdMolTransforms
from rdkit.Chem.rdMolTransforms import CanonicalizeConformer
from rdkit.Chem.rdDistGeom import EmbedMolecule

def vol_from_mol(mol):
    AllChem.EmbedMolecule(mol)
    volume = AllChem.ComputeMolVolume(mol)
    return volume

def gen_solvatedmol4packmol(box_size, molecule_name, solvent_name, num_solvent_molecules, path_for_final_structures):
    
    pdb_filepath = "/scratch/scw1977/dan/polymer_md/pdb_files/molecules"
    
    solvent_pdb_path = pdb_filepath + "/" + solvent_name
    molecule_pdb_path = pdb_filepath + "/" + molecule_name

    box_size2write = str(int(float(box_size)))
    
    input_content = f"""tolerance 2.0
output {path_for_final_structures}/{molecule_name}_solvated_{solvent_name}_{box_size2write}.pdb
filetype pdb
structure {molecule_pdb_path}.pdb
    number 1
    inside cube 0. 0. 0. {box_size2write}.
    fixed 10. 10. 10. 0. 0. 0.
end structure
structure {solvent_pdb_path}.pdb
    number {num_solvent_molecules}
    inside cube 0. 0. 0. {box_size2write}.
end structure
"""
    input_filename = f"{molecule_name}_{solvent_name}_{box_size2write}.inp"
    packmol_input_directory = "/scratch/scw1977/dan/polymer_md/simulation_systems/packmol_inputs"
    input_filepath = packmol_input_directory + "/" + input_filename
    with open(input_filepath, "w") as input_file:
        input_file.write(input_content)

def calculate_and_construct_input(box_size, molecule_name, solvent_name, path_for_final_structures):
    # Step 1: Determine Volume of the Box
    box_volume = float(box_size)**3

    # Step 2: Determine Volume of the Molecule and solvent
    pdb_master_dir = "/scratch/scw1977/dan/polymer_md/pdb_files/molecules"
    mol_filename = f'{molecule_name}.pdb'  # Replace with the actual MOL file name
    mol_pdb_filepath = pdb_master_dir + "/" + mol_filename
    molecule = Chem.MolFromPDBFile(mol_pdb_filepath)
    
    solvent_filename = f'{solvent_name}.pdb'
    solvent_pdb_filepath = pdb_master_dir + "/" + solvent_filename
    solvent = Chem.MolFromPDBFile(solvent_pdb_filepath)

    if molecule is not None and solvent is not None:        
        molecule = Chem.AddHs(molecule)
        AllChem.Compute2DCoords(molecule, nFlipsPerSample = 10, nSample = 500, bondLength = 3.0)
        EmbedMolecule(molecule, useRandomCoords = True, forceTol = 0.000075)
        for conf in molecule.GetConformers():
            rdMolTransforms.CanonicalizeConformer(conf)
        molecule_volume = vol_from_mol(molecule)
        
        solvent = Chem.AddHs(solvent)
        AllChem.Compute2DCoords(solvent, nFlipsPerSample = 10, nSample = 500, bondLength = 3.0)
        EmbedMolecule(solvent, useRandomCoords = True, forceTol = 0.000075)
        for conf in molecule.GetConformers():
            rdMolTransforms.CanonicalizeConformer(conf)
        solvent_volume = vol_from_mol(solvent)

        # Step 3: Calculate Remaining Space and Solvent Molecules
        remaining_space = box_volume - molecule_volume
        num_solvent_molecules = int(remaining_space / solvent_volume)

        # Step 4: Construct Input File
        gen_solvatedmol4packmol(box_size, molecule_name, solvent_name, num_solvent_molecules, path_for_final_structures)

        print(f"Box Volume: {box_volume} A^3")
        print(f"Molecule Volume: {molecule_volume} A^3")
        print(f"Remaining Space: {remaining_space} A^3")
        print(f"Number of Solvent Molecules: {num_solvent_molecules}")
    else:
        print("Failed to load the molecule.")
'''
# Example usage:
box_size = 20.0  # Adjust the box size as needed
molecule_name = "3HB"  # These are searching in the pdb_files/molecules directory for the molecule and solvent
solvent_name = "water"  # As above

calculate_and_construct_input(box_size, molecule_name, solvent_name)
'''
