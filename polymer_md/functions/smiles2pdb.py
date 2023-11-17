import csv
import itertools
from rdkit import Chem
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdMolTransforms
from rdkit.Chem.rdMolTransforms import CanonicalizeConformer
from rdkit.Chem.rdDistGeom import EmbedMolecule
import sys

def SmilesToPDB(smiles, name, directory):
    forbidden_codes = ["AAA", "BBB", "CCC", "UNL"]
    # Load existing residue codes from the CSV file
    residue_codes = load_residue_codes()

    # Check if the name or smiles is already in the database
    existing_entry = find_existing_entry(residue_codes, name, smiles)

    if existing_entry:
        # Use existing residue code
        residue_code = existing_entry[2]  # Assuming code is the third column
    else:
        # Generate a unique 3-letter residue code excluding forbidden codes
        residue_code = generate_unique_residue_code(residue_codes, forbidden_codes)

        # Update the CSV file with the new entry
        update_residue_codes_csv(name, smiles, residue_code)

    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol, explicitOnly=False)
    EmbedMolecule(mol, useExpTorsionAnglePrefs=True, useBasicKnowledge=True)
    for conf in mol.GetConformers():
        rdMolTransforms.CanonicalizeConformer(conf)
    Chem.AssignAtomChiralTagsFromStructure(mol)

    # Replace default "UNL" codes in the PDB content
    mol_block = Chem.MolToPDBBlock(mol)
    pdb_content = mol_block.replace(" UNL", f" {residue_code}")

    # Add residue information to the PDB file
    pdb_filepath = directory + "/" + name + ".pdb"
    with open(pdb_filepath, "w") as pdb_file:
        pdb_file.write(pdb_content)

def load_residue_codes():
    # Load existing residue codes from the CSV file
    residue_codes = []
    with open("/scratch/scw1977/dan/polymer_md/pdb_files/residue_codes.csv", "r") as csv_file:
        reader = csv.reader(csv_file)
        residue_codes = [row for row in reader]
    return residue_codes

def generate_unique_residue_code(residue_codes, forbidden_codes=None):
    # Generate a unique 3-letter residue code not already in the database and not in the forbidden codes list
    existing_codes = set(entry[2] for entry in residue_codes)  # Assuming code is the third column
    forbidden_codes = set(forbidden_codes) if forbidden_codes else set()

    # Generate all possible combinations of three letters
    all_combinations = [''.join(combination) for combination in itertools.product('ABCDEFGHIJKLMNOPQRSTUVWXYZ', repeat=3)]

    # Find the first unused code not in the forbidden codes list
    new_code = next((code for code in all_combinations if code not in existing_codes and code not in forbidden_codes), None)

    if new_code is None:
        raise ValueError("Unable to generate a unique residue code.")

    return new_code

def update_residue_codes_csv(name, smiles, residue_code):
    # Update the CSV file with the new entry if it doesn't exist
    residue_codes = load_residue_codes()
    existing_entry = find_existing_entry(residue_codes, name, smiles)

    if not existing_entry:
        with open("/scratch/scw1977/dan/polymer_md/pdb_files/residue_codes.csv", "a", newline="") as csv_file:
            writer = csv.writer(csv_file)

            # Add a new entry to the CSV file
            writer.writerow([name, smiles, residue_code])

def find_existing_entry(residue_codes, name, smiles):
    for entry in residue_codes:
        if len(entry) >= 2 and (entry[0] == name or entry[1] == smiles):
            return entry
    return None