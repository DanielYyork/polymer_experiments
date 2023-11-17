The python script here "pdsbsfromcsvs" generates starting structures for molecular dynamics simulation of solvated molecules (individual molecules)

To run the script in Supercomputing wales there are a series of steps to follow:
  1) run the command: module load mamba/4.14.0-0
  2) run the command: source activate /scratch/scw1977/openmm-7.7/
  - The 2 commands above activate the environment as the python scripts in "polymer_md/functions" require RDkit and a future version will require openmm (both of which are in this environment)]
  3) run the script: python3 spdsbsfromcsvs.py <csv_file> <box_size> <solvent_type> <input_directory>

  - Where each argument is as follows:
       <csv_file>: CSV file with names and smilestrings of each compound to solvate.
       <box_size>: Size of the box for simulations in the form 20.0 (MUST BE A FLOAT).
       <solvent_type>: Type of the solvent for solvation, e.g., water/ethanol.
       <input_directory>: Filepath to the solvated monomers input directory. i.e. /scratch/scw1977/dan/polymer_md/simulation_systems/packmol_inputs/ 

  - Recommended: python3 spdsbsfromcsvs.py <csv_file> <box_size> <solvent_type> <input_directory> > experiment_1.txt
    - "experiment_1.txt" will contain any outputs from the script

  -NOTE: A series of files are made with this script:
    pdbfiles of each molecule/solvent in the csv in: "polymer_md/pdb_files/molecules/"
    packmol input files for solvated systems of monomers in: "polymer_md/simulation_systems/packmol_inputs"
    final PDB files of MD starting systems in: "polymer_md/simulation_systems/solvated_monomers"

- This code is very new and the altering of some paths may be neccesarily so it is recommended to run this script interactively and address any errors in pathing.

DEPENDANCIES:
python enviroment containing openmm and rdkit
packmol software - a path to the executable command
