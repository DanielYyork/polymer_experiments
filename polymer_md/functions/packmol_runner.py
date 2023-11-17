import os
import subprocess as subprocess

def find_matching_files(directory, names, solvent):
    matching_files = []

    for filename in os.listdir(directory):
        for name in names:
            if name in filename and solvent in filename:
                matching_files.append(os.path.join(directory, filename))

    return matching_files

def submit_packmol_input(packmol_input_filepath, path_to_packmol):     
    try:
        subprocess.run(f'{path_to_packmol} < {packmol_input_filepath}', shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running packmol: {e}")
  
    
