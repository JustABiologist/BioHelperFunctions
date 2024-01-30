from Bio import PDB
import numpy as np
import os
from multiprocessing import Pool

def process_pdb_file(pdb_file):
    parser = PDB.PDBParser()
    structure = parser.get_structure('X', pdb_file)
    dihedrals = []

    for model in structure:
        for chain in model:
            polypeptides = PDB.PPBuilder().build_peptides(chain)
            for poly_index, poly in enumerate(polypeptides):
                phi_psi = poly.get_phi_psi_list()
                for res_index, dihedral in enumerate(phi_psi):
                    if None not in dihedral:  # Ignore residues missing an angle
                        phi, psi = dihedral
                        # Convert radians to degrees and adjust range from [0, 360) to [-180, 180)
                        phi = (phi * 180 / np.pi) % 360 - 180
                        psi = (psi * 180 / np.pi) % 360 - 180
                        dihedrals.append((phi, psi))
    return dihedrals

def compile_frequencies(dihedral_list):
    frequency_matrix = np.zeros((360, 360), dtype=int)
    for phi, psi in dihedral_list:
        # Update the frequency matrix, note that we need to shift the range to index properly
        frequency_matrix[int(phi) + 180, int(psi) + 180] += 1
    return frequency_matrix

# List all PDB files
pdb_directory = '/media/florian/dc434ad3-38cd-442f-b7af-41802cfa5baa/top2018'
pdb_files = [os.path.join(pdb_directory, f) for f in os.listdir(pdb_directory) if f.endswith('.pdb')]

# Define the number of processes to use
num_processes = 4  # Adjust based on your system

# Create a pool of processes and map the function over the PDB files
with Pool(processes=num_processes) as pool:
    results = pool.map(process_pdb_file, pdb_files)

# Initialize the master frequency matrix
master_frequency_matrix = np.zeros((360, 360), dtype=int)

# Compile results into the master frequency matrix
for frequencies in results:
    master_frequency_matrix += compile_frequencies(frequencies)

# Normalize the master frequency matrix to turn it into a probability map
probability_matrix = master_frequency_matrix / np.sum(master_frequency_matrix)

# Save the master frequency matrix to a file for future plotting
np.savetxt('/home/florian/Schreibtisch/Bachelorarbeit_Koch/Ramachandran/test.txt', probability_matrix)
