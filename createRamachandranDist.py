from Bio import PDB
import numpy as np
import os
from multiprocessing import Pool
from tqdm import tqdm
import pathlib
import matplotlib.pyplot as plt

def process_pdb_file(pdb_file):
    parser = PDB.PDBParser(QUIET=True)  # QUIET suppresses warnings
    structure = parser.get_structure('X', pdb_file)
    
    dihedrals_by_residue = {residue: [] for residue in PDB.Polypeptide.standard_aa_names}

    for model in structure:
        for chain in model:
            polypeptides = PDB.PPBuilder().build_peptides(chain)
            for poly in polypeptides:
                phi_psi = poly.get_phi_psi_list()
                for res_index, dihedral in enumerate(phi_psi):
                    residue = poly[res_index].get_resname()
                    if residue in dihedrals_by_residue and None not in dihedral:
                        phi, psi = (angle * 180 / np.pi + 180 for angle in dihedral)
                        dihedrals_by_residue[residue].append((phi % 360, psi % 360))

    return dihedrals_by_residue

def compile_frequencies(dihedral_data):
    frequencies = {residue: np.zeros((360, 360), dtype=int) for residue in PDB.Polypeptide.standard_aa_names}

    for residue, dihedrals in dihedral_data.items():
        for phi, psi in dihedrals:
            frequencies[residue][int(phi), int(psi)] += 1

    return frequencies

def plot_ramachandran(residue, frequency_matrix, output_path):
    plt.figure(figsize=(8, 6))
    plt.title(f'Ramachandran Plot for {residue}')
    plt.xlabel('Phi angles (degrees)')
    plt.ylabel('Psi angles (degrees)')
    plt.contourf(range(-180, 180), range(-180, 180), frequency_matrix, levels=50, cmap='viridis')
    plt.colorbar(label='Frequency')
    plt.xlim(-180, 180)
    plt.ylim(-180, 180)
    plt.grid(True)
    plt.savefig(os.path.join(output_path, f'Ramachandran_{residue}.png'))
    plt.close()

def main():
    input_path = '/media/florian/dc434ad3-38cd-442f-b7af-41802cfa5baa/top2018'
    output_path = '~/Schreibtisch/Bachelorarbeit_Koch/Ramachandran'

    pdb_files = [os.path.join(input_path, f) for f in os.listdir(input_path) if f.endswith('.pdb')]

    num_processes = 14  # Adjust based on your system

    with Pool(processes=num_processes) as pool:
        all_results = list(tqdm(pool.imap(process_pdb_file, pdb_files), total=len(pdb_files)))

    master_frequencies = {residue: np.zeros((360, 360), dtype=int) for residue in PDB.Polypeptide.standard_aa_names}

    for result in all_results:
        for residue in master_frequencies:
            residue_data = result.get(residue, [])
            if residue_data:
                # Convert list to numpy array with dtype=int and update frequency matrix
                for phi, psi in residue_data:
                    master_frequencies[residue][int(phi), int(psi)] += 1
            else:
                # If no data for residue, skip
                continue
    expanded_output_path = os.path.expanduser(output_path)
    pathlib.Path(os.path.expanduser(output_path)).mkdir(parents=True, exist_ok=True)

    for residue, matrix in master_frequencies.items():
        normalized_matrix = matrix / np.sum(matrix)
        file_name = os.path.join(os.path.expanduser(output_path), f'{residue}_frequency_matrix.txt')
        np.savetxt(file_name, normalized_matrix)
        plot_ramachandran(residue, normalized_matrix, expanded_output_path)

if __name__ == "__main__":
    main()

