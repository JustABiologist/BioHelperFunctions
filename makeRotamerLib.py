from Bio import PDB
import os
import numpy as np
import csv
from concurrent.futures import ProcessPoolExecutor
import argparse

def calculate_angles(residue):
    """Calculate phi, psi, and chi angles for a given residue."""
    angles = []

    # Phi angle
    if residue.get_id()[1] > 1:
        phi = PDB.Polypeptide.phi(residue)
        angles.append(np.degrees(phi) if phi is not None else None)

    # Psi angle
    psi = PDB.Polypeptide.psi(residue)
    angles.append(np.degrees(psi) if psi is not None else None)

    # Chi angles (only for side chain)
    for i in range(1, 5):  # Chi1 to Chi4
        try:
            chi = getattr(PDB.Polypeptide, f'Chi{i}')(residue).degrees
            angles.append(chi if chi is not None else None)
        except:
            angles.append(None)

    return angles

def process_single_pdb(filename):
    pdb_parser = PDB.PDBParser(QUIET=True)
    data = []

    structure = pdb_parser.get_structure(filename, filename)
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname() == 'CYS':
                    angles = calculate_angles(residue)
                    data.append([os.path.basename(filename), *angles])

    return data

def process_pdb_files(input_directory, output_file, cores):
    pdb_files = [os.path.join(input_directory, f) for f in os.listdir(input_directory) if f.endswith(".pdb")]
    all_data = []

    with ProcessPoolExecutor(max_workers=cores) as executor:
        results = executor.map(process_single_pdb, pdb_files)

        for result in results:
            all_data.extend(result)

    # Write data to CSV
    with open(output_file, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Filename', 'Phi', 'Psi', 'Chi1', 'Chi2', 'Chi3', 'Chi4'])
        for data in all_data:
            writer.writerows(data)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process PDB files and generate a cysteine rotamer CSV.")
    parser.add_argument('-i', '--input', required=True, help="Input directory containing PDB files")
    parser.add_argument('-o', '--output', required=True, help="Output CSV file path")
    parser.add_argument('-c', '--cores', type=int, default=4, help="Number of cores to use (default: 4)")

    args = parser.parse_args()
    process_pdb_files(args.input, args.output, args.cores)