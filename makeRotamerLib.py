from Bio import PDB
import os
import numpy as np
import csv
from concurrent.futures import ProcessPoolExecutor

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

def process_pdb_files(directory):
    pdb_files = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith(".pdb")]
    all_data = []

    with ProcessPoolExecutor(max_workers=16) as executor:
        results = executor.map(process_single_pdb, pdb_files)

        for result in results:
            all_data.extend(result)

    # Write data to CSV
    with open('cysteine_rotamers.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Filename', 'Phi', 'Psi', 'Chi1', 'Chi2', 'Chi3', 'Chi4'])
        for data in all_data:
            writer.writerows(data)

# Directory containing PDB files
pdb_directory = '/path/to/your/pdb/files'

# Process the PDB files and output to CSV
process_pdb_files(pdb_directory)