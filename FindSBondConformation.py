import numpy as np
import os

def read_pdb(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    cysteines = [line for line in lines if line.startswith("ATOM") and "CB" in line and "CYS" in line]
    return cysteines

def get_coordinates(atom_line):
    x = float(atom_line[30:38].strip())
    y = float(atom_line[38:46].strip())
    z = float(atom_line[46:54].strip())
    return np.array([x, y, z])

def calculate_distances(cysteines):
    distances = np.zeros((len(cysteines), len(cysteines)))
    for i, cysteine_i in enumerate(cysteines):
        for j, cysteine_j in enumerate(cysteines):
            if i != j:
                coord_i = get_coordinates(cysteine_i)
                coord_j = get_coordinates(cysteine_j)
                distances[i, j] = np.linalg.norm(coord_i - coord_j)
    return distances

def find_pairs(distances, threshold=5.0):
    pairs = []
    used_indices = set()
    for i in range(len(distances)):
        if i in used_indices:
            continue
        for j in range(i + 1, len(distances)):
            if j in used_indices:
                continue
            if distances[i, j] <= threshold:
                pairs.append((i, j))
                used_indices.add(i)
                used_indices.add(j)
                break
    return pairs

def process_pdb_file(pdb_file):
    cysteines = read_pdb(pdb_file)
    distances = calculate_distances(cysteines)
    pairs = find_pairs(distances)
    return len(pairs) * 2 == len(cysteines)

def process_directory(directory:str):
    results = []
    for filename in os.listdir(directory):
        if filename.endswith(".pdb"):
            result = process_pdb_file(os.path.join(directory, filename))
            results.append(result)
    return results

def frame_picker(results, directory, output_path):
    counter = 0
    if not os.path.isdir(output_path):
        os.system(f"mkdir {output_path}")
        print(f"Created directory for output at {output_path}")
    l_dir = os.listdir(directory)
    for result, file in zip(results, l_dir):
        if result:
            os.system(f"cp {directory}/{file} {output_path}/{file}")
            counter += 1
    print(f"There were {counter} frames in the trajectory that satisfied the set conditions")
        


def main():
    pdb_file = '/Users/floriangrun/Downloads/all-2/trajectories/' # Replace with your PDB file name
    output_directory = "/Users/floriangrun/Downloads/all-2/relevant_frames"
    results = process_directory(pdb_file)
    frame_picker(results, pdb_file, output_directory)
    
if __name__ == "__main__":
    main()
