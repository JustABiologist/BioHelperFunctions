def create_topology_file(input_filename, output_topology_filename):
    with open(input_filename, 'r') as file:
        lines = file.readlines()

    with open(output_topology_filename, 'w') as topo_file:
        for line in lines:
            if line.startswith("ATOM") or line.startswith("SHEET") or line.startswith("CONECT") or line.startswith("TER"):
                topo_file.write(line)
            elif line.startswith("ENDMDL"):
                topo_file.write("END\n")
                break

def create_trajectory_files(input_filename, trajectory_folder):
    with open(input_filename, 'r') as file:
        lines = file.readlines()

    frame_lines = []
    frame_count = 0
    for line in lines:
        if line.startswith("REMARK") and frame_lines:
            with open(f"{trajectory_folder}/frame_{frame_count}.pdb", 'w') as frame_file:
                frame_file.writelines(frame_lines)
                frame_file.write("END\n")
            frame_lines = []
            frame_count += 1
        if line.startswith("ATOM"):
            frame_lines.append(line)
        elif line.startswith("ENDMDL"):
            frame_lines.append("END\n")

    # Write the last frame
    if frame_lines:
        with open(f"{trajectory_folder}/frame_{frame_count}.pdb", 'w') as frame_file:
            frame_file.writelines(frame_lines)

input_filename = '/Users/floriangrun/Downloads/all-2/file_MD000.pdb'  # Replace with your file path
output_topology_filename = '/Users/floriangrun/Downloads/all-2/topology.pdb'
trajectory_folder = '/Users/floriangrun/Downloads/all-2/trajectories'

create_topology_file(input_filename, output_topology_filename)
create_trajectory_files(input_filename, trajectory_folder)