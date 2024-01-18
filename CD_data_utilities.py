import pandas as pd 
import numpy as np
import re
import os
from pathlib import Path

def parse_file_v4(filepath, keytuple: tuple):
    """Takes a path to a .pcd file. Will output a dataframe for the Data in it. Additionally, it will output
       a directory with relevant metadata. The relevant metadata is specified with the keytuple."""
    meta_dict = {}
    data_section = False
    data_entries = []
    pattern = r'\s{6,}([^ ]+(?: [^ ]+)*)'
    
    with open(filepath) as file:
        for line in file.readlines():
            # Check if the line starts with any of the keys in the keytuple
            if any(line.startswith(key) for key in keytuple):
                key = next(key for key in keytuple if line.startswith(key))
                value = re.search(pattern, line).group(1)
                meta_dict[key] = value
            
            # Data section begins
            if line.startswith("DATA"):
                data_section = True
                continue
            
            # Process data lines
            if data_section and not line.startswith("CALIBRATION"):
                values = line.strip().split()
                if values:
                    data_entries.append(values)
            # Stop parsing if 'CALIBRATION' is encountered
            if line.startswith("CALIBRATION"):
                break

    # Convert data entries to a DataFrame for easier manipulation
    if data_entries:
        data_df = pd.DataFrame(data_entries, columns=['Wavelength', 'Final', 'HT', 'Smoothed', 'Avg. Sample', 'Avg. Baseline'])
        data_df = data_df.apply(pd.to_numeric, errors='coerce')  # Convert all columns to numeric
    else:
        data_df = pd.DataFrame()

    return meta_dict, data_df

def do_for_directory(dirname, keytuple):
    l_data = []
    dir_path = Path(dirname)
    for path in os.listdir(dirname):
        joined = str(dir_path.joinpath(path))
        l_data.append(parse_file_v4(joined, keytuple))
    return l_data

def split_data_pdb(l_data):
    pdb_id = [entry for entry in l_data if entry[0]["PDB ID"] != "No data provided"]
    non_pdb = list(set(l_data) - set(pdb_id))
    return pdb_id, non_pdb






def main():
    # Test the function with a file path and keytuple
    keytuple = ("Experimental Temperature (C)", "PCDDBID", "PDB ID")
    metadata_v4, data_df_v4 = parse_file_v4("/Users/floriangrun/Desktop/Bachelorarbeit_Koch/CD_sim/cd_bsp.pcd", keytuple)
    print(metadata_v4)
    print(data_df_v4.head())

if __name__ == "__main__":
    main()