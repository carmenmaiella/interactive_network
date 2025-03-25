#script where I will call the single script for the COMPELEX WITH MERGED NODE (less nodes) and NOT MERGED NODE(more nodes)
#only the compute button function that will call everything (try to understand this part)
#they will produce as an output the interactive network
import complex_not_merged as no_merg
import complex_merged as merg
import protein_network as protein_net
import argparse
import sys 
import os
import pandas as pd
import igraph as ig
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import numpy as np
import json
from pyvis.network import Network

def parse_args():

    parser = argparse.ArgumentParser()
    
    if len(sys.argv)==1:
        parser.print_help()
        # parser.print_usage() # for just the usage line
        parser.exit()

    parser.add_argument('-i', dest='in_dir', default='',
                        help='path to the input directory where the ouput of AlphBridge are stored')
    
    parser.add_argument('-o', dest='o_dir', default='',
                        help='path to a directory where output folder are stored')
    
    parser.add_argument('-t', dest= 'threshold', choices=['0.5', '0.75', '0.9'] , default='0.9',
                        help='output from different threshold. Options: 0.5, 0.7, 0.75, 0.9')
    
    #parser.add_argument('-tab', dest='table', type=lambda x: x.lower() == 'true', default=False,
    #                    help='Save or not the comparison dataframe. Options: True or False')
    
    args = parser.parse_args()

    return args

def from_json_to_df(in_dir, threshold):
    #CONVERTING THE JSON FILE IN A DF
    all_data = []
        
    # Extract folder name from `in_dir`, ignore the first part 
    folder_name = os.path.basename(in_dir).split('_', 1)[1] 

    # Extract protein names from the folder, and assign each protein a letter
    proteins = folder_name.split('_')  # Split the folder name to get protein names
    protein_to_letter = {protein: chr(65 + i) for i, protein in enumerate(proteins)}

    # Traverse through the input directory to check all subdirectories
    for subdir in os.listdir(in_dir):
        subdir_path = os.path.join(in_dir, subdir)

        # Ensure we are processing only valid directories
        if not os.path.isdir(subdir_path):
            continue

        # Check if this directory itself is "AlphaBridge" or contains it
        if "AlphaBridge" in subdir:
            # Case 1: "AlphaBridge" is an immediate subdirectory of `in_dir`
            alphabridge_path = subdir_path
        else:
            # Case 2: "AlphaBridge" is inside another subdirectory
            alphabridge_path = os.path.join(subdir_path, "AlphaBridge")

        # Process if the "AlphaBridge" directory exists
        if os.path.isdir(alphabridge_path):
            for file in os.listdir(alphabridge_path):
                if file.endswith(".json"):
                    json_file_path = os.path.join(alphabridge_path, file)
                    with open(json_file_path, 'r') as json_file:
                        json_data = json.load(json_file)

                    # Process JSON data and store in a list
                    rows = []
                    if "interactions" in json_data:
                        for interaction in json_data["interactions"]:
                            if interaction.get("cut-off") == threshold:
                                for interface in interaction["interfaces"]:
                                    interface_name = interface["interface_name"]
                                    for link in interface["links"]:
                                        row = {
                                            "folder_name": folder_name,
                                            "interface": interface_name,
                                            "prot_1": link["first"]["asym_id"],
                                            "start_1": link["first"]["link_range"]["start"],
                                            "end_1": link["first"]["link_range"]["end"],
                                            "prot_2": link["second"]["asym_id"],
                                            "start_2": link["second"]["link_range"]["start"],
                                            "end_2": link["second"]["link_range"]["end"],
                                        }
                                        rows.append(row)

                    # Add rows to the main list
                    all_data.extend(rows)

    # Convert the list to a DataFrame, THIS IS THE ONE THAT IS NEEDED FOR BOTH
    df = pd.DataFrame(all_data)
    #print(df)
    # Replace letters with protein names in 'prot_1' and 'prot_2' columns
    for column in ['prot_1', 'prot_2']:
        df[column] = df[column].apply(
            lambda protein_code: next((protein for protein, letter in protein_to_letter.items() if letter == protein_code), protein_code)
        )
    print(df)
    return df



'''MAIN FUNCTION FOR OBTANING BOTH NETWORK: WITHOUT MERGING (MORE NODES) AND WITH MERGIN (LESS NODES)'''

def main():
    args = parse_args()
    in_dir = args.in_dir
    threshold = float(args.threshold)
    # funtion that creates the df that will be used as an inpt for the next fuction
    df = from_json_to_df(in_dir, threshold)
    #function --> INFORMATION FOR THE NETWORK WITH MORE NODES(NO MERGED)
    j_not_merged = no_merg.get_protein_network_no_merging(df)
    print("FINISHED THE NOT MERGED NETWORK")
    #function --> INFORMATION FOR THE NETWORK WITH LESS NODES(MERGED)
    j_merged = merg.get_protein_network_merging(df)
    print("FINISHED THE MERGED NETWORK")
    #function --> INFORMATION FOR THE NETWORK WHERE NOES == PROTEIN
    j_proteins = protein_net.get_protein_network(df)
    print("FINISHED THE WHOLE PROTEIN NETWORK")
    # Combine them into a single dictionary
    combined_networks = {
        "network_not_merged": j_not_merged,
        "network_merged": j_merged,
        "protin network": j_proteins 
    }

    # Save to a JSON file
    with open(f"{args.o_dir}/combined_networks.json", "w") as f:
        json.dump(combined_networks, f, indent=4)

if __name__ == '__main__':
    main()
