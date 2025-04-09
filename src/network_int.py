'''MAIN SCRIPT
input: confident interface interaction between proteins predictd by alphabridge
ouput: json file containing the 3 different networks:
    1) network that shows confident contact interface between protein (not merged)
    2) network that shows confident contact interface between protein (merged)
    3) network that shows a general overview of the confident physical interaction between proteins'''

#import required packages

import complex_not_merged as no_merg
import complex_merged as merg
import protein_network as protein_net
import argparse
import sys 
import os
import pandas as pd
import json

#running the script form the command line

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
    
    parser.add_argument('-t', dest='threshold', nargs='+', choices=['0.5', '0.75', '0.9'],
    default=['0.9'],help='Output from different thresholds. Options: 0.5, 0.75, 0.9. You can pass multiple values separated by space.'
)

    
    #parser.add_argument('-tab', dest='table', type=lambda x: x.lower() == 'true', default=False,
    #                    help='Save or not the comparison dataframe. Options: True or False')
    
    args = parser.parse_args()

    return args


#CONVERTING THE JSON FILE (output of alphabridge) IN A DF

def from_json_to_df(in_dir, threshold):
    all_data = []
        
    # handle situatio where the folder name has complex at the beginning or just the name of the protein
    #result --> just the protein name are kept 
    #folder_name = os.path.basename(in_dir)
    #if "omplex" not in folder_name: #so the c can be either C or c + every name that you prefer
    #    folder_name = folder_name
    #else:
    #    folder_name = folder_name.split('_', 1)[1]

    # Extract protein names from the folder, and assign each protein a letter
    #proteins = folder_name.split('_')  # Split the folder name to get protein names
    #protein_to_letter = {protein: chr(65 + i) for i, protein in enumerate(proteins)}

    # Traverse through the input directory to check all subdirectories
    for subdir in os.listdir(in_dir):
        subdir_path = os.path.join(in_dir, subdir)

        # Ensure we are processing only valid directories
        if not os.path.isdir(subdir_path):
            continue

        # Check where is the alphabridge output subdiectory
        if "AlphaBridge" in subdir:
            # Case 1: "AlphaBridge" is an immediate subdirectory of `in_dir`
            alphabridge_path = subdir_path
        #THIS WAS FOR THE BINARY INTERACTION, FOR NOW I PUT IT COMMENTED    
        #else:
            # Case 2: "AlphaBridge" is inside another subdirectory
            #alphabridge_path = os.path.join(subdir_path, "AlphaBridge")

        # Process if the "AlphaBridge" directory exists
        if os.path.isdir(alphabridge_path):
            for file in os.listdir(alphabridge_path):
                if file.endswith(".json"):
                    json_file_path = os.path.join(alphabridge_path, file)
                    with open(json_file_path, 'r') as json_file:
                        json_data = json.load(json_file)

                    # Process JSON data and store in a list
                    #from auth to label
                    auth2label = {}
                    chains = json_data['structure'][0]['chains']
                    for entity_type, rec_list in chains.items():
                        for rec in rec_list:
                            #print(rec)
                            #print(rec['auth_asym_id'], rec['label_asym_id'])
                            auth_asym_id = rec['auth_asym_id']
                            label_asym_id = rec['label_asym_id']
                            if not auth_asym_id in auth2label:
                                auth2label[auth_asym_id] = str()
                            auth2label[auth_asym_id] = label_asym_id
                    #from label to auth
                    label2auth = {value:key for key, value in auth2label.items()}
                    rows = []
                    if "interactions" in json_data:
                        for interaction in json_data["interactions"]:
                            if interaction.get("cut-off") == threshold:
                                for interface in interaction["interfaces"]:
                                    interface_name = interface["interface_name"]
                                    for link in interface["links"]:
                                        row = {
                                            #folder_name": folder_name,
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
                    #create the df for storing the score
                    pairwise_interaction = json_data['structure'][0]['pairwise_interaction']
                    #print(pairwise_interaction)
                    df_pairwise_interaction = pd.DataFrame(pairwise_interaction)
                    #print(df_pairwise_interaction)
 
    
    #print(auth2label)
    #print(label2auth)


    # Convert the list to a DataFrame, THIS IS THE ONE THAT IS NEEDED FOR BOTH
    df = pd.DataFrame(all_data)
    #print(df)


    # form auth to label for whatever is not a protein
    #for column in ['prot_1', 'prot_2']:
    #    df[column] = df[column].replace(auth2label)
    #print(df)
    return df, label2auth, df_pairwise_interaction



'''MAIN FUNCTION FOR OBTANING BOTH NETWORK: WITHOUT MERGING (MORE NODES) AND WITH MERGIN (LESS NODES)'''

def main():
    args = parse_args()
    in_dir = args.in_dir
    threshold = [float(t) for t in args.threshold]
    if not os.path.exists(args.o_dir):
        os.makedirs(args.o_dir)
    all_threshold = {}
    #iterate for every possible threshold
    for th in threshold:
        #print(f"th:{th}")
        # funtion that creates the df that will be used as an inpt for the next fuction
        df, dit, df_pairwise_interaction = from_json_to_df(in_dir, th)
        #function --> INFORMATION FOR THE NETWORK WITH MORE NODES(NO MERGED)
        j_not_merged = no_merg.get_protein_network_no_merging(df,dit)
        print("FINISHED THE NOT MERGED NETWORK")
        #function --> INFORMATION FOR THE NETWORK WITH LESS NODES(MERGED)
        #j_merged = merg.get_protein_network_merging(df)
        #print("FINISHED THE MERGED NETWORK")
        #function --> INFORMATION FOR THE NETWORK WHERE NOES == PROTEIN
        j_proteins = protein_net.get_protein_network(df,dit,df_pairwise_interaction,th)
        print("FINISHED THE WHOLE PROTEIN NETWORK")
        # Combine them into a single dictionary
        combined_networks = {"cut-off": th,
            "network_not_merged": j_not_merged,
        # "network_merged": j_merged,
            "protein_network": j_proteins 
        }
        #combine everything
        all_threshold[f"network at {th}"] = combined_networks

    # Save to a JSON file
    with open(f"{args.o_dir}/combined_networks.json", "w") as f:
        json.dump(all_threshold, f, indent=4)

if __name__ == '__main__':
    main()
