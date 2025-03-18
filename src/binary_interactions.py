#script where the svg binary interaction netwrok is computed
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
from networkx.readwrite import json_graph;
def parse_args():

    parser = argparse.ArgumentParser()
    
    if len(sys.argv)==1:
        parser.print_help()
        # parser.print_usage() # for just the usage line
        parser.exit()

    parser.add_argument('-i', dest='in_dir', default='',
                        help='main directory where all trhe subdirectory of each binary interaction are stored')
    
    parser.add_argument('-o', dest='o_dir', default='',
                        help='path to a directory where output folder are stored')
    
    parser.add_argument('-t', dest= 'threshold', choices=['0.5', '0.75', '0.9'] , default='0.9',
                        help='output from different threshold. Options: 0.5, 0.7, 0.75, 0.9')
    
    #parser.add_argument('-tab', dest='table', type=lambda x: x.lower() == 'true', default=False,
    #                    help='Save or not the comparison dataframe. Options: True or False')
    
    args = parser.parse_args()

    return args

#FUNCTION FOR:

#1 MERGING INTERVALS (START,END)

def defining_intervals(interval1, interval2, min_value_for_merging):
    #case 1 cheeking if they are equal
    if interval1 == interval2:
        return interval1
    start_1, end_1 = interval1
    start_2, end_2 = interval2
    overlap =min(end_1,end_2)-max(start_1,start_2)+1
    #case 2 check if one intervals contains another
    #1st example wjere the interval1 contains the interval 2
    if start_1 <=start_2 and end_1 >= end_2:
        return interval1
    #2nd example where the interval2 contains the int 1
    if start_2 <=start_1 and end_2 >= end_1:
        return interval2
    #case 3 check if the two interval overlap with a length > of the min_value_for_merging value
    if overlap > min_value_for_merging:
        if not (end_1<start_2 or end_2<start_1): 
            return[min(start_1,start_2), max(end_1,end_2)]
    else:
        return interval1, interval2
    

#2 MERGING INTERVALS --> USING THE defining_intervals FUNCTION IN THE BIG DATAFRAME

def get_merging_intervals(df, min_value_for_merging):
    for i in range(len(df['protein'])):
        interval_1 = df.loc[i, 'start_merged'], df.loc[i, 'end_merged']
        for j in range(len(df['protein'])-1):
            interval_2 = df.loc[j+1, 'start_merged'], df.loc[j+1, 'end_merged']
            new_intervals = defining_intervals(interval_1, interval_2, min_value_for_merging)
            if isinstance(new_intervals, tuple) and all(isinstance(i,tuple)for i in new_intervals):
               pass
            else:
                df.at[i,'start_merged'] = new_intervals[0]
                df.at[i,'end_merged'] = new_intervals[1]
                df.at[j+1,'start_merged'] = new_intervals[0]
                df.at[j+1,'end_merged'] = new_intervals[1]
    return df


#3 MERGING INTERFACES

#function that put  the intrval of the same interface in a list of list
def create_new_column_interface_intervals(df):
    #grouped = df.groupby(["num_binary_combination", "interface"])
    grouped = df.groupby(["bin_comb_with", "interface"])
    #grouped = df.groupby(["numn_int","interface"])
    #list for intervas for each interfaces
    interface_intervals = [None] * len(df)
    for _, group in grouped:
        intervals = group[["start_merged", "end_merged"]].drop_duplicates().values.tolist()
        for idx in group.index:
            interface_intervals[idx] = intervals
    df["interface_intervals"] = interface_intervals
    #create also a column for updating later the merged ones (now they have the same value)
    df["interface_intervals_merged"] = interface_intervals
    return df

#function that check if there are overlipping interfaces
def check_overlapping_interfaces(inte_1, inte_2):
    # check if every int in inte_1 are in inte_2
    if all(i in inte_2 for i in inte_1):
        return inte_2
    # check if every int in inte_2 are in inte_1
    if all(i in inte_1 for i in inte_2):
        return inte_1
    # check overlap
    if any(i in inte_2 for i in inte_1) or any(i in inte_1 for i in inte_2):
        #union of intervals (removing duplicates)
        union = inte_1 + [i for i in inte_2 if i not in inte_1]
        return sorted(union)
    #otherwise
    return None


#4 OBTAINING THE SAME OUTPUT OF ALPHABRIDGE BUT WITH MERGED INTERVALS/INTERFACES

def factorial(n):    
    if n == 0 or n ==1:        
        return 1    
    else:    
        return n * factorial(n-1)

def each_combination(start, building_combination, k, l, all_combi):          
    if len(building_combination) == k:
              all_combi.append(building_combination)           
              return   
    for i in range(start, len(l)):        
        nuova_comb = building_combination + [l[i]]            
        each_combination(i +1, nuova_comb, k, l, all_combi)

def combinations(k,l):    
    all_combi = []        
    each_combination(0,[],k,l, all_combi)    
    return all_combi


#5 FUNCTION TO GENERATE INTRA-PROTEIN EDGES

def generate_intraprotein_edges(protein_dict):
    edges = []
    for protein, intervals in protein_dict.items():
        protein_edges = []
        for i in range(len(intervals)):
            for j in range(i + 1, len(intervals)):
                # Create a connection between every pair of intervals
                protein_edges.append((str(intervals[i]), str(intervals[j])))
        edges.extend(protein_edges)
    return edges


#6 OBTAINING THE PROTEIN NETWORK AND THE NEW DATFRAME MODIFIED COMPARED WITH THE ALPHABRIDGE DATAFRAME USED AS AN INPUT

def get_binary_interaction_network(in_dir, threshold):
    # STEP 1--> putting all the information of binding interfaces in an unique dataframe    
    alpha_bridge_dict = {}
    all_data = []
    count = 1  # Counter for num_binary_combination

    # Extract the base name of the input directory (without the "fold_" prefix)
    in_dir_name = os.path.basename(in_dir).lstrip("fold_")

    # Traverse the input directory to check all subdirectories
    for subdir in os.listdir(in_dir):
        subdir_path = os.path.join(in_dir, subdir)

        # Ensure we are processing only valid directories
        if not os.path.isdir(subdir_path):
            continue

        # Check if this directory itself is "AlphaBridge" or contains it
        if "AlphaBridge" in subdir:
            # Case 1: "AlphaBridge" is an immediate subdirectory of `in_dir`
            alphabridge_path = subdir_path
            folder_name = in_dir_name  # Use the name of `in_dir` without "fold_"
        else:
            # Case 2: "AlphaBridge" is within another subdirectory
            alphabridge_path = os.path.join(subdir_path, "AlphaBridge")
            folder_name = subdir.lstrip("fold_")  # Name of the subdirectory without "fold_"
            

        # Process if the "AlphaBridge" directory exists
        if os.path.isdir(alphabridge_path):
            for file in os.listdir(alphabridge_path):
                if file.endswith(".json"):
                    json_file_path = os.path.join(alphabridge_path, file)
                    with open(json_file_path, 'r') as json_file:
                        json_data = json.load(json_file)

                    # Process JSON data and store in DataFrame
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

                    # Add additional columns for prot_1, prot_2, and binary combination
                    if rows:
                        df = pd.DataFrame(rows)

                        # Add information based on pathways
                        parts = folder_name.split('_')
                        #if len(parts) >= 3:
                        prot_1 = parts[0]
                        prot_2 = parts[1]

                        # Add additional columns to the DataFrame
                        df['prot_1'] = prot_1
                        df['prot_2'] = prot_2
                        df['num_binary_combination'] = count
                        count += 1

                        # Add the DataFrame to the dictionary and list
                        key = f"{prot_1}_{prot_2}"
                        alpha_bridge_dict[key] = df
                        all_data.append(df)

    # Combine all DataFrames into a single DataFrame for further processing
    if all_data:
        combined_df = pd.concat(all_data, ignore_index=True)
        # Add a column to track the original index
        combined_df['original_index'] = combined_df.index
        print("Successfully combined all data!")
    else:
        print("No data found.")

    #STEP 2 --> from the big datafr0	ame create more small dayaframe for a sigle protein
    prot_names = set(combined_df['prot_1']).union(set(combined_df['prot_2']))

    #dictionary for storing them
    protein_df_dict = {}
    
    for prot in prot_names:
        #whent the slected prtoein is in the first column
        prot1_rows = combined_df[combined_df['prot_1'] == prot][['interface','prot_1', 'start_1', 'end_1','num_binary_combination','prot_2','original_index']].rename(
            columns={'prot_1': 'protein', 'start_1': 'start', 'end_1': 'end', 'prot_2':'bin_comb_with'}
        )
        # when the selected protein is in the sec column
        prot2_rows = combined_df[combined_df['prot_2'] == prot][['interface','prot_2', 'start_2', 'end_2','num_binary_combination','prot_1','original_index']].rename(
            columns={'prot_2': 'protein', 'start_2': 'start', 'end_2': 'end','prot_1':'bin_comb_with'}
        )
        #combine the two rows information
        combined_rows = pd.concat([prot1_rows, prot2_rows], ignore_index=True)
        #add the combined one as the value of the dictionary
        protein_df_dict[prot] = combined_rows

    #adding new columns for updating the intervals later (so we can keep track of the orginal ones)
    for prot, df in protein_df_dict.items():
            df['start_merged'] = df['start']
            df['end_merged'] = df['end']

    print(f"protein_df_dict AFTER STEP 2 (from the big to the small df)-->{protein_df_dict}")

    # STEP 3 --> function for mergin intervals in the big dataframe

    for prot, df in protein_df_dict.items():
        df = get_merging_intervals(df,3)

    print(f"protein_df_dict AFTER STEP 3(mergin intervals in the big dataframe)-->{protein_df_dict}")
    #STEP 4 --> merge intervals in interfaces intervals, here you keep track of the single interface for a single binary interaction
    for prot, df in protein_df_dict.items():
        df = create_new_column_interface_intervals(df)

    print(f"protein_df_dict AFTER STEP 4(merge intervals in interfaces intervals)-->{protein_df_dict}")
    
    #STEP 5 --> check if there are overlapping interfaces in the dataframe of the signle protein (here you check all togheter bc it could be)
    #that "different" interfaces in different binary interaction are actually overlapping
    #apply the function in the dataframe
    for prot, df in protein_df_dict.items():
        for i in range(len(df['protein'])):
            interval_1 = df.loc[i, 'interface_intervals_merged']
            for j in range(len(df['protein'])-1):
                interval_2 = df.loc[j+1, 'interface_intervals_merged']
                new_interfaces = check_overlapping_interfaces(interval_1, interval_2)
                if new_interfaces == None:
                    pass
                else:
                    df.at[i,'interface_intervals_merged'] = new_interfaces
                    df.at[j+1,'interface_intervals_merged'] = new_interfaces

    print(f"protein_df_dict AFTER STEP 5(check if there are overlapping interfaces)-->{protein_df_dict}")


    #STEP 6 --> obtaining again the output of alphabridge but merged
    
    list_of_possible_combination = []
    for key, value in alpha_bridge_dict.items():
        list_of_possible_combination.append(key)
    possible_combinations = [key.split("_") for key in list_of_possible_combination]

    merged_protein_dataframe_dict = {}

    for prot1, prot2 in possible_combinations:
        if prot1 in protein_df_dict and prot2 in protein_df_dict:
            df1 = protein_df_dict[prot1]
            df2 = protein_df_dict[prot2]

            merged_df = pd.merge(df1, df2, how="inner", on = 'original_index', suffixes=('_1', '_2'))

            key = f"{prot1}_{prot2}"
            merged_protein_dataframe_dict[key] = merged_df
            for key, df in merged_protein_dataframe_dict.items():
                df["interface_intervals_mergerd_1_labels"] = df["interface_intervals_merged_1"].apply(lambda x: "_".join(map(str, x)) if isinstance(x, list) else str(x)) + "_" + df["protein_1"].astype(str)
                df["interface_intervals_mergerd_2_labels"] = df["interface_intervals_merged_2"].apply(lambda x: "_".join(map(str, x)) if isinstance(x, list) else str(x)) + "_" + df["protein_2"].astype(str)
            
                merged_protein_dataframe_dict[key] = df
 
        else:
            print(f"missing datframe: {prot1}, {prot2}")

    

    print(f"merged_protein_dataframe_dict  STEP 6(obtaining again the output of alphabridge but merged-->{merged_protein_dataframe_dict}")


    #STEP 7 --> comparing the dataframe before and after merging
    comparison_dataframe_dict = {}

    # Merge the dataframes from alpha_bridge_dict and merged_dataframe by matching keys
    for key in alpha_bridge_dict:
        if key in merged_protein_dataframe_dict:
            # Concatenate the dataframes horizontally
            concatenated_df = pd.concat([alpha_bridge_dict[key], merged_protein_dataframe_dict[key]], axis=1)
            comparison_dataframe_dict[key] = concatenated_df

    #STEP 8--> START TO TAKE INFORMATION FOR PLOTTING

    # Initialize the dictionary to store proteins and their unique nodes
    protein_nodes = {}

    # Iterating through each row of the DataFrame
    for protein, df in merged_protein_dataframe_dict.items():
        for _, row in df.iterrows():
            # Extract protein IDs and their corresponding residue intervals
            prot_1_id = row['protein_1']
            prot_2_id = row['protein_2']
            prot_1_interface = row["interface_intervals_mergerd_1_labels"] 
            prot_2_interface = row["interface_intervals_mergerd_2_labels"] 

            # Initialize the protein in the dictionary if not already present
            if prot_1_id not in protein_nodes:
                protein_nodes[prot_1_id] = []
            if prot_2_id not in protein_nodes:
                protein_nodes[prot_2_id] = []

            # Add unique intervals for prot_1_id
            if prot_1_interface not in protein_nodes[prot_1_id]:  # Check if the interval is unique
                protein_nodes[prot_1_id].append(prot_1_interface)

            # Add unique intervals for prot_2_id
            if prot_2_interface not in protein_nodes[prot_2_id]:  # Check if the interval is unique
                protein_nodes[prot_2_id].append(prot_2_interface)

    print(f'protein nodes {protein_nodes}')

    
    
    # Calculate the number of unique nodes (intervals) for each protein
    nodes_for_each_protein = [len(value) for value in protein_nodes.values()]

    # INTERACTION BETWEEN DIFFERENT PROTEINS
    inter_protein_interactions = []

    # Iterate over each row in the dataframe
    for protein, df in merged_protein_dataframe_dict.items():
        for _, row in df.iterrows():
            # Extract data for prot_1 and prot_2
            #prot_1_id = row["protein_1"]
            aa_prot_1 = row["interface_intervals_mergerd_1_labels"]
            #prot_2_id = row["protein_2"]
            aa_prot_2 = row["interface_intervals_mergerd_2_labels"]

            # Create a list of interactions between protein pairs
            interaction = (aa_prot_1, aa_prot_2)
            inter_protein_interactions.append(interaction)

    # Print all protein interactions
    print(f"All protein interactions: {inter_protein_interactions}")

    # To get unique interactions, use a set
    unique_inter_protein_interactions = []

    for interaction in inter_protein_interactions:
        if interaction not in unique_inter_protein_interactions:
            unique_inter_protein_interactions.append(interaction)

    # Print unique protein interactions
    print(f"Unique protein interactions: {unique_inter_protein_interactions}")


    #INTERACTION WHITIN THE SAME PROTEIN
    edges_intraprpt = generate_intraprotein_edges(protein_nodes)
    print(f'edges_intraprt {edges_intraprpt}')

    
    #STEP 9 --> PLOTTING
    g = ig.Graph()
    number_of_nodes = 0
    labels = []

    
    #PROTEINS for the legend
    proteins_leg = []
    for key, value in protein_nodes.items():
        proteins_leg.append(key)

    
    #NODES/LABELS
    #calculate the number of nodes that will be plotted
    for nodes in nodes_for_each_protein:
        number_of_nodes = number_of_nodes + nodes
        
    #labels--> each node has intervals of aa
    for key, value in protein_nodes.items():
        for interval_aa_str in value:
            labels.append(interval_aa_str)

    print(f"labels: {labels}")
    print(f"number of nodes:{number_of_nodes}")
    
    #adding nodes
    g.add_vertices(number_of_nodes)
    
    #assign labels to each node (the order follow the dictionary protein_nodes)
    g.vs['label'] = labels
    for vertex in g.vs:     
        print(f"ID: {vertex.index}, Label: {vertex['label']}")


    #COLORS 
    # Generate as many unique colors as the number of protein groups
    num_proteins = len(nodes_for_each_protein)
    cmap = plt.get_cmap("tab10")  # 'tab10' colormap has good distinct colors, can handle up to 10 easily
    colors = [mcolors.rgb2hex(cmap(i)) for i in np.linspace(0, 1, num_proteins)]
    
    # Assign colors to nodes based on their protein group
    node_colors = []
    for protein_index, size in enumerate(nodes_for_each_protein):
        node_colors.extend([colors[protein_index]] * size) 
        
    # Assign colors to graph nodes
    g.vs["color"] = node_colors

    
    ## EDGES
    # different interfaces
    edges_diff = []
    for source, target in unique_inter_protein_interactions:
        # Convert source and target to match the label format
        source_index = g.vs.find(label=f"{source}").index
        target_index = g.vs.find(label=f"{target}").index
        edges_diff.append((source_index, target_index))
    print(f"edges diff:{edges_diff}")
    g.add_edges(edges_diff)
    
    # same protein edges
    edges_same = []
    for key, value in edges_intraprpt:
        # Convert source and target to match the label format
        source_index_1 = g.vs.find(label=f"{key}").index
        target_index_1 = g.vs.find(label=f"{value}").index
        edges_same.append((source_index_1, target_index_1))
    print(f"edges same:{edges_same}")
    g.add_edges(edges_same)
    g.es['weight'] =[500.0] *g.ecount()
    g.es['color'] = "black"
    g.es["widht"] = 100
    for edge in edges_same:
            edges_id = g.get_eid(edge[0],edge[1])
            g.es[edges_id]["weight"]= 1.0
            g.es[edges_id]["color"]= "light gray"
            g.es[edges_id]["widht"]= 0.5

   #LAYAOUT

    layout = g.layout("fruchterman_reingold")

    
   #PLOT 
    
    #fig, ax = plt.subplots(figsize=(10, 10))
    #ig.plot(g,target = ax, layout = layout, bbox =(800,800), margin = 20, vertex_color = g.vs['color'],vertex_size=10, edge_weight = g.es["weight"], vertex_label = labels, color =g.es["color"], widht = g.es["widht"])


    #LEGEND 

    #nodes_patches = [mpatches.Patch(color=colors[i], label=proteins_leg[i]) for i in range(len(proteins_leg))]
    #edge_patches = [
    #    mpatches.Patch(color='lightgray', label='Same Protein Connections'),
    #    mpatches.Patch(color='black', label='Interface Connections')
    #]
    #plt.legend(handles=nodes_patches + edge_patches, loc="upper right", title="Proteins and Interfaces")
    
    #network in network x
    g_networkx = g.to_networkx()
    
    #informaton for the json file
    jobs = json_graph.node_link_data(g_networkx)

    return jobs


def main():
    args = parse_args()
    in_dir = args.in_dir
    threshold = float(args.threshold)
    
    # Function that creates the network and the comparison dataframe
    jobs = get_binary_interaction_network(in_dir, threshold)

    # Save to a JSON file
    with open(f"{args.o_dir}/combined_networks.json", "w") as f:
        json.dump(jobs, f, indent=4)


if __name__ == '__main__':
    main()

