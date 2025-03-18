#script for obtaining ONLY THE INTERACTIVE NETWORK WITH NO MERGE(more nodes)
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
#FUCNTION USED IN BOTH NETWORK WITH MERGING NODES AND NOT MERGIN NODES
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

#function for modify the labels
def formatting_labels(s):
    parts = s.split(" ")
    result = [parts[0]]  
    temp_list = []

    for part in parts[1:]:  
        temp_list.append(part[:-1])
    for i in range(len(temp_list)-1):
        if i%2 == 0:
            temp_list[i]=temp_list[i].strip("(")
            if temp_list[i] == temp_list[i+1]:
                result.append(f"({temp_list[i] })")
            else:
                result.append(f"({temp_list[i]}-{temp_list[i+1]})")
    return " ".join(result)

#FUCTION THAT ARE NEEDED IN THE MAIN FUNCTION FOR THE NETWORK WITH NO MERGIN NODES
def create_new_column_interface_intervals_no_merge(df):
    grouped = df.groupby(["interface"])
    
    # List for intervals for each interface (start_1, end_1)
    interface_intervals_1 = [None] * len(df)
    for _, group in grouped:
        # Convert to a list of tuples
        intervals = [tuple(x) for x in group[["start_1", "end_1"]].drop_duplicates().values]
        for idx in group.index:
            interface_intervals_1[idx] = intervals
    df["interface_intervals_1"] = interface_intervals_1

    grouped = df.groupby(["interface"])
    
    # List for intervals for each interface (start_2, end_2)
    interface_intervals_2 = [None] * len(df)
    for _, group in grouped:
        # Convert to a list of tuples
        intervals = [tuple(x) for x in group[["start_2", "end_2"]].drop_duplicates().values]
        for idx in group.index:
            interface_intervals_2[idx] = intervals
    df["interface_intervals_2"] = interface_intervals_2

    return df

'''MAIN FUNCTION'''
def get_protein_network_no_merging(df):
    #df_interfaces_nomerg --> DF THAT IS USED FOR OBTAINING THE NO MERGED ETWORK (==MORE NODES)

    df_interfaces = create_new_column_interface_intervals_no_merge(df)
    # Save the DataFrame with the updated protein names
    df_interfaces["interface_intervals_1_labels"] = df["prot_1"]+" "+ df_interfaces["interface_intervals_1"].apply(lambda x: " ".join(map(str, x)) if isinstance(x, list) else str(x)).astype(str)
    df_interfaces["interface_intervals_1_labels"] = df_interfaces["interface_intervals_1_labels"].apply(formatting_labels)
    df_interfaces["interface_intervals_2_labels"] = df["prot_2"]+ " "+ df_interfaces["interface_intervals_2"].apply(lambda x: " ".join(map(str, x)) if isinstance(x, list) else str(x)).astype(str)
    df_interfaces["interface_intervals_2_labels"] = df_interfaces["interface_intervals_2_labels"].apply(formatting_labels)

    #df_interfaces.to_csv(output_csv, index=False)

    #print(df_interfaces)
    #START TO TAKE INFORMATION FOR PLOTTING

    # Initialize the dictionary to store proteins and their unique nodes
    protein_nodes_old = {}

    # Iterating through each row of the DataFrame
    
    for _, row in df_interfaces.iterrows():
        # Extract protein IDs and their corresponding residue intervals
        prot_1_id = row['prot_1']
        prot_2_id = row['prot_2']
        prot_1_interface = row["interface_intervals_1_labels"] 
        prot_2_interface = row["interface_intervals_2_labels"] 

        # Initialize the protein in the dictionary if not already present
        if prot_1_id not in protein_nodes_old:
            protein_nodes_old[prot_1_id] = []
        if prot_2_id not in protein_nodes_old:
            protein_nodes_old[prot_2_id] = []

        # Add unique intervals for prot_1_id
        if prot_1_interface not in protein_nodes_old[prot_1_id]:  # Check if the interval is unique
            protein_nodes_old[prot_1_id].append(prot_1_interface)

        # Add unique intervals for prot_2_id
        if prot_2_interface not in protein_nodes_old[prot_2_id]:  # Check if the interval is unique
            protein_nodes_old[prot_2_id].append(prot_2_interface)

    print(f'protein nodes old {protein_nodes_old}')

    protein_nodes = dict(sorted(protein_nodes_old.items()))
    print(f'protein nodes  {protein_nodes}')

    

    # Calculate the number of unique nodes (intervals) for each protein
    nodes_for_each_protein = [len(value) for value in protein_nodes.values()]

    # INTERACTION BETWEEN DIFFERENT PROTEINS
    inter_protein_interactions = []

    # Iterate over each row in the dataframe
    for _, row in df_interfaces.iterrows():
        # Extract data for prot_1 and prot_2
        #prot_1_id = row["protein_1"]
        aa_prot_1 = row["interface_intervals_1_labels"]
        #prot_2_id = row["protein_2"]
        aa_prot_2 = row["interface_intervals_2_labels"]

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
    
    #assign labels to each node (the order follow the dictionary protein_nodes)order
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
    print(f"node colors{node_colors}")

    
    #EDGES
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
    #print(f"edges sapython3 define_interfaces.py -i path-to-AlphaFold3-folderme:{edges_same}")
    g.add_edges(edges_same)
    #g.es['weight'] =[500.0] *g.ecount()
    g.es['weight'] = 500
    g.es['color'] = "black"
    g.es["width"] = 1.5
    #iterating all ossible edges
    for edge in edges_same:
        #for every tuple (=edge inside the list) we take the FROM node and the TO node
        source_index, target_index = edge
        #get ID f the edge
        edges_id = g.get_eid(source_index, target_index)
        
        # take the actual color of the node
        node_color = g.vs[source_index]["color"]
        
        g.es[edges_id]["weight"] = 3
        g.es[edges_id]["color"] = node_color  #we give for each edge the color of his node
        g.es[edges_id]["width"] = 0.5


   #LAYAOUT

    layout = g.layout("fruchterman_reingold")

    
   #PLOT 
    
    #fig, ax = plt.subplots(figsize=(10, 10))
    #ig.plot(g,target = ax, layout = layout, bbox =(800,800), margin = 20, vertex_color = g.vs['color'],vertex_size=10, edge_weight = g.es["weight"], vertex_label = None, color =g.es["color"], widht = g.es["widht"])


    #LEGEND 

    #PROTEINS for the legend
    #proteins_leg = []
    #for key, value in protein_nodes.items():
    #    proteins_leg.append(key)

        #creating the content of the first legend
    #nodes_patches = [mlines.Line2D([0], [0], marker='o', color='w',
    #                            markerfacecolor=colors[i], markersize=12,
    #                            label=proteins_leg[i]) for i in range(len(proteins_leg))]
    
    #creating the content of the second legend
    #edge_patches = [
    #    mpatches.Patch(color='lightgray', label='Same Protein Connections'),
    #    mpatches.Patch(color='black', label='Interface Connections')
    #]
    
    #first legend
    #first_legend = ax.legend(handles=nodes_patches, loc="upper right")
    #ax.add_artist(first_legend)
    
    #second legend
    #ax.legend(handles=edge_patches, loc="upper left", title=r"$\bf{Interfaces}$")

    #network in network x
    g_networkx = g.to_networkx()
    
    #informaton for the json file
    jobs = json_graph.node_link_data(g_networkx)

    return jobs

#g = get_protein_network_no_merging('/local_data/carmen/complexes_trial_old/complex1_pfdn1_vbp1_pfdn4_pfdn5', 0.75)
#print(g)