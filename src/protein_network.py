#script for obtaining the protein interaction network
#EVERY NODE IS A PROTEIN
#EVERY EDGE IS A PHYSICAL INTERECTION BETWEEN THE PROTEINS

#import required packages

import igraph as ig
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from networkx.readwrite import json_graph;

#FUCNTION USED IN NETWORK WITH NOT MERGIN NODES AND THE PROTEIN==NODE NETWORK

#formatting labels --> display the label nodes as "protein name (form aa,to aa) (form aa,to aa)""

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

#all the intervals that belong to the same interface --> togheter

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

'''MAIN FUNCTION!!!'''

 #STEP1--> HAVING THE DF WITH ALL POSSIBLE INFORMATION
def get_protein_network(df):
    df_interfaces = create_new_column_interface_intervals_no_merge(df)
    df_interfaces["interface_intervals_1_labels"] = df["prot_1"]+" "+ df_interfaces["interface_intervals_1"].apply(lambda x: " ".join(map(str, x)) if isinstance(x, list) else str(x)).astype(str)
    df_interfaces["interface_intervals_1_labels"] = df_interfaces["interface_intervals_1_labels"].apply(formatting_labels)
    df_interfaces["interface_intervals_2_labels"] = df["prot_2"]+ " "+ df_interfaces["interface_intervals_2"].apply(lambda x: " ".join(map(str, x)) if isinstance(x, list) else str(x)).astype(str)
    df_interfaces["interface_intervals_2_labels"] = df_interfaces["interface_intervals_2_labels"].apply(formatting_labels)

    # Initialize the dictionary to store proteins and their unique nodes
    protein_labels = {}

    # Iterating through each row of the DataFrame
    for _, row in df_interfaces.iterrows():
        # Extract protein IDs and their corresponding residue intervals
        prot_1_id = row['prot_1']
        prot_2_id = row['prot_2']
        prot_1_interface = row["interface_intervals_1_labels"] 
        prot_2_interface = row["interface_intervals_2_labels"] 

        # Initialize the protein in the dictionary if not already present
        if prot_1_id not in protein_labels:
            protein_labels[prot_1_id] = []
        if prot_2_id not in protein_labels:
            protein_labels[prot_2_id] = []

        # Add unique intervals for prot_1_id
        if prot_1_interface not in protein_labels[prot_1_id]:  
            protein_labels[prot_1_id].append(prot_1_interface)

        # Add unique intervals for prot_2_id
        if prot_2_interface not in protein_labels[prot_2_id]:  
            protein_labels[prot_2_id].append(prot_2_interface)

    #print(f'Protein labels before formatting: {protein_labels}')

    # Format labels correctly
    for protein in protein_labels:
        intervals = protein_labels[protein]
        #print(protein)
        
        # take out the name of the protein in the other string (I only want it once)
        formatted_intervals = [intervals[0]]  # only the first interval keep he name
        #print(formatted_intervals)
        #for the rest of them:
        for interval in intervals[1:]:
            #replace protein == name of the proten with nothing
            formatted_intervals.append(interval.replace(f"{protein} ", "")) 
            #print(formatted_intervals)
        
        # formatted labels --> sring separated by , --> joing them in an unique one
        protein_labels[protein] = ", ".join(formatted_intervals)

    #print(f'Protein labels after formatting: {protein_labels}')

    #orering them --> for having the same color as the other network
    protein_labels_sorted = dict(sorted(protein_labels.items()))

    #print(f'Protein labels after sorting: {protein_labels_sorted}')

    # INTERACTIONS with name of the protein
    protein_pairs = []
    protein_pairs_unique = []

    for prot1, prot2 in zip(df_interfaces["prot_1"], df_interfaces["prot_2"]):
        pair = [prot1, prot2]
        protein_pairs.append(pair)

    for protein in protein_pairs:
        if protein not in protein_pairs_unique:
            protein_pairs_unique.append(protein)

    #print(protein_pairs)
    #print(protein_pairs_unique)

    #print("All interactions:", protein_pairs)
    #print("Unique interactions:", protein_pairs_unique)

    # INTERACTIONS with the all label of the protein
    protein_pairs_unique_labels = []

    for prot1, prot2 in protein_pairs_unique:
        labels_1 = protein_labels.get(prot1, "No Label")
        labels_2 = protein_labels.get(prot2, "No Label")
        protein_pairs_unique_labels.append([labels_1, labels_2])

    #print("Unique interactions with labels:", protein_pairs_unique_labels)
    g = ig.Graph()
    number_of_nodes = (len(protein_labels_sorted))
    labels = []
    #adding nodes
    g.add_vertices(number_of_nodes)

    #assign labels
    for key, value in protein_labels_sorted.items():
        labels.append(value)

    g.vs['label'] = labels

    #for vertex in g.vs:     
    #    print(f"ID: {vertex.index}, Label: {vertex['label']}")

    #COLORS 
    # Generate as many unique colors as the number of protein groups
    cmap = plt.get_cmap("rainbow") 
    colors = [mcolors.rgb2hex(cmap(i)) for i in np.linspace(0, 1, number_of_nodes)]
        
    # Assign colors to graph nodes
    g.vs["color"] = colors
    #print(f"node colors{colors}")

    #edges
    edges = []
    for source, target in protein_pairs_unique_labels:
        # Convert source and target to match the label format
        source_index = g.vs.find(label=f"{source}").index
        target_index = g.vs.find(label=f"{target}").index
        edges.append((source_index, target_index))
    #print(f"edges :{edges}")
    g.add_edges(edges)
    #network in network x
    g_networkx = g.to_networkx()
    
    #informaton for the json file
    jobs = json_graph.node_link_data(g_networkx)

    return jobs

