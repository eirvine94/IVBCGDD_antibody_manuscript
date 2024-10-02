#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Edward B. Irvine

Description:
This script reads in adjacency list data from CSV files, generates network graphs using the `networkx` library, 
and visualizes them based on correlation values. The `net_from_adj_list` function is defined to create and plot network graphs with specified layouts, node sizes, and font sizes. 
The script processes three datasets ('BAL_pre', 'plasma_pre', and 'plasma_post') to produce circular network graphs, saving each visualization as a PDF file.
    

Created: 06 January 2021

Modified: 02 October 2024
"""

################
### Housekeeping --------------------------------------------------------------------------------------------------------------------------------------------------------------
################

# Import required libraries 
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

# Read in data
BAL_pre = pd.read_csv("eth_BAL_pre_correlates_threshold.csv")
plasma_pre = pd.read_csv("eth_plasma_pre_correlates_threshold.csv")
plasma_post = pd.read_csv("eth_plasma_post_correlates_threshold.csv")





######################
### Define function(s) ------------------------------------------------------------------------------------------------------------------------------------------------------------
######################

def net_from_adj_list(adj_list, layout, node_size, font_size, output_file):
    """
    Generates and visualizes a network graph from an adjacency list with edge weights and colors 
    based on correlation values. The function also saves the graph to a PDF file.
    
    Parameters:
    - adj_list (DataFrame): A pandas DataFrame containing an adjacency list with columns:
        * 'Var1': Source nodes
        * 'Var2': Target nodes
        * A correlation value for each edge
    - layout (str): The layout for the graph visualization ('kamada_kawai', 'spring', or 'circular').
    - node_size (int): Size of the nodes in the graph.
    - font_size (int): Size of the labels in the graph.
    - output_file (str): Name of the output file (without extension) where the graph will be saved.
    
    Returns:
    - G (networkx.Graph): The generated network graph object.
    """
    # Initialize graph
    G = nx.Graph()
    
    # Add nodes
    nodes = set(adj_list.Var1.unique()).union(set(adj_list.Var2.unique()))
    G.add_nodes_from(nodes)
    
    # Loop through each line
    for row in adj_list.itertuples():
        
        # Store parameters
        source_node = row[2]
        target_node = row[3]
        spearman = float(row[4])
       
        # Add edge (red if positive correlation, blue if negative) 
        if spearman > 0:
            G.add_edge(source_node, target_node, weight = spearman*10, color = "#AAAAAA")
        else: 
            G.add_edge(source_node, target_node, weight = spearman*10, color = "#005493")
    
    # Creates lists for edges and weights
    edges, weights = zip(*nx.get_edge_attributes(G, 'weight').items())
    edges, colors = zip(*nx.get_edge_attributes(G, 'color').items())
    
    # Set graph layout
    if layout == "kamada_kawai":
        positions = nx.kamada_kawai_layout(G, weight = weights, scale=100)
    elif layout == "spring":
        positions = nx.spring_layout(G)
    elif layout == "circular":
        positions = nx.circular_layout(G)
    else:
        positions = nx.circular_layout(G)

    # Figure size
    plt.figure(figsize=(20,20))

    # Plot nodes
    nx.draw_networkx_nodes(G, 
                           pos = positions, 
                           node_size = node_size,
                           node_color = "#AAAAAA") 
    
    # Style labels
    nx.draw_networkx_labels(G, pos = positions, font_size = font_size, font_family = 'sans-serif')
    
    # Draw edges
    nx.draw_networkx_edges(G, 
                           pos = positions, 
                           style='solid',
                           width = weights, 
                           edge_color = colors, 
                           edge_vmin = min(weights), 
                           edge_vmax = max(weights))

    # Display graph without axis
    plt.axis('off')
    
    # Save image
    plt.savefig(output_file + ".pdf", format = "PDF")
    return(G)
       




########
### Main ------------------------------------------------------------------------------------------------------------------------------------------------------------
########
net_BAL_pre_circle = net_from_adj_list(adj_list = BAL_pre, layout = "circular", node_size = 9000, font_size = 25, output_file = "eth_net_BAL_pre_circle")
net_plasma_pre_circle = net_from_adj_list(adj_list = plasma_pre, layout = "circular", node_size = 9000, font_size = 25, output_file = "eth_net_plasma_pre_circle")
net_plasma_post_circle = net_from_adj_list(adj_list = plasma_post, layout = "circular", node_size = 9000, font_size = 25, output_file = "eth_net_plasma_post_circle")
