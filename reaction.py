import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

def draw_reaction_graph(filename):
    # Read the network data from an Excel file
    nw = pd.read_excel(filename)
    num_of_nodes = len(nw['Nodes'])
    node_names = nw['Nodes'].tolist()

    # Initialize a directed graph
    G = nx.DiGraph()

    # Add nodes
    for node in node_names:
        G.add_node(node)

    # Add edges for activation and inhibition
    for i in range(num_of_nodes):
        activators = str(nw['Activators'][i]).split(',')
        inhibitors = str(nw['Inhibitors'][i]).split(',')

        # Clean the node names from spaces and check if they are not NaN
        activators = [act.strip() for act in activators if act.strip() in node_names]
        inhibitors = [inh.strip() for inh in inhibitors if inh.strip() in node_names]

        current_node = node_names[i]

        # Add edges for each activator
        for act in activators:
            if act:  # Check to make sure 'act' is not an empty string
                G.add_edge(act, current_node, color='green', relationship='activates')

        # Add edges for each inhibitor
        for inh in inhibitors:
            if inh:  # Check to make sure 'inh' is not an empty string
                G.add_edge(inh, current_node, color='red', relationship='inhibits')

    # Use a different layout to spread nodes more effectively
    pos = nx.kamada_kawai_layout(G)  # You can also try nx.circular_layout or nx.fruchterman_reingold_layout

    # Draw the graph
    plt.figure(figsize=(24, 24))  # Increase figure size
    edges = G.edges(data=True)
    colors = [data['color'] for u, v, data in edges]
    nx.draw(G, pos, edge_color=colors, with_labels=True, node_color='lightblue', node_size=3000, font_size=12)  # Adjust node size and font size
    edge_labels = dict([((u, v,), d['relationship'])
                        for u, v, d in G.edges(data=True)])
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_color='blue')  # Adjust font color for visibility

    plt.title('Reaction Graph from Network Data')
    plt.savefig('reaction.png')  # Save as a larger image
    plt.show()

# Example usage:
draw_reaction_graph('CRT.xlsx')
