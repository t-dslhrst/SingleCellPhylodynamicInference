import networkx as nx 
import numpy as np
from Bio import Phylo

tree = Phylo.read("example_data.tree", "newick")

branch_lengths = [clade.branch_length for clade in tree.find_clades()][1:]

# Initialize graph
G = nx.Graph()

# Recursive function to add clades and edges to the graph
def add_clades_to_graph(clade, parent=None):
    # Add the clade as a node
    G.add_node(clade)

    if parent is not None:
        # Add an edge from parent to this clade with the branch length + 1 as weight
        branch_length = clade.branch_length if clade.branch_length is not None else 0
        G.add_edge(parent, clade, weight=int(branch_length+1))

    # Recursively add child clades
    for child in clade.clades:
        add_clades_to_graph(child, clade)

# Start adding clades from the root
add_clades_to_graph(tree.root)
weights_p1 = [data['weight'] for _, _, data in G.edges(data=True)]
adj = nx.to_numpy_array(G, dtype=int)


np.savetxt("example_data_adj_matrix_p1.csv", adj, delimiter=",", fmt='%i')