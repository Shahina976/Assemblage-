#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
import textwrap
import matplotlib.pyplot as plt
matplotlib.use("Agg")

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"

def isfile(path): # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file
    
    :raises ArgumentTypeError: If file doesn't exist
    
    :return: (str) Path 
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="k-mer size (default 22)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file (default contigs.fasta)")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as an image (png)")
    return parser.parse_args()




def read_fastq(fastq_file):
    """Extract reads from fastq files.

    :param fastq_file: (str) Path to the fastq file.
    :return: A generator object that iterate the read sequences. 
    """

    with open(fastq_file, "r") as file:
        for line in file:
            yield(next(file).strip())
            next(file)
            next(file)



def cut_kmer(read, kmer_size):
    """Cut read into kmers of size kmer_size.
    
    :param read: (str) Sequence of a read.
    :return: A generator object that iterate the kmers of of size kmer_size.
    """
    for i in range(0, len(read) - kmer_size + 1):
        kmer = read[i:i + kmer_size]
        yield kmer



def build_kmer_dict(fastq_file, kmer_size):
    """Build a dictionary object of all kmer occurrences in the fastq file.

    :param fastq_file: (str) Path to the fastq file.
    :param kmer_size: (int) Length of the k-mers.
    :return: A dictionary object that identifies all kmer occurrences.
    """
    kmer_dict = {}  # Dictionary to store k-mers and their occurrences

    for read in read_fastq(fastq_file):
        for kmer in cut_kmer(read, kmer_size):
            if kmer in kmer_dict:
                kmer_dict[kmer] += 1  # Increment count for existing k-mer
            else:
                kmer_dict[kmer] = 1  # Add new k-mer to the dictionary

    return kmer_dict



def build_graph(kmer_dict):
    """Build the debruijn graph

    :param kmer_dict: A dictionnary object that identify all kmer occurrences.
    :return: A directed graph (nx) of all kmer substring and weight (occurrence).
    """
    # Créer un graphe dirigé
    graph = nx.DiGraph()

    # Parcourir chaque k-mer et son occurrence dans le dictionnaire
    for kmer, occurrence in kmer_dict.items():
        # Extraire le préfixe et le suffixe du k-mer
        prefix = kmer[:-1]
        suffix = kmer[1:]

        # Vérifier si le préfixe et le suffixe existent déjà dans le graphe
        if not graph.has_node(prefix):
            graph.add_node(prefix)
        if not graph.has_node(suffix):
            graph.add_node(suffix)

        # Ajouter un arc dirigé du préfixe au suffixe avec le poids (occurrence)
        graph.add_edge(prefix, suffix, weight=occurrence)

    return graph



def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """Remove a list of path in a graph. A path is set of connected node in
    the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of path
    :param delete_entry_node: (boolean) True->We remove the first node of a path 
    :param delete_sink_node: (boolean) True->We remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
    pass


def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    """Select the best path between different paths

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of path
    :param path_length_list: (list) A list of length of each path
    :param weight_avg_list: (list) A list of average weight of each path
    :param delete_entry_node: (boolean) True->We remove the first node of a path 
    :param delete_sink_node: (boolean) True->We remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
    pass

def path_average_weight(graph, path):
    """Compute the weight of a path

    :param graph: (nx.DiGraph) A directed graph object
    :param path: (list) A path consist of a list of nodes
    :return: (float) The average weight of a path
    """
    return statistics.mean([d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)])

def solve_bubble(graph, ancestor_node, descendant_node):
    """Explore and solve bubble issue

    :param graph: (nx.DiGraph) A directed graph object
    :param ancestor_node: (str) An upstream node in the graph 
    :param descendant_node: (str) A downstream node in the graph
    :return: (nx.DiGraph) A directed graph object
    """
    pass

def simplify_bubbles(graph):
    """Detect and explode bubbles

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    pass

def solve_entry_tips(graph, starting_nodes):
    """Remove entry tips

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    pass

def solve_out_tips(graph, ending_nodes):
    """Remove out tips

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    pass

def get_starting_nodes(graph):
    """Obtient les nœuds sans prédécesseurs (nœuds d'entrée).

    :param graph: Un objet de graphe dirigé (DiGraph).
    
    :return: Une liste de nœuds sans prédécesseurs.
    """
    starting_nodes = [node for node in graph.nodes() if not any(graph.predecessors(node))]
    return starting_nodes


def get_sink_nodes(graph):
    """Obtient les nœuds de sortie (nœuds sans successeurs).

    :param graph: Un objet de graphe dirigé (DiGraph).
    
    :return: Une liste de nœuds sans successeurs.
    """
    sink_nodes = []
    for node in graph.nodes():
        if not any(graph.successors(node)):
            sink_nodes.append(node)
    return sink_nodes


def get_contigs(graph, starting_nodes, ending_nodes):
    """Extract the contigs from the graph.

    :param graph: A directed graph object (nx.DiGraph).
    :param starting_nodes: A list of nodes without predecessors.
    :param ending_nodes: A list of nodes without successors.
    :return: List of tuples, each containing a contiguous sequence and its length.
    """
    contigs = []

    for start_node in starting_nodes:
        for end_node in ending_nodes:
            if nx.has_path(graph, start_node, end_node):
                paths = nx.all_simple_paths(graph, source=start_node, target=end_node)
                for path in paths:
                    contig = path[0]
                    for i in range(1, len(path)):
                        contig += path[i][-1]
                    contig_length = len(contig)
                    contigs.append((contig, contig_length))
    return contigs



def save_contigs(contigs_list, output_file):
    """Write all contigs in fasta format

    :param contig_list: (list) List of [contiguous sequence and their length]
    :param output_file: (str) Path to the output file
    """
    pass


def draw_graph(graph, graphimg_file): # pragma: no cover
    """Draw the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param graphimg_file: (str) Path to the output file
    """                                   
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


#==============================================================
# Main program
#==============================================================
def main(): # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit 
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)

        # 1) fonction read_fastq()
    kmer_size = 7 
    fastq = read_fastq(args.fastq_file)
    print(fastq)
    print("")

    # 2) Lire toutes les séquences du fichier FASTQ et couper avec cut_kmer()
    for read in fastq:
        kmer_coupe = list(cut_kmer(read, kmer_size))
        print(kmer_coupe)
        print("")

    # 3) fonction dict_kmer
    dict_kmer = build_kmer_dict(args.fastq_file, kmer_size)
    print(dict_kmer)
    print("")

    # 4) fonction igraph
    # fonction pour construire le graphe de De Bruijn
    debruijn_graph = build_graph(dict_kmer)

    graph_length = len(debruijn_graph)
    num_nodes = len(debruijn_graph.nodes)
    num_edges = len(debruijn_graph.edges)

    print(f"Length of the De Bruijn Graph: {graph_length}")
    print(f"Number of nodes: {num_nodes}")
    print(f"Number of edges: {num_edges}")
    print("")
    
    # # affichage du graphe 
    # pos = nx.spring_layout(debruijn_graph)
    # nx.draw(debruijn_graph, pos, with_labels = False, node_size = 10, node_color='b', font_size=8)
    # plt.savefig("igraph.png")
    # plt.close

    # 5) et 6) fonction get_starting_nodes et get_sink_nodes
    starting_nodes = get_starting_nodes(debruijn_graph)
    print("entrance nodes : ", starting_nodes)
    ending_nodes = get_sink_nodes(debruijn_graph)
    print("output nodes : ", ending_nodes)
    print("")

    # 7) fonction get_contigs
    contigs = get_contigs(debruijn_graph, starting_nodes, ending_nodes)
    for contig, length in contigs : 
        print(f"Contig : {contig}, Longueur : {length}")
        print("")


if __name__ == '__main__': # pragma: no cover
    main()









