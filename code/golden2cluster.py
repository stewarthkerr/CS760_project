#!/usr/bin/env python3
#Author: Heloise Letissier

import argparse
import csv
import json
import networkx as nx
import os
import pandas as pd
import sys
import urllib
import urllib.request as urequest

def parse_args():
    parser = argparse.ArgumentParser(prog=sys.argv[0], description="Converts a biogrid-derived graph of gene interactions to one of 'clustered' modules.")
    parser.add_argument('-i', '--input', action='store', required=True, help="An input file with list of genes in experiment that have gene expression data, in a tab-delimited format.")
    parser.add_argument('--cluster', action='store', required=True, help="Input cluster CSV-formatted file containing gene names and their module category.")
    parser.add_argument('-c', '--gene_column', action='store', default="Gene ID", help="Name of column header containing gene ids for lookup in BioGRID repo.")
    parser.add_argument('--delimiter', action='store', default="\t", help="Default delimiter separating fields from input file.")
    parser.add_argument('--pickle', action='store', required=True, help="If specified, input python pickled networkx graph previously generated.")
    parser.add_argument('--opickle', action='store', required=True, help="If specified, input python pickled networkx module graph previously generated.")
    parser.add_argument('-o', '--output', action='store', required=True, help="An output file as adjacency list format with lists of genes and their associated parents.")
    parsed = parser.parse_args(sys.argv[1:])
    return(parsed)

if __name__ == '__main__':
    parsed = parse_args()
    try:
        G = nx.read_gpickle(parsed.pickle)
    except Exception as e:
        sys.stderr.write("Failed to read networkx graph from pickle file %s.  Exception: %s\n" % (parsed.pickle, e))

    mG = nx.Graph() #Module graph

    cluster_df = pd.read_csv(parsed.cluster, header=None, names=["id","module"])
    modules = list(cluster_df.groupby('module'))
    num_modules = len(modules)-1
    G.graph['modules'] = [[] for i in range(num_modules)]
    G.graph['edges'] = [[] for i in range(num_modules)]
    G.graph['size'] = [0] * num_modules
    for i in range(0,num_modules):
        module_num = int(modules[i][0])
        ids        = modules[i][1]['id']
        for n in ids:
            if n in G:
                G.node[n]['module'] = module_num
                G.graph['modules'][module_num].append(n)
        G.graph['size'][module_num] = len(G.graph['modules'][module_num])

    for i in range(0,num_modules):
        mG.add_node(i,size=len(G.graph['modules'][i]),genes=G.graph['modules'][module_num])

    for i in range(0,num_modules):
        for j in range(i+1, num_modules):
            cut_size = nx.algorithms.cuts.cut_size(G,G.graph['modules'][i],G.graph['modules'][j])
            if cut_size > 0:
                cut_size = cut_size / ((G.graph['size'][i] + G.graph['size'][j])/2)
                mG.add_edge(i, j, weight=cut_size)

    nx.write_gpickle(mG, parsed.opickle)
