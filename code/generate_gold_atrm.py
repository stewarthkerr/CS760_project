#!/usr/bin/env python3
#Author: Heloise Letissier

import argparse
import csv
import networkx as nx
import os
import sys

def parse_args():
    parser = argparse.ArgumentParser(prog=sys.argv[0], description="Utility to generate the gold standard graph of gene interactions by pulling interactions from ATRM file.")
    parser.add_argument('-i', '--input', action='store', required=True, help="An input file with list of transcription factor gene and target gene pairs.")
    parser.add_argument('-s', '--source_column', action='store', default="TF ID", help="Name of column header containing the source transcription factor gene name.")
    parser.add_argument('-d', '--dest_column', action='store', default="Target ID", help="Name of column header containing the destination target gene name.")
    parser.add_argument('--delimiter', action='store', default="\t", help="Default delimiter separating fields from input file.")
    parser.add_argument('-o', '--output', action='store', required=True, help="An output file as adjacency list format with lists of genes and their associated parents.")
    parser.add_argument('--pickle', action='store', help="If specified, is the file path for the pickle file representing the NetworkX graph that is being built.")
    parsed = parser.parse_args(sys.argv[1:])
    return(parsed)

if __name__ == '__main__':
    parsed = parse_args()
    G = nx.Graph()
    with open(parsed.input,'r') as input_fh:
        input_csv   = csv.DictReader(input_fh, delimiter=parsed.delimiter)
        for r in input_csv:
            tf = r[parsed.source_column]
            target = r[parsed.dest_column]
            G.add_edge(tf, target, effect=r['Activate/Repress'])
    nx.write_gpickle(G, parsed.pickle)
    df = nx.to_pandas_adjacency(G)
    df.to_csv(parsed.output, header=df.columns)
