#!/usr/bin/env python3
#Author: Ongo Gablogian

import argparse
import networkx as nx
import numpy as np
import pandas as pd
import sys

def parse_args():
    parser = argparse.ArgumentParser(prog=sys.argv[0], description="Generate golden network adjacency matrix, using specified weight threshold to generate binary edges.")
    parser.add_argument('--pickle', action='store', required=True, help="If specified, input python pickled networkx module graph previously generated.")
    parser.add_argument('-o', '--output', action='store', required=True, help="An flattened tab-delimited adjacency list with source module, destination module, 1 for edge, 0 for no-edge.")
    parser.add_argument('-t', '--threshold', action='store', type=float, default=0.0, help="Threshold")
    parsed = parser.parse_args(sys.argv[1:])
    return(parsed)

if __name__ == '__main__':
    parsed = parse_args()
    try:
        G = nx.read_gpickle(parsed.pickle)
    except Exception as e:
        sys.stderr.write("Failed to read networkx graph from pickle file %s.  Exception: %s\n" % (parsed.pickle, e))

    df = nx.to_pandas_adjacency(G)
    df = df > parsed.threshold
    df = df.astype('uint8')
    #This is a n x n array. Flatten it.
    out_fh = open(parsed.output, 'w')    
    for index,val in np.ndenumerate(df):
        rvalue = index[0]
        cvalue = index[1]
        out_fh.write("ME%d\tME%d\t%d\n" % (rvalue,cvalue,val))
    out_fh.close()

