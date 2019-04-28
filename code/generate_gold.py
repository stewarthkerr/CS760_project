#!/usr/bin/env python3
#Author: Heloise Letissier

import argparse
import csv
import json
import networkx as nx
import os
import sys
import urllib
import urllib.request as urequest

def parse_args():
    parser = argparse.ArgumentParser(prog=sys.argv[0], description="Utility to generate the gold standard graph of gene interactions by pulling interactions from BioGRID repo.")
    parser.add_argument('-i', '--input', action='store', required=True, help="An input file with list of genes in experiment that have gene expression data, in a tab-delimited format.")
    parser.add_argument('-c', '--gene_column', action='store', default="Gene ID", help="Name of column header containing gene ids for lookup in BioGRID repo.")
    parser.add_argument('--delimiter', action='store', default=',', help="Default delimiter separating fields from input file.")
    parser.add_argument('-a', '--access_key', action='store', default="ce5f7787eb9eed0fe6f1461c4d568162", help="BioGRID access key for using their REST API service.")
    parser.add_argument('-t', '--tax_id', action='store', type=int, default=3702, help="BioGRID taxonomy id of species to search genes in.")
    parser.add_argument('-o', '--output', action='store', required=True, help="An output file as adjacency list format with lists of genes and their associated parents.")
    parser.add_argument('--interval', action='store', type=int, default=250, help="At what interval to save graph and print status upon buildup of graph.")
    parser.add_argument('--pickle', action='store', help="If specified, is the file path for the pickle file representing the NetworkX graph that is being built.")
    parsed = parser.parse_args(sys.argv[1:])
    return(parsed)

edge_metas=['BIOGRID_INTERACTION_ID','OFFICIAL_SYMBOL_A','OFFICIAL_SYMBOL_B','EXPERIMENTAL_SYSTEM','ONTOLOGY_TERMS']
#All potential evidence codes for literature-based evidence of interaction
#["affinity capture-luminescence","affinity capture-ms","affinity capture-rna","affinity capture-western","biochemical activity","co-crystal structure","co-fractionation","co-localization","co-purification","dosage growth defect","dosage lethality","dosage rescue","far western","fret","negative genetic","pca","phenotypic enhancement","phenotypic suppression","positive genetic","protein-peptide","protein-rna","proximity label-ms","reconstituted complex","synthetic growth defect","synthetic haploinsufficiency","synthetic lethality","synthetic rescue","two-hybrid"]
#Remove 'pca' as valid experimental evidence codes
evidence_codes=["affinity capture-luminescence","affinity capture-ms","affinity capture-rna","affinity capture-western","biochemical activity","co-crystal structure","co-fractionation","co-localization","co-purification","dosage growth defect","dosage lethality","dosage rescue","far western","fret","negative genetic","phenotypic enhancement","phenotypic suppression","positive genetic","protein-peptide","protein-rna","proximity label-ms","reconstituted complex","synthetic growth defect","synthetic haploinsufficiency","synthetic lethality","synthetic rescue","two-hybrid"]
evidenceList = '|'.join(evidence_codes)
biogrid_url="https://webservice.thebiogrid.org/interactions?searchIds=true&geneList=%s&includeInteractors=true&taxId=%d&interSpeciesExcluded=true&evidenceList=%s&includeEvidence=true&format=json&accesskey=%s"

if __name__ == '__main__':
    parsed = parse_args()
    if( parsed.pickle and os.path.exists(parsed.pickle) ):
        G = nx.read_gpickle(parsed.pickle)
    else:
        G = nx.Graph()
    with open(parsed.input,'r') as input_fh:
        input_csv   = csv.DictReader(input_fh, delimiter=parsed.delimiter)
        geneids     = [e[parsed.gene_column] for e in input_csv]
        #geneList    = '|'.join(geneids)
        nodes       = [e[0] for e in G.edges()]
        geneids     = list(set(geneids) - set(nodes)) #Remove anything that has already been searched
        count = 0
        for gid in geneids:
            while 1:
                try:
                    url_fh = urequest.urlopen(biogrid_url % (gid,parsed.tax_id,evidenceList,parsed.access_key))
                    break
                except (TimeoutError,urllib.error.URLError) as e:
                    print(e)
                    pass
            data   = json.loads(url_fh.read())
            url_fh.close()
            if isinstance(data,dict): #Required because empty strings are cast as lists by json.loads()
                for interaction in data.values():
                    geneA = interaction["SYSTEMATIC_NAME_A"]
                    geneB = interaction["SYSTEMATIC_NAME_B"]
                    G.add_edge(geneA,geneB)
                    for m in edge_metas:
                        if m in interaction:
                            G[geneA][geneB][m] = interaction[m]
            if (count % parsed.interval) == 0:
                print("Extracting gene number %d" % (count))
                if( parsed.pickle):
                    nx.write_gpickle(G, parsed.pickle)
            count = count + 1
        df = nx.to_pandas_adjacency(G, nodelist=geneids)
        df.to_csv(parsed.output, header=df.columns)
