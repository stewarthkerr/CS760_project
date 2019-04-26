#!/usr/bin/env python3
#Author: Heloise Letissier

import argparse
import csv
import json
import sys
import urllib.request as urequest

def parse_args():
    parser = argparse.ArgumentParser(prog=sys.argv[0], description="Utility to generate the gold standard graph of gene interactions by pulling interactions from BioGRID repo.")
    parser.add_argument('-i', '--input', action='store', required=True, help="An input file with list of genes in experiment that have gene expression data, in a tab-delimited format.")
    parser.add_argument('-c', '--gene_column', action='store', default="Gene ID", help="Name of column header containing gene ids for lookup in BioGRID repo.")
    parser.add_argument('-a', '--access_key', action='store', default="ce5f7787eb9eed0fe6f1461c4d568162", help="BioGRID access key for using their REST API service.")
    parser.add_argument('-t', '--tax_id', action='store', type=int, default=3702, help="BioGRID taxonomy id of species to search genes in.")
    parser.add_argument('-o', '--output', action='store', required=True, help="An output file in JSON format with lists of genes and their associated parents.")
    parsed = parser.parse_args(sys.argv[1:])
    return(parsed)

biogrid_url="https://webservice.thebiogrid.org/interactions?searchIds=true&geneList=%s&includeInteractors=true&taxId=%d&format=json&accesskey=%s"

if __name__ == '__main__':
    parsed = parse_args()

    with open(parsed.input,'r') as input_fh:
        input_csv   = csv.DictReader(input_fh, delimiter="\t")
        geneids     = [e[parsed.gene_column] for e in input_csv]
        for gid in geneids:
            url_fh = urequest.urlopen(biogrid_url % (gid,parsed.tax_id,parsed.access_key))
            data   = json.loads(url_fh.read())
