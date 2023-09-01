#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 12:34:43 2021

@author: ruth


"""

import pandas as pd
from pathlib import Path
import argparse



def run(raw_data_folder, data_folder):
    data = pd.read_csv(raw_data_folder / 'taxidlineage.dmp', sep = '\t')
    seq_org = pd.read_csv(data_folder / 'seq_org.tsv', sep = '\t', header=None)
    required_orgs = set(seq_org[1])


    results = {}
    for i, row in data.iterrows():

        leaf = row[0]
        if leaf not in required_orgs:
            continue
        
        if row[2] == row[2]:
            parents = row[2].split()[::-1]
            if parents[0] != str(leaf):
                parents.insert(0, str(leaf))
        
        results[leaf]= parents


    lost = required_orgs - set(results.keys())
    print(len(results), 'covered', len(lost), 'lost', lost)

    with open(data_folder / 'org_lineage.csv', 'w') as f:
        for x in results.values():
            f.write( ','.join(x) +'\n')





def arguments(args=None):
    parser = argparse.ArgumentParser(description='SeqFind script for Selenzy')
    parser.add_argument('data_folder', 
                        help='specify data directory for new files, please end with slash')
    parser.add_argument('raw_data_folder',
                        help='specify data directory for raw databases files, please end with slash')

    arg = parser.parse_args(args=args)
    return arg


if __name__ == '__main__':
    arg = arguments()
    raw_data_folder = Path(arg.raw_data_folder)
    data_folder = Path(arg.data_folder)
    run(raw_data_folder, data_folder)




