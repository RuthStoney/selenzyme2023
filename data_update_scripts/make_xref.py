#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 18 11:44:46 2023

@author: ruth
"""

import pandas as pd
from pathlib import Path
import argparse


def run(raw_data_folder, data_folder):
    data_head = []
    data = []
    with open(raw_data_folder / 'reac_xref.tsv') as f:
        for line in f:
            if '#' in line:
                data_head.append(line.strip())
            elif 'MNXR' in line:
                source = line.strip().split('\t')[0]
                if 'R:' in source:
                    source = source.replace('R:', ':')
                if '.reaction:' in source:
                    source = source.replace('.reaction:', ':')
                
                mnxr = line.strip().split('\t')[1]
                data.append([source, mnxr])
                if 'MNXR' not in mnxr:
                    print('check', mnxr)
    
    # data = list(set(data))
    res = list(set(map(lambda i: tuple(i), data)))
    data_sorted = sorted(res, key = lambda x: x[1])
    
    with open(data_folder / 'reac_xref2.tsv', 'w') as f:
        for x in data_head:
            f.write(x + '\n')
        for x in data_sorted:
            f.write('\t'. join(x) + '\n')
                
                
            
        
    
    # data = pd.read_csv(raw_data_folder / 'taxidlineage.dmp', sep = '\t')

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