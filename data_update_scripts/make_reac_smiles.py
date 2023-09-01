#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 21 19:12:38 2023

@author: ruth
"""


from pathlib import Path
import pandas as pd


def arguments(args=None):
    parser = argparse.ArgumentParser(description='SeqFind script for Selenzy')
    parser.add_argument('data_folder', 
                        help='specify data directory for new files, please end with slash')
    parser.add_argument('raw_data_folder',
                        help='specify data directory for raw databases files, please end with slash')

    arg = parser.parse_args(args=args)
    return arg


def run(raw_data_folder, data_folder):
    # raw_data_folder = raw_data_folder = Path('/home/ruth/code/update_selenzyme/run_folder_min_dec/raw_data_update')
    # data_folder = Path('/home/ruth/code/update_selenzyme/selenzyme_aug2023/data_2023_t/')

    reac_prop = pd.read_csv(raw_data_folder / 'reac_prop.tsv', skiprows=351, sep='\t')
    chem_prop = pd.read_csv(raw_data_folder / 'chem_prop.tsv', skiprows=351, sep='\t')


    # get the substrates and products for the reactions
    reaction_compounds = {}
    compounds = set()
    for i, row in reac_prop.iterrows():
        if row['#ID'] == "EMPTY" or type(row.mnx_equation) != str: continue
        sub, prod = row.mnx_equation.split(' = ')
        sub = set([x.split('@')[0] for x in sub.split(' ') if 'MNXM' in x])
        prod = set([x.split('@')[0] for x in prod.split(' ') if 'MNXM' in x])
        reaction_compounds[row['#ID']] = [sub, prod]
        compounds = compounds  | set(sub) | set(prod)

    # get the compound smiles
    comp_smiles = {}
    for n, row in chem_prop.iterrows(): 
        if row['#ID'] not in compounds: 
            continue

        if row.SMILES != row.SMILES: continue
        smiles = min(row.SMILES.split('.'), key=len)
        comp_smiles[row['#ID']] = smiles.strip()
    
    # make the reaction smiles 
    reac_smi = []
    for reac, comps in reaction_compounds.items():
        sub = '.'.join([comp_smiles[x] for x in comps[0] if x in comp_smiles])
        prod = '.'.join([comp_smiles[x] for x in comps[1] if x in comp_smiles])
        
        if sub and prod:
            reac_smi.append([reac, '>>'.join([sub, prod])])
    
    # write to file 
    with open(data_folder/'reac_smi.csv', 'w+') as f:
        f.write('RID,SMILES\n')
        for x in reac_smi:
            f.write(','.join(x) + '\n')
 
    
    
if __name__ == '__main__':
    arg = arguments()
    raw_data_folder = Path(arg.raw_data_folder)
    data_folder = Path(arg.data_folder)
    run(raw_data_folder, data_folder)