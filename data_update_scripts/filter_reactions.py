#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 25 13:41:01 2023

@author: ruth

filter unnessesary reactions and compounds and get stats

"""

import pandas as pd
from pathlib import Path
from collections import Counter
from statistics import median
import argparse

# data_folder     = Path('/home/ruth/code/update_selenzyme/selenzyme_update/data_google_cloud_edited/data_untouched2/')
# raw_data_folder = Path('/home/ruth/code/update_selenzyme/selenzyme_update/data_google_cloud_edited/data_untouched2/')

# chem_prop = pd.read_csv(raw_data_folder / 'chem_prop.tsv', comment='#', sep='\t')
# chem_prop.columns = ['#ID', 'name', 'formula', 'charge', 'mass', 'InChI', 'InChIKey', 'SMILES']
# reac_prop = pd.read_csv(raw_data_folder / 'reac_prop.tsv',  comment='#', sep='\t')
# reac_prop.columns = ['#ID', 'mnx_equation', 'eq', 'is_balanced', 'classifs',  'reference']



def run(raw_data_folder, data_folder):

    chem_prop = pd.read_csv(raw_data_folder / 'chem_prop.tsv', comment='#', sep='\t')
    chem_prop.columns = ['#ID', 'name', 'reference', 'formula', 'charge', 'mass', 'InChI', 'InChIKey', 'SMILES']
    reac_prop = pd.read_csv(raw_data_folder / 'reac_prop.tsv',  comment='#', sep='\t')
    reac_prop.columns = ['#ID', 'mnx_equation', 'reference', 'classifs', 'is_balanced', 'is_transport']



    reac_seqs = data_folder / 'reac_seqs.tsv'
    seq_org = data_folder / 'seq_org.tsv'
    reac_smi = data_folder / 'reac_smi.csv'

    # get the substrates and products 
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
    chem_prop = chem_prop[chem_prop['#ID'].isin(compounds)]
    chem_prop = chem_prop.dropna(subset = ['SMILES'])
    comp_w_smiles = set(chem_prop['#ID'])

    # filter out any reactions where the substrates/ products dont have smiles 
    reaction_compounds2 = {k: [v[0].intersection(comp_w_smiles), v[1].intersection(comp_w_smiles)] for k, v in reaction_compounds.items() if len(v[0].intersection(comp_w_smiles))>0 and  len(v[1].intersection(comp_w_smiles))>0 }
    compounds2 = set.union(*[x for k, v in reaction_compounds2.items() for x in v])



    # make file to filter the input for make_fingerprints.py 
    if reac_seqs.exists():
        
        # filter react_prop by the reactions with enzymes 
        reac_seqs = pd.read_csv(reac_seqs, sep='\t', header=None, names=['mnxr', 'up', 'uniprot', 'ref', 'ec'])
        reactions_enzymes = set(reac_seqs['mnxr']).intersection(reac_prop['#ID'])

        # filter out any reactions where the substrates/ products dont have enzymes
        reaction_compounds3 = {k: v for k, v in reaction_compounds2.items() if k in reactions_enzymes }
        compounds3 = set.union(*[x for k, v in reaction_compounds3.items() for x in v])

        print('\n filter reactions by enzymes')
        print('no. reactions before filtering by smiles', len(reaction_compounds))
        print('no. reactions after filtering by smiles', len(reaction_compounds2))  
        print('no. reactions after filtering by enzymes', len(reaction_compounds3))
        print('')
        print('no. compounds before filtering by reac_prop', len(compounds))
        print('no. compounds after filtering by smiles', len(compounds2))    
        print('no. compounds after filtering by enzymes', len(compounds3))
        print('no. enzymes', len(set(reac_seqs.uniprot)))
        
        with open(raw_data_folder / 'reaction_smiles_enz_filter.tsv', 'w') as f:
            for k, v in reaction_compounds3.items():
                f.write(k + '\t' + ','.join(list(v[0] | v[1]))+ '\n')
                
                
        if reac_smi.exists():
            reac_smi = pd.read_csv(reac_smi, skiprows=0, header=None, names=['mnxr', 'smiles'])
            reac_smi = set(reac_smi.mnxr)
            reaction_compounds3 = {k: v for k, v in reaction_compounds3.items() if k in reac_smi }
            compounds3 = set.union(*[x for k, v in reaction_compounds3.items() for x in v])
            reactions_enzymes = set(reaction_compounds3.keys()).intersection(reac_prop['#ID'])
            
            print('\n filter by processed reactions')
            print('no. reactions after processing', len(reaction_compounds3))
            print('no. compounds after processing', len(compounds3))
            print('no. enzymes', len(set(reac_seqs.uniprot).intersection(set(reac_seqs.uniprot))))
            
        
        if seq_org.exists():
            seq_org = pd.read_csv(seq_org, sep='\t', header=None, names = ['uniprot', 'org', 'org_name'])

            enz_org = set(seq_org.uniprot)
            reac_seqs2 = reac_seqs[reac_seqs.uniprot.isin(enz_org)]
            reactions_enzymes2 = set(reac_seqs2['mnxr']).intersection(reactions_enzymes)

            reaction_compounds4 = {k: v for k, v in reaction_compounds3.items() if k in reactions_enzymes2 }
            compounds4 = set.union(*[x for k, v in reaction_compounds4.items() for x in v])

            print('\n filter by exnzymes with organisms')
            print('no. reactions after filtering by enzymes with org', len(reaction_compounds4))
            print('no. compounds after filtering by enzymes with org ', len(compounds4))       
            
            print('no. enzymes with org', len(set(seq_org.uniprot)))   
            print('no. orgs', len(set(seq_org.org)))  

            
            # get extra data 
            reac_prop2 = reac_prop[reac_prop['#ID'].isin(reaction_compounds4)]
            reac_prop2.loc[:,'source'] = [x.split(':')[0] for x in reac_prop2.reference]
            print(Counter(reac_prop2.source))
            
            reac_seqs2 = reac_seqs.loc[reac_seqs.uniprot.isin(enz_org)]
            reac_seqs2 = reac_seqs2.loc[reac_seqs2.mnxr.isin(reaction_compounds4)]
            enz_count = Counter(reac_seqs2.mnxr)
            hist = Counter(enz_count.values())
            
            print('median enzymes per reaction', median(list(enz_count.values())))

            
            # range_list = [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], 
            #               [6, 11], [11, 21], [21, 31], [31, 51], [51, 101], [101, max(hist)]]

            range_list = [[1, 2], [2, 6], [6, 11], [11, 21], [21, 51], [51, max(hist)+1]]
            
            print()
            for x in range_list:
                print(x[0], '-', x[1]-1, '\t', sum([hist[x] for x in list(range(x[0], x[1]))]), round(sum([hist[y] for y in list(range(x[0], x[1]))])/len(reaction_compounds4), 3) )
     
          
                
    # make file to filter the input for make_reac_seqs.py     
    else:

        print('no. reactions before filtering by smiles', len(reaction_compounds))
        print('no. reactions after filtering by smiles', len(reaction_compounds2))
        print('')
        print('no. compounds before filtering by smiles', len(compounds))
        print('no. compounds after filtering by smiles', len(compounds2))
        print('\n made reaction_smiles_filter.txt')
        with open(raw_data_folder / 'reaction_smiles_filter.txt', 'w') as f:
            for k in reaction_compounds2.keys():
                f.write(k + '\n' )
            



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