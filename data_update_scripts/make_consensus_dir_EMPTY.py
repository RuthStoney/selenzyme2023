#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 22 11:48:52 2023

@author: ruth
"""


import pandas as pd
from pathlib import Path
from collections import Counter

data_folder = Path('/home/ruth/code/update_selenzyme/run_folder_min_dec/data_2023/')
raw_data_folder = Path('/home/ruth/code/update_selenzyme/run_folder_min_dec/raw_data_update/')


reac_prop = pd.read_csv(raw_data_folder / 'reac_prop.tsv', skiprows=351, sep='\t')

cofactors = set(['WATER', 'MNXM13', 'MNXM735438', 'MNXM3', 'MNXM40333', 'MNXM64096', 'MNXM10'])

mnxr = []
subs = []
prods = []
subs_co = []
prods_co = []
ecs = []


for i, row in reac_prop.iterrows():
    if row['#ID'] == "EMPTY" or type(row.mnx_equation) != str: continue
    mnxr.append(row['#ID'])

    sub, prod = row.mnx_equation.split(' = ')
    sub = set([x.split('@')[0] for x in sub.split(' ') if 'MNXM' in x])
    prod = set([x.split('@')[0] for x in prod.split(' ') if 'MNXM' in x])
    
    sub_co = sub.intersection(cofactors)
    prod_co = prod.intersection(cofactors)
    
    sub = sub - cofactors
    prod = prod - cofactors
    
    subs.append('|'.join(sub))
    prods.append('|'.join(prod))
    
    if len(sub_co)>0:
        subs_co.append('|'.join(sub_co))
    else:
        subs_co.append('')
    if len(prod_co)>0:
        prods_co.append('|'.join(sub_co))
    else:
        prods_co.append('')
    
    ec= row.classifs
    if type(ec) == str:
        ec = ec.split(';')[0]
    else:
        ec = 'NOEC'
    
    ecs.append(ec)

df = pd.DataFrame({'ec':ecs, 'id': mnxr, 'dir':[0]*len(mnxr), 's': subs, 'p':prods, 'sc': subs_co, 'pc': prods_co})
df.to_csv(data_folder / 'rxn_consensus_20160612.txt', header = False, index=False, sep='\t')






