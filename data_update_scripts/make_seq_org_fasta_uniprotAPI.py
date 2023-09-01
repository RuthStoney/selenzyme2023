#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 17:23:54 2023

@author: ruth
making seq_org
"""
from pathlib import Path
import pandas as pd
import re
import requests
import time 
import argparse


def create_taxonomy_dict(fasta_file):
    uniprot_dict = {}
    taxonomy_dict = {}
    taxonomy_code_dict = {}
    tax_names = {}

    with open(fasta_file, 'r') as file:
        for line in file:
            if '>' in line:
                sp, unip, info = line.split('|')

                tax_code = info.split()[0].split('_')[1]
                
                # name = re.findall('\w+_\w+\s(.*?)\sOS=', info)[0]
                
                tax = re.findall('OX=\w+', info)[0].replace('OX=', '')
                
                tax_n = re.findall('OS=(.*?) OX=', info)[0]
                
                uniprot_dict[unip] = tax
                taxonomy_dict[tax_n] = tax
                taxonomy_code_dict[tax_code] = tax
                tax_names[tax] = tax_n

    return uniprot_dict, taxonomy_dict, taxonomy_code_dict, tax_names

def write_seq_org1(enzymes, uniprot_dict, tax_names):
    l = []
    for enz in enzymes:
        tax_id = uniprot_dict[enz]
        l.append([enz, tax_id, tax_names[tax_id] ])
    return l

def write_seq_org2(brenda, enzymes, taxonomy_dict):
    l = []
    for enz in enzymes:
        org_name = brenda['org'][brenda['enz'] == enz].values[0]
        l.append([enz, taxonomy_dict[org_name], org_name] )
    return l
    

def write_seq_rest(enzymes):
    print('sending ', + len(enzymes), 'to uniprot api. Estimated time in hours', (len(enzymes)*2 /60) /60)
    
    l=[]
    retrieved_enz = set()
    
    for enz in enzymes:
        time.sleep(2)
        response = requests.get('https://www.uniprot.org/uniprot/'+ enz +'.json')
        
        # Check if the request was successful (status code 200)
        if response.ok:
            # Extract the taxonomic code from the response JSON
            data = response.json()
            if 'organism' not in data:
                print(data)
                continue
            tax_code = data['organism']['taxonId']
            if 'commonName' in data['organism']:
                name = data['organism']['commonName']
            elif 'scientificName' in data['organism']:
                name = data['organism']['scientificName']
            else:
                name = ''
                
            l.append([enz, tax_code, name])
            retrieved_enz.add(enz)
         
            # print("Taxonomic code for "+ enz +":", tax_code)
        else:
            print("Uniprot query failed:", response.status_code, enz)
        
    return l, retrieved_enz

def names_dmp(taxlin, brenda_lost):
    taxlin = pd.read_csv(taxlin, sep = '|', header=None)
    taxlin.columns = ['taxid', 'name', 'lin', 'x', 'y']
    
    n = [x.lower().strip() for x in taxlin['name']]
    d = dict(zip(n, taxlin['taxid'].astype(str).str.strip()))  #taxlin['name'].str.strip()
    
    l = []
    covered = set()
    for i, row in brenda_lost.iterrows():
        if row.org.lower() in d:
            l.append([row.enz, d[row.org.lower()], row.org])
            covered.add(row.enz)
    return l, covered, d




def run(raw_data_folder, data_folder, legacy_folder):
    reac_seqs = pd.read_csv(data_folder / 'reac_seqs.tsv', sep='\t', header=None, names=['mnxr', 'up', 'uniprot', 'ref', 'ec'])
    required_enz = set(reac_seqs['uniprot'])


    ### Get most of the data from the uniprot fasta file
    fasta_file = str(raw_data_folder) + '/uniprot_sprot.fasta'
    uniprot_dict, taxonomy_dict, taxonomy_code_dict, tax_names = create_taxonomy_dict(fasta_file)

    # make reac_seq for the enzymes in the uniprot file 
    covered = required_enz.intersection(set(uniprot_dict.keys()))
    dataset1 = write_seq_org1(covered, uniprot_dict, tax_names)

    #### make reac_seq for the enzymes not in the uniprot file
    lost1 = required_enz - covered

    ### recover missing enzymes from the brenda file (no lost enzymes were found in expasy)
    # expasy = pd.read_csv(output_folder/ 'expasy_data.tsv', sep='\t')
    # expasy_recovery = lost1.intersection(set( expasy['enz']))

    brenda = pd.read_csv(raw_data_folder / 'brenda_data.tsv', sep='\t').drop_duplicates()
    brenda_recovery = brenda[ (brenda['enz'].isin(lost1)) & (brenda['org'].isin(set(taxonomy_dict.keys()))) ]
    covered = set(brenda_recovery.enz)
    dataset2 = write_seq_org2(brenda_recovery, covered, taxonomy_dict)

    brenda_lost = brenda[ (brenda['enz'].isin(lost1)) & ~(brenda['enz'].isin(covered)) & ~(brenda['org'].isin(set(taxonomy_dict.keys()))) ]
    datasetx, covered2, d = names_dmp(raw_data_folder / 'names.dmp', brenda_lost)

    lost2 = lost1 - (covered | covered2)
    lost_orgs = brenda_lost.org[brenda_lost['enz'].isin(lost2)]

    print('Brenda - couldnt find taxonomy name in fasta file or names.dmp to recover taxID:', len(lost2), 'found taxonmy name ', len(covered | covered2) )



    ### recover missing enzymes from the previous file seq_org file
    seq_org_old = pd.read_csv(legacy_folder / 'seq_org.tsv', header=None, sep='\t', names = ['unip', 'tax', 'tax_name'])

    recovered_previous_file = lost2.intersection(set(seq_org_old['unip']))
    dataset3 = seq_org_old[ seq_org_old['unip'].isin(recovered_previous_file) ]
    lost3 = lost2 - recovered_previous_file


    ### get the remaining enzymes from the uniprot REST API 
    dataset4, covered = write_seq_rest(list(lost3))
    lost4 = lost3 - covered


    seq_org_new = pd.concat( [dataset3, pd.DataFrame(data = dataset1 + dataset2 + dataset4, columns = dataset3.columns )]) 
    print('\nEnzymes input:', len(required_enz))
    print('Enzymes covered:', len(set(seq_org_new.unip)))
    print('Enzymes lost:', len(lost4))

    seq_org_new.to_csv(data_folder / 'seq_org.tsv', sep = '\t', index=False, header=False)




def arguments(args=None):
    parser = argparse.ArgumentParser(description='SeqFind script for Selenzy')
    parser.add_argument('data_folder', 
                        help='specify data directory for new files, please end with slash')
    parser.add_argument('raw_data_folder',
                        help='specify data directory for raw databases files, please end with slash')
    parser.add_argument('legacy_folder',
                        help='specify data directory for raw databases files, please end with slash')

    arg = parser.parse_args(args=args)
    return arg


if __name__ == '__main__':
    arg = arguments()
    raw_data_folder = Path(arg.raw_data_folder)
    data_folder = Path(arg.data_folder)
    legacy_folder = Path(arg.legacy_folder)
    run(raw_data_folder, data_folder, legacy_folder)



