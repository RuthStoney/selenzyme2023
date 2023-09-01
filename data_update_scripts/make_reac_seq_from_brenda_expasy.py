#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 11 19:45:10 2023

@author: ruth
Update reac_seqs from Brenda and Expasy, then compare it to the old files
use Brenda and Expasy to link ec numbers to enzymes
use the EC numbers to link enzymes to MextaNetX reactions

"""


import re
import pandas as pd
from pathlib import Path
import argparse


class DataSet():
    
    def __init__(self):
        self.data = None
        self.ecs = set()
        self.reactions = set()
        self.enzymes = set()

    # get the sets of values 
    def get_ecs(self):
        data_types = set([type(value).__name__ for value in self.data['ec']])
        if data_types == {'str'}:
            self.ecs = set(self.data['ec'])
        else:
            self.ecs = set([y for x in self.data['ec'] for y in x])

    def get_reactions(self):
        self.reactions =  set(self.data['mnxr'])    

    def get_enzymes(self):
        self.enzymes =  set(self.data['enz'])
    
    #join datasets then get the unqiue values
    def join_data(self, df_list):
        data = pd.concat(df_list)
        self.data = data.drop_duplicates()
        self.get_ecs()
        self.get_enzymes()
    
    def print_data(self, data_set, report):
        print('\n' + data_set)
        x = []
        for k, v in report.items():
            x = x + [k+':'] + [str(len(v))]
        print(' '.join(x))
    
    def write_file(self, out_file, header = True):
        self.data.to_csv(out_file, sep='\t', index=False)
        


class MataNetxDataSet(DataSet):
    
    def __init__(self, file_path):
        super().__init__()
        if 'reac_prop' in str(file_path):
            self.read_reac_prop_tsv(file_path)
        else:
            self.read_reac_seqs_tsv(file_path)            
        self.get_ecs()         
        self.get_reactions()

    
    def read_reac_seqs_tsv(self, file_path):
        self.data = pd.read_csv(file_path, sep='\t', header=None, names=['mnxr', 'up', 'uniprot', 'ref', 'ec'])
        self.data['ec'] = [[y for y in re.split(r'\||;',  str(x))] for x in self.data['ec'] ]

    def read_reac_prop_tsv(self, file_path):
        # the metanetx reactions
        self.data = pd.read_csv(file_path, skiprows=351, sep='\t')
        self.data = self.data[self.data.is_transport != 'T']
        self.data.columns = ['mnxr', 'mnx_equation', 'reference', 'classifs', 'is_balanced', 'is_transport']
        self.data.reference =[i.split('#')[0] for i in self.data['reference']]
        
        self.data['db'] = self.data['reference'].str.split('#').str[0] 
        self.data['db_ref'] =   self.data['reference'].str.split('#').str[1] 
        self.data['ec'] = self.data['classifs'].fillna('').apply(lambda x: [y for y in x.split(';') if re.search(r'\d+\.\d+\.\d+\.\d+$', y)])
               


class Expasy(DataSet):
    
    def __init__(self, file_path):
        super().__init__()
        self.read_file(file_path)
        self.get_ecs()
        self.get_enzymes()
        
    def read_file(self, file_path):
        data_raw = open(file_path,'r').read().split('\n//\nID')[1:]
        
        data = []
        for x in data_raw:    
            lines = x.split('\n')
            ec = lines[0].strip()
            for line in lines:
                if 'DR   ' in line:
                    line = line.replace('DR   ', '').split(';')
                    for y in line:
                        if y:
                            data.append([ec, y.split(',')[0].strip(), y.split(',')[1].split('_')[1].strip()])       
        self.data = pd.DataFrame(data, columns=['ec', 'enz', 'org']).drop_duplicates()
        


class Brenda(DataSet):
    
    def __init__(self, file_path):
        super().__init__()
        self.read_file(file_path)
        self.get_ecs()
        self.get_enzymes()
    
    def read_file(self, file_path):
        data_raw = open(file_path,'r').read().split('\n///\nID\t')        
        data = {}
        moved = {}
            
        for ec_no in data_raw:
            
            lines = ec_no.split(r'\n(?!\t)')
            lines = re.split(r'\n(?!\t)', ec_no)
            ec = lines[0]
            if '(' in ec:
                start_ec = re.search(r'^\d+.\d+.\d+.\d+', ec)[0]
                other_ecs = set(re.findall(r'\d+.\d+.\d+.\d+', ec)) - set([start_ec])
                moved[start_ec] = other_ecs
                continue
            data[ec] = []
            
            for line in lines:
                if 'PR\t' in line:
                    unip=''
                    # remove brackets because they correpond to different organisms
                    line_mod = re.sub(r'\([^()]*\)', '', line).replace('\t', '').replace('\n', '')
                    unip_search = re.findall(r'(\s(\d|[A-Z]){6}\s|\s(\d|[A-Z]){10}\s)', line_mod)
                    unip = [x[0].strip() for x in unip_search]
                            
                    if unip:  
                        org_search =  re.search(r'#\d+#(.*?)'+ unip[0], line_mod).group(1).strip()
                        if org_search:
                            # org = ' '.join(org_search[0].split()[-2:])
                            data[ec].append(['|'.join(unip), org_search])
                        else:
                            data[ec].append(['|'.join(unip)])
                            
        # remove missing ecs
        missing=[]
        for k, v in data.items():
            if not v:
                missing.append(k)
        
        for x in missing: data.pop(x)
        
        # fill in moved ecs
        for k, v in moved.items():
            for x in v:
                if x in data:
                    data[k] = data[x]
        
        
        for k, v in moved.items():
            for x in v:
                if x in data:
                    data[k] = data[x]
        
        
        d2 = [[k, z, y[1]] if len(y)==2 else [k, z, ''] for k, x in data.items() for y in x for z in y[0].split('|')]
        self.data = pd.DataFrame(d2, columns=['ec', 'enz', 'org']).drop_duplicates()


class ReacSeqs(DataSet):
    
    def __init__(self, file_path):
        self.read_file(file_path)
        self.get_ecs()
        self.get_enzymes()
        self.get_reactions()
    
    def read_file(self, file_path):
        data = pd.read_csv(file_path, sep='\t', 
                        header=None, names=['mnxr', 'up', 'enz', 'ref', 'ec'])
        self.data = data.fillna('')



def run(raw_data_folder, data_folder, legacy_folder):
    ### Read in the data 
    # Get the EC - enzyme data 
    brenda = Brenda(raw_data_folder / 'brenda_2023_1.txt')
    brenda.print_data('brenda', {'ecs' : brenda.ecs, 'enzymes': brenda.enzymes})
    brenda.write_file(raw_data_folder/ 'brenda_data.tsv')


    expasy = Expasy(raw_data_folder / 'expasy_dat.txt')
    expasy.print_data('expasy', {'ecs' : expasy.ecs, 'enzymes': expasy.enzymes})
    expasy.write_file(raw_data_folder/ 'expasy_data.tsv')


    # Combine the Brenda and Expasy datasets
    combi = DataSet()
    combi.join_data([brenda.data, expasy.data])
    combi.print_data('combi', {'ecs' : combi.ecs, 'enzymes': combi.enzymes})

    # Read in the MetaNetX reactions 
    reac_prop = MataNetxDataSet(raw_data_folder / 'reac_prop.tsv')
    reac_prop.print_data('reac_prop', {'ecs' : reac_prop.ecs, 'reactions': reac_prop.reactions})

    # Get the previous reac_seq for comparison 
    reac_seqs_old = ReacSeqs(legacy_folder / 'reac_seqs.tsv')
    reac_seqs_old.print_data('reac_seqs_old', {'ecs' : reac_seqs_old.ecs, 'enzymes': reac_seqs_old.enzymes})


    ######  rebuild reac_seqs
    reac_seqs = []
    lost_ecs =set()
    valid_ecs = set(combi.data['ec'])

    for i, row in reac_prop.data.iterrows():
        h = 0
        ecs = row.ec

        # Check if any EC in ecs is a valid ECS
        if any(ec in valid_ecs for ec in ecs):
            h = 1
            d = combi.data.loc[combi.data['ec'].isin(ecs), 'enz']
            for x in d:
                reac_seqs.append([row['mnxr'], 'uniprot', x, row.reference, ';'.join(ecs)])
        
        if h == 0:
            lost_ecs.update(ecs)
       
                
    reac_seqs_new = pd.DataFrame(data=reac_seqs).drop_duplicates()


    ###### Compare old and new reac_seqs files
    reac_prop.print_data('total input reactions (reac_prop)', {'ecs' : reac_prop.ecs, 'reactions': reac_prop.reactions})

    print('\nnew reac_seq')
    print('reactions:', len(set(reac_seqs_new[0])), 'enzymes:',  len(set(reac_seqs_new[2])), 'ecs:', len(set(reac_seqs_new[4])))
    reac_seqs_old.print_data('old reac_seqs', 
                         {'reactions': reac_seqs_old.reactions, 'enzymes' : reac_seqs_old.enzymes , 'ecs' : reac_seqs_old.ecs })

    combi.reactions = set(reac_seqs_new[0])
    missing_reactions = reac_prop.data[~reac_prop.data['mnxr'].isin(combi.reactions)]

    # write file 
    reac_seqs_new.to_csv(data_folder / 'reac_seqs.tsv', sep='\t', header=False, index=False)




def arguments(args=None):
    parser = argparse.ArgumentParser(description='SeqFind script for Selenzy')
    parser.add_argument('data_folder', 
                        help='specify data directory for new files, please end with slash')
    parser.add_argument('raw_data_folder',
                        help='specify data directory for raw databases files, please end with slash')
    parser.add_argument('legacy_folder',
                        help='specify data directory for previous daata files, please end with slash')
    arg = parser.parse_args(args=args)
    return arg


if __name__ == '__main__':
    arg = arguments()
    raw_data_folder = Path(arg.raw_data_folder)
    data_folder = Path(arg.data_folder)
    legacy_folder = Path(arg.legacy_folder)


    run(raw_data_folder, data_folder, legacy_folder)