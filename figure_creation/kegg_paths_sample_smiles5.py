#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  2 12:45:39 2023

@author: ruth
"""

import os
import requests
import json
import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
from matplotlib import pyplot as plt, colors
import matplotlib
from collections import Counter
import time
#import shutil
from statistics import median
# import requests

        
def ec_comp(refs, tests):
    #refs = refs.split(';')
    # tests = tests.split(';')
    # refs, tests = exact_match_ec, r['EC Number']
    if type(tests) == str:
        tests = tests.split(';')
    if type(refs) == str:
        refs = refs.split(';')
    
    tophit =0
    for ref in refs:
        rf = ref.split('.')
        #if len(rf)<4: continue
    
        for test in tests:
            t= test.split('.')
            #if len(t)<4: continue
            
            for i in list(range(0, min(len(rf), len(t)))):
                if rf[i] == t[i]:
                    i1 = i+1
                    if i1> tophit: tophit=i1
                else:
                    break
        if tophit ==4: break
    return tophit


def sort_extras(extras, easymax, n=2):
    extras = extras[extras['Rxn. ID - additional'].notna()]
    if n == 2:
        sort1, sort2 = 'Rxn Sim RF.', 'Rxn Sim.' 
    else:
        sort1, sort2 = 'Rxn Sim.', 'Rxn Sim RF.'  
    
    hits = []
    for i, x in extras.iterrows():
        for y in x['Rxn. ID - additional'].split('|'):
            y = y.split('_')
            if y[n] == 'nan': continue
            if float(y[1]) == 1.0: continue
            if float(y[n])>=easymax:
                #data[['Rxn. ID', 'EC Number','Rxn Sim.', 'Rxn Sim RF.']]
                hits.append([y[0], x['EC Number'], float(y[1]), float(y[2]), x['Rxn. ID']])

    if len(hits) ==0:
        return([[], -1])
    
    hits = pd.DataFrame(hits, columns=['Rxn. ID', 'EC Number','Rxn Sim.', 'Rxn Sim RF.', 'in row'] )
    hits = hits.drop_duplicates()
    hits = hits.sort_values([sort1, sort2], ascending = [False, False])
    return([hits, hits[sort1][0]])
           

def make_plot(Results, outfile_scatter, fail_file):
    print('reactions', Results.shape[0], 'set', len(set(Results.react_input)))
    
    outfile_scatter = outfile_scatter.replace('.png', '_'+str(Results.shape[0])+'.png')
    fail_file = fail_file.replace('.csv', '_'+str(Results.shape[0])+'.csv')
    data_file = fail_file.replace('fails', 'data_summary')
        
    scatter_data = {}
    for i, row in Results.iterrows():
        k = (row.score_rf, row.score_2018)
        if k not in scatter_data:
            scatter_data[k] = 0
        scatter_data[k]+=1
    
    sl = list(scatter_data.keys())
    l1 = [x[0] for x in sl]
    l2 = [x[1] for x in sl]
    l3 = [scatter_data[x] for x in sl]
    l3_scaled = [(scatter_data[x]*50)+5 for x in sl]
    
    s2018 = {0:0, 1:0, 2:0, 3:0, 4:0}
    for x in sl:
        s2018[x[1]]+= scatter_data[x]/sum(scatter_data.values())
    s2018 ={k: round(v, 3) for k, v in s2018.items()}
    
    if -1 in l1:
        srf = {-1:0, 0:0, 1:0, 2:0, 3:0, 4:0}
    else:
        srf = {0:0, 1:0, 2:0, 3:0, 4:0}
        
    for x in sl:
        srf[x[0]]+= scatter_data[x]/sum(scatter_data.values())
    srf ={k: round(v, 3) for k, v in srf.items()}
    
    # cmap = plt.cm.get_cmap('winter_r')  # winter_r   gist_earth_r
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ['#243c4c', '#5f9b8c', '#a0c382', "#fac846", "#ff7d2d", "#cc6424" ])
    #cmap = mpl.colors.LinearSegmentedColormap.from_list('winter', cmaplist, cmap.N)
    norm = colors.Normalize(vmax=max(l3), vmin=min(l3))
    norm = colors.BoundaryNorm(np.arange(0, max(l3)+1, 1), cmap.N)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [12, 1]})
    fig.tight_layout(pad=3.0)
    
    ax1.plot(list(range(0, 5)), color='grey', ls="--")

    ax1.yaxis.grid(color='lightgray')
    ax1.xaxis.grid(color='lightgray')
    ax1.scatter(l2, l1, l3_scaled, color=cmap(norm(l3)), edgecolors='dimgrey')
    # ax1.yaxis = ['Inavaliable', '1', '22', '3', '4']
    ax1.set_title('EC of Closest Non-Identical Reactions')
    if -1 in l1:
        ax1.set_yticklabels(['', 'UA', '0', '1', '2', '3', '4'])
    
    # for i, x in enumerate(l3):
    #     ax1.annotate(x, (l2[i], l1[i]), color='red')
    
    ax1.set_xlabel("Correct EC digits sim_2018")
    ax1.set_ylabel("Correct EC digits sim_RF")
    
    ax3 = ax1.twiny()
    ax3.set_xlim(ax1.get_xlim())
    ax3.set_xticklabels(['']+[str(round(s2018[x]*100, 1))+'%' for x in list(range(0, 5))] )
    
    ax4 = ax1.twinx()
    ax4.set_ylim(ax1.get_ylim())
    if -1 in l1:  
        ax4.set_yticklabels([''] + [str(round(srf[x]*100, 1))+'%' for x in list(range(-1, 5))])
    else:
        ax4.set_yticklabels([''] + [str(round(srf[x]*100, 1))+'%' for x in list(range(0, 5))])
    
    ax3.tick_params(axis='x', colors='dimgrey', top=False) 
    ax4.tick_params(axis='y', colors='dimgrey', right=False)
    
    
    ax3.set_xlabel("% Correct EC digits sim_2018" , color='dimgrey')
    ax4.set_ylabel("% Correct EC digits sim_RF" , color='dimgrey')
    
    cbar=fig.colorbar(matplotlib.cm.ScalarMappable(norm = norm,
                   cmap = cmap), ticks= list(range(min(l3), max(l3)+1)),
                   cax = ax2, orientation ='vertical')
                  # label ='Correct EC digits')
    cbar.ax.set_yticklabels(list(range(min(l3), max(l3)+1))) 
    plt.savefig(outfile_scatter, dpi=300, bbox_inches='tight')
    plt.show()   
  
    fail_dot = Results[(Results.score_rf<Results.score_2018) & (Results.score_rf> -1)]
    print('worse than 2018', fail_dot.shape[0], fail_dot.shape[0]/Results.shape[0])
    
    fail_dot = Results[Results.score_rf == -1]
    print('AAM failed', fail_dot.shape[0], fail_dot.shape[0]/Results.shape[0])   

    fail_dot = Results[(Results.score_rf<Results.score_2018) | (Results.score_rf == 0)]
    fail_dot.to_csv(fail_file, index=False)
    
    data_summary = Results[['path', 'source_reaction', 'score_rf', 'score_2018', 'react_input', 'ec_input', 'react_rf', 'react_2018']]
    data_summary.to_csv(data_file, index = False)
    
    # fail_dot['path1'] = [ Data.pathway[Data.reaction == x].values[0] for x in fail_dot.source_reaction ]   
    
    # fail_dot['smile'] = [ Data.smiles[Data.reaction == x].values[0] for x in fail_dot.source_reaction ]
     
    return (fail_dot)

 


      
chem_prop = "/home/ruth/code/update_selenzyme/run_folder_min_dec/raw_data_update/chem_prop.tsv"
outfile_data = "/home/ruth/code/update_selenzyme/analysis_ouput/kegg_paths2.csv"
upload_path = '/home/ruth/code/update_selenzyme/selenzyme_2023/selenzyme2/selenzyPro/uploads/'


#chem_prop = "/home/ruth/code/update_selenzyme/run_folder/selenzyme/selenzyPro/data/chem_prop.tsv"
chem_prop= pd.read_csv(chem_prop, sep='\t', skiprows=351)
chem_prop['name'] = chem_prop['name'].str.lower()
chem_prop['source'] = [x.split(':')[0] for x in chem_prop.reference]
chem_prop['source_id'] = [x.split(':')[1] for x in chem_prop.reference]
reac_prop_kegg = set(chem_prop['source_id'][chem_prop['source'] == 'keggC'])



# paths_O = sorted(['M00375', 'M00029', 'M00123', 'M00927', 'M00003', 'M00009', 'M00013',  'M00082', 'M00049', 'M00019', 
#                 'M00016', 'M00127', 'M00137', 'M00672', 'M00101', 'M00048', 'M00741', 'M00056', 
#          'M00097', 'M00915', 'M00962', 'M00793', 'M00027', 'M00166', 'M00804', 'M00872', 'M00941','M00736'])



# paths = set(['M00001', 'M00002', 'M00003', 'M00004', 'M00005', 'M00006', 'M00007', 'M00008', 'M00009',  
#          'M000010', 'M00011', 'M00012', 'M00013', 'M00014', 'M00015', 'M00016', 'M00017', 'M00018', 'M00019',
#          'M00375', 'M00029', 'M00123', 'M00927', 'M00003', 'M00009', 'M00013',  'M00082', 'M00049', 'M00019', 
#                 'M00016', 'M00127', 'M00137', 'M00672', 'M00101', 'M00048', 'M00741', 'M00056', 
#          'M00097', 'M00915', 'M00962', 'M00793', 'M00027', 'M00166', 'M00804', 'M00872', 'M00941','M00736'])


paths = ['M00003', 'M00009', 'M00012', 'M00022', 'M00101' ]
paths = ['M00001', 'M00004', 'M00009', 'M00018', 'M00022', 'M00101', 'M00110', 'M00125', 'M00137', 'M00138',  'M00372', 'M00944']  #  'M00002', 'M00003', 
# paths = ['M00868', 'M00121', 'M00125', 'M00116', 'M00372', 'M00889', 'M00012', 'M00044',  ]
# paths = ['M00001']

# M00001 R01070 
# M00004 R01827 
# M00004 R01830 
# M00001_M00004 R02740 
# M00001 R04779

# M00048 pretty fail 
# M00973 Nitrogen metabolism

pathways = {}
pathways2 = {}
reactions = {}

reactions['R01070'] = [[{'C05378'}, {'C00111', 'C00118'}]]
reactions['R04779'] = [[{'C05345'}, {'C05378'}]]
reactions['R02740'] = [[{'C00668'}, {'C05345'}]]
reactions['R01827'] = [[{'C05382', 'C00118'}, {'C00279', 'C05345'}]]
reactions['R01830'] = [[{'C05345', 'C00118'}, {'C00279', 'C00231'}]]

pathways['R01070'] = ['M00001']
pathways['R04779'] = ['M00001']
pathways['R02740'] = ['M00001_M00004']
pathways['R01827'] = ['M00004']
pathways['R01830'] = ['M00004']

ecs = {}
ecs['R01070'] = ['4.1.2.13']
ecs['R04779'] = ['2.7.1.11']
ecs['R02740'] = ['5.3.1.9']
ecs['R01827'] = ['2.2.1.2']
ecs['R01830'] = ['2.2.1.1']

compounds = ['C00002', 'C00008',  'C05345',  'C05378'  ]
names = ['ATP', 'ADP', 'beta-D-Fructose 6-phosphate', 'beta-D-Fructose 1,6-bisphosphate']

for path in paths:
    # get reactions for pathway
    url ='https://rest.kegg.jp/get/'+path
    r = requests.post( url).text
    lines = r.splitlines()
    pathways2[path] = []
    
    if path == 'M00001':
        pathways2[path] = pathways2[path] + ['R04779', 'R02740', 'R01070']
    if path == 'M00004':
        pathways2[path] = pathways2[path] + ['R02740', 'R01827', 'R01830'] 

    for x in lines: 
        if re.search("^\s+R\d{5}", x) or 'REACTION' in x:
            s = set(re.findall("C\d{5}", x.split('->')[0]))
            p = set(re.findall("C\d{5}", x.split('->')[1]))
            r = re.findall("R\d{5}", x)
            if r[0] not in pathways:
                pathways[r[0]] = []
                reactions[r[0]] = []
                
            pathways[r[0]].append(path)
            if [s, p] not in reactions[r[0]] and [p, s] not in reactions[r[0]] :
                reactions[r[0]].append([s, p])
            pathways2[path].append(r[0])
            
            url ='https://rest.kegg.jp/get/'+r[0]
            r1 = requests.post( url).text
            lines2 = r1.splitlines()
            for y in lines2:
                if y[0:6] == 'ENZYME':
                    ecs[r[0] ] = y.split()[1:]
    

        if re.search("^\s+C\d{5}", x) or 'COMPOUND' in x:
            c = re.search("C\d{5}", x).group()
            compounds.append(c)
            names.append(x.split(c)[1].strip())

print('\nstarting reaction count', len(reactions))
print('starting compound count', len(set(compounds)))
  


man = {'14-Demethyllanosterol':set(['MNXM728840']),
       '4alpha-Methylzymosterol-4-carboxylate':set(['MNXM732410']),
       '4alpha-Methylzymosterol': set(['MNXM1104216']),
       'Lathosterol': set(['MNXM810']),
       'Dihydrokaempferol':set(['MNXM1125']),
       'cis-3,4-Leucopelargonidin':set(['MNXM2221']),
       '2-Dehydro-3-deoxy-D-arabino-heptonate 7-phosphate': set(['MNXM1219']),
       'Shikimate 3-phosphate':set(['MNXM1265']),
       '5,10-Methenyltetrahydromethanopterin':set(['MNXM740011']),
       '2-(Methylthio)ethanesulfonate':set(['MNXM671']),
       'Coenzyme M 7-mercaptoheptanoylthreonine-phosphate heterodisulfide':set(['MNXM1145']),
       'D-Glucono-1,5-lactone 6-phosphate':set(['MNXM1104491']),
       '2-Phospho-D-glycerate': set(['MNXM275']),
       '3-Phospho-D-glyceroyl phosphate':set(['MNXM1108073']),
       'Glycerone phosphate':set(['MNXM77']),
       '(S)-2,3-Epoxysqualene':set(['MNXM727928']),
       '3-Phospho-D-glycerate':set(['MNXM727604']),
       'D-Ribose 5-phosphate': set(['MNXM1104056']),
       'acetyl-CoA': set(['MNXM732449']),
       '(2R)-2-Hydroxy-3-(phosphonooxy)-propanal': set(['MNXM739946']),
       '3-Dehydroquinate': set(['MNXM736874']),
       '7-Dehydrocholesterol': set(['MNXM730454']),
       'Cholesterol': set(['MNXM726122']),
       'Lanosterol': set(['MNXM736942']),
       'Malonyl-CoA': set(['MNXM1106093']),
       'Sedoheptulose 7-phosphate': set(['MNXM733314']),
       'Succinyl-CoA': set(['MNXM1104774']),
       '3-Keto-4-methylzymosterol': set(['MNXM730540']),
       'Acetyl-CoA':set(['MNXM1104266']),
       'Isocitrate':set(['MNXM89661']),
       '4alpha-Carboxy-5alpha-cholesta-8,24-dien-3beta-ol':set(['MNXM739474']),
       '5-Phosphoribosylamine': set(['MNXM1104549']),
       "5'-Phosphoribosylglycinamide": set(['MNXM1104522']),
       "5'-Phosphoribosyl-N-formylglycinamide":set(['MNXM1094915']),
       "2-(Formamido)-N1-(5'-phosphoribosyl)acetamidine": set(['MNXM1094916']),
       'Aminoimidazole ribotide': set(['MNXM1105983']),
       '1-(5-Phospho-D-ribosyl)-5-amino-4-imidazolecarboxylate': set(['MNXM737371']),
       "1-(5'-Phosphoribosyl)-5-amino-4-(N-succinocarboxamide)-imidazole": set(['MNXM1104647']),
       "1-(5'-Phosphoribosyl)-5-amino-4-imidazolecarboxamide": set(['MNXM1103365']),
       "1-(5'-Phosphoribosyl)-5-formamido-4-imidazolecarboxamide": set(['MNXM1104518']),
       'Dehydroepiandrosterone':set(['MNXM731293']),
       'Androstenedione':set(['MNXM730402']),
       'trans-Cinnamate':set(['MNXM736095']),
       '4-Coumarate': set(['MNXM505']),
        'p-Coumaroyl-CoA': set(['MNXM1102789']),
        'Naringenin chalcone': set(['MNXM727324']),
        'Salutaridinol': set(['MNXM732800']),
        '7-O-Acetylsalutaridinol': set(['MNXM733592']),
        'Oripavine': set(['MNXM737665']),  
        "5-Amino-6-(5'-phosphoribosylamino)uracil": set(['MNXM1532']),
        "5-Amino-6-(5'-phosphoribitylamino)uracil": set(['MNXM1178']),
        '5-Amino-6-(1-D-ribitylamino)uracil': set(['MNXM737103']),
        '3,4-Dihydroxy-2-butanone 4-phosphate': set(['MNXM1102294']),
        '6,7-Dimethyl-8-(1-D-ribityl)lumazine': set(['MNXM1107755']),
        'FMN': set(['MNXM1105927']),
        'FAD': set(['MNXM1105936']),
        '2-Succinyl-5-enolpyruvyl-6-hydroxy-3-cyclohexene-1-carboxylate': set(['MNXM730931']),
        '(1R,6R)-6-Hydroxy-2-succinylcyclohexa-2,4-diene-1-carboxylate': set(['MNXM727528']),
        'Menaquinol': set(['MNXM736183']),
        'beta-Cryptoxanthin': set(['MNXM729605']),
        'Zeaxanthin': set(['MNXM731637']),
        'Antheraxanthin': set(['MNXM586']),
        'Xanthoxin': set(['MNXM726697']),
        'Abscisic aldehyde': set(['MNXM1998']),
        'Abscisate': set(['MNXM725931']),
        "3'-Keto-3'-deoxy-ATP": set(['MNXM4978']),
        "N6,N6,O-Tridemethylpuromycin-5'-phosphate": set(['MNXM6338'])                                    
       }

       # '3-Keto-4-methylzymosterol': set(['MNXM730540', 'MNXM730541']),


# get smiles
smiles={}
mnxr = {}
lost_c = set()
lc = set()
for i in list(range(0, len(names))):
    x = names[i]
    smi=set()

    if (x not in man) and x.lower() in chem_prop.name.values:
        smi_ds = chem_prop[['#ID', 'SMILES']][chem_prop.name == x.lower()].dropna(subset=['SMILES'])
        smi = set(smi_ds.SMILES.values)
        if len(smi) == 0:
            print('no smiles', compounds[i], x)
            lc.add(compounds[i])
            #crash
        smi2 = [x for x in smi if '*' not in x]
        if len(smi2) == 0:
            print('no star free smiles', compounds[i], x)
            lc.add(compounds[i])
            continue
        mnxr[compounds[i]] = set(smi_ds['#ID'].values)
        smiles[compounds[i]]  = smi2
        if len(set(smi_ds['#ID'].values))>1: print('multiples', x, compounds[i])

    elif x in man or len(smi) == 0:
        #smi = chem_prop.SMILES[chem_prop['#ID'].isin(man[x])].dropna()
        try:
            smi_ds = chem_prop[['#ID', 'SMILES']][chem_prop['#ID'].isin(man[x])].dropna(subset=['SMILES'])
        except:
            print('go get', compounds[i], x)
            continue
        smi = set(smi_ds.SMILES.values)
        if len(smi) == 0:
            print('no smiles', compounds[i], x)
            lc.add(compounds[i])
            if x in man:
                print('man', x, man[x])
                crash
                
        smi2 = [x for x in smi if '*' not in x]
        if len(smi2) == 0:
            print('no star free smiles', compounds[i], x)
            lc.add(compounds[i])
            continue
        mnxr[compounds[i]] = set(smi_ds['#ID'].values)
        smiles[compounds[i]]  = smi2
        
        if len(set(smi_ds['#ID'].values))>1: 
            if x in man:
                print('man', x, man[x])
            print('multiples', x, compounds[i])
    else:
        print('wtf', compounds[i], x)


ambiguous_reactions = []
for k, v in reactions.items():
    if len(v)>1:
        print('ambiguous reaction', k)
        ambiguous_reactions.append(k)
    else:
        reactions[k] = v[0]

reactions = {k:v for k, v in reactions.items() if k not in ambiguous_reactions}  

reac_lost = []
for reac, v in reactions.items():
    if v[0].intersection(lc) or v[1].intersection(lc):
        reac_lost.append(reac)
        print('lost', reac, pathways[reac])

reactions = {k:v for k, v in reactions.items() if k not in reac_lost}
print('\nreaction count - missing compounds removed', len(reactions), '\t lost', len(reac_lost), reac_lost)
print('compound count - missing compounds removed', len(smiles), '\t lost', len(lc))


reactions = {k:v for k, v in reactions.items() if k in ecs}
print('\nreaction count - missing ecs removed', len(reactions))


# reaction smiles
smiles_r = []
sm_react = []
reac_lost2 = []
missing_smi = set()
for r, comps in reactions.items():

    if r not in ecs or r not in pathways:
        print('lost1', r)
        reac_lost2.append(r)
        continue
    for s in comps[0]:
        if s not in smiles: 
            print('lost2', r, s)
            reac_lost2.append(r)
            missing_smi.add(s)
            continue
        for p in comps[1]:
            if p not in smiles: 
                print('lost3', r, p)
                reac_lost2.append(r)
                missing_smi.add(p)
                continue
            for x in smiles[s]:
                for y in smiles[p]:
                    smiles_r.append(x+'>>'+y)
                    sm_react.append(r)

                    

print('\nreaction count - reactions with smiles',  'total smiles:', len(set(smiles_r)), 'reactions covered:', len(set(sm_react)))
multi_smiles = {k: v for k, v in dict(Counter(sm_react)).items() if v>1}

# s, p = reactions['R07215']
# for x in s: print(mnxr[x])
# for x in p: print(mnxr[x])

                                                  
# make reactions     
Data = pd.DataFrame({'smiles': smiles_r, 'reaction': sm_react,  
                     'ec': ['_'.join(ecs[x]) for x in sm_react], 'pathway':['_'.join(pathways[x]) for x in sm_react] })        

#Data.to_csv(outfile_data, index = False) 




Data = pd.read_csv(outfile_data)
Data = Data.drop_duplicates(subset = ['smiles'])
Data['ec'] = [[y for y in x.split("_")] for x in Data.ec ]


print('\nreactions', Data.shape[0], len(set(Data.reaction)), '\n')
issues = ["R01826", "R01518",   "R02291", "R03698", "R05640", "R07509", "R07215", "R01826", "R01826", "R06953", "R08517", "R00351", "R00549", "R01518", "R01837", "R02035", "R04457"]
#Data = Data[Data.reaction.isin(issues)]


###### SELENZYME
# url='http://0.0.0.0:32784/REST'
url='http://0.0.0.0:5001/REST'
hits= 5000

# Data_copy = Data.copy()
# Data = Data_copy[:50]



# change inflection point
#for fs in list(range(0, 9)):  0, 1, 3, 
for fs in [5]:
    print('')
    print('FS', fs)
    
    out_folder =      "/home/ruth/code/update_selenzyme/analysis_ouput/final/"
    outfile_results = "/home/ruth/code/update_selenzyme/analysis_ouput/final/"+ str(fs) + "_kegg_paths_results_bb.csv"
    outfile_scatter = "/home/ruth/code/update_selenzyme/analysis_ouput/final/"+ str(fs) + "_kegg_scatter_update_bb.png"
    fail_file =       "/home/ruth/code/update_selenzyme/analysis_ouput/final/"+ str(fs) + "_fails_bb.csv"
    
    sim_2018 = []
    sim_rf = []
    fails = []
    fails2 = []
    # reactions without an exact match 
    fails3={}
    fail_total = []
    fail_rf = {}
    results = []
    source_reaction = []
    rf_match_min =[]
    rf_match_max = []
    
    #Data_safe = Data.copy()
    Data = Data_safe
    
    
    for i, row in Data.iterrows():  
        print('')
        print(i, row.reaction, row.pathway)
        # shutil.rmtree(upload_path, ignore_errors=True)
        #if row.smiles == '[NH3+][C@@H](CC([O-])=O)C([O-])=O>>N[C@@H](CCCNC(=N)N[C@@H](CC(O)=O)C(O)=O)C(O)=O': crash
    
        smarts = row.smiles
       # smarts = 'O=C(CO)[C@@H](O)[C@H](O)[C@H](O)[C@H](O)COP(=O)([O-])[O-]>>O=C[C@H](O)[C@H](O)COP(=O)([O-])[O-]'

        # if smarts == 'C=C(OP(=O)([O-])[O-])C(=O)[O-]>>O=C([O-])C(=O)C[C@@H](O)[C@H](O)[C@H](O)COP(=O)([O-])[O-]': crash
        
        try:        
            r = requests.post( os.path.join(url, 'Query') , json={'smarts': smarts, 'targets': 5000, 'frag_size': fs} )  
            time.sleep(5)
            res = json.loads( r.content.decode('utf-8') )
        
            val1 = json.loads( res['data'] )
            data = pd.DataFrame(val1)
            data['Rxn Sim RF.'].fillna(value =np.nan, inplace=True)
            d2 = data.drop_duplicates(subset=['Rxn. ID'])
        except:
            print('skipping', i, row.reaction)
            fail_total.append(row.reaction)
            continue
        
        split_enz = []
        for i2, row2 in data.iterrows():
            if row2['Rxn. ID - additional']:
                ex = row2['Rxn. ID - additional'].split('|')
                for y in ex:
                    if y == '-': 
                        continue
                    mxr, s18, s23 = y.split('_')
                    row2['Rxn. ID'] = mxr
                    row2['Rxn Sim.'] = float(s18)
                    row2['Rxn Sim RF.'] = float(s23)
                    split_enz.append(row2)
                    
        split_enz = pd.DataFrame(split_enz, columns = data.columns)
        split_enz = split_enz.drop_duplicates()
        #split_enz.to_csv(out_folder + row.reaction + '.csv')
        print(data.shape[0], split_enz.shape[0], '=', pd.concat([data, split_enz]).shape[0])
        data = pd.concat([data, split_enz])
    
        # find the input values
        # exact_matches = data[['Rxn. ID', 'EC Number', 'Rxn Sim.', 'Rxn Sim RF.']][data['Rxn Sim.']==1].dropna(subset = ['EC Number']) 
        exact_matches = data[['Rxn. ID', 'EC Number', 'Rxn Sim.', 'Rxn Sim RF.']][data['Rxn Sim.']>0.99].dropna(subset = ['EC Number']) 
        # exact_matches = exact_matches.dropna(subset = ['Rxn Sim RF.']) 
        if exact_matches.shape[0] == 0:
            print('COULDNT MATCH', row.reaction)
            fails3[row.reaction] = data[['Rxn Sim.', 'EC Number', 'Query']][data['Rxn Sim.']==max(data['Rxn Sim.'])].drop_duplicates()
            continue
        
        exact_matches = exact_matches.drop_duplicates()
        exact_matches['hits']=0
        exact_matches['ecs'] = 0
        exact_match_ec = set()
        for i1, x in exact_matches.iterrows():
            ecs_retrieved = set(x['EC Number'].split(';')) 
            exact_matches['hits'][i1] = len(set(row.ec).intersection(ecs_retrieved))
            exact_match_ec = exact_match_ec | ecs_retrieved
        
        # if max(exact_matches.hits) == 0:
        if len(exact_match_ec) == 0:
            print('COULDNT MATCH EC', row.reaction)
            fails3[row.reaction] = data[['Rxn Sim.', 'EC Number']][data['Rxn Sim.']==max(data['Rxn Sim.'])].drop_duplicates()
            continue
        exact_match = '|'.join(exact_matches['Rxn. ID'])
        
        data = data[-data['Rxn. ID'].isin(exact_matches['Rxn. ID'])].dropna(subset = ['EC Number']) 
        data['Rxn Sim RF.'] = data['Rxn Sim RF.'].round(4)
        data['Rxn Sim.'] = data['Rxn Sim.'].round(4)
    #try:
        if exact_matches.dropna(subset = ['Rxn Sim RF.']).shape[0]>0:        
           
            d =  data.dropna(subset = ['Rxn Sim RF.']) 
            if d.shape[0] ==0: 
                print('RF missing for', row.reaction)
                fail_rf[row.reaction] = max(data['Rxn Sim.'][(data['Rxn Sim.'] <1)])
                continue
    
            easymax_rf = max(d['Rxn Sim RF.'])
    
            tophits = d[['Rxn. ID', 'EC Number','Rxn Sim.', 'Rxn Sim RF.']] [d['Rxn Sim RF.'] == easymax_rf]
            tophits_rf = tophits.drop_duplicates()
            print(tophits_rf.shape[0])
    
            if tophits_rf.shape[0] == 0: 
                crash2
                continue
        
            # store the results
            s=0
            for i1, r in tophits_rf.iterrows(): 
                ec_score = ec_comp(exact_match_ec, r['EC Number'])
                if ec_score >= s: 
                    s =ec_score 
                    score = [s, exact_match, r]
        else:
            score = [-1, exact_match, pd.Series({'Rxn. ID': None, 'EC Number': None, 'Rxn Sim.': None, 'Rxn Sim RF.': None})]
            print('no RF for ', row.reaction)
            
        print('RF', score[0], exact_match_ec, 'recomendation',  score[2][1] )
    
    
        
        ########################################## sim_2018      
        # get the easymax
        easymax = max(data['Rxn Sim.'])
        
        # # get the top scoring reactions
        tophits = data[['Rxn. ID', 'EC Number','Rxn Sim.', 'Rxn Sim RF.']] [data['Rxn Sim.'] == easymax]
        tophits_2018 = tophits.drop_duplicates()
        print(tophits_2018.shape[0])
        
        if tophits_2018.shape[0] == 0: 
            crash4
            continue
        
        ########################################## store the results
        s=0
        for x, r in tophits_2018.iterrows():   
            ec_score = ec_comp(exact_match_ec, r['EC Number'] )
            if ec_score >= s: 
                s =ec_score 
                score_2018 = [s, exact_match, r]             
        print('2018', score_2018[0], exact_match_ec, 'recomendation',  score_2018[2][1] )
               
    
        if score[0]<score_2018[0]: 
            print(row.reaction)
            fails.append([exact_match, row])
            #crash
        if max([score[0], score_2018[0]])==0: 
            print(row.smiles)
            print('new', score[0], 'old', score_2018[0])
    
        if score[0] >2 and score_2018[0]<2: 
            print('win example', exact_match, row.smiles)
        
        
        results.append([exact_match, score[0], score_2018[0], 
                        exact_match_ec,  score[2]['EC Number'], score_2018[2]['EC Number'],
                        score[2]['Rxn. ID'], score[2]['Rxn Sim RF.'], score[2]['Rxn Sim.'], 
                        score_2018[2]['Rxn. ID'], score_2018[2]['Rxn Sim RF.'], score_2018[2]['Rxn Sim.'], 
                        row.reaction, min(exact_matches['Rxn Sim RF.']), 
                        max(exact_matches['Rxn Sim RF.']), row.smiles])
    
    
    
    print('total couldnt match ec', len(fails3), '/', Data.shape[0])     
    
    Results = pd.DataFrame(results, columns=['react_input', 'score_rf', 'score_2018',
                                             'ec_input', 'ec_rf', 'ec_2018', 
                                             'react_rf', 'react_rf_rf_score', 'react_rf_2018_score', 
                                             'react_2018', 'react_2018_rf_score', 'react_2018_2018_score', 
                                             'source_reaction', 'rf_match_min', 'rf_match_max', 'smile'])  
    
    
    Results['path'] = [Data.pathway[Data.reaction == x].values[0] for x in Results.source_reaction ]
    Results = Results.drop_duplicates(subset=['smile'])
    
    outfile_results = outfile_results.replace('bb', str(round(sum(Results.score_rf)/Results.shape[0]*100)))
    outfile_scatter = outfile_scatter.replace('bb', str(round(sum(Results.score_rf)/Results.shape[0]*100)))
    
    Results.to_csv(outfile_results, index=False)
    Results = pd.read_csv(outfile_results)

    
    
    print('\nreactions in dataset', Results.shape[0], len(set(Results.path)))
    fail_dot_all = make_plot(Results, outfile_scatter, fail_file)
    
    print('rf median', median(Results.score_rf), '2018 median', median(Results.score_2018))
    
    lost_rfs = fail_dot_all[fail_dot_all.score_rf == -1]
    worse_tests = fail_dot_all[fail_dot_all.score_rf < fail_dot_all.score_2018]
    
    print('RF fails ', lost_rfs.shape[0]/Results.shape[0])
    print('RF fails ', worse_tests.shape[0]/Results.shape[0])
    
    #Counter(Results_filter.path)
    
    Results2 = Results[Results.path.isin(['M00001','M00003', 'M00009', 'M00012', 'M00022', 'M00101'])]
    Counter(Results2.path)
    
    print('\nreactions in dataset', Results2.shape[0], len(set(Results2.path)))
    fail_dot_all = make_plot(Results, outfile_scatter, fail_file)
    
    print('rf median', median(Results.score_rf), '2018 median', median(Results.score_2018))
    
    lost_rfs = fail_dot_all[fail_dot_all.score_rf == -1]
    valid_rfs =  fail_dot_all[fail_dot_all.score_rf != -1]
    worse_tests = valid_rfs[valid_rfs.score_rf < valid_rfs.score_2018]
    
    print('RF lost ', lost_rfs.shape[0]/Results.shape[0], lost_rfs.shape[0])
    print('RF fails ', worse_tests.shape[0]/Results.shape[0], worse_tests.shape[0])
    

    
    


x = Results[['path', 'score_rf', 'score_2018', 'source_reaction']].sort_values(['path', 'score_rf','score_2018'])

bp = []
cc = Counter(Results.path)
fc= Counter(fail_dot_all[fail_dot_all.score_rf<2].path)
for k, v in cc.items():
    if k in fc:
        bp.append( [round(fc[k]/v, 3), v, fc[k], k])
    else:
        bp.append( [0, v, 0, k])
bp = sorted(bp)





#Results = pd.read_csv('/home/ruth/Downloads/data_summary_bb_103.csv' )
fail_dot_all = make_plot(Results, outfile_scatter, fail_file)



# Results2 = Results.drop_duplicates(subset = ['source_reaction'])
# print('number of reactions', Results2.shape[0], len(set(Results2.source_reaction)))
# fail_dot_drop_dup = make_plot(Results2, outfile_scatter, fail_file)

# Results_filter = Results[~Results.path.isin(['M00137', 'M00012', 'M00138'])]
# Results_filter = Results[~Results.path.isin(['M00044', 'M00012'])]
# print('number of reactions', Results_filter.shape[0], len(set(Results_filter.path)), len(set(Results_filter.source_reaction)))
# fail_dot_50 = make_plot(Results_filter, outfile_scatter, fail_file)



# Results['path'] = [Data.pathway[Data.reaction == x].values[0] for x in Results.source_reaction ]
# Results2 = Results[Results.path.isin(['M00003', 'M00004', 'M00007', 'M00016', 'M00018', 'M00027', 'M00048','M00127', 'M00029', 'M00736', 'M00941'])]
# test_sample = ['M00003', 'M00009', 'M00012', 'M00022', 'M00101' ]
# Results = Results[Results.path.isin(test_sample)]


# r_safe = Results


# Counter(Data.pathway)
# print('total couldnt match ec', len(fails3), '/', Data.shape[0], 'sets ', len(set(fails3)), len(set(Data.reaction)), '\n')  

# print('number of reactions', Results.shape[0], len(set(Results.source_reaction)))
# Counter(Results.path)


# Counter(Results2.path)
# Results = Results2

# # Results = Results[~Results.path.isin(['M00138', 'M00044'])] #, 'M00007', 'M00016', 'M00018', 'M00027', 'M00048','M00127', 'M00029', 'M00736', 'M00927', 'M00941'])]
# #Res = Results[0:50]

# fail_dot = Results[Results.score_rf == 0]
# fs =Data.smiles[Data.reaction.isin(fail_dot.source_reaction)].values


# fail_dot_drop_dup = make_plot(Results2, outfile_scatter)



    # easymax = max(data['Rxn Sim.'][data['Rxn Sim.'] !=1])
    # extradata, extramax = sort_extras(data[data['Rxn Sim.']>=easymax], easymax, 1)
    # #if extramax == 1 or easymax ==1: crash
    
    # if easymax > extramax:
    #     tophits = data[['Rxn. ID', 'EC Number','Rxn Sim.', 'Rxn Sim RF.']] [data['Rxn Sim.'] == easymax]
    #     tophits = tophits.drop_duplicates()
    # else:
    #     tophits = extradata
    # tophits = tophits[tophits['EC Number'].astype(bool)]
    # if tophits.shape[0] == 0: continue
    
    # score_2018=[0]
    # for x, r in tophits.iterrows():   
    #     s=ec_comp(row.ec, r['EC Number'])
    #     if s>=score_2018[0]: score_2018 = [s, exact_match, r]
    # sim_2018.append(score_2018) 
    
    
    
    
    
    # altmax = None   
    
    # for i2, row2 in data.iterrows():
    #     if row2['Rxn Sim RF.'] > easymax:
    #         break
    #     if row2['Rxn. ID - additional'] != None and len(row2['Rxn. ID - additional'])>0:  
    #         extras = row2['Rxn. ID - additional'].split('|')
    #         if float(row2['Rxn. ID - additional'].split('_')[1]) > easymax:
    #             altmax = float(row2['Rxn. ID - additional'].split('_')[1]) 
    #             tophits = 
            
    
    # if altmax == None:
    #     tophits = data[['Rxn. ID', 'EC Number','Rxn Sim.', 'Rxn Sim RF.']] [data['Rxn Sim.'] == easymax]
    #     tophits = tophits.drop_duplicates()
    #     for x, r in tophits.iterrows():   
    #         s=ec_comp(row.ec, r['EC Number'])



# url = 'https://rest.kegg.jp/conv/chebi/'+'+'.join(compounds)
# r = requests.post( url).text.split('\n')
# chebi = { x.split('\t')[0].split()[0].split(':')[1] : x.split('\t')[1] for x in r if len(x)>0}
# for kegg, che in chebi.items():
#     if che in chem_prop.reference.values:
#         print('win')
#         smiles[kegg] = chem_prop.SMILES[chem_prop.reference == che].values[0]
#     else:
#         print(kegg, che)
