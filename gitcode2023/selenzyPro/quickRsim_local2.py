'''
quickRsim (c) University of Manchester 2017

quickRsim is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  Pablo Carbonell, SYNBIOCHEM
@description: Compute a fast reaction similarity with the database
@examples: 
1. Compute similarity to one reaction in the database (requires -chem parameter)
python quickRsim.py data/reac_prop.tsv data/metanetx.fs -rid MNXR3215 -chem data/chem_prop.tsv
2. Compute similarity to a given reaction file for a threshold above 0.9
python quickRsim.py data/reac_prop.tsv data/metanetx.fs -rxn rhea15870.rxn -th 0.9
'''

from __future__ import print_function
import argparse
import subprocess
import re
#from os import path
import os
from rdkit.Chem import rdChemReactions
import math
import numpy as np
from rdkit import DataStructs, Chem
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprint #, GetAtomPairFingerprint, GetTopologicalTorsionFingerprint
from rdkit.Chem.rdmolops import PatternFingerprint, RDKFingerprint
from rdkit.Chem import AllChem
#from mcs_functions2 import get_MCS2
from mcs_functions2 import rxnMapper_fun, reactFragDists, process_compF, getReactingFrags, get_MCS2, process_comp, neutralize_atoms, mol_with_atom_index,  getAtoms, getAtomFragments, get_common_atoms_frags, dict_to_sparse
from func_timeout import func_timeout, FunctionTimedOut
#from scoringAlgOct import generate_RFscore2, simpleWeight, dict_to_sparse
import pandas as pd
from rxnmapper import RXNMapper
from rdkit.Chem import Draw
import mcs_functions2
from statistics import mean 


#from rdkit.Chem import AllChem
#from rdkit import Chem
#from rdkit.Chem import rdMolDescriptors, Draw
#from rdkit.Chem import rdChemReactions
#from rdkit.Chem import Draw

#from rdkit.Chem import rdChemReactions
#from rdkit.Chem.Draw import IPythonConsole
#IPythonConsole.ipython_useSVG = True
#IPythonConsole.molSize = 800,400



def fingerprint():
    """ I keep the three fingerprints that give better results in the tests """
    # fpd =  {'Pattern': ('ptfp', None, True, PatternFingerprint, 2), 'RDK': ('rdkfp', None, True, RDKFingerprint, 1),
    #         'Morgan' : ('mgfp5', 5, False, GetMorganFingerprint, 3)}
    fpd =  {'Pattern': ('ptfp', None, True, PatternFingerprint, 2), 'RDK': ('rdkfp', None, True, RDKFingerprint, 1),
            'Morgan' : ('FP_Morg', 5, False, GetMorganFingerprint, 3)}
    fpd =  {'Morgan' : ('FP_Morg', 5, False, GetMorganFingerprint, 3)}
    return fpd


# def fingerprintRA(fp):
#     """ I keep the three fingerprints that give better results in the tests """
#     fpd =  {'Pattern': ('FP_MorgRF', 'FP_MorgRF1'), 
#             'RDK': ('FP_MorgRF', 'FP_MorgRF1'), 
#             'Morgan' : ('FP_MorgRF', 'FP_MorgRF1')}
#     return fpd[fp]

def loadFingerprint(datadir, fpid):
    fpi = fingerprint()[fpid]
    fpfile = os.path.join(datadir, fpi[0]+'.npz')
    # data = np.load(fpfile) RS change
    data = np.load(fpfile, allow_pickle=True)
    fp = data['x']
    fpn = data['y']
    fpparam = fpi[1]
    # Some fingerprints are stored as bit strings
    if fpi[2] == True:
        fp = [DataStructs.CreateFromBitString(z) for z in fp]
    fpfun = fpi[3]
    data.close()
    return fp, fpn, fpparam, fpfun

def loadFingerprint2(datadir, fpid):
    fpi = fpid
    fpfile = os.path.join(datadir, fpi+'.npz')
    # data = np.load(fpfile) RS change
    data = np.load(fpfile, allow_pickle=True)
    fp = data['x']
    fpn = data['y']
    fpr = data['z']
    data.close()

    # make it into a dictionary
    fpDict={}
    for i, r in enumerate(fpr):
        n = fpn[i]
        d = fp[i]

        if r not in fpDict.keys():
            fpDict[r]={}

        if n not in fpDict[r].keys():
            fpDict[r][n]={}

        fpDict[r][n] =d
    
    return fp, fpn, fpr, fpDict


def loadFingerprint3(datadir, fpid):
    # datadir, fpid = arg.datadir, 'FP_MorgRF'
    fpi = fpid
    fpfile = os.path.join(datadir, fpi+'.npz')
    # data = np.load(fpfile) RS change
    data = np.load(fpfile, allow_pickle=True)
    fp = data['x']
    fpn = data['y']
    fpr = data['z']
    fpd = data['d']
    data.close()

    # make it into a dictionary
    fpDict={}
    rfDist={}
    for i, r in enumerate(fpr):
        n = fpn[i]
        d = fp[i]
        dist = fpd[i]

        if r not in fpDict.keys():
            fpDict[r]={}
            rfDist[r]={}

        if n not in fpDict[r].keys():
            fpDict[r][n]={}
            rfDist[r][n]={}

        fpDict[r][n] = d
        rfDist[r][n]=  {x.split("=")[0]: x.split("=")[1] for x in dist }
    
    return fp, fpn, fpr, fpDict, rfDist


def storeReaction(smi, rfile):
    left, right = smi.split('>>')
    subs = left.split('.')
    prods = right.split('.')
    sd, pd = ({},{})
    for s in subs:
        if s not in sd:
            sd[s] = 0
        sd[s] += 1
    for p in prods:
        if p not in pd:
            pd[p] = 0
        pd[p] += 1
    rsp = {rfile: (sd, pd)}
    return rsp

def getReaction(rfile):
    rxn = rdChemReactions.ReactionFromRxnFile(rfile)
    smi = rdChemReactions.ReactionToSmiles(rxn)
    return storeReaction(smi, rfile), smi

def getReactionFromSmiles(smi, rxnfile):
    smi = '.'.join([x for x in smi.split('>>')[0].split('.') if x !='*' and x!= '[*]']) + '>>' + '.'.join([x for x in smi.split('>>')[1].split('.') if x !='*' and x!= '[*]'])
    rxn = rdChemReactions.ReactionFromSmarts(smi)
    mdl = rdChemReactions.ReactionToRxnBlock(rxn)
    with open(rxnfile, 'w') as handler:
        handler.write(mdl)
    return storeReaction(smi, rxnfile), smi

def getReactionFromSmilesFile(smartsfile, rxnfile):
    with open(smartsfile) as handler:
        smarts = handler.readline()
    rxn = rdChemReactions.ReactionFromSmarts(smarts)
    smi = rdChemReactions.ReactionToSmiles(rxn)
    mdl = rdChemReactions.ReactionToRxnBlock(rxn)
    with open(rxnfile, 'w') as handler:
        handler.write(mdl)
    return storeReaction(smi, rxnfile), smi

def getClosest(smi, fpfile, th=0.8, fp=None, fpn=None, fpp=None, fpfun=None, marvin=False):
    dist = {}
    if fp is None:
        print('Reading fingerprints')
        data = np.load(fpfile)
        fp = data['x'] 
        fpn = data['y']
        fpp = 8
        fpfun = GetMorganFingerprint
        data.close()

    targetMol = Chem.MolFromSmiles(smi)
    # If RDkit fails, we sanitize first using molconvert from ChemAxon, which is more robust
    if targetMol is None and marvin:
        try:
            cmd = ['molconvert', 'mol', smi]
            cmd2 = ['molconvert', 'smiles']
            job = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            job2 = subprocess.Popen(cmd2, stdin=job.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            job.stdout.close()
            out, err = job2.communicate()
            targetMol = Chem.MolFromSmiles(out)
        except:
            pass
    if fpp is not None:
        # targetFp = fpfun(targetMol, fpp) # targetMol = Chem.MolFromSmiles(smi), fpp=5
        info1 = {}
        targetFp = AllChem.GetMorganFingerprint(targetMol, 8, bitInfo=info1, invariants=AllChem.GetConnectivityInvariants(targetMol, includeRingMembership=False))
    else:
        # targetFp = fpfun(targetMol)
        info1 = {}
        targetFp = AllChem.GetMorganFingerprint(targetMol, 8, bitInfo=info1, invariants=AllChem.GetConnectivityInvariants(targetMol, includeRingMembership=False))
        
    tn = DataStructs.BulkTanimotoSimilarity(targetFp, list(fp))
    for i in sorted(range(0, len(tn))):
        dist[fpn[i]] = tn[i]
    return dist, fp, fpn

def bulkTani(targetFp, fp, fpn):
    tn = DataStructs.BulkTanimotoSimilarity(targetFp, list(fp))
    dist={}
    for i in sorted(range(0, len(tn))):
        dist[fpn[i]] = tn[i]
    return dist

def getReactants(equation):
    reactants = {}
    for x in equation.split(' + '):
        n, c = x.split(' ')
        try:
            n = int(n)
        except:
            pass
        reactants[c] = n
    return reactants

def reacSubsProds(dbfile):
    # dbfile = '/home/ruth/code/update_selenzyme/run_folder_min_dec/data_2023/data2/reac_prop.tsv'
    rsp = {}
    for l in open(dbfile):
        if l.startswith('#'):
            continue
        m = l.rstrip().split('\t')
        rid = m[0]
        subs = {}
        prods = {}
        m = l.rstrip().split('\t')
        # print(m[1])
        left, right = m[1].split(' = ')
        subs = getReactants(left)
        prods = getReactants(right)
        # ec = m[4]
        ec = m[3]
        if len(subs) > 0 and len(prods) > 0:
            rsp[rid] = (subs, prods, ec)
    return rsp

def getStructs(dbfile):
    structs = {}
    for l in open(dbfile):
        if l.startswith('#'):
            continue
        m = l.rstrip().split('\t')
        cid = m[0]
        smiles = m[6]
        if len(smiles) > 0:
            structs[cid] = smiles
    return structs

def getRSim(s1, p1, s2, p2, sim):
    cl = {'s1': s1, 'p1': p1, 's2': s2, 'p2':p2}
    ss = {} 
    simm = {}
    pairs = [('s1','s2'), ('s1', 'p2'), ('p1', 's2'), ('p1', 'p2')]
    compPairs = {}
    for p in pairs:
        #print('\n', p)
        pairings = set()
        simm[p] = {}
        compPairs[p]=[]

        for x in cl[p[0]]:
            simm[p][x] = (0.0, x, None)
            if x in sim:
                for y in cl[p[1]]:
                    if y in sim[x]:
                        pairings.add( (sim[x][y], x, y) )
                        #print(sim[x][y], y, x) 
                           # import pdb
                            #pdb.set_trace()
        found = {'left': set(), 'right': set()}
        for v in sorted(pairings, key = lambda h: -h[0]):
            if v[1] not in found['left'] and v[2] not in found['right']:
                # if similarity is greater that zero
                if v[0] > simm[p][v[1]][0]:
                    simm[p][v[1]] = v
                    found['left'].add(v[1])
                    found['right'].add(v[2])
                    compPairs[p].append([v[1], v[2]])
        s = []
        for x in simm[p]:
            s.append(simm[p][x][0])
        if len(s) > 0:
            ss[p] = sum(s)/len(s)
        else:
            ss[p] = 0.0
    S1 = math.sqrt(ss[pairs[0]]**2 + ss[pairs[3]]**2)/math.sqrt(2)
    S2 = math.sqrt(ss[pairs[1]]**2 + ss[pairs[2]]**2)/math.sqrt(2)
    #for k, v in simm.items(): print(k, [x[0] for x in v.values()])
    return(S1, S2, compPairs)



def generate_score( subsRF, prodsRF, subsRF_DB, prodsRF_DB, subProdPairs):         
    # compare the query RF for the database RF
    subsRF_DB = [x for x in subsRF_DB if x !=[]]
    prodsRF_DB = [x for x in prodsRF_DB if x !=[]]
    
    if len(subsRF_DB) == 0 or len(prodsRF_DB) == 0:
        return 0, 0
    
    maxSubSim= [max( DataStructs.BulkTverskySimilarity(x, subsRF_DB, 1, 0  ) ) for x in subsRF]
    maxProdSim= [max( DataStructs.BulkTverskySimilarity(x, prodsRF_DB, 1, 0  ) ) for x in prodsRF]    
    
    meanSub = sum(maxSubSim)/len(maxSubSim)
    meanProd = sum(maxProdSim)/len(maxProdSim)
    
    return meanSub, meanProd          



def generate_score2( queryRF, databaseRF, subProdPairs ):         
    # compare the query RF for the database RF
    databaseRF = {k:v for k, v in databaseRF.items() if v !=[]}
    
    if len(databaseRF) == 0:
        return None, None
    
    scores={}
    for p, x in subProdPairs.items():
        scores[p]=[]
        for comp1, comp2 in x:
            if comp1 in queryRF.keys() and comp2 in databaseRF.keys():
                #ts = DataStructs.TanimotoSimilarity(queryRF[comp1], databaseRF[comp2])
                # Tversky because the database has RFs from all the comparisons
                ts = DataStructs.TverskySimilarity(queryRF[comp1], databaseRF[comp2], 1, 0  )
                scores[p].append(ts)
                
        if len(scores[p])>0:
            if sum(scores[p])>0:        
                scores[p] = sum(scores[p])/len(scores[p])
            else:
                scores[p] = 0

    if scores[('s1', 's2')] != [] and scores[('p1','p2')] !=[]:
        S1 = math.sqrt(scores[('s1', 's2')]**2 + scores[('p1','p2')]**2) /math.sqrt(2)
    else:
        S1 = None
    if scores[('s1', 'p2')] !=[] and scores[('p1','s2')] != []:
        S2 = math.sqrt(scores[('s1', 'p2')]**2 + scores[('p1','s2')]**2) /math.sqrt(2)   
    else:
        S2 = None
    return S1, S2


def generate_RFscore( subsRF, prodsRF, subsRF_DB,prodsRF_DB):         
    # compare the query RF for the database RF
    subsRF_DB = [x for x in subsRF_DB if x !=[]]
    prodsRF_DB = [x for x in prodsRF_DB if x !=[]]
    
    if len(subsRF_DB) == 0 or len(prodsRF_DB) == 0:
        return 0, 0
    
    maxSubSim= [max( DataStructs.BulkTverskySimilarity(x, subsRF_DB, 1, 0  ) ) for x in subsRF]
    maxProdSim= [max( DataStructs.BulkTverskySimilarity(x, prodsRF_DB, 1, 0  ) ) for x in prodsRF]    
    
    meanSub = sum(maxSubSim)/len(maxSubSim)
    meanProd = sum(maxProdSim)/len(maxProdSim)
    
    return meanSub, meanProd 

#################################################
                        
# def get_MCS2(comps1, comps2, fragSize):
    
#     # process the comps
#     Comps1 = []
#     sucesses = []
#     for c1 in comps1:
#         comp, atomMap, fpM, info, fragAtoms = process_comp(c1, fragSize)
#         Comps1.append([comp, atomMap, fpM, info, fragAtoms])
#         if comp is not None:
#             sucesses.append(c1)

#     Comps2 = []
#     for c2 in comps2:
#         comp, atomMap, fpM, info, fragAtoms = process_comp(c2, fragSize)
#         Comps2.append([comp, atomMap, fpM, info, fragAtoms])

#     # for each compound in comps1 get the biggest MCS for the compounds in comp2, then use this largest MCS to get the reacting frags
#     rfrags =[]
#     adj = []
#     fps=[]
#     noOverlap=0
#     totalOverlap=0
#     cantMap=0
#     noComp2=0
#     GetSubstructFail=0
    
#     for x in Comps1:
#         comp1, atomMap1, fpM1, info1, fragAtoms1 = x
        
#         if comp1 is None:
#             print('no comp1')
#             rfrags.append([])
#             adj.append([])  
#             fps.append([])
#             continue
#         fps.append(fpM1)
#         comp1RFs={}
#         comp1RFAs={}
        
#         # for each comp2 get the biggest mcs
#         m=0 # number of atoms
#         # biggestMcs = [] # mcs smarts
#         # i = [] # index
#         for i1, y in enumerate(Comps2):
#             comp2, atomMap2, fpM2, info2, fragAtoms2 = y
#             if comp2 is None:
#                 noComp2+=1
#                 continue


#             res = Chem.rdFMCS.FindMCS([comp1, comp2], timeout=1)
#             m = res.numAtoms  
#             # print('m', m)
#             # this comp2 might not have RFs but others will
#             if m==0:
#                 noOverlap+=1
#                 continue
#             bMcs = res.smartsString
            
#             # get the reacting frags
#             if m < comp1.GetNumAtoms():

#                 bMcs = Chem.MolFromSmarts(bMcs)
#                 # outerAtoms = set([y  for x in comp1.GetSubstructMatches(bMcs) for y in set(range(0,comp1.GetNumAtoms()))- set(x)])

#                 # mcsAtoms=set([y for x in comp1.GetSubstructMatches(bMcs) for y in x - ])
#                 # reactAtoms = [set([k]) | v-outerAtoms for k, v in atomMap1.items() if k in outerAtoms and len(v-outerAtoms)>0]
#                 reactAtoms = [set([x]) | atomMap1[x]-set(substruct)  for substruct in comp1.GetSubstructMatches(bMcs) for x in substruct if len(atomMap1[x]-set(substruct))>0 ]
#                 if len(reactAtoms)==0 and len(comp1.GetSubstructMatches(bMcs))==0:
#                     print("GetSubstructMatches failed")
#                     GetSubstructFail +=1
#                     continue
                
#                 if len(reactAtoms) == comp1.GetNumAtoms():
#                     cantMap+=1
#                     continue
                
#                 # get all the overlapping fragments
#                 rf = {k:v for k, v in fragAtoms1.items() for reacAtomSet in reactAtoms if len(v.intersection(reacAtomSet))>0  }
                
#                 # # separate the smallest fragments from the larger ones
#                 rf1, rfAdj1 = get_common_atoms_frags(rf, fragAtoms1, reactAtoms)
#                 if len(rf1)==0:
#                     print('problem overlapping frags not substrate')
#                     crash
                    
#                 for k, v in rf1.items(): 
#                     if k not in comp1RFs.keys():
#                         comp1RFs[k] =v
                        
#                 for k, v in rfAdj1.items(): 
#                     if k not in comp1RFAs.keys():
#                         comp1RFAs[k] =v 
                        
                        
#             elif comp2.GetSubstructMatches(comp1):
#                 ### if comp is a substructure of comp 1
#                 if comp1.GetNumAtoms() == comp2.GetNumAtoms() == m:
#                     print('compound unchanged')
#                     totalOverlap+=1
#                     continue
                    
#                 # print('the whole compound is substructure')
#                 # get the atoms in comp2 that are in comp1
#                 mapAt= comp2.GetSubstructMatches(comp1) 
#                 # map the atoms comp2:comp1
#                 mapAtoms = {}
#                 for mapping in mapAt:
#                     for i, x in enumerate(mapping):
#                         if x not in mapAtoms.keys():
#                             mapAtoms[x]=set()
#                         mapAtoms[x].add(i)

#                 # get the atoms on the reacting end of comp2
#                 mcsAtoms2=set(mapAtoms.keys())
#                 #reactAtoms2 = [set([k]) | v-mcsAtoms2 for k, v in atomMap2.items() if k in mcsAtoms2 and len(v-mcsAtoms2)>0]
                
#                 outerAtoms2 = set([y  for x in comp2.GetSubstructMatches(comp1) for y in set(range(0,comp2.GetNumAtoms()))- set(x)])
#                 if len(outerAtoms2) == comp1.GetNumAtoms():
#                     cantMap+=1
#                     continue

#                 # reactAtoms2 = [set([k]) | v-outerAtoms2 for k, v in atomMap2.items() if k in outerAtoms2 and len(v-outerAtoms2)>0]
#                 reactAtoms2 = [set([x]) | atomMap2[x]-set(substruct)  for substruct in comp2.GetSubstructMatches(comp1) for x in substruct if len(atomMap2[x]-set(substruct))>0 ]
                         
#                 # map reacting atoms onto compound1
#                 reactAtoms = []
#                 for x in reactAtoms2:
#                     for y in x:
#                         if y in mapAtoms.keys():
#                             reactAtoms.append(mapAtoms[y])
                
#                 # if the reacting Atoms aren't touching then separate them out
#                 removeI = []
#                 addI = []
#                 for x in reactAtoms:
#                     if len(x) == 1: continue
#                     linked=set()
#                     for i in x:
#                         linked= linked|atomMap1[i]
                        
#                     if x == (x - linked):
#                         removeI.append(x)
#                         for y in x:
#                             addI.append(set([y]))
                 
#                 reactAtoms = reactAtoms + addI
#                 reactAtoms = [x for x in reactAtoms if x not in removeI]
                
#                 # get all the overlapping fragments
#                 rf = {k:v for k, v in fragAtoms1.items() for reacAtomSet in reactAtoms if len(v.intersection(reacAtomSet))>0  }
                
#                 # # separate the smallest fragments from the larger ones
#                 rf1, rfAdj1 = get_common_atoms_frags(rf, fragAtoms1, reactAtoms)

#                 if len(rf1)==0:
#                     print('problem overlapping frags subset')
#                     crash
               
#                 for k, v in rf1.items(): 
#                     if k not in comp1RFs.keys():
#                         comp1RFs[k] =v
                        
#                 for k, v in rfAdj1.items(): 
#                     if k not in comp1RFAs.keys():
#                         comp1RFAs[k] =v
                        
#             else:
#                 #changed bonds same atoms
#                 # bMcs = Chem.MolFromSmarts(bMcs)
#                 # compMCS, atomMapMCS, fpMMCS, infoMCS, fragAtomsMCS = process_comp(Chem.MolToSmiles(bMcs))
#                 bMcs = Chem.MolFromSmarts(bMcs)
                
#                 if bMcs.GetNumBonds() == comp1.GetNumBonds():
#                     # print('compound unchanged')
#                     totalOverlap+=1
#                     continue             
                
#                 # map the atoms
#                 mapAt= comp1.GetSubstructMatches(bMcs) 
#                 mapAtoms = {}
#                 for mapping in mapAt:
#                     for i, x in enumerate(mapping):
#                         if x not in mapAtoms.keys():
#                             mapAtoms[x]=set()
#                         mapAtoms[x].add(i)
#                         # mapAtoms[i]=set()
#                         # mapAtoms[i].add(x)
                
#                 # get the bonds in MCS
#                 bondsMCS = [ [x.GetBeginAtomIdx(), x.GetEndAtomIdx(), str(x.GetBondType())] for x in bMcs.GetBonds()]
#                 bondsMCSR = [[x[1], x[0], x[2]] for x in bondsMCS]
#                 # get the bonds for comp1
#                 bondsComp1 = [ [x.GetBeginAtomIdx(), x.GetEndAtomIdx(), str(x.GetBondType())] for x in comp1.GetBonds()] 

#                 reactAtoms = []
#                 for x in bondsComp1:
#                     if len(list(mapAtoms[x[0]]))>1 or len(list(mapAtoms[x[1]]))>1:
#                         print('death')
#                         crash
#                     a1, a2, t = list(mapAtoms[x[0]])[0], list(mapAtoms[x[1]])[0], x[2]
#                     if [a1, a2, t] not in bondsMCS and [a1, a2, t] not in bondsMCSR:# and [a1, a2, t] not in atomMapList2 :
#                         print('comp1', x, "\n",'mcs', [a1, a2, t], "\n\n")
#                         reactAtoms.append(set([x[0], x[1]]))
                
#                 # get all the overlapping fragments
#                 rf = {k:v for k, v in fragAtoms1.items() for reacAtomSet in reactAtoms if len(v.intersection(reacAtomSet))>0  }
                
#                 # # separate the smallest fragments from the larger ones
#                 rf1, rfAdj1 = get_common_atoms_frags(rf, fragAtoms1, reactAtoms)

#                 if len(rf1)==0:
#                     print('problem changed bonds/same atoms')
#                     crash
               
#                 for k, v in rf1.items(): 
#                     if k not in comp1RFs.keys():
#                         comp1RFs[k] =v
                        
#                 for k, v in rfAdj1.items(): 
#                     if k not in comp1RFAs.keys():
#                         comp1RFAs[k] =v                    
                    
                        
#             if len(rf1)==0 and totalOverlap==0 and cantMap == 0:
#                 print('no RFs for a comp1')
#                 crash
        
#         # store all the fragments from mapping comp1 onto all comp2s
#         if len(comp1RFs)>0:
#             rfrags.append(dict_to_sparse(comp1RFs)[0])
#         else:
#             rfrags.append([])
#             if len(comps2)-noComp2 == totalOverlap:
#                 print('compounds are all the same')
#             elif len(comps2)-noComp2==noOverlap:
#                 print('no overlap between comps1 and comps2')
#             elif len(comps2)-noComp2==cantMap:
#                 print('couldnt map comps1 and comps2')
#             elif len(comps2)-noComp2==GetSubstructFail:
#                 print('GetSubstructure failed')
#             else:
#                 print('combined problem', totalOverlap,noOverlap, cantMap, GetSubstructFail )
#                 # crash
                       
#         if len(comp1RFAs)>0:
#             adj.append(dict_to_sparse(comp1RFAs)[0])
#         else:
#             adj.append([])
            
#         # print('comp1 end rfrags', comp1RFs.keys())

#     if len(rfrags)==0: 
#         print('no rfs')
#         crash
                                                     
#     return  rfrags, adj, fps, sucesses

# def process_comp(comp1, fragSize):
#     # get all the details for the fingerprints
    
#     if type(comp1)==str:
#         try:
#             comp1 = Chem.MolFromSmiles(comp1)
#         except:
#             print('failed MolFromSmiles')
#             return(None, None, None, None, None)
            
#     try:
#         # Chem.SanitizeMol(comp1)
#         # rdMolStandardize.Uncharger().uncharge(comp1)
#         comp1 = neutralize_atoms(comp1)
#         comp1 = mol_with_atom_index(comp1)
 
#         # Draw.DrawMorganBits(prints, molsPerRow=4, legends=[str(x) for x in list(fpM1.GetNonzeroElements().keys()) ])
               
#         # get each atoms adjacent neighbours
#         atomMap1 = getAtoms(comp1)
        
#         # Morgan fingerprints
#         info1={}
#         fpM1 = AllChem.GetMorganFingerprint(comp1, fragSize, bitInfo=info1, invariants=AllChem.GetConnectivityInvariants(comp1, includeRingMembership=False))
        
#         # get all the atoms in each fragment, returns {fragname:[{atoms}, startnode, radius]}           
#         fragAtoms1 = getAtomFragments(fpM1, info1, atomMap1)   
    
#         return(comp1, atomMap1, fpM1, info1, fragAtoms1)
#     except:
#         return(None, None, None, None, None)
     

def neutralize_atoms(mol):
    if type(mol)==str:
        mol = Chem.MolFromSmiles(mol)
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
    return mol

def mol_with_atom_index( mol ):
    atoms = mol.GetNumAtoms()
    for idx in range( atoms ):
        mol.GetAtomWithIdx( idx ).SetProp( 'molAtomMapNumber', str( mol.GetAtomWithIdx( idx ).GetIdx() ) )
    return mol

def getAtoms(comp1):
    atomMap = {}
    # loop through the bonds to get each atoms adjacent atoms
    for b in comp1.GetBonds():
        start = b.GetBeginAtomIdx()
        end = b.GetEndAtomIdx()
        
        if start not in atomMap.keys():
            atomMap[start] = set()
        atomMap[start].add(end)
 
        if end not in atomMap.keys():
            atomMap[end] = set()
        atomMap[end].add(start)     
    return atomMap


def getAtomFragments(fp1, info1, atomMap):
    # get all the atoms in each fragment
    fragAtoms = {}

    for bitName in fp1.GetNonzeroElements().keys():
        bitInfo = info1[bitName]
    
        for x in bitInfo:
            
            startNode, rad = x
            scope = set([startNode])
            # if radius is zero, the start node is the only node in the scope
            if rad ==0:
                # print(bitName, scope)
                fragAtoms[str(bitName)+'_'+ str(startNode)] = scope
                
            # if not iterate through atomMap to get the other nodes
            else:
                adjNodes = atomMap[startNode]
                scope = scope | adjNodes
                i=0
                while i < rad:
                    adjNodes2 = set()
                    for x in adjNodes:
                        adjNodes2 = adjNodes2 | atomMap[x]
                        scope = scope | adjNodes2
                        i+=1
                fragAtoms[str(bitName)+'_'+ str(startNode)] = scope
                
    return fragAtoms

def get_common_atoms_frags(rf, fragAtoms, reactAtoms):
    # reaction at one site:

    minRF = {}
    maxRF = {}
    
    for Aset in reactAtoms:
        trueRF = {k:v for k,v in rf.items() if len(v.intersection(Aset)) == len(Aset)}        
        minsize = min([len(x) for x in trueRF.values()])
        minRF.update({k:v for k, v in trueRF.items() if len(v) == minsize})
        maxRF.update({k:v for k, v in fragAtoms.items() if k not in minRF.keys() and len(v.intersection(Aset))>0 })
    
    return minRF, maxRF


def generate_RFscore2(subSmiles, queryRF, queryDists, s2, r2, rfDict, rfdist, subProdPairs):  
    # subSmiles, queryRF, queryDists, s2, r2, rfDict, rfdist, subProdPairs = list(s1.keys()), queryRF, queryDists, s2, r2, rfDict, rfdist, subProdPairs[('s1', 's2')] 
    # subSmiles, queryRF, queryDists, s2, r2, rfDict, rfdist, subProdPairs = list(p1.keys()), queryRF, queryDists, p2, r2, rfDict, rfdist, subProdPairs[('p1', 'p2')]
    
    subList = list(subSmiles)
    # # deal with empty queries with RFs
    # queryNoRFs = [x for x in subList if x not in queryRF.keys() or len(queryRF[x])==0]
    # unweightedScores = {}
    # for x in queryNoRFs:
    #     y =  [y for x1, y in subProdPairs if x1==x]
    #     if y not in rfDict[r2] or len(rfDict[r2]) ==0:
    #         unweightedScores[x, y][0] = 1
            
        
    
    querySubsRF = {x:queryRF[x] for x in subList if x in queryRF.keys() and len(queryRF[x])>0}
    querySubsDist = {x:queryDists[x] for x in subList if x in queryRF.keys() and len(queryDists[x])>0}
    
    dbList = list(s2.keys())
    dbSubsRF = {x:rfDict[r2][x] for x in dbList if x in rfDict[r2]}
    dbSubsDist = {x:rfdist[r2][x] for x in dbList if x in rfDict[r2]}
    
    spPairs = ["_".join([x[0], x[1]]) for x in subProdPairs]
    unweightedScores = {}
    
    for sb1, qSubRF in querySubsRF.items() :
        qSubsDist = querySubsDist[sb1]
        qSortDists = {}
        
        # iterate through the query fragments
        for k, v in qSubsDist.items():
            dists1 = [int(x[0]) for x in v]
            for d in dists1:
                if d not in qSortDists.keys():
                    qSortDists[d]=[]
                qSortDists[d].append(k)


        # iterate through the db fragments
        for sb2, dSubRF in dbSubsRF.items():
            if dSubRF =={}: 
                print('empty', sb2)
                continue
            if '_'.join([sb1, sb2]) not in spPairs:
                continue
            #print(sb1, sb2)


            dSubsDist = dbSubsDist[sb2]
            dSortDists = {}
            for fragNo, k in enumerate(dSubRF.GetNonzeroElements().keys()):
                dists2 =  [int(y.split('_')[1]) for y in dSubsDist[str(k)].split('|')  ]
                for d in dists2:            
                    if d not in dSortDists.keys():
                        dSortDists[d]=[]
                    dSortDists[d].append(k)
           
            # get jaccard sim
            unweightedScores[sb1, sb2] =  {}
            
            for fragDist in set(qSortDists.keys()).union(set(qSortDists.keys())):
                if fragDist not in set(qSortDists).intersection(set(dSortDists)):
                    unweightedScores[sb1, sb2][fragDist] =0
                    
                else:
                    intersect = {frag: min([dSortDists[fragDist].count(frag), qSortDists[fragDist].count(frag)]) for frag in set(dSortDists[fragDist] + qSortDists[fragDist]) if frag in dSortDists[fragDist] and frag in qSortDists[fragDist]}
                    intersectScore = sum(intersect.values())
                    # unionScore = (len(qSortDists[fragDist]) + len(dSortDists[fragDist])) /2
                    # if intersectScore/unionScore > 0.95:
                    #     print(sb1, sb2, intersectScore/unionScore)
                    # unweightedScores[sb1, sb2][fragDist] = intersectScore/unionScore
                    unweightedScores[sb1, sb2][fragDist] = intersectScore/len(qSortDists[fragDist])
    
    
    weightedScores={}
    for pair, scores in unweightedScores.items():
        # weights = simpleWeight(max(scores.keys()))
        weights = log_fun_weight(max(scores.keys()))
        # if len(weights) != len(scores): print(len(weights), len(scores))
        weightedScores[pair] = sum([ scores[x]*weights[x] if x in scores and scores[x]>0 else 0 for x in list(range(0, max(scores.keys())+1)) ])


    return(weightedScores)



def log_fun_weight(c):   
    # log function 
    weights1 = []
    k = 0.75
    i = 3
    for x in list(range(0, c+1)):
        weights1.append((1 / (1 + math.e ** (k*(x - i)) ) )+0)
    return([ x/sum(weights1) for x in weights1])

def simpleWeight(c):
    n=1
    m=[]
    for x in list(range(0, c)):
        m.append(n)
        n = n*1.2
    if len(m)>0:
        m.append(m[-1])
    else:
        m.append(1)
    return([x/sum(m) for x in sorted(m, reverse=True)])




def dict_to_sparse(rfs, MorganRad):
    rfList1 = [x for k, v in rfs.items() for x in [k]*len(v) ]
    # make an empty sparse int vector
    SparseIntVect1 = AllChem.GetMorganFingerprint(Chem.MolFromSmiles(''), MorganRad)
    # update the empty sparse int vector with the reacting fragments
    SparseIntVect1.UpdateFromSequence(rfList1)    
    return(SparseIntVect1)


# def dict_to_sparse(rfs):
#     # make sure the input is a list
#     if type(rfs) != list:
#         rfs = [rfs]
    
#     # make into sparse int vectors
#     rf_list = []
#     for rf in rfs:
#         RF1 =  [int(x.split('_')[0]) for x in rf.keys()]
    
#         # Make empty UIntSparseIntVect, then update them with the sparse int vectors       
#         rf1 = AllChem.GetMorganFingerprint(Chem.MolFromSmiles(''), 5)
#         rf1.UpdateFromSequence(RF1)
#         rf_list.append(rf1)
#     return(rf_list)

def tidy(queryRF, prodsRF, hits1, queryAdj, prodsAdj, hits2):
    queryRF = dict(zip(hits1, queryRF))
    queryAdj = dict(zip(hits1, queryAdj))
    prodsRF = dict(zip(hits2, prodsRF))
    prodsAdj = dict(zip(hits2, prodsAdj))
    queryRF.update(prodsRF)
    queryAdj.update(prodsAdj)
    return(queryRF, queryAdj)




# def drawReactions(smi1, smi2):

#     rxn1 = rdChemReactions.ReactionFromSmarts(smi1, useSmiles=True)
#     rxn2 = rdChemReactions.ReactionFromSmarts(smi2, useSmiles=True)
#     return [Draw.ReactionToImage(rxn1), Draw.ReactionToImage(rxn2)]
                        
def test():
    #import random
    #import Selenzy_local
    chem_prop = "/home/ruth/code/update_selenzyme/new_data/chem_prop.tsv"
    chem_prop= pd.read_csv(chem_prop, sep='\t', skiprows=347)
    
    reac_prop= "/home/ruth/code/update_selenzyme/new_data/reac_prop.tsv"
    reac_prop= pd.read_csv(reac_prop, sep='\t', skiprows=347)  
    
    reac_smiles = "/home/ruth/code/update_selenzyme/run_folder/selenzyme2/selenzyPro/data/reac_smi.csv"
    reac_smi = pd.read_csv(reac_smiles)
    
    pc = Selenzy_local.readData('/home/ruth/code/update_selenzyme/run_folder/selenzyme2/selenzyPro/data/')
    
    sample = random.sample(range(reac_prop.shape[0]), 10000)

    
    HighWholeCompScores, HighRFScores, HighW_lowRF, HighRF_lowW = [], [], [], []
    i=0   

    for x in sample: # problems
        testReaction = reac_prop.iloc[x]['#ID']
        print('\n\nrebuilding from chem prob', testReaction)
        eq = reac_prop['mnx_equation'][reac_prop['#ID']==testReaction].values[0].split('=')
        subs, prods = re.findall(r'MNXM\d+', eq[0]), re.findall(r'MNXM\d+', eq[1])
        
        # filter for those without smiles or only 
        subs = set([x for x in subs for y in chem_prop['SMILES'][chem_prop['#ID']== x].values if y ==y])
        prods = set([x for x in prods for y in chem_prop['SMILES'][chem_prop['#ID']== x].values if y ==y]) 

        if set(subs) == set(prods):
            print('no transformations', testReaction)
            continue
        
        cofactors = set(['MNXM3', 'MNXM7', 'MNXM8', 'MNXM9', 'MNXM10', 'MNXM13', 'WATER'])
        if not subs - cofactors or not prods - cofactors:
            print('only cofactors')
            continue
        
        if testReaction not in fprRF:
            continue
               
        subsC = {x:int(re.findall(r'\d+ '+x, eq[0])[0].split(' ')[0]) for x in subs}
        prodsC = {x:int(re.findall(r'\d+ '+x, eq[1])[0].split(' ')[0]) for x in prods}
        subsmi = [y for x in subs for y in list(chem_prop['SMILES'][chem_prop['#ID']== x].values)*subsC[x] if y ==y]
        prosmi = [y for x in prods for y in list(chem_prop['SMILES'][chem_prop['#ID']== x].values)*prodsC[x] if y ==y]
        if len(set(subsC.values()) | set(prodsC.values()))>1: print('multi-values')
        
        if set(subsmi) == set(prosmi):
            print('no transformations', testReaction)
            continue     
        
        smi = '.'.join(subsmi) +'>>'+'.'.join(prosmi)           
        arg = arguments(['/home/ruth/code/update_selenzyme/run_folder/selenzyme2/selenzyPro/data/', 
            'Morgan',
            '-smarts', smi, 
            '-out', '/home/ruth/code/update_selenzyme/RSquickRsim_new.txt'] )
        
        ##### run #####
        try:
            highWholeCompScores, highRFScores, highW_lowRF, highRF_lowW = run(arg, pc)
        except: 
            pass
        ###############
        highWholeCompScores.sort(key=lambda x:x[2] == testReaction, reverse=True)
        HighWholeCompScores.extend([[testReaction] + x +[smi] for x in highWholeCompScores])
        highRFScores.sort(key=lambda x:x[2] == testReaction, reverse=True)
        HighRFScores.extend([[testReaction] + x +[smi] for x in highRFScores])
        HighW_lowRF.append([testReaction] + highW_lowRF +[smi])
        HighRF_lowW.append([testReaction] +highRF_lowW +[smi] )
        
        i+=1
        print('finished run ', i)
        if i>50: break
    # ([max([S1, S2]), RF_score, r2,  srf, prf])
    
    hw = pd.DataFrame(HighWholeCompScores, columns = ['QueryMNXR', 'wholeCompScore', 'RFScore', 'DBMNXR', 'DBsubs', 'DBprods', 'Qsmiles'])
    hw.to_csv('highWhole.csv', index=False)
    hits1 = [1 if x in hw.DBMNXR[hw.QueryMNXR == x].values else 0  for x in set(hw.QueryMNXR)]

    hrf = pd.DataFrame(HighRFScores, columns = ['QueryMNXR', 'wholeCompScore', 'RFScore', 'DBMNXR', 'DBsubs', 'DBprods', 'Qsmiles'])
    hrf.to_csv('highRF.csv', index=False)
    hits2 = [1 if x in hrf.DBMNXR[hrf.QueryMNXR == x].values else 0  for x in set(hrf.QueryMNXR)]
    
    hwlr = pd.DataFrame(HighW_lowRF, columns = ['QueryMNXR','ScoreDisparity', 'wholeCompScore', 'RFScore', 'DBMNXR', 'DBsubs', 'DBprods', 'Qsmiles'])
    hwlr.to_csv('highWhole_lowRF.csv', index=False)
   
    
    hrlw = pd.DataFrame(HighRF_lowW, columns = ['QueryMNXR','ScoreDisparity', 'wholeCompScore', 'RFScore', 'DBMNXR', 'DBsubs', 'DBprods', 'Qsmiles'])
    hrlw.to_csv('highRF_lowWhole.csv', index=False)
    maxRF = hrlw[hrlw.RFScore>0.9]
    
    # get the similarity of the high RFS and low whole compounds
    hrlw= hrlw[hrlw['RFScore']>0.8]
    for i, row in hrlw.iterrows():
        print(reac_prop.classifs[reac_prop['#ID'] == row.QueryMNXR].values, reac_prop.classifs[reac_prop['#ID'] == row.DBMNXR].values, row.QueryMNXR, row.DBMNXR)


def draw_rfs(comp, fp, bi, rfs):
    prints = [(comp, x, bi) for x in fp.GetNonzeroElements().keys() if x in rfs] #fp.GetOnBits()]
    return(Draw.DrawMorganBits(prints, molsPerRow=4))
        



########################################################################



# def drawFrags(mol, fp, info):
    
#     #mol, fp, info = [comp[0], comp[2], comp[3]]
#     tpls = [(mol, x, info) for x in fp.GetNonzeroElements()]
#     p = Draw.DrawMorganBits(tpls, molsPerRow=5, legends=[str(x) for x in fp.GetNonzeroElements()])
#     return p


def run(arg, pc, rxn_mapper):

    print('in run')
    # fragSize=8

    # read in the data
    if arg.out:
        fileObj = open(arg.out, 'w')    
    if arg.high:
        fileObj = open(arg.high, 'w')
    rsp = reacSubsProds(os.path.join(arg.datadir, 'reac_prop.tsv'))

    smiles = ''
    if arg.rxn is not None:
        rTarget, smiles = getReaction(arg.rxn)
        
    elif arg.smarts is not None:
        rxnfile = os.path.join(os.path.dirname(arg.out), 'reaction.rxn')
        try:
            rTarget, smiles = getReactionFromSmiles(arg.smarts, rxnfile)
        except:
            print('smile could not be processed')
            #return
    elif arg.smartsfile is not None:
        rxnfile = os.path.join(os.path.dirname(arg.out), 'reaction.rxn')
        rTarget, smiles = getReactionFromSmilesFile(arg.smartsfile, rxnfile)        
    elif arg.rid is not None:
        struct = getStructs(arg.chem)
        rTarget = {arg.rid: [{},{}]}
        for side in (0,1):
            for s in rsp[arg.rid][side]:
                if s in struct:
                    rTarget[arg.rid][side][struct[s]] = rsp[arg.rid][side][s]
    else:
        raise Exception('No target')
    
        
    print('reading fp data')
    # Read fingerprint info from preload data if available
    if pc is not None: 
        if arg.fp in pc.fp:
            fp, fpn, fpp, fpfun = pc.fp[arg.fp]
            fpRF, fpnRF, fprRF, rfDict, rfdist = pc.fpr[arg.fp]
        else:
            fp, fpn, fpp, fpfun = loadFingerprint(arg.datadir, arg.fp)
            fpRF, fpnRF, fprRF, rfDict, rfdist = loadFingerprint3(arg.datadir, 'FP_MorgRF')     
    else:       
        fp, fpn, fpp, fpfun = loadFingerprint(arg.datadir, arg.fp)
        fpRF, fpnRF, fprRF, rfDict, rfdist = loadFingerprint3(arg.datadir, 'FP_MorgRF')
  
    
    print('\nGet reacting fragments for the Query reaction') 
    queryRF = {}
    queryDists = {}
    sim={}
    
    for r in rTarget:
        subSmiles = rTarget[r][0].keys()
        prodSmiles = rTarget[r][1].keys()
        if set(subSmiles) == set(prodSmiles):
            print('no transformation')
            
        #comp, atomMap, fpM, info, fragAtoms2, smile, count
        subMols = [process_compF(k, 8) +(k, v,) for k, v in rTarget[r][0].items()]
        prodMols = [process_compF(k, 8) +(k, v,) for k, v in rTarget[r][1].items()]
        #images = []
    
        # AAM using RXNmapper
        RFavail = 0
        try:
            # for unbalances reactions there need to be more atoms on the substrate side
            if sum([x[0].GetNumAtoms() for x in subMols]) >= sum([x[0].GetNumAtoms() for x in prodMols]):
                smile = '.'.join(subSmiles) +'>>'+'.'.join(prodSmiles)
                reactingAtoms, conf = mcs_functions2.rxnMapper_fun(smile, subMols, prodMols, rxn_mapper)
            else:
                smile = '.'.join(prodSmiles) +'>>'+'.'.join(subSmiles)
                reactingAtoms, conf = mcs_functions2.rxnMapper_fun(smile, prodMols ,subMols, rxn_mapper)  
            RFavail = 1
        except RuntimeError: 
            print('AAM failed query - compound too large')

        except ValueError: 
            print('AAM failed query - issue with input smiles')
            if '*' in smile: print('please use smiles without *s')

        except: 
            print('AAM failed query')
 
        
        # measure compound similarities
        for comp in subMols + prodMols:
            smile = comp[5]
            # measure distances
            sim[smile] = bulkTani(comp[2], fp, fpn) 
        
        # get the reacting fragments
        if RFavail == 1:
            for comp in subMols + prodMols:                   
                smile = comp[5]
                if smile not in reactingAtoms: continue
                # rfs are fragNumber:(startAtom, dist), dists are fragNumber: distFromRA, startAtom, RA                
                rfs, dists = reactFragDists(comp, reactingAtoms[smile])  
                # draw_rfs(comp[0], comp[2], comp[3], rfs)
                if len(rfs)>0:
                    queryRF[smile]=rfs
                    queryDists[smile]=dists  


    # sim_comp = []
    # for x in sim: 
    #     t=[]
    #     for k, v in sim[x].items(): 
    #         if v>0.5:
    #             t.append([v, x[0:30], k])
    #     sim_comp.append(sorted(t, reverse=True)[0:5])
             
    # print('similar compounds')
    # for x in sim_comp: 
        # for x2 in x: print(x2)
                
    print("\nget scores")            
    # highWholeCompScores = []
    # highRFScores = []
    # highW_lowRF = []
    # highRF_lowW = []
    
    for r1 in rTarget:
        s1, p1 = rTarget[r1]
        
        for r2 in rsp:   
                        
            s2, p2, ec2 = rsp[r2]   
            if s2 == p2: continue
            S1, S2, subProdPairs = getRSim(s1, p1, s2, p2, sim)
            # if r2=='MNXR94784': crash
            #if max(S1, S2)>0.5: crash
            
            # if aam failed continue loop
            if RFavail ==0:
                if arg.out:
                   print(r1, r2, S1, S2, smiles,  float("Nan"), ec2, file = fileObj)
                if arg.high:
                    print(r1, r2, max([S1, S2]),  smiles,  float("Nan"), ec2, file = fileObj) 
                continue
            
            # if RFs missing for that reaction
            if r2 not in rfDict.keys():
                if arg.out:
                   print(r1, r2, S1, S2, smiles,  float("Nan"), ec2, file = fileObj)
                if arg.high:
                    print(r1, r2, max([S1, S2]),  smiles,  float("Nan"), ec2, file = fileObj) 
                continue                
            
            
            RF_score=float("Nan")
            if r2 in rfDict.keys():
                if S1 > 0 and S2 > 0:
                    subProdPairs = {k:[x for x in v if x[0] in queryRF.keys()]  for k, v in subProdPairs.items()}
                    try:
                        # calculate the scores for the RFs
                        # forward score
                        if S1 >= S2:
                            S_RF = generate_RFscore2(list(s1.keys()), queryRF, queryDists, s2, r2, rfDict, rfdist, subProdPairs[('s1', 's2')] )
                            P_RF = generate_RFscore2(list(p1.keys()), queryRF, queryDists, p2, r2, rfDict, rfdist, subProdPairs[('p1', 'p2')] )
                            if (sum(S_RF.values())+sum(P_RF.values()))>0: 
                                RF_score = (sum(S_RF.values())+sum(P_RF.values())) / len(queryRF) #(len(S_RF) + len(P_RF))
                            else: RF_score=0
                        else:
                            # reverse score
                            S_RF = generate_RFscore2(list(s1.keys()), queryRF, queryDists, p2, r2, rfDict, rfdist, subProdPairs[('s1', 'p2')])
                            P_RF = generate_RFscore2(list(p1.keys()), queryRF, queryDists, s2, r2, rfDict, rfdist, subProdPairs[('p1', 's2')])
                            if (sum(S_RF.values())+sum(P_RF.values()))>0: 
                                # geometric mean 
                                RF_score = math.sqrt(mean(S_RF.values()) * mean(P_RF.values()))
                                # RF_score = (sum(S_RF.values())+sum(P_RF.values())) / len(queryRF) #(len(S_RF) + len(P_RF))
                            else: 
                                RF_score=0
                    except:
                        # print('rf score generation failed', r2)
                        RF_score = float("Nan")


                    # if RF_score > 0.6: print('RF score', RF_score, S1, S2)
                        
                    # if S1>0.9 or S2>0.9 or RF_score>0.8 : 
                    #     print(r2, 'S1', S1, 'S2', S2, 'RF', RF_score)
                    
                    if arg.out:
                        print(r1, r2, S1, S2, smiles, RF_score, ec2,   file = fileObj) 

                    if arg.high:
                        if S1 >= S2:
                            print(r1, r2, S1, smiles,  RF_score, ec2,  file = fileObj)
                        else:
                            print(r1, r2, S2, smiles,  RF_score, ec2,   file = fileObj)                

    
    # highWholeCompScores = sorted(highWholeCompScores, reverse=True)[0:10]  #sorted(highWholeCompScores, reverse=True)
    # highRFScores.sort(key=lambda x:x[1], reverse=True)
    # highRFScores= highRFScores[0:10]
    # return [highWholeCompScores, highRFScores, highW_lowRF, highRF_lowW]




 
def arguments(args=None):
    parser = argparse.ArgumentParser(description='quickRSim Pablo Carbonell, SYNBIOCHEM, 2016')
    parser.add_argument('datadir', help='Data folder')
    parser.add_argument('fp', help='Fingerprint for reactants')
    parser.add_argument('-fpr', help='Fingerprint for reactants', default='FP_MorgRF')
    parser.add_argument('-rxn', 
                        help='Input reaction rxn file')
    parser.add_argument('-rid', 
                        help='Input reaction id')
    parser.add_argument('-smarts', 
                        help='Input reaction SMARTS')
    parser.add_argument('-smartsfile', 
                        help='Input reaction SMARTS file')
    parser.add_argument('-chem', 
                        help='Metanetx chemical structures (if input is reaction id)')
    parser.add_argument('-th', type=float, default=0.8, 
                        help='Similarity threshold [default=0.8]')
    parser.add_argument('-out', 
                        help='Output results in .txt file, please specify file name')
    parser.add_argument('-high', 
                        help='Output results in .txt file with highest similarity score from both forwards and backwards reactions, please specify file name')
    parser.add_argument('-marvin', 
                        help='Call marvin if needed (skip if fails)')

    if args is not None:
        arg = parser.parse_args(args=args)
    else:
        arg = parser.parse_args()
    return arg      



if __name__ == '__main__':
    arg = arguments()

    # smi = 'BrCCOCCBr.CCN(C(C)C)C(C)C.CCOC(C)=O.CN(C)C=O.Cl.NCC(F)(F)CO>>OCC(F)(F)CN1CCOCC1'
    # MNXR189106 / MNXR106859 2-oxoglutarate MNXM20 + 1 CoA ->  succinyl-CoA MNXM741428
    # subProdPairs[('s1', 'p2')]
    # [['CC(C)(COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1OP(=O)([O-])[O-])[C@@H](O)C(=O)NCCC(=O)NCCS',  'MNXM727276'],
    # ['O=C([O-])CCC(=O)C(=O)[O-]', 'MNXM20']]
    # smi = 'O=C([O-])CCC(=O)C(=O)[O-].CC(C)(COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1OP(=O)([O-])[O-])[C@@H](O)C(=O)NCCC(=O)NCCS>>CC(C)(COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1OP(=O)([O-])[O-])[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)CCC(=O)[O-]'

    # #  (E)-coniferol -> coniferyl acetate|2.3.1.224
    # # reaction hits {'MNXR113611'}
    # # smi = 'COc1cc(\C=C\CO)ccc1O>>COc1cc(\C=C\COC(C)=O)ccc1O'
    
    
    # smi = '[NH3+][C@@H](CCC(=O)N[C@@H](CS)C([O-])=O)C([O-])=O.Nc1ncnc2n(cnc12)[C@@H]1O[C@H](COP([O-])(=O)OP([O-])([O-])=O)[C@@H](O)[C@H]1O.OP([O-])([O-])=O>>Nc1ncnc2n(cnc12)[C@@H]1O[C@H](COP([O-])(=O)OP([O-])(=O)OP([O-])([O-])=O)[C@@H](O)[C@H]1O.[NH3+][C@@H](CS)C([O-])=O'
    # smi = 'CC1=C2C3OC3C3=C(C=CC=C3)C2=C(C)C2=CC=CC=C12.[H]O[H]>>CC1=C2[C@H](O)[C@@H](O)C3=C(C=CC=C3)C2=C(C)C2=CC=CC=C12'
    # # smi = 'O=C([O-])CCC(=O)C(=O)[O-].NC(CC(=O)[O-])C(=O)O>>O=C([O-])CC(=O)C(=O)[O-].NC(CCC(=O)[O-])C(=O)O'
    #smi = '[*].CC=C(C)C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)N1C=NC2=C(N)N=CN=C12>>CCC(C)C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)N1C=NC2=C(N)N=CN=C12'
    # smi = 'CC=C(C)C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)N1C=NC2=C(N)N=CN=C12>>CCC(C)C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)N1C=NC2=C(N)N=CN=C12'
    # MNXR106362
    
    
    # # smi = 'O=C1CC(c2ccc(O)c(O)c2)Oc2cc(O)cc(O)c21>>O=C1c2c(O)cc(O)cc2OC(c2ccc(O)c(O)c2)C1O'
    # smi = 'O=C([O-])CCC(=O)[O-].CC(=O)CC(=O)SCCN=C([O-])CCN=C([O-])C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)([O-])[O-]>>CC(=O)CC(=O)[O-].CC(C)(COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)([O-])[O-])C(O)C([O-])=NCCC([O-])=NCCSC(=O)CCC(=O)[O-]'
    # # smi = 'O[C@H](COP([O-])([O-])=O)[C@@H](O)C=O>>O[C@H](COP([O-])([O-])=O)[C@@H](O)[C@H](O)CC(=O)C([O-])=O'
    # smi = 'Cn1cnc2[nH]c(=O)[nH]c(=O)c21>>O=c1[nH]c(=O)c2nc[nH]c2[nH]1'
    
    # # # fp = '/home/ruth/code/update_selenzyme/run_folder_min_dec/data/'

    # smi = 'C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@@]43C)[C@@H]1CC[C@]2(O)C(=O)CO>>C[C@]12C[C@H](O)[C@H]3[C@@H](CCC4=CC(=O)CC[C@@]43C)[C@@H]1CC[C@]2(O)C(=O)CO'    
    # fp = '/home/ruth/code/update_selenzyme/run_folder_min_dec/data_2023/data2/'
    # arg = arguments([fp, 
    #     'Morgan',
    #     '-smarts', smi, 
    #     '-out', '/home/ruth/code/update_selenzyme/RSquickRsim_new.txt'] )
    # pc=None
    
    
    rxn_mapper = RXNMapper()     
    run(arg, pc, rxn_mapper)

 
 

                       
#python quickRsim_local2.py \
#      /home/ruth/code/update_selenzyme/run_folder/selenzyme/selenzyPro/data/ \
#      Morgan \
#      -smarts 'CC(=O)CC(=O)CC(=O)[O-].[H+]>>CC1=CC(O)=CC(=O)O1.[H]O[H]' \
#      -chem /home/ruth/code/update_selenzyme/selenzyme_update/data_google_cloud_edited/data_test/chem_prop.tsv \
#      -out /home/ruth/code/update_selenzyme/run_folder/RSquickRsim_ring_MTrue.txt 
    
    # datadir='/home/ruth/code/update_selenzyme/run_folder/selenzyme/selenzyPro/data/'
    # fp='Morgan'
    # fpr='FP_MorgRF'
    # out='/home/ruth/code/update_selenzyme/run_folder/selenzyme/selenzyPro/uploads/results_quickRsim.txt'
    # smarts='[H]Oc1c(OC([H])([H])[H])c([H])c(C([H])=C([H])C(=O)OC2([H])OC([H])(C([H])([H])O[H])C([H])(O[H])C([H])(O[H])C2([H])O[H])c([H])c1OC([H])([H])[H]>>[H]OC(=O)C([H])([H])C([H])(OC(=O)C([H])=C([H])c1c([H])c(OC([H])([H])[H])c(O[H])c(OC([H])([H])[H])c1[H])C(=O)O[H].[H]OC([H])([H])C([H])(O[H])C([H])(O[H])C([H])(O[H])C([H])(O[H])C([H])=O'
    # th=0.8

