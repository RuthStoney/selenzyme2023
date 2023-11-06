'''
quickRsim (c) University of Manchester 2017

quickRsim is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  Pablo Carbonell, Ruth Stoney SYNBIOCHEM
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
import os
import math
import numpy as np

from rdkit import DataStructs, Chem
from rdkit.Chem import rdChemReactions, AllChem
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprint 
from rdkit.Chem.rdmolops import PatternFingerprint, RDKFingerprint
# from rdkit.Chem import Draw
from statistics import mean
 
from rxnmapper import RXNMapper
from rf_functions import Compounds, ReactingAtoms, ReactingFragments


def fingerprint():
    """ I keep the three fingerprints that give better results in the tests """
    fpd = {'Morgan': ('FP_Morg', 5, False, GetMorganFingerprint, 3)}
    fpd = {'Morgan': ('FP_Morg', 5, False, GetMorganFingerprint, 3)}
    return fpd



def loadFingerprint(datadir, fpid):
    fpi = fingerprint()[fpid]
    fpfile = os.path.join(datadir, fpi[0]+'.npz')
    # data = np.load(fpfile)
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


def loadFingerprintRF(datadir, fpid):
    fpi = fpid
    fpfile = os.path.join(datadir, fpi + '.npz')
    data = np.load(fpfile, allow_pickle=True)
    fp = data['x']
    fpn = data['y']
    fpr = data['z']
    fpd = data['d']
    data.close()

    # make it into a dictionary
    fpDict = {}
    rfDist = {}
    for i, r in enumerate(fpr):
        n = fpn[i]
        d = fp[i]
        dist = fpd[i]

        if r not in fpDict.keys():
            fpDict[r] = {}
            rfDist[r] = {}

        if n not in fpDict[r].keys():
            fpDict[r][n] = {}
            rfDist[r][n] = {}

        fpDict[r][n] = d
        rfDist[r][n] = {x.split("=")[0]: x.split("=")[1] for x in dist}
    
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

# def getClosest(smi, fpfile, th=0.8, fp=None, fpn=None, fpp=None, fpfun=None, marvin=False):
#     dist = {}
#     if fp is None:
#         print('Reading fingerprints')
#         data = np.load(fpfile)
#         fp = data['x'] 
#         fpn = data['y']
#         fpp = 8
#         fpfun = GetMorganFingerprint
#         data.close()

#     targetMol = Chem.MolFromSmiles(smi)
#     # If RDkit fails, we sanitize first using molconvert from ChemAxon, which is more robust
#     if targetMol is None and marvin:
#         try:
#             cmd = ['molconvert', 'mol', smi]
#             cmd2 = ['molconvert', 'smiles']
#             job = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#             job2 = subprocess.Popen(cmd2, stdin=job.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#             job.stdout.close()
#             out, err = job2.communicate()
#             targetMol = Chem.MolFromSmiles(out)
#         except:
#             pass
#     if fpp is not None:

#         info1 = {}
#         targetFp = AllChem.GetMorganFingerprint(targetMol, 8, bitInfo=info1, invariants=AllChem.GetConnectivityInvariants(targetMol, includeRingMembership=False))
#     else:
#         info1 = {}
#         targetFp = AllChem.GetMorganFingerprint(targetMol, 8, bitInfo=info1, invariants=AllChem.GetConnectivityInvariants(targetMol, includeRingMembership=False))
        
#     tn = DataStructs.BulkTanimotoSimilarity(targetFp, list(fp))
#     for i in sorted(range(0, len(tn))):
#         dist[fpn[i]] = tn[i]
#     return dist, fp, fpn


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
    rsp = {}
    for l in open(dbfile):
        if l.startswith('#'):
            continue
        m = l.rstrip().split('\t')
        rid = m[0]
        subs = {}
        prods = {}
        m = l.rstrip().split('\t')
        left, right = m[1].split(' = ')
        subs = getReactants(left)
        prods = getReactants(right)
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
        pairings = set()
        simm[p] = {}
        compPairs[p]=[]

        for x in cl[p[0]]:
            simm[p][x] = (0.0, x, None)
            if x in sim:
                for y in cl[p[1]]:
                    if y in sim[x]:
                        pairings.add( (sim[x][y], x, y) )

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

    return(S1, S2, compPairs)


def bulkTani(targetFp, fp, fpn):
    tn = DataStructs.BulkTanimotoSimilarity(targetFp, list(fp))
    dist={}
    for i in sorted(range(0, len(tn))):
        dist[fpn[i]] = tn[i]
    return dist


def run(arg, pc, rxn_mapper):

    ### Data processing 
    # read in the data
    if arg.out:
        fileObj = open(arg.out, 'w')    
    if arg.high:
        fileObj = open(arg.high, 'w')
    rsp = reacSubsProds(os.path.join(arg.datadir, 'reac_prop.tsv'))

    # read in the query reaction
    smiles = ''
    if arg.rxn is not None:
        rTarget, smiles = getReaction(arg.rxn)
    elif arg.smarts is not None:
        rxnfile = os.path.join(os.path.dirname(arg.out), 'reaction.rxn')
        try:
            rTarget, smiles = getReactionFromSmiles(arg.smarts, rxnfile)
        except:
            print('smile could not be processed')
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
      
    # Read fingerprint info from preload data if available
    if pc is not None: 
        if arg.fp in pc.fp:
            fp, fpn, fpp, fpfun = pc.fp[arg.fp]
            fpRF, fpnRF, fprRF, rfDict, rfdist = pc.fpr[arg.fp]
        else:
            fp, fpn, fpp, fpfun = loadFingerprint(arg.datadir, arg.fp)
            fpRF, fpnRF, fprRF, rfDict, rfdist = loadFingerprintRF(arg.datadir, 'FP_MorgRF')     
    else:      
        print('reading fp data')
        fp, fpn, fpp, fpfun = loadFingerprint(arg.datadir, arg.fp)
        fpRF, fpnRF, fprRF, rfDict, rfdist = loadFingerprintRF(arg.datadir, 'FP_MorgRF')
  
    ### Compound Calculations
    reacting_atoms = ReactingAtoms()
    reacting_fragments = ReactingFragments()
    queryRF = {}
    queryDists = {}
    sim = {}
    
    # Break the input query into chemical components then calculate
    #     1. the reacting atoms for each compound (reactingAtoms)
    #     2. the distance between compounds and the compounds in the database (sim)
    #     3. the reacting fragments and  distances (queryRF, queryDists)    
    for r in rTarget:        
        print('\nAAM of the Query reaction') 
        subSmiles = rTarget[r][0].keys()
        prodSmiles = rTarget[r][1].keys()
        if set(subSmiles) == set(prodSmiles):
            print('no transformation')
            
        process_compounds = Compounds()
        subMols = [process_compounds.process_compF(k, 8) +(k, v,) for k, v in rTarget[r][0].items()]
        prodMols = [process_compounds.process_compF(k, 8) +(k, v,) for k, v in rTarget[r][1].items()]
    
        # AAM using RXNmapper 
        reactingAtoms, conf = reacting_atoms.get_reacting_atoms(subMols, prodMols, subSmiles, prodSmiles, rxn_mapper)  
 
        # measure compound similarities
        print('\nCompound-compound distances')
        for comp in subMols + prodMols:
            smile = comp[5]
            # measure distances for whole compounds - equivalent to getClosest
            sim[smile] = bulkTani(comp[2], fp, fpn) 
            
            # use reacting atoms to generate reacting fragments
            if conf>0 and smile in reactingAtoms:
                # rfs are fragNumber:(startAtom, dist), dists are fragNumber: distFromRA, startAtom, RA                
                rfs, dists = reacting_fragments.get_distances(comp, reactingAtoms[smile])  
                
                if len(rfs)>0:
                    queryRF[smile]=rfs
                    queryDists[smile]=dists  
                    # draw_rfs(comp[0], comp[2], comp[3], rfs)

                           
    ### Reaction similarities
    for r1 in rTarget:
        s1, p1 = rTarget[r1]      
        for r2 in rsp: 
            ### whole compound similarity                        
            s2, p2, ec2 = rsp[r2]   
            if s2 == p2: continue
            S1, S2, subProdPairs = getRSim(s1, p1, s2, p2, sim)
            if S1 > 0 and S2 > 0:
            
                # if aam failed, record scores replacing these values with Nan
                if conf==0 or r2 not in rfDict.keys():
                    if arg.out:
                       print(r1, r2, S1, S2, smiles,  float("Nan"), ec2, file = fileObj)
                       
                    if arg.high:
                        print(r1, r2, max([S1, S2]),  smiles,  float("Nan"), ec2, file = fileObj) 
                    continue
                

                ### RF similarity
                RF_score=float("Nan")

                subProdPairs = {k:[x for x in v if x[0] in queryRF.keys()]  for k, v in subProdPairs.items()}
                try:
                    # calculate the scores for the RFs
                    # forward score
                    if S1 >= S2:
                        S_RF = reacting_fragments.generate_RFscore2(list(s1.keys()), queryRF, queryDists, s2, r2, rfDict, rfdist, subProdPairs[('s1', 's2')] )
                        P_RF = reacting_fragments.generate_RFscore2(list(p1.keys()), queryRF, queryDists, p2, r2, rfDict, rfdist, subProdPairs[('p1', 'p2')] )                       
                        if (sum(S_RF.values())+sum(P_RF.values()))>0: 
                            RF_score = (sum(S_RF.values())+sum(P_RF.values())) / len(queryRF) #(len(S_RF) + len(P_RF))
                        else: RF_score=0
                    else:
                        # reverse score
                        S_RF = reacting_fragments.generate_RFscore2(list(s1.keys()), queryRF, queryDists, p2, r2, rfDict, rfdist, subProdPairs[('s1', 'p2')])
                        P_RF = reacting_fragments.generate_RFscore2(list(p1.keys()), queryRF, queryDists, s2, r2, rfDict, rfdist, subProdPairs[('p1', 's2')])
                        if (sum(S_RF.values())+sum(P_RF.values()))>0: 
                            # geometric mean 
                            RF_score = math.sqrt(mean(S_RF.values()) * mean(P_RF.values()))
                        else: RF_score=0
                except:
                    RF_score = float("Nan")
                
                if arg.out:
                    print(r1, r2, S1, S2, smiles, RF_score, ec2,   file = fileObj) 

                if arg.high:
                    if S1 >= S2:
                        print(r1, r2, S1, smiles,  RF_score, ec2,  file = fileObj)
                    else:
                        print(r1, r2, S2, smiles,  RF_score, ec2,   file = fileObj)                

 
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
    # arg = arguments()

    smi = 'O=C[C@H](O)COP(=O)([O-])[O-]>>O=C[C@H](O)[C@H](O)COP(=O)([O-])[O-]'
    arg = arguments(['/home/ruth/code/update_selenzyme/run_folder_min_dec/data/', 
        'Morgan',
        '-smarts', smi, 
        '-out', '/home/ruth/code/update_selenzyme/RSquickRsim_new.txt',
        '-frag_size', '1'] )
    pc=None


    rxn_mapper = RXNMapper()   
    run(arg, pc, rxn_mapper)