#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  8 10:13:22 2021

@author: ruth

Get the reacting fragments for all the reactions in reac_prop

reaction data   - reac_prop.tsv
smiles          - chem_prop.tsv

make into fingerprints using updated parameters
save as mgfp5.npz 
    where the type is numpy.lib.npyio.NpzFile, 
    x is UIntSparseIntVect objects and y is a numpy array of MNXM Id 

get the reacting fragments 
save as mgfp5_rf.npz


"""

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import os
from pathlib import Path
from rdkit.Chem import Draw
from rxnmapper import RXNMapper
import argparse


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
                i=1
                while i < rad:
                    # get the nodes for the next iteration
                    adjNodes2 = set()
                    for x in adjNodes:
                        adjNodes2 = adjNodes2 | atomMap[x]
                        scope = scope | adjNodes2
                    adjNodes = adjNodes2
                    i+=1
                fragAtoms[str(bitName)+'_'+ str(startNode)] = scope   
    return fragAtoms

def mergeDicts(d2, d1):
    for k, v in d1.items():
        if k not in d2.keys():
            d2[k]=[]
        d2[k].append(v)
    return(d2)  
    
def reactingFragDists(subDistReact, hitFrags1, fragAtoms1, bi1):
    # find the closest instances and measure distances
    fragDist ={}
    bi1_filtered = {}
    for ra, dists in subDistReact.items():
        # get the closest instance for each fragment, by finding the smallest radius
        fd={}
        bi={}
        
        for frag in hitFrags1:          
            # get the fragment instance details
            fragName, startAtom = list(map(int, frag.split('_')))
                  
            # get the fragment atoms
            frag_atoms = fragAtoms1[frag]
            if ra not in frag_atoms: continue
            
            # look for the furthest atom for the reacting atoms
            for d, a in dists:           
                if a in frag_atoms: 
                    # see if its the smallest radius for the reacting atom
                    if fragName in fd.keys() and d> fd[fragName][0]: 
                        break
                    else:
                        fd[fragName]=[d, startAtom, ra]
                        
                        # bi dictionary contains all instances, filter for this o
                        bi_instance = [x for x in bi1[fragName] if x[0] == startAtom]
                        bi[fragName]= bi_instance[0]    
                        break

	# merge lists while maintaining list structure of values
        fragDist = mergeDicts(fragDist, fd)
        bi1_filtered = mergeDicts(bi1_filtered, bi)

    return(bi1_filtered, fragDist)

def reactFragDists(comp1, reactAtoms):
    ### set up the requirements  
    molN = comp1[0]
    atomMap1 = comp1[1]
    fp1 = comp1[2]
    bi1 = comp1[3]
     
    ### Get the distances between the reacting atoms and the furthest atom in each RF
    distM = Chem.rdmolops.GetDistanceMatrix(molN)
    distReact = {ra:sorted([[x, i] for i, x in enumerate(distM[ra])], reverse=True) for ra in reactAtoms }
    
    # get fragments that contain the reacting atom
    fragAtoms1 = getAtomFragments(fp1, bi1, atomMap1) 
    hitFrags1 = set([k for k, v in fragAtoms1.items() if len(v.intersection(reactAtoms))>0])
    
    ### Filter frags to those containsing reacting atoms, then measure the distance to the furthest atoms
    bi1_filtered, fragDist1 = reactingFragDists(distReact, hitFrags1, fragAtoms1, bi1)
    bi1_filtered = {k:tuple(set(v)) for k, v in bi1_filtered.items()}
    
    return(bi1_filtered, fragDist1)




def clean_up_data(sub_ReactingFragsM):
    sub_RF = [[k, v] for k, v in sub_ReactingFragsM.items() if type(v) !=int and type(v) !=np.int64 ]
    
    sub_RF2=[]
    for k, v in sub_RF:
        if type(v) == pd.core.series.Series:
            for x in v:
                sub_RF2.append([k,x])
        else:
            sub_RF2.append([k,v])
    return(sub_RF2)








def rxnMapper_fun(subsmiles, prodsmiles, subs_fp, prods_fp, rxn_mapper):   
    
    ### generate a reaction smile and map it - this will rearrange the order of the compounds and ordewr of the atoms
    react_smile = '.'.join(subsmiles) +'>>'+'.'.join(prodsmiles)
        
    # get_attention_guided_atom_maps_for_reactions
    results = rxn_mapper.get_attention_guided_atom_maps([react_smile])[0]  
    smileM, conf = results['mapped_rxn'], results['confidence'] 
    rxn1 = AllChem.ReactionFromSmarts(smileM, useSmiles=True)
    rxn2 = AllChem.ReactionFromSmarts(smileM.split('>>')[1] + '>>' + smileM.split('>>')[0], useSmiles=True)

    ### get the reacting atoms within the mapped reactions 
    rxn1.Initialize()
    rxn2.Initialize()
    react_atoms_subs = rxn1.GetReactingAtoms(mappedAtomsOnly=True)     
    react_atoms_prods = rxn2.GetReactingAtoms(mappedAtomsOnly=True) 

    ### get the reacting fragments - substrates
    reacting_fragments = {}
    # figlist=[]
    for i, m in enumerate(rxn1.GetReactants()):
        # rxnmapper rearranges the atoms/compound, so they need to be fingerprinted again
        fp, smile, size, info,atomMap1 = get_morg(m, 8)
        # get the reacting fragments for the new fingerprints 
        rfs, dists = reactFragDists([m, atomMap1, fp, info], react_atoms_subs[i])
        
   
        
        # map to the input compounds by matching their fingerprints
        for k, v in subs_fp.items():
            if v == fp:
                reacting_fragments[k] = [rfs, dists] 


    if len(reacting_fragments) != len(subs_fp):
        raise Exception("MappingFailure")    

    ### get the reacting fragments - products
    for i, m in enumerate(rxn2.GetReactants()):
        fp, smile, size, info,atomMap1 = get_morg(m, 8)
        rfs, dists = reactFragDists([m, atomMap1, fp, info], react_atoms_prods[i])

        for k, v in prods_fp.items():
            if v == fp:
                reacting_fragments[k] = [rfs, dists] 
                
    if len(reacting_fragments) != len(subs_fp) + len(prods_fp):
        raise Exception("MappingFailure")
        
    return reacting_fragments, conf, react_smile


def get_morg(mol, inchi=0):   

    Chem.SanitizeMol(mol)
    info1 = {}
    fpM = AllChem.GetMorganFingerprint(mol, 8, bitInfo=info1, invariants=AllChem.GetConnectivityInvariants(mol, includeRingMembership=False))
    atomMap1 = get_atoms(mol)
    return fpM, Chem.MolToSmiles(mol), mol.GetNumAtoms(), info1, atomMap1

def get_atoms(mol):
    atomMap = {}
    # loop through the bonds to get each atoms adjacent atoms
    for b in mol.GetBonds():
        start = b.GetBeginAtomIdx()
        end = b.GetEndAtomIdx()
        
        if start not in atomMap.keys():
            atomMap[start] = set()
        atomMap[start].add(end)
 
        if end not in atomMap.keys():
            atomMap[end] = set()
        atomMap[end].add(start)     
    return atomMap

def get_inchi(smiles):
    i= Chem.MolToInchi(Chem.MolFromSmiles(smiles)).split('/')
    return '/'.join(i[0: min(6, len(i)-1)])


def run(raw_data_folder, data_folder):

    reac_prop = pd.read_csv(raw_data_folder / 'reac_prop.tsv', skiprows=351, sep='\t')
    chem_prop = pd.read_csv(raw_data_folder / 'chem_prop.tsv', skiprows=351, sep='\t')
    filter_reactions = pd.read_csv(raw_data_folder / 'reaction_smiles_enz_filter.tsv', sep='\t', header=None)
    compounds_in_reactions = set([y for x in filter_reactions[1] for y in str(x).split(',')])
    filter_reactions = set(filter_reactions[0])


    #### Get fingerprints for chemicals in the reactions file
    mnxmCovered=set()
    FingerprintsM=[]
    fpd = {}

    MNXM=[]
    not_in_reactions = set()
    fail = dict()
    comp_smiles = dict()
    comp_size = dict()

    morgan_lost = 0
    for n, row in chem_prop.iterrows(): 
        if row['#ID'] not in compounds_in_reactions: 
            continue

        if row.SMILES != row.SMILES: continue
        smiles = row.SMILES.split('.')
        inchi = row.InChI

        
        # if the smile contains multiple compounds then the fingerprint is the sum 
        try:
            fps = []
            smiles_store = []
            sizes=[]
            for smile in smiles:
                fp, smile, size, info, atomMap = get_morg(Chem.MolFromSmiles(smile))
                fps.append(fp)
                smiles_store.append(smile) 
                sizes.append(size)
        except:
            try:
                fps = []
                smiles_store = []
                
                if len(smiles) ==1:
                    fp, smile, size, info, atomMap = get_morg(Chem.MolFromInchi(inchi), 1)
                    fps.append(fp)
                    smiles_store.append(smile)
                    sizes.append(size)
                else:
                    smiles_store = Chem.MolToSmiles(Chem.MolFromInchi(inchi)).split('.')
                    for smile in smiles:
                        fp, smile, size, info, atomMap = get_morg(Chem.MolFromSmiles(smile))
                        fps.append(fp)  
                        smiles_store.append(smile)
                        sizes.append(size)
            except:
                fail[row['#ID']] = smile

        if len(fps)==0:
            continue
        
        fp = fps[0]
        if len(fps)>1:
            for i in list(range(1, len(fps))):
                fp = fp + fps[i]

        FingerprintsM.append(fp)    
        MNXM.append(row['#ID'])

        comp_smiles[row['#ID']] = '.'.join(smiles)
        comp_size[row['#ID']] = sum(sizes)

    print('\ncompounds', len(MNXM), 'out of', len(compounds_in_reactions), 'fail', len(fail))


    fpd = dict(zip(MNXM, FingerprintsM))

    ### Get the reaction fragments

    reac_prop = reac_prop[reac_prop['#ID'].isin(filter_reactions)].reset_index()


    rxn_mapper = RXNMapper()

    # reaction_smiles = {}
    aam_issues = {'tooBig':[], 'starSmiles' :[], 'unknown' : [], 'mappingFailure': []}
    reaction_issues = {'same_sub_prod' :set(), 'emptyReactions': set(), 'emptyReactions_fp': set(), 'emptyReactions_stars': set(),  'missingRFs': {}}
    compound_issues = {}

    MNXM_RF = []
    MNXR_RF = []
    FP_react = []
    Dists = []

    # get the chemical components from reac_prop and reconstruct the smile compounds
    for rowNo, row in reac_prop.iterrows():

        if row.is_transport == 'T': 
            reaction_issues['same_sub_prod'].add(row['#ID'])
            continue

        # read in the row data
        subs, prods = row.mnx_equation.split(' = ')
        if len(subs)==0 or len(prods) ==0: 
            reaction_issues['emptyReactions'].add(row['#ID'])
            continue
        
        subs_count = {x1.split(' ')[1].split('@')[0] : int(x1.split(' ')[0]) for x1 in subs.split(' + ') }
        prods_count ={x1.split(' ')[1].split('@')[0] : int(x1.split(' ')[0]) for x1 in prods.split(' + ') }
        reaction = row['#ID']
     
        # only process compounds with fingerprints
        subs =  set([x for x in subs_count.keys()  if x in MNXM and x in comp_smiles.keys()])
        prods = set([x for x in prods_count.keys()  if x in MNXM and x in comp_smiles.keys()])

        if len(subs)==0 or len(prods) ==0: 
            reaction_issues['emptyReactions_fp'].add(reaction)
            continue
        
        # get the smiles for reac_smi 
        # reaction_smiles[reaction] = '.'.join([comp_smiles[x] for x in subs]) + '>>' + '.'.join([comp_smiles[x] for x in prods])

        # the AAM can't process smiles with stars
        subs = set([x for x in subs if '*' not in comp_smiles[x]])
        prods = set([x for x in prods if '*' not in comp_smiles[x]])    

        
        # check that we have viable substrates + products
        if len(subs)==0 or len(prods)==0:
            reaction_issues['emptyReactions_stars'].add(reaction)
            continue


       
        subsmiles =  [x1 for x in subs for x1 in  [comp_smiles[x]]*subs_count[x]]
        prodsmiles = [x1 for x in prods for x1 in [comp_smiles[x]]*prods_count[x]] 
        
        subs_inchi =  set([ get_inchi(comp_smiles[x]) for x in subs])
        prods_inchi = set([  get_inchi(comp_smiles[x]) for x in prods])
        if set(subs) == set(prods) or subsmiles == prodsmiles or subs_inchi == prods_inchi:
            reaction_issues['same_sub_prod'].add(reaction)
            continue
            
        
        ### AAM
        try:
            # for unbalances reactions there need to be more atoms on the substrate side
            if sum([comp_size[x]* subs_count[x] for x in subs]) >= sum([comp_size[x]* prods_count[x] for x in prods]):
                s = '.'.join(subsmiles) +'>>'+'.'.join(prodsmiles)
                reactingAtoms, conf, react_smile = rxnMapper_fun(subsmiles, prodsmiles, {x: fpd[x] for x in subs}, {x: fpd[x] for x in prods}, rxn_mapper)            
            else:
                s = '.'.join(prodsmiles) +'>>'+'.'.join(subsmiles)
                reactingAtoms, conf, react_smile = rxnMapper_fun( prodsmiles, subsmiles, {x: fpd[x] for x in prods},  {x: fpd[x] for x in subs}, rxn_mapper)  
          
        except RuntimeError: 
            aam_issues['tooBig'].append([reaction, rowNo,  s])
            continue
        except ValueError: 
            aam_issues['starSmiles'].append([reaction, rowNo, s])
            print('star smile',rowNo)
            continue
        except Exception:
            aam_issues['mappingFailure'].append([reaction, rowNo, s])
            continue       
        except: 
            print('other error', rowNo)
            aam_issues['unknown'].append([reaction, s])
            continue

        # check that molecules arent missing rfs 
        lostRAs = set([k for k,v in reactingAtoms.items() if len(v[0])==0])
        single_atom = set([x for x in lostRAs if comp_size[x] == 1])
        if lostRAs:
            if len(lostRAs)> len(single_atom):
                for x in lostRAs - single_atom:
                    if x not in compound_issues:
                        compound_issues[x] = set([reaction])
                    else:
                        compound_issues[x].add(reaction)
                        
            if subs.intersection(lostRAs) == subs or prods.intersection(lostRAs) == prods:
                if reaction not in reaction_issues['missingRFs']:
                    reaction_issues['missingRFs'][reaction] = []
                    
                if subs.intersection(lostRAs) == subs:
                    reaction_issues['missingRFs'][reaction].append(subs)
                if prods.intersection(lostRAs) == prods:
                    reaction_issues['missingRFs'][reaction].append(prods)                           
                # print('missingRFS', lostRAs - single_atom, rowNo)
                continue


         ### get the reacting fragments for every molecule
        for comp in subs | prods:

            ### get the reacting fragments into a sparse int vector
            # get the reacting frags in a list
            rfs, dists = reactingAtoms[comp]
            rfList1 = [x for k, v in rfs.items() for x in [k]*len(v) ]
            if len(rfList1) == 0:
                continue
            
            # make an empty sparse int vector
            SparseIntVect1 = AllChem.GetMorganFingerprint(Chem.MolFromSmiles(''), 8)
            # update the empty sparse int vector with the reacting fragments
            SparseIntVect1.UpdateFromSequence(rfList1)
            
            # prepare list of distancces
            distList = []
            for k in rfs.keys():
                v = dists[k]
                distStr = str(k) + '=' + '|'.join([str(x[2]) +'_' + str(int(x[0])) for x in v])
                distList.append( distStr)
            
            # save the data
            MNXM_RF.append(comp)
            MNXR_RF.append(reaction)
            FP_react.append(SparseIntVect1)
            Dists.append(distList)


    total_reactions = len(set(reac_prop['#ID']))

    print('\n\nsucessful reactions', len(set(MNXR_RF)), 'sucessful compounds', len(set(MNXM_RF)))
    print('\ntotal reactions', total_reactions, 'missing reactions', total_reactions - len(set(MNXR_RF)), 
          '\t', round(  (len(set(MNXR_RF)) / total_reactions)*100, 3), '%')

    print('\nreaction issues', sum([len(x) for x in  reaction_issues.values()]), '\t', round( (sum([len(x) for x in  reaction_issues.values()])/total_reactions)*100 ,3) , '%'  )
    for k, v in reaction_issues.items(): print(k, '\t', len(v), '\t', round( (len(v)/total_reactions)*100 ,3) , '%' )

    print('\nmapping issues', sum([len(x) for x in  aam_issues.values()]), '\t', round( ( sum([len(x) for x in  aam_issues.values()]) /total_reactions)*100 ,3) , '%' )
    for k, v in aam_issues.items(): print(k,  '\t',len(v), '\t', round( (len(v)/total_reactions)*100 ,3) , '%' )

    # save to npz file - full compounds
    outfolderM = data_folder / 'Morgan/'
    outfolderM_RF = data_folder / 'Morgan/RF/'

    if not os.path.exists(outfolderM): os.makedirs(outfolderM)
    if not os.path.exists(outfolderM_RF): os.makedirs(outfolderM_RF)


    # save to npz file 
    #  Morgan data
    np.savez_compressed(outfolderM / 'FP_Morg.npz', x=FingerprintsM , y=MNXM)
    np.savez_compressed(outfolderM / 'RF/FP_MorgRF.npz', x=FP_react, y=MNXM_RF, z=MNXR_RF, d=Dists )



    # save the reac_smi.csv 
    # with open(data_folder / 'reac_smi.csv', 'w') as f:
    #     f.write('RID,SMILES')
    #     for k, v in  reaction_smiles.items():
    #         f.write('\n' + k+','+v)

    # save cleaned reac_prop 
    reac_prop2 = reac_prop.drop('index', axis=1)
    reac_prop2.mnx_equation = [' '.join([y.split('@')[0] for y in x.split()]) for x in reac_prop.mnx_equation]
    reac_prop2.to_csv(data_folder / 'reac_prop.tsv', sep='\t', header=None, index=False)


def arguments(args=None):
    parser = argparse.ArgumentParser(description='Generate fingerprints')
    parser.add_argument('data_folder', 
                        help='specify data directory for new files, please end with slash')
    parser.add_argument('raw_data_folder',
                        help='specify data directory for raw databases files, please end with slash')

    arg = parser.parse_args(args=args)
    return arg


if __name__ == '__main__':
    arg = arguments()
    
    arg = arguments(['/home/ruth/code/update_selenzyme/selenzyme_aug2023/data_2023_t/',
        '/home/ruth/code/update_selenzyme/run_folder_min_dec/raw_data_update/'] )
    
    
    
    raw_data_folder = Path(arg.raw_data_folder)
    data_folder = Path(arg.data_folder)

    run(raw_data_folder, data_folder)

