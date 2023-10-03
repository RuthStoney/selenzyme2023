#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  3 12:56:49 2022

@author: Ruth Stoney
"""
from rdkit import Chem
from rdkit.Chem import AllChem
import re
from rxnmapper import RXNMapper
import math


class Compounds:

    # mol_with_atom_index 
    def _index_atoms(self, mol):
        atoms = mol.GetNumAtoms()
        for idx in range(atoms):
            mol.GetAtomWithIdx(idx).SetProp('molAtomMapNumber', str(mol.GetAtomWithIdx(idx).GetIdx()))
        return mol

    # _getAtoms 
    def _map_atoms(self, comp1):
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

    def _get_fragment_atoms(self, fp1, info1, atomMap):
        # get all the atoms in each fragment
        fragAtoms = {}
    
        for bitName in fp1.GetNonzeroElements().keys():
            bitInfo = info1[bitName]
    
            for x in bitInfo:
                startNode, rad = x
                scope = set([startNode])
                # if radius is zero, the start node is the only node in the scope
                if rad == 0:
                    fragAtoms[str(bitName) + '_' + str(startNode)] = scope
                    
                # if not iterate through atomMap to get the other nodes
                else:
                    adjNodes = atomMap[startNode]
                    scope = scope | adjNodes
                    i = 1
                    while i < rad:
                        # get the nodes for the next iteration
                        adjNodes2 = set()
                        for x in adjNodes:
                            adjNodes2 = adjNodes2 | atomMap[x]
                            scope = scope | adjNodes2
                        adjNodes = adjNodes2
                        i += 1
                    fragAtoms[str(bitName) + '_' + str(startNode)] = scope
        return fragAtoms
    
    def process_compF(self, comp1, size):
        # get all the details for the fingerprints
        
        if type(comp1) == str:
            try:
                comp1 = Chem.MolFromSmiles(comp1)
            except:
                print('failed MolFromSmiles')
                
        try:
            Chem.SanitizeMol(comp1)
            comp1 = self._index_atoms(comp1)
    
            # get each atoms adjacent neighbours
            atomMap1 = self._map_atoms(comp1)
            
            # Morgan fingerprints
            info1 = {}
            fpM1 = AllChem.GetMorganFingerprint(
                comp1, size, bitInfo=info1,
                invariants=AllChem.GetConnectivityInvariants(comp1, includeRingMembership=False))
            
            # get all the atoms in each fragment, returns {fragname:[{atoms}, startnode, radius]}   
            fragAtoms1 = self._get_fragment_atoms(fpM1, info1, atomMap1)
        
            return(comp1, atomMap1, fpM1, info1, fragAtoms1)
        except:
            print('unable to process compound')    


class ReactingAtoms:

    def get_reacting_atoms(self, subMols, prodMols, subSmiles, prodSmiles, rxn_mapper):
        # AAM using RXNmapper
        try:
            # for unbalances reactions there need to be more atoms on the substrate side
            if sum([x[0].GetNumAtoms() for x in subMols]) >= sum([x[0].GetNumAtoms() for x in prodMols]):
                smile = '.'.join(subSmiles) + '>>' + '.'.join(prodSmiles)
                reactingAtoms, conf = self._run_rxnMapper(smile, subMols, prodMols, rxn_mapper)
            else:
                smile = '.'.join(prodSmiles) + '>>' + '.'.join(subSmiles)
                reactingAtoms, conf = self._run_rxnMapper(smile, prodMols, subMols, rxn_mapper)
            
            return reactingAtoms, conf
    
        except RuntimeError:
            print('AAM failed query - compound too large')
            return ([None, 0])
    
        except ValueError:
            print('AAM failed query - issue with input smiles', smile)
            if '*' in smile:
                print('please use smiles without *s')
            return ([None, 0])
    
        except:
            print('AAM failed query')
            return ([None, 0])
        
    # _rxnMapper_fun
    def _run_rxnMapper(self, smile, subMols, prodMols, rxn_mapper):
        print('in RXNMapper_fun processing', smile)
      
        # generate a SMILE with atoms mapped between substrate and products
        results = rxn_mapper.get_attention_guided_atom_maps([smile])[0]
        smileM, conf = results['mapped_rxn'], results['confidence']
    
        # get the reacting atoms 
        rxn1 = AllChem.ReactionFromSmarts(smileM, useSmiles=True)
        rxn2 = AllChem.ReactionFromSmarts(smileM.split('>>')[1] + '>>' + smileM.split('>>')[0], useSmiles=True)
        rxn = AllChem.ReactionFromSmarts(smile, useSmiles=True)
        rxn1.Initialize()
        rxn2.Initialize()
        react_atoms_subs = rxn1.GetReactingAtoms(mappedAtomsOnly=True)  
        react_atoms_prods = rxn2.GetReactingAtoms(mappedAtomsOnly=True)
        
        # map the atom numbers from the RXNMapper index to the RDkit index 
        subReact, subMap = self._get_RAs(rxn.GetReactants(), rxn1.GetReactants(), react_atoms_subs)
        prodReact, prodMap = self._get_RAs(rxn.GetProducts(), rxn2.GetReactants(), react_atoms_prods)
     
        # make string: reacting atoms dict
        if len(re.split('\\.|>>', smile)) != len(subMap + prodMap):
            print('mapping error')
        
        reactAtoms = {}
        for i, x in enumerate(smile.split('>>')[0].split('.')):
            if x not in reactAtoms:
                reactAtoms[x] = set()
            reactAtoms[x] = reactAtoms[x] | set(subReact[i])
     
        for i, x in enumerate(smile.split('>>')[1].split('.')):
            if x not in reactAtoms:
                reactAtoms[x] = set()
            reactAtoms[x] = reactAtoms[x] | set(prodReact[i])
        
        print('reacting atoms:', reactAtoms.values())
        return(reactAtoms, conf)
    
    def _get_RAs(self, rxn, rxn1, react_atoms_subs):
        # map the compounds and atoms (RXNMapper can change the order of the compounds)
        subMap = []
        subReact = []
    
        # get the inchis - not matching compounds can match inchi
        for x in rxn:
            x.SetProp('inchi', str(Chem.inchi.MolToInchi(x)))
        for x in rxn1:
            x.SetProp('inchi', str(Chem.inchi.MolToInchi(x)))
    
        already_mapped = []
        for i1, comp1 in enumerate(rxn):
            Chem.SanitizeMol(comp1)
            comp1 = Chem.RemoveHs(comp1)
            # check if the unmapped compound got paired with a mapped compound
            # show_atom_number(comp1, 'atomLabel')
            mpd = 0
            for i2, comp2 in enumerate(rxn1):
                if i2 in already_mapped:
                    continue
                comp2 = Chem.RemoveHs(comp2)
                
                if comp1.GetProp('inchi') == comp2.GetProp('inchi'):
                    already_mapped.append(i2)
    
                    # get the reacting atoms
                    react_a = react_atoms_subs[i2]
                    
                    # map the atom indices of the mapped onto the unmapped
                    a_map = comp1.GetSubstructMatch(comp2)
                    if comp1.GetNumAtoms() - 1 < max(a_map) or comp2.GetNumAtoms() - 1 < max(a_map):
                        print('mapping failed')
                    subMap.append(a_map)
                    
                    # get the reacting atoms
                    subReact.append([a_map[x] for x in react_a])
                    # mark this compound as taken      
                    mpd = 1
                    break
                
            # if the compound wasnt mapped, add an empty list
            if mpd == 0:
                print('e')
                subMap.append([])
                subReact.append([])
                
        return(subReact, subMap) 
    

class ReactingFragments:
    # reactFragDists
    def get_distances(self, comp1, reactAtoms):
        # set up the requeirements
        molN = comp1[0]
        atomMap1 = comp1[1]
        fp1 = comp1[2]
        bi1 = comp1[3]
         
        # Get the distances between the reacting atoms and the furthest atom in each RF
        distM = Chem.rdmolops.GetDistanceMatrix(molN)
        distReact = {ra: sorted([[x, i] for i, x in enumerate(distM[ra])], reverse=True) 
                     for ra in reactAtoms}
        
        # get fragments that contain the reacting atom
        fragAtoms1 = self._get_fragment_atoms(fp1, bi1, atomMap1)
        hitFrags1 = set([k for k, v in fragAtoms1.items() if len(v.intersection(reactAtoms)) > 0])
        
        # Filter frags to those containing reacting atoms, then measure the distance to the furthest atoms
        bi1_filtered, fragDist1 = self._calculate_distances(distReact, hitFrags1, fragAtoms1, bi1)
        bi1_filtered = {k: tuple(set(v)) for k, v in bi1_filtered.items()}
        
        return(bi1_filtered, fragDist1)
   

    def generate_RFscore2(self, subSmiles, queryRF, queryDists, s2, r2, rfDict, rfdist, subProdPairs):  
        # subSmiles, queryRF, queryDists, s2, r2, rfDict, rfdist, subProdPairs = list(s1.keys()), queryRF, queryDists, s2, r2, rfDict, rfdist, subProdPairs[('s1', 's2')] 
        
        subList = list(subSmiles)
        
        querySubsRF = {x: queryRF[x] 
                       for x in subList if x in queryRF.keys() and len(queryRF[x]) > 0}
        querySubsDist = {x: queryDists[x] 
                         for x in subList if x in queryRF.keys() and len(queryDists[x]) > 0}
        
        dbList = list(s2.keys())
        dbSubsRF = {x: rfDict[r2][x] for x in dbList if x in rfDict[r2]}
        dbSubsDist = {x: rfdist[r2][x] for x in dbList if x in rfDict[r2]}
        
        spPairs = ["_".join([x[0], x[1]]) for x in subProdPairs]
        unweightedScores = {}
        
        for sb1, qSubRF in querySubsRF.items():
            qSubsDist = querySubsDist[sb1]
            qSortDists = {}
            
            # iterate through the query fragments
            for k, v in qSubsDist.items():
                dists1 = [int(x[0]) for x in v]
                for d in dists1:
                    if d not in qSortDists.keys():
                        qSortDists[d] = []
                    qSortDists[d].append(k)
     
            # iterate through the db fragments
            for sb2, dSubRF in dbSubsRF.items():
                if dSubRF == {}: 
                    print('empty', sb2)
                    continue
                if '_'.join([sb1, sb2]) not in spPairs:
                    continue
    
                dSubsDist = dbSubsDist[sb2]
                dSortDists = {}
                for fragNo, k in enumerate(dSubRF.GetNonzeroElements().keys()):
                    dists2 = [int(y.split('_')[1]) for y in dSubsDist[str(k)].split('|')]
                    for d in dists2:            
                        if d not in dSortDists.keys():
                            dSortDists[d] = []
                        dSortDists[d].append(k)
               
                # get jaccard sim
                unweightedScores[sb1, sb2] = {}
                
                for fragDist in set(qSortDists.keys()).union(set(qSortDists.keys())):
                    if fragDist not in set(qSortDists).intersection(set(dSortDists)):
                        unweightedScores[sb1, sb2][fragDist] = 0
                        
                    else:
                        intersect = {frag: min([dSortDists[fragDist].count(frag), qSortDists[fragDist].count(frag)]) for frag in set(dSortDists[fragDist] + qSortDists[fragDist]) if frag in dSortDists[fragDist] and frag in qSortDists[fragDist]}
                        intersectScore = sum(intersect.values())
                        unweightedScores[sb1, sb2][fragDist] = intersectScore / len(qSortDists[fragDist])
           
        weightedScores = {}
        for pair, scores in unweightedScores.items():
            weights = self._log_fun_weight(max(scores.keys()))
            weightedScores[pair] = sum([scores[x] * weights[x] 
                                        if x in scores and scores[x] > 0 else 0 
                                        for x in list(range(0, max(scores.keys()) + 1))])
    
        return(weightedScores)

    # getAtomFragments
    def _get_fragment_atoms(self, fp1, info1, atomMap):
        # get all the atoms in each fragment
        fragAtoms = {}
    
        for bitName in fp1.GetNonzeroElements().keys():
            bitInfo = info1[bitName]
    
            for x in bitInfo:
                startNode, rad = x
                scope = set([startNode])
                # if radius is zero, the start node is the only node in the scope
                if rad == 0:
                    fragAtoms[str(bitName) + '_' + str(startNode)] = scope
                    
                # if not iterate through atomMap to get the other nodes
                else:
                    adjNodes = atomMap[startNode]
                    scope = scope | adjNodes
                    i = 1
                    while i < rad:
                        # get the nodes for the next iteration
                        adjNodes2 = set()
                        for x in adjNodes:
                            adjNodes2 = adjNodes2 | atomMap[x]
                            scope = scope | adjNodes2
                        adjNodes = adjNodes2
                        i += 1
                    fragAtoms[str(bitName) + '_' + str(startNode)] = scope
        return fragAtoms

    # reactingFragDists 
    def _calculate_distances(self, subDistReact, hitFrags1, fragAtoms1, bi1):
        # find the closest instances and measure distances
        fragDist = {}
        bi1_filtered = {}
        for ra, dists in subDistReact.items():
            # get the closest instance for each fragment, by finding the smallest radius
            fd = {}
            bi = {}
            
            for frag in hitFrags1:
                # get the fragment instance details
                fragName, startAtom = list(map(int, frag.split('_')))
                      
                # get the fragment atoms
                frag_atoms = fragAtoms1[frag]
                if ra not in frag_atoms:
                    continue
                
                # look for the furthest atom for the reacting atoms
                for d, a in dists:
                    if a in frag_atoms:
                        # see if its the smallest radius for the reacting atom
                        if fragName in fd.keys() and d > fd[fragName][0]:
                            break
                        else:
                            fd[fragName] = [d, startAtom, ra]
                            
                            # bi dictionary contains all instances, filter for this o
                            bi_instance = [x for x in bi1[fragName] if x[0] == startAtom]
                            bi[fragName] = bi_instance[0]
                            break
                        
            fragDist = self._dict_update(fragDist, fd)
            bi1_filtered = self._dict_update(bi1_filtered, bi)
    
        return(bi1_filtered, fragDist)

    def _dict_update(self, d2, d1):
        for key, value in d1.items():
            d2.setdefault(key, [])
            d2[key].append(value)
        return d2

    def _log_fun_weight(self, c):   
        # log function 
        weights1 = []
        k = 0.75
        i = 3
        for x in list(range(0, c + 1)):
            weights1.append((1 / (1 + math.e ** (k * (x - i)))) + 0)
        return ([x / sum(weights1) for x in weights1])


if __name__ == "__main___":
    print('test functions')

    subsmiles = ['O=C([O-])C(=O)Cc1ccc(O)cc1', '[NH3+][C@@H](CCC(=O)[O-])C(=O)[O-]']
    prodsmiles = ['O=C([O-])CCC(=O)C(=O)[O-]', '[NH3+][C@@H](Cc1ccc(O)cc1)C(=O)[O-]']
    
    # set the RDKit atom indices
    process_compounds = Compounds()
    subMols = [process_compounds.process_compF(smile, 8) + (smile,) for smile in subsmiles]
    prodMols = [process_compounds.process_compF(smile, 8) + (smile,) for smile in prodsmiles]
    subSmiles = [Chem.MolToSmiles(x[5]) for x in subMols]
    prodSmiles = [Chem.MolToSmiles(x[5]) for x in prodMols]
    
    print('')
    reacting_atoms = ReactingAtoms()
    rxn_mapper = RXNMapper() 
    reactingAtoms, conf = reacting_atoms.get_reacting_atoms(subMols, prodMols, 
                                                            subsmiles, prodsmiles, rxn_mapper)
    
    print('')
    reacting_fragments = ReactingFragments()
    queryRF = {}
    queryDists = {}
    for comp in subMols + prodMols:
        rfs, dists = reacting_fragments.get_distances(comp, reactingAtoms[comp[5]])
        queryRF[comp[5]] = rfs
        queryDists[comp[5]] = dists