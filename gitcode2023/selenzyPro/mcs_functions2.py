#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  3 12:56:49 2022

@author: ruth
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rxnmapper import RXNMapper
from rdkit.Chem import AllChem
from collections import Counter
import re
from rdkit.Chem import Draw

# import sys
# import os
#from rdkit.Chem import rdFMCS

#from rdkit.Chem import Draw



# https://gist.github.com/greglandrum/61c1e751b453c623838759609dc41ef1
def draw_chemical_reaction(smiles, highlightByReactant=False, font_scale=1.5):
    from rdkit import Chem
    from rdkit.Chem import Draw
    from rdkit.Chem import rdChemReactions
    from rdkit.Chem.Draw import rdMolDraw2D
    from rdkit.Chem.Draw import IPythonConsole
    from IPython.display import SVG, display
    rxn = rdChemReactions.ReactionFromSmarts(smiles,useSmiles=True)
    trxn = rdChemReactions.ChemicalReaction(rxn)
    # move atom maps to be annotations:
    for m in trxn.GetReactants():
        moveAtomMapsToNotes(m)
    for m in trxn.GetProducts():
        moveAtomMapsToNotes(m)
    d2d = rdMolDraw2D.MolDraw2DSVG(800,300)
    d2d.drawOptions().annotationFontScale=font_scale
    d2d.DrawReaction(trxn,highlightByReactant=highlightByReactant)

    d2d.FinishDrawing()

    return d2d.GetDrawingText()

def moveAtomMapsToNotes(m):
    for at in m.GetAtoms():
        if at.GetAtomMapNum():
            at.SetProp("atomNote",str(at.GetAtomMapNum()))

def get_RAs(rxn, rxn1, react_atoms_subs):
    #rxn, rxn1 = rxn.GetReactants(), rxn1.GetReactants()
    #rxn, rxn1, react_atoms_subs = rxn.GetProducts(), rxn2.GetReactants(), react_atoms_prods
    # map the compounds and atoms
    subMap=[]
    subReact=[]

    # get the inchis - not matching compounds can match inchi
    for x in rxn: x.SetProp('inchi', str(Chem.inchi.MolToInchi(x)))
    for x in rxn1: x.SetProp('inchi', str(Chem.inchi.MolToInchi(x)))


    #         # map the mapped compound to the unmapped compound
    #         for x in react_mapped_subs[i2].GetSubstructMatches(react_subs[i2]):
    
    #             atom_ind_ind = dict(zip(x, list(range(0, comp1.GetNumAtoms()))))
    #             reactingAtoms[s]= set([atom_ind_ind[x] for x in reacting_atoms])

    p = []
    for i1, comp1 in enumerate(rxn):
        Chem.SanitizeMol(comp1)
        comp1=Chem.RemoveHs(comp1)
        # check iif the unmapped compound got paired with a mapped compound
        # show_atom_number(comp1, 'atomLabel')
        mpd = 0
        for i2, comp2 in enumerate(rxn1):
            if i2 in p: continue
            comp2=Chem.RemoveHs(comp2)
            
            if comp1.GetProp('inchi') == comp2.GetProp('inchi'):
            #if comp1.HasSubstructMatch(comp2) and comp2.HasSubstructMatch(comp1):
                # print(i1, i2)
                p.append(i2)


                # get the reacting atoms
                react_a = react_atoms_subs[i2]   
                
                # map the atom indices of the mapped onto the unmapped
                a_map = comp1.GetSubstructMatch(comp2)
                if comp1.GetNumAtoms() -1 < max(a_map) or comp2.GetNumAtoms() -1 < max(a_map): print('mapping failed')
                subMap.append(a_map)
                
                # get the reacting atoms
                subReact.append( [a_map[x] for x in react_a] )
                # mark this compound as taken
                
                mpd = 1
                break
            
        # if the compound wasnt mapped, add an empty list
        if mpd == 0:
            print('e')
            subMap.append([])
            subReact.append([])
            
    return(subReact, subMap)

def show_atom_number(mol, label):
    for atom in mol.GetAtoms():
        atom.SetProp(label, str(atom.GetIdx()))
    return mol


def moveAtomMapsToNotes(m):
    for at in m.GetAtoms():
        if at.GetAtomMapNum():
            at.SetProp("atomNote",str(at.GetAtomMapNum()))
            print(at.GetIdx(),at.GetAtomMapNum() )



# def rxnMapper_fun(smile, subMols, prodMols, rxn_mapper):
#     print('in RXNMapper_fun')
#     from rdkit.Chem import Draw

#     # subMols, prodMols = [compDetails[comp] for x in subs], [compDetails[comp] for x in prods]
#     # subMols, prodMols =     [compDetails[x] for x in prods], [compDetails[x] for x in subs]
        
#     # get_attention_guided_atom_maps_for_reactions - never shown in examples
#     results = rxn_mapper.get_attention_guided_atom_maps([smile])[0]  
#     smileM, conf = results['mapped_rxn'], results['confidence'] 
#     #display(SVG(draw_chemical_reaction(smileM,  highlightByReactant=True)))
#     rxn1 = AllChem.ReactionFromSmarts(smileM, useSmiles=True)
#     rxn2 = AllChem.ReactionFromSmarts(smileM.split('>>')[1] + '>>' + smileM.split('>>')[0], useSmiles=True)
#     rxn = AllChem.ReactionFromSmarts(smile, useSmiles=True)
#     rxn1.Initialize()
#     rxn2.Initialize()
#     react_atoms_subs = rxn1.GetReactingAtoms(mappedAtomsOnly=True)     
#     react_atoms_prods = rxn2.GetReactingAtoms(mappedAtomsOnly=True) 
#     # show_atom_number(rxn1.GetReactants()[0], 'atomLabel')


#     # move atom maps to be annotations:
#     for m in rxn1.GetReactants():
#         moveAtomMapsToNotes(m)
#         comp, atomMap, fpM1, info, fragAtoms1 = process_comp(m, 8)
#     for m in rxn1.GetProducts():
#         moveAtomMapsToNotes(m)


#     subReact, subMap = get_RAs(rxn.GetReactants(), rxn1.GetReactants(), react_atoms_subs)
#     prodReact, prodMap = get_RAs(rxn.GetProducts(), rxn2.GetReactants(), react_atoms_prods)
 
#     # make string: reacting atoms dict
#     if len(re.split('\\.|>>', smile)) != len(subMap + prodMap): print('mapping error')
    
#     reactAtoms={}
#     for i, x in enumerate(smile.split('>>')[0].split('.')):
#         if x not in reactAtoms: 
#             reactAtoms[x] =set()
#         reactAtoms[x] = reactAtoms[x] | set(subReact[i])
 
#     for i, x in enumerate(smile.split('>>')[1].split('.')):
#         if x not in reactAtoms: 
#             reactAtoms[x] =set()
#         reactAtoms[x] = reactAtoms[x] | set(prodReact[i])
        
#     return(reactAtoms,  conf)
            






def rxnMapper_fun(smile, subMols, prodMols, rxn_mapper):
    print('in RXNMapper_fun')
    from rdkit.Chem import Draw

    # subMols, prodMols = [compDetails[comp] for x in subs], [compDetails[comp] for x in prods]
    # subMols, prodMols =     [compDetails[x] for x in prods], [compDetails[x] for x in subs]
        
    # get_attention_guided_atom_maps_for_reactions - never shown in examples
    results = rxn_mapper.get_attention_guided_atom_maps([smile])[0]  
    smileM, conf = results['mapped_rxn'], results['confidence'] 
    #display(SVG(draw_chemical_reaction(smileM,  highlightByReactant=True)))
    rxn1 = AllChem.ReactionFromSmarts(smileM, useSmiles=True)
    rxn2 = AllChem.ReactionFromSmarts(smileM.split('>>')[1] + '>>' + smileM.split('>>')[0], useSmiles=True)
    rxn = AllChem.ReactionFromSmarts(smile, useSmiles=True)
    rxn1.Initialize()
    rxn2.Initialize()
    react_atoms_subs = rxn1.GetReactingAtoms(mappedAtomsOnly=True)     
    react_atoms_prods = rxn2.GetReactingAtoms(mappedAtomsOnly=True) 
    # show_atom_number(rxn1.GetReactants()[0], 'atomLabel')

    
    subReact, subMap = get_RAs(rxn.GetReactants(), rxn1.GetReactants(), react_atoms_subs)
    prodReact, prodMap = get_RAs(rxn.GetProducts(), rxn2.GetReactants(), react_atoms_prods)
 
    # make string: reacting atoms dict
    if len(re.split('\\.|>>', smile)) != len(subMap + prodMap): print('mapping error')
    
    reactAtoms={}
    for i, x in enumerate(smile.split('>>')[0].split('.')):
        if x not in reactAtoms: 
            reactAtoms[x] =set()
        reactAtoms[x] = reactAtoms[x] | set(subReact[i])
 
    for i, x in enumerate(smile.split('>>')[1].split('.')):
        if x not in reactAtoms: 
            reactAtoms[x] =set()
        reactAtoms[x] = reactAtoms[x] | set(prodReact[i])
        
    return(reactAtoms,  conf)
            


# def rxnMapper_fun(smile, subMols, prodMols, rxn_mapper):
    # smile = 'CC(C)S.CN(C)C=O.Fc1cccnc1F.O=C([O-])[O-].[K+].[K+]>>CC(C)Sc1ncccc1F'
    # smile = 'BrCCOCCBr.CCN(C(C)C)C(C)C.CCOC(C)=O.CN(C)C=O.Cl.NCC(F)(F)CO>>OCC(F)(F)CN1CCOCC1'
    # smile = 'CCOCC.C[Mg+].O=Cc1ccc(F)cc1Cl.[Br-]>>CC(O)c1ccc(F)cc1Cl'
    # smile = 'BrCCOCCBr.CCN(C(C)C)C(C)C.CCOC(C)=O.CN(C)C=O.Cl.NCC(F)(F)CO>>OCC(F)(F)CN1CCOCC1'
    # smile ='C1=CCC=C1.O=C1C=CC(=O)C1>>O=C1CC(=O)C2C3C=CC(C3)C12'
    # smile = 'OC1C=CC=CC=1C(C=C)(C)C>>C1=CC=C(OCC=C(C)C)C=C1'
    # react_subs =  [process_comp(x, 8)[0] for x in smile.split('>>')[0].split('.')]
    # react_prods = [process_comp(x, 8)[0] for x in smile.split('>>')[1].split('.')]
    # subMols, prodMols = [compDetails[x] for x in subs], [compDetails[x] for x in prods]
    
    # mapp the atoms between substrates/products
    #sys.stdout = NullIO()
    #rxn_mapper = RXNMapper()  
    #sys.stdout = sys.__stdout__   


    # # get_attention_guided_atom_maps_for_reactions - never shown in examples
    # mapped_smiles = rxn_mapper.get_attention_guided_atom_maps([smile])[0]  
    # smile_m, conf = mapped_smiles['mapped_rxn'], mapped_smiles['confidence']
    
    # # reorder the compounds 
    # subs1 =  [ Counter(re.sub(r"[^A-Za-z]+", '', x)) for x in smile.split('>>')[0].split('.')]
    # subs2 =  [ Counter(re.sub(r"[^A-Za-z]+", '', x)) for x in smile_m.split('>>')[0].split('.')]
    # subs1 =  [{k:v for k, v in x.items() if k !='H'}  for x in subs1 ]
    # subs2 =  [{k:v for k, v in x.items() if k !='H'}  for x in subs2 ]
    # sub_order = []
    # for i1, x1 in enumerate(subs1):
    #     for i2, x2 in enumerate(subs2):
    #         if x1 == x2:
    #             sub_order.append(i2)
    #             break
                
    # prods1 =  [ Counter(re.sub(r"[^A-Za-z]+", '', x)) for x in smile.split('>>')[1].split('.')]
    # prods2 =  [ Counter(re.sub(r"[^A-Za-z]+", '', x)) for x in smile_m.split('>>')[1].split('.')]
    # prods1 =  [{k:v for k, v in x.items() if k !='H'}  for x in prods1 ]
    # prods2 =  [{k:v for k, v in x.items() if k !='H'}  for x in prods2 ] 
    # prod_order = []
    # for i1, x1 in enumerate(prods1):
    #     for i2, x2 in enumerate(prods2):
    #         if x1 == x2:
    #             prod_order.append(i2)
    #             break
            
    # smile_m = '.'.join([smile_m.split('>>')[0].split('.')[i] for i in sub_order]) + '>>' + '.'.join([smile_m.split('>>')[1].split('.')[i] for i in prod_order])                
    # smile_m2 = smile_m.split('>>')[1] + '>>'+ smile_m.split('>>')[0]
    
    # # get reacting atoms by index (not mapped number) 
    # rxn = AllChem.ReactionFromSmarts(smile, useSmiles=True)
    # rxn1 = AllChem.ReactionFromSmarts(smile_m, useSmiles=True)
    # rxn1.Initialize()
    # rxn2 = AllChem.ReactionFromSmarts(smile_m2, useSmiles=True)
    # rxn2.Initialize()
    
    # react_atoms_subs = rxn1.GetReactingAtoms(mappedAtomsOnly=True)     
    # react_atoms_prods = rxn2.GetReactingAtoms(mappedAtomsOnly=True)
    
    # # get the substrates / products for the reactions    
    # react_subs_d =  {x[5]: x[0] for x in subMols}
    # react_prods_d = {x[5]: x[0] for x in prodMols}
    # react_subs = [react_subs_d[x] for x in smile.split('>>')[0].split('.')]
    # react_prods = [react_prods_d[x] for x in smile.split('>>')[1].split('.')]    
    # react_mapped_subs =  rxn1.GetReactants()
    # react_mapped_prods = rxn1.GetProducts()
    
    # # map over the reacting compounds to the unmapped atom indexes
    # reactingAtoms ={}  
    # for i2, reacting_atoms in enumerate(react_atoms_subs):
    #     s = smile.split('>>')[0].split('.')[i2]
    #     reactingAtoms[s]= set() 
    #     if len(reacting_atoms) > 0:
    #         # map the mapped compound to the unmapped compound
    #         for x in react_mapped_subs[i2].GetSubstructMatches(react_subs[i2]):
    #             atom_ind_ind = dict(zip(x, list(range(0, react_subs[i2].GetNumAtoms()))))
    #             reactingAtoms[s]= set([atom_ind_ind[x] for x in reacting_atoms])
   
    # for i2, reacting_atoms in enumerate(react_atoms_prods):
    #     s = smile.split('>>')[1].split('.')[i2]
    #     reactingAtoms[s]= set() 
    #     if len(reacting_atoms) > 0:
    #         # map the mapped compound to the unmapped compound
    #         for x in react_mapped_prods[i2].GetSubstructMatches(react_prods[i2]):
    #             atom_ind_ind = dict(zip(x, list(range(0, react_prods[i2].GetNumAtoms()))))
    #             reactingAtoms[s]= set([atom_ind_ind[x] for x in reacting_atoms])   
    # # print(reactingAtoms)
    # # react_subs[0]
    # # react_prods[0] 
    # return(reactingAtoms, conf)



def test_multi(x):
    return(x*2)


def get_MCS2(comps1, comps2, fragSize):
    
    # process the comps
    Comps1 = []
    sucesses = []
    for c1 in comps1:
        comp, atomMap, fpM, info, fragAtoms = process_comp(c1, fragSize)
        Comps1.append([comp, atomMap, fpM, info, fragAtoms])
        if comp is not None:
            sucesses.append(c1)

    Comps2 = []
    for c2 in comps2:
        comp, atomMap, fpM, info, fragAtoms = process_comp(c2, fragSize)
        Comps2.append([comp, atomMap, fpM, info, fragAtoms])

    # for each compound in comps1 get the biggest MCS for the compounds in comp2, then use this largest MCS to get the reacting frags
    rfrags =[]
    adj = []
    fps=[]
    noOverlap=0
    totalOverlap=0
    cantMap=0
    noComp2=0
    GetSubstructFail=0
    
    for x in Comps1:
        comp1, atomMap1, fpM1, info1, fragAtoms1 = x
        
        if comp1 is None:
            print('no comp1')
            rfrags.append([])
            adj.append([])  
            fps.append([])
            continue
        fps.append(fpM1)
        comp1RFs={}
        comp1RFAs={}
        
        # for each comp2 get the biggest mcs
        m=0 # number of atoms
        # biggestMcs = [] # mcs smarts
        # i = [] # index
        for i1, y in enumerate(Comps2):
            comp2, atomMap2, fpM2, info2, fragAtoms2 = y
            if comp2 is None:
                noComp2+=1
                continue


            res = Chem.rdFMCS.FindMCS([comp1, comp2], timeout=1)
            m = res.numAtoms  
            # this comp2 might not have RFs but others will
            if m==0:
                noOverlap+=1
                continue
            bMcs = res.smartsString
            
            # get the reacting frags
            if m < comp1.GetNumAtoms():

                bMcs = Chem.MolFromSmarts(bMcs)
                # outerAtoms = set([y  for x in comp1.GetSubstructMatches(bMcs) for y in set(range(0,comp1.GetNumAtoms()))- set(x)])

                # mcsAtoms=set([y for x in comp1.GetSubstructMatches(bMcs) for y in x - ])
                # reactAtoms = [set([k]) | v-outerAtoms for k, v in atomMap1.items() if k in outerAtoms and len(v-outerAtoms)>0]
                reactAtoms = [set([x]) | atomMap1[x]-set(substruct)  for substruct in comp1.GetSubstructMatches(bMcs) for x in substruct if len(atomMap1[x]-set(substruct))>0 ]
                if len(reactAtoms)==0 and len(comp1.GetSubstructMatches(bMcs))==0:
                    print("GetSubstructMatches failed")
                    GetSubstructFail +=1
                    continue
                
                if len(reactAtoms) == comp1.GetNumAtoms():
                    cantMap+=1
                    continue
                
                # get all the overlapping fragments
                rf = {k:v for k, v in fragAtoms1.items() for reacAtomSet in reactAtoms if len(v.intersection(reacAtomSet))>0  }
                
                # # separate the smallest fragments from the larger ones
                rf1, rfAdj1 = get_common_atoms_frags(rf, fragAtoms1, reactAtoms)
                if len(rf1)==0:
                    print('problem overlapping frags not substrate')
                    continue
                    
                for k, v in rf1.items(): 
                    if k not in comp1RFs.keys():
                        comp1RFs[k] =v
                        
                for k, v in rfAdj1.items(): 
                    if k not in comp1RFAs.keys():
                        comp1RFAs[k] =v 
                        
                        
            elif comp2.GetSubstructMatches(comp1):
                ### if comp is a substructure of comp 1
                if comp1.GetNumAtoms() == comp2.GetNumAtoms() == m:
                    print('compound unchanged')
                    totalOverlap+=1
                    continue
                    
                # print('the whole compound is substructure')
                # get the atoms in comp2 that are in comp1
                mapAt= comp2.GetSubstructMatches(comp1) 
                # map the atoms comp2:comp1
                mapAtoms = {}
                for mapping in mapAt:
                    for i, x in enumerate(mapping):
                        if x not in mapAtoms.keys():
                            mapAtoms[x]=set()
                        mapAtoms[x].add(i)

                # get the atoms on the reacting end of comp2
                mcsAtoms2=set(mapAtoms.keys())
                #reactAtoms2 = [set([k]) | v-mcsAtoms2 for k, v in atomMap2.items() if k in mcsAtoms2 and len(v-mcsAtoms2)>0]
                
                outerAtoms2 = set([y  for x in comp2.GetSubstructMatches(comp1) for y in set(range(0,comp2.GetNumAtoms()))- set(x)])
                if len(outerAtoms2) == comp1.GetNumAtoms():
                    cantMap+=1
                    continue

                # reactAtoms2 = [set([k]) | v-outerAtoms2 for k, v in atomMap2.items() if k in outerAtoms2 and len(v-outerAtoms2)>0]
                reactAtoms2 = [set([x]) | atomMap2[x]-set(substruct)  for substruct in comp2.GetSubstructMatches(comp1) for x in substruct if len(atomMap2[x]-set(substruct))>0 ]
            
                # map reacting atoms onto compound1
                reactAtoms = []
                for x in reactAtoms2:
                    for y in x:
                        if y in mapAtoms.keys():
                            reactAtoms.append(mapAtoms[y])
                
                # if the reacting Atoms aren't touching then separate them out
                removeI = []
                addI = []
                for x in reactAtoms:
                    if len(x) == 1: continue
                    linked=set()
                    for i in x:
                        linked= linked|atomMap1[i]
                        
                    if x == (x - linked):
                        removeI.append(x)
                        for y in x:
                            addI.append(set([y]))
                 
                reactAtoms = reactAtoms + addI
                reactAtoms = [x for x in reactAtoms if x not in removeI]
                
                # get all the overlapping fragments
                rf = {k:v for k, v in fragAtoms1.items() for reacAtomSet in reactAtoms if len(v.intersection(reacAtomSet))>0  }
                
                # # separate the smallest fragments from the larger ones
                rf1, rfAdj1 = get_common_atoms_frags(rf, fragAtoms1, reactAtoms)

                if len(rf1)==0:
                    print('problem overlapping frags subset')
                    crash
               
                for k, v in rf1.items(): 
                    if k not in comp1RFs.keys():
                        comp1RFs[k] =v
                        
                for k, v in rfAdj1.items(): 
                    if k not in comp1RFAs.keys():
                        comp1RFAs[k] =v
                        
            else:
                #changed bonds same atoms
                # bMcs = Chem.MolFromSmarts(bMcs)
                # compMCS, atomMapMCS, fpMMCS, infoMCS, fragAtomsMCS = process_comp(Chem.MolToSmiles(bMcs))
                bMcs = Chem.MolFromSmarts(bMcs)
                
                if bMcs.GetNumBonds() == comp1.GetNumBonds():
                    # print('compound unchanged')
                    totalOverlap+=1
                    continue             
                
                # map the atoms
                mapAt= comp1.GetSubstructMatches(bMcs) 
                mapAtoms = {}
                for mapping in mapAt:
                    for i, x in enumerate(mapping):
                        if x not in mapAtoms.keys():
                            mapAtoms[x]=set()
                        mapAtoms[x].add(i)
                        # mapAtoms[i]=set()
                        # mapAtoms[i].add(x)
                
                # get the bonds in MCS
                bondsMCS = [ [x.GetBeginAtomIdx(), x.GetEndAtomIdx(), str(x.GetBondType())] for x in bMcs.GetBonds()]
                bondsMCSR = [[x[1], x[0], x[2]] for x in bondsMCS]
                # get the bonds for comp1
                bondsComp1 = [ [x.GetBeginAtomIdx(), x.GetEndAtomIdx(), str(x.GetBondType())] for x in comp1.GetBonds()] 

                reactAtoms = []
                for x in bondsComp1:
                    if len(list(mapAtoms[x[0]]))>1 or len(list(mapAtoms[x[1]]))>1:
                        print('death')
                        crash
                    a1, a2, t = list(mapAtoms[x[0]])[0], list(mapAtoms[x[1]])[0], x[2]
                    if [a1, a2, t] not in bondsMCS and [a1, a2, t] not in bondsMCSR:# and [a1, a2, t] not in atomMapList2 :
                        print('comp1', x, "\n",'mcs', [a1, a2, t], "\n\n")
                        reactAtoms.append(set([x[0], x[1]]))
                
                # get all the overlapping fragments
                rf = {k:v for k, v in fragAtoms1.items() for reacAtomSet in reactAtoms if len(v.intersection(reacAtomSet))>0  }
                
                # # separate the smallest fragments from the larger ones
                rf1, rfAdj1 = get_common_atoms_frags(rf, fragAtoms1, reactAtoms)

                if len(rf1)==0:
                    print('problem changed bonds/same atoms')
                    crash
               
                for k, v in rf1.items(): 
                    if k not in comp1RFs.keys():
                        comp1RFs[k] =v
                        
                for k, v in rfAdj1.items(): 
                    if k not in comp1RFAs.keys():
                        comp1RFAs[k] =v                    
                    
                        
            if len(rf1)==0 and totalOverlap==0 and cantMap == 0:
                print('no RFs for a comp1')
                crash
        
        # store all the fragments from mapping comp1 onto all comp2s
        if len(comp1RFs)>0:
            rfrags.append(dict_to_sparse(comp1RFs)[0])
        else:
            rfrags.append([])
            if len(comps2)-noComp2 == totalOverlap:
                print('compounds are all the same')
            elif len(comps2)-noComp2==noOverlap:
                print('no overlap between comps1 and comps2')
            elif len(comps2)-noComp2==cantMap:
                print('couldnt map comps1 and comps2')
            elif len(comps2)-noComp2==GetSubstructFail:
                print('GetSubstructure failed')
            else:
                print('combined problem', totalOverlap,noOverlap, cantMap, GetSubstructFail )
                # crash
                       
        if len(comp1RFAs)>0:
            adj.append(dict_to_sparse(comp1RFAs)[0])
        else:
            adj.append([])
            
        # print('comp1 end rfrags', comp1RFs.keys())

    if len(rfrags)==0: 
        print('no rfs')
        crash
                                                     
    return  rfrags, adj, fps, sucesses

def process_comp(comp1,fragSize):
    # get all the details for the fingerprints
    
    if type(comp1)==str:
        try:
            comp1 = Chem.MolFromSmiles(comp1)
        except:
            print('failed MolFromSmiles')
            return(None, None, None, None, None)
            
    try:
        # Chem.SanitizeMol(comp1)
        # rdMolStandardize.Uncharger().uncharge(comp1)
        comp1 = neutralize_atoms(comp1)
        Chem.rdmolops.RemoveStereochemistry(comp1)
        comp1 = mol_with_atom_index(comp1)
 
        # Draw.DrawMorganBits(prints, molsPerRow=4, legends=[str(x) for x in list(fpM1.GetNonzeroElements().keys()) ])
               
        # get each atoms adjacent neighbours
        atomMap1 = getAtoms(comp1)
        
        # Morgan fingerprints
        info1={}
        fpM1 = AllChem.GetMorganFingerprint(comp1, fragSize, bitInfo=info1, invariants=AllChem.GetConnectivityInvariants(comp1, includeRingMembership=False))
        
        # get all the atoms in each fragment, returns {fragname:[{atoms}, startnode, radius]}           
        fragAtoms1 = getAtomFragments(fpM1, info1, atomMap1)   
    
        return(comp1, atomMap1, fpM1, info1, fragAtoms1)
    except:
        print('process_comp failed')

        return(None, None, None, None, None)
    

# def get_morg(mol, inchi=0):   
#     # rdMolStandardize.Uncharger().uncharge(mol)  
#     # mol = Chem.MolFromSmiles(neutralize_atoms(mol)) 
#     # Chem.rdmolops.RemoveStereochemistry(mol)  #  turns it into a string!!
#     Chem.SanitizeMol(mol)
#     info1 = {}
#     fpM = AllChem.GetMorganFingerprint(mol, 8, bitInfo=info1, invariants=AllChem.GetConnectivityInvariants(mol, includeRingMembership=False))
#     atomMap1 = get_atoms(mol)
#     return fpM, Chem.MolToSmiles(mol), mol.GetNumAtoms(), info1, atomMap1

def process_compF(comp1, size):
    # get all the details for the fingerprints
    
    if type(comp1)==str:
        try:
            comp1 = Chem.MolFromSmiles(comp1)
        except:
            print('failed MolFromSmiles')
            #return(None, None, None, None, None)
            
    try:
        Chem.SanitizeMol(comp1)
        # rdMolStandardize.Uncharger().uncharge(comp1)
        # comp1 = neutralize_atoms(comp1)
        # Chem.rdmolops.RemoveStereochemistry(comp1)
        
        comp1 = mol_with_atom_index(comp1)
 
        # Draw.DrawMorganBits(prints, molsPerRow=4, legends=[str(x) for x in list(fpM1.GetNonzeroElements().keys()) ])
               
        # get each atoms adjacent neighbours
        atomMap1 = getAtoms(comp1)
        
        # Morgan fingerprints
        info1={}
        fpM1 = AllChem.GetMorganFingerprint(comp1, size, bitInfo=info1, invariants=AllChem.GetConnectivityInvariants(comp1, includeRingMembership=False))
        
        # get all the atoms in each fragment, returns {fragname:[{atoms}, startnode, radius]}           
        fragAtoms1 = getAtomFragments(fpM1, info1, atomMap1)   
    
        return(comp1, atomMap1, fpM1, info1, fragAtoms1)
    except:
        print('unable to process compound')    

def process_compF_safe(comp1, size):
    # get all the details for the fingerprints
    
    if type(comp1)==str:
        try:
            comp1 = Chem.MolFromSmiles(comp1)
        except:
            print('failed MolFromSmiles')
            #return(None, None, None, None, None)
            
    try:
        # Chem.SanitizeMol(comp1)
        # rdMolStandardize.Uncharger().uncharge(comp1)
        comp1 = neutralize_atoms(comp1)
        Chem.rdmolops.RemoveStereochemistry(comp1)
        comp1 = mol_with_atom_index(comp1)
 
        # Draw.DrawMorganBits(prints, molsPerRow=4, legends=[str(x) for x in list(fpM1.GetNonzeroElements().keys()) ])
               
        # get each atoms adjacent neighbours
        atomMap1 = getAtoms(comp1)
        
        # Morgan fingerprints
        info1={}
        fpM1 = AllChem.GetMorganFingerprint(comp1, size, bitInfo=info1, invariants=AllChem.GetConnectivityInvariants(comp1, includeRingMembership=False))
        
        # get all the atoms in each fragment, returns {fragname:[{atoms}, startnode, radius]}           
        fragAtoms1 = getAtomFragments(fpM1, info1, atomMap1)   
    
        return(comp1, atomMap1, fpM1, info1, fragAtoms1)
    except:
        print('unable to process compound')    
    

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




def get_common_atoms_frags(rf, fragAtoms, reactAtoms):
    # reaction at one site:

    minRF = {}
    maxRF = {}
    
    for Aset in reactAtoms:
        try:
            #trueRF = {k:v for k,v in rf.items() if len(v.intersection(Aset)) == len(Aset)}        
            minsize = min([len(x) for x in rf.values()])
            minRF.update({k:v for k, v in trueRF.items() if len(v) == minsize})
            maxRF.update({k:v for k, v in fragAtoms.items() if k not in minRF.keys() and len(v.intersection(Aset))>0 })
        except:
            continue  
    return minRF, maxRF


def sort_rf(rf, reactAtoms):
    # min_rfs contain the fewest possable non-reacting atoms
    # all fragments containing no non-reacting atoms are min-rfs
    # if not all the reacting atoms can be covered by reacting fragment only fragments, choose those with the highest proportion reacting fragments
    # score the rfs for reacting atoms vs none reacting atoms, prefect score is 1
    # plus tiny score based on size

    scores=[]
    min_rf= {}
    covered = set()
    for k, fragset in rf.items():
        hits = fragset.intersection(reactAtoms)
        len_bonus = len(hits)/1000
        score = (len(hits)/len(fragset)) + len_bonus
        if score >1:
            min_rf[k]= fragset
            covered = covered | fragset
        else:
            scores.append([score, k, fragset])
    
    # all reacting atoms covered by prefect fragments
    if len(covered) == len(reactAtoms):
        rfa = {k:v for k, v in rf.items() if k not in min_rf.keys()}
        return(min_rf, rfa )

    # find the fragments with the highest proportion reacting atoms to plug the gaps
    else:
        while len(covered)> len(reactAtoms):
            print(reactAtoms-covered)
            for score, k, fragset in scores:
                if len(fragset.intersection(covered)) == len(fragset):
                    continue
                min_rf[k]= fragset
                covered = covered | fragset
    rfa = {k:v for k, v in rf.items() if k not in min_rf.keys()}
    return(min_rf, rfa )




def dict_to_sparse(rfs):
    # make sure the input is a list
    if type(rfs) != list:
        rfs = [rfs]
    
    # make into sparse int vectors
    rf_list = []
    for rf in rfs:
        RF1 =  [int(x.split('_')[0]) for x in rf.keys()]
    
        # Make empty UIntSparseIntVect, then update them with the sparse int vectors       
        rf1 = AllChem.GetMorganFingerprint(Chem.MolFromSmiles(''), 5)
        rf1.UpdateFromSequence(RF1)
        rf_list.append(rf1)
    return(rf_list)



def reactFragDists(comp1, reactAtoms):
    ### set up the requeirements  
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


def duplicateMultipleComps(Comps1):
    Comps1multi = []
    for x in Comps1: 
        for i, y in enumerate([x]*int(x[6])):
            newStr = y[5]+'_'+str(i+1)
            newProt = list(y)
            newProt[5] = newStr
            Comps1multi.append(tuple(newProt))
    return Comps1multi


### available atoms 
    
#def MCS_avilAtoms       
def getReactingFrags(comps1, comps2):
    # comps1, comps2 = list(rTarget[r][0].keys()), list(rTarget[r][1].keys())
    # process comps failing to remove N+
    # 'C1=CC(=C[N+](=C1)C2C(C(C(O2)COP(=O)([O-])OP(=O)(O)OCC3C(C(C(O3)N4C=NC5=C(N=CN=C54)N)OP(=O)(O)O)O)O)O)C(=O)N'
    # comps1, comps2 = [comp_smiles[x] for x in subs], [comp_smiles[x] for x in prods]
    # comps1, comps2 = [['C1=CC(=CC=C1CC(C(=O)O)N)O'], ['C1=CC(=C(C=C1CC(C(=O)O)N)O)O']]
    # comps1, comps2 = [compDetails[x]+[subs_count[x]] for x in subs], [compDetails[x]+[prods_count[x]] for x in prods]
    
    #print('getReactingFrags fun')
    # process the comps
    if type(comps1[0])==str:
        Comps1 = []
        for c1 in comps1:
            if type(c1) == str:
                smile1 = c1
                comp, atomMap, fpM, info, fragAtoms1 = process_comp(c1, 8)
            else:
                comp, atomMap, fpM, info, fragAtoms1, smile1 = c1
            Comps1.append([comp, atomMap, fpM, info, fragAtoms1, smile1])
    
        Comps2 = []
        for c2 in comps2:
            if type(c2) == str:
                smile2 = c2
                comp, atomMap, fpM, info, fragAtoms2 = process_comp(c2, 8)
            else:
                comp, atomMap, fpM, info, fragAtoms2, smile2 = c2
            Comps2.append([comp, atomMap, fpM, info, fragAtoms2, smile2])
    else:
        Comps1 = comps1
        Comps2 = comps2
        
    Comps1, Comps2 = duplicateMultipleComps(Comps1), duplicateMultipleComps(Comps2)
    availableAtoms = {x[5]: set(range(0, x[0].GetNumAtoms()))  for x in Comps1+ Comps2}
    
    reactingAtoms =  {key: [] for key in availableAtoms.keys()}
    reactingBonds =  {key: [] for key in availableAtoms.keys()}
    BestPairs = []

    
    # loop through the MCS matches until all compounds have been mached
    while 1==1:
        bestPair, rbonds_rAtoms_cAtoms1, rbonds_rAtoms_cAtoms2 = getReactingAtomsBonds(Comps1, Comps2, availableAtoms)

        # if the best score is decent
        if bestPair==None or bestPair[0]<=0.5: 
            break
    
        # print('\ncomp 1:', bestPair[3], ', comp 2:', bestPair[4])
        # print('reacting Atoms', rbonds_rAtoms_cAtoms1[1] )
        # print('reacting Atoms', rbonds_rAtoms_cAtoms2[1] )
        
        # save the results
        reactingBonds[bestPair[1]].append(rbonds_rAtoms_cAtoms1[0]) 
        reactingAtoms[bestPair[1]].append(rbonds_rAtoms_cAtoms1[1])
        reactingBonds[bestPair[2]].append(rbonds_rAtoms_cAtoms2[0]) 
        reactingAtoms[bestPair[2]].append(rbonds_rAtoms_cAtoms2[1])   
        
        # remove the atoms covered by this match from the available atoms 
        availableAtoms[bestPair[1]] = availableAtoms[bestPair[1]] - rbonds_rAtoms_cAtoms1[2]
        availableAtoms[bestPair[2]] = availableAtoms[bestPair[2]] - rbonds_rAtoms_cAtoms2[2]
        #for k, v in availableAtoms.items(): print(k[0:10], len(v))
        BestPairs.append(bestPair[1:3])
 
    # clean up reacting atoms dict
    reactingAtoms = {k:[y for x in v for y in x ] for k, v in reactingAtoms.items() }
    for k, v in reactingAtoms.items():
        flat = set()
        for x in v:
            if type(x) == set:
                flat = flat | x
            else:
                flat.add(x)
        reactingAtoms[k]=flat 

    reactingAtoms2 = {}
    for k, v in reactingAtoms.items():
        s, i = k.split('_')
        if s not in reactingAtoms.keys():
            reactingAtoms2[s]=v
        else:
            reactingAtoms2[s]= reactingAtoms[s] | v
    
    # RF fragments
    RFs1, RFAs1 = get_RFs(reactingAtoms, Comps1)
    RFs2, RFAs2 = get_RFs(reactingAtoms, Comps2)
    # frags1 = [x[2] for x in Comps1]
    # frags2 = [x[2] for x in Comps2]
    # info1 = [x[3] for x in Comps1]
    # info2 = [x[3] for x in Comps2]

    return(reactingAtoms2, BestPairs )


#def MCS_avilAtoms       
def getReactingFrags2(comps1, comps2):
    # comps1, comps2 = list(rTarget[r][0].keys()), list(rTarget[r][1].keys())
    # process comps failing to remove N+
    # 'C1=CC(=C[N+](=C1)C2C(C(C(O2)COP(=O)([O-])OP(=O)(O)OCC3C(C(C(O3)N4C=NC5=C(N=CN=C54)N)OP(=O)(O)O)O)O)O)C(=O)N'
    # comps1, comps2 = [comp_smiles[x] for x in subs], [comp_smiles[x] for x in prods]
    # comps1, comps2 = [['C1=CC(=CC=C1CC(C(=O)O)N)O'], ['C1=CC(=C(C=C1CC(C(=O)O)N)O)O']]
    
    # process the comps
    if type(comps1[0])==str:
        Comps1 = []
        for c1 in comps1:
            if type(c1) == str:
                smile1 = c1
                comp, atomMap, fpM, info, fragAtoms1 = process_comp(c1, 8)
            else:
                comp, atomMap, fpM, info, fragAtoms1, smile1 = c1
            Comps1.append([comp, atomMap, fpM, info, fragAtoms1, smile1])
    
        Comps2 = []
        for c2 in comps2:
            if type(c2) == str:
                smile2 = c2
                comp, atomMap, fpM, info, fragAtoms2 = process_comp(c2, 8)
            else:
                comp, atomMap, fpM, info, fragAtoms2, smile2 = c2
            Comps2.append([comp, atomMap, fpM, info, fragAtoms2, smile2])
    else:
        Comps1 = comps1
        Comps2 = comps2
    
    
    availableAtoms = {x[5]: set(range(0, x[0].GetNumAtoms()))  for x in Comps1+ Comps2}
    #for k, v in availableAtoms.items(): print(k[0:10], len(v))
    
    reactingAtoms =  {key: [] for key in availableAtoms.keys()}
    reactingBonds =  {key: [] for key in availableAtoms.keys()}
    BestPairs = []

    # loop through the MCS matches until all compounds have been mached
    while 1==1:
        bestPair, rbonds_rAtoms_cAtoms1, rbonds_rAtoms_cAtoms2 = getReactingAtomsBonds(Comps1, Comps2, availableAtoms)
        # if the best score is decent
        if bestPair==None or bestPair[0]<=0.5: 
            #print('breaking')
            break
    
        # print('\ncomp 1:', bestPair[3], ', comp 2:', bestPair[4])
        # print('reacting Atoms', rbonds_rAtoms_cAtoms1[1] )
        # print('reacting Atoms', rbonds_rAtoms_cAtoms2[1] )
        
        # save the results
        reactingBonds[bestPair[1]].append(rbonds_rAtoms_cAtoms1[0]) 
        reactingAtoms[bestPair[1]].append(rbonds_rAtoms_cAtoms1[1])
        reactingBonds[bestPair[2]].append(rbonds_rAtoms_cAtoms2[0]) 
        reactingAtoms[bestPair[2]].append(rbonds_rAtoms_cAtoms2[1])   
        
        # remove the atoms covered by this match from the available atoms 
        availableAtoms[bestPair[1]] = availableAtoms[bestPair[1]] - rbonds_rAtoms_cAtoms1[2]
        availableAtoms[bestPair[2]] = availableAtoms[bestPair[2]] - rbonds_rAtoms_cAtoms2[2]
        #for k, v in availableAtoms.items(): print(k[0:10], len(v))
        BestPairs.append(bestPair[1:3])
 
    # clean up reacting atoms dict
    reactingAtoms = {k:[y for x in v for y in x ] for k, v in reactingAtoms.items() }
    for k, v in reactingAtoms.items():
        flat = set()
        for x in v:
            if type(x) == set:
                flat = flat | x
            else:
                flat.add(x)
        reactingAtoms[k]=flat 
    
    # RF fragments
    RFs1, RFAs1 = get_RFs(reactingAtoms, Comps1)
    RFs2, RFAs2 = get_RFs(reactingAtoms, Comps2)
    frags1 = [x[2] for x in Comps1]
    frags2 = [x[2] for x in Comps2]
    info1 = [x[3] for x in Comps1]
    info2 = [x[3] for x in Comps2]

    return(reactingAtoms, BestPairs )





def get_RFs(reactingAtoms, Comps1):
    # get the reacting fragments
    comp1RFs={}
    comp1RFAs={}
    for c1 in Comps1:
        comp, atomMap, fpM, info, fragAtoms1, smile, molCount = c1

        reactAtoms = reactingAtoms[smile]
        rf = {k:v for k, v in fragAtoms1.items() if len(v.intersection(reactAtoms))>0  }

        # # separate the smallest fragments from the larger ones
        rf1, rfAdj1 = get_common_atoms_frags(rf, fragAtoms1, reactAtoms)
        
        if len(rf1)==0:
            #print('problem overlapping frags not substrate')
            return(None, None)

        for k, v in rf1.items(): 
            if k not in comp1RFs.keys():
                comp1RFs[k] =v
                
        for k, v in rfAdj1.items(): 
            if k not in comp1RFAs.keys():
                comp1RFAs[k] =v
    return(comp1RFs, comp1RFAs)
   

        

    
    
def getReactingAtomsBonds(Comps1, Comps2, availableAtoms):
    # Comps1, Comps2 = [compDetails[x]+[subs_count[x]] for x in subs], [compDetails[x]+[prods_count[x]] for x in prods]
    # get the best scoring mcs
    bestScore=0
    bestPair = []

    for i1, x in enumerate(Comps1):
        comp1, atomMap1, fpM1, info1, fragAtoms1, smile1, molCount = x
        available1 = availableAtoms[smile1]
        
        if comp1 is None or len(available1) ==0:
            continue
        
        # for each comp2 get the biggest mcs
        for i2, y in enumerate(Comps2):
            comp2, atomMap2, fpM2, info2, fragAtoms2, smile2, molCount = y
            #if '*' in smile1 and '*' in smile2: crash
            available2 = availableAtoms[smile2]
            
            if comp2 is None or len(available2) ==0:
                continue
           
            # if not all the atoms are available then take them out of the MCS:
            if len(available1) < comp1.GetNumAtoms() or len(available2) < comp2.GetNumAtoms():
                comp1Av = remove_unavailable_atoms(comp1, available1)
                comp2Av = remove_unavailable_atoms(comp2, available2)
                # so get the MCS for the available atoms
                res = Chem.rdFMCS.FindMCS([comp1Av, comp2Av], timeout=10)
                resAtomsBonds = Chem.rdFMCS.FindMCS([comp1Av, comp2Av],  bondCompare= Chem.rdFMCS.BondCompare.CompareOrderExact, timeout=10)

            else:
                res = Chem.rdFMCS.FindMCS([comp1, comp2], timeout=10) 
                resAtomsBonds = Chem.rdFMCS.FindMCS([comp2, comp1],  bondCompare= Chem.rdFMCS.BondCompare.CompareOrderExact, timeout=10)
            
            if res.numAtoms==0: continue

            # gets the matches for the MCS
            bMcs = res.smartsString
            bMcs = Chem.MolFromSmarts(bMcs)
            matches1 = get_mcs_atom_numbers( comp1, bMcs) 
            matches2 = get_mcs_atom_numbers( comp2, bMcs)
            if sum([len(x) for x in matches1])==0 or sum([len(x) for x in matches2])==0: continue
        
            # check for changing bonds 
            if res.numBonds > resAtomsBonds.numBonds:
                # print('inner bonds')
                bondsMatch = 0
                diffBonds, bondTypeChanged1, bondTypeChanged2, rAtoms1, rAtoms2 = changedBondsToRA(comp1, comp2, bMcs, atomMap1, atomMap2)
                
                if res.numAtoms - len(rAtoms1) > resAtomsBonds.numAtoms:
                    res = resAtomsBonds
                    bondsMatch = 1
                
                # these bonds that have changed make up the reacting atoms
                # reactingAtoms2 = {}
                # reactingAtoms2[smile1]= set([x[0] for x in diffAtoms])
                # reactingAtoms2[smile2]= set([x[1] for x in diffAtoms])  
                # print('inner bonds', reactingAtoms2)


            else:
                # print('mcs match up')
                bondsMatch = 0
                bondTypeChanged1 ={} 
                bondTypeChanged2 = {}
        
       
            # get the bonds
            bondsMCS =   [ [x.GetBeginAtomIdx(), x.GetEndAtomIdx(), str(x.GetBondType())] for x in bMcs.GetBonds() ]
            bondsComp1 = [ [x.GetBeginAtomIdx(), x.GetEndAtomIdx(), str(x.GetBondType())] for x in comp1.GetBonds()]
            bondsComp2 = [ [x.GetBeginAtomIdx(), x.GetEndAtomIdx(), str(x.GetBondType())] for x in comp2.GetBonds()]
            
            # for each atom get all adjacent bonds
            bondsM = make_bonds_dict(bondsMCS)[1]
            bonds1 = make_bonds_dict(bondsComp1)[1]
            bonds2 = make_bonds_dict(bondsComp2)[1]
            
            # compare each compound to the MCS, score is atoms covered in mcs/ bonds broken +1
            score1, reactingBonds1, reactingAtoms1, coveredAtoms1, match1 = best_msc_RAs(matches1, available1, bonds1, bondsM, bondTypeChanged1)
            score2, reactingBonds2, reactingAtoms2, coveredAtoms2, match2 = best_msc_RAs(matches2, available2, bonds2, bondsM, bondTypeChanged2)
            if score1 == None or score2 == None: continue

            
            # get the missing reacting atoms, so that both compounds have the same number of reacting atoms
            missing_ra1, missing_ra2 = get_boundary_atoms(bondsM, bonds1, bonds2, match1, match2)
            reactingAtoms1 = reactingAtoms1 + list(missing_ra1)
            reactingAtoms2 = reactingAtoms2 + list(missing_ra2)
        
            if min([score1, score2])>bestScore:
                bestScore = min([score1, score2])
                rbonds_rAtoms_cAtoms1 = [reactingBonds1, reactingAtoms1, coveredAtoms1]
                rbonds_rAtoms_cAtoms2 = [reactingBonds2, reactingAtoms2, coveredAtoms2]
                bestPair = [bestScore, smile1, smile2, i1, i2,  comp1, comp2]
       
         
    if bestScore >0:    
        return(bestPair, rbonds_rAtoms_cAtoms1, rbonds_rAtoms_cAtoms2 )
    else:
        return(None, None, None)
                


def remove_unavailable_atoms(comp, available):
    # remove unavailable atoms
    
    if len(available) == comp.GetNumAtoms():
        return (comp)
                   
    removeAtoms = [ x for x in list(range(0, comp.GetNumAtoms())) if x not in available]
    edit_mol = Chem.EditableMol(comp)
    for x in sorted(removeAtoms, reverse=True):
        edit_mol.RemoveAtom(x)
        
    return ( edit_mol.GetMol() )


def get_boundary_atoms(bondsM, bonds1, bonds2, match1, match2):
    
    subset1 = set()
    subset2 = set()
    for m, cp1 in match1.items():
        cp2 = match2[m]
        
        # if this atom in compound1 has more bonds than MCS, but compound 2 doesn't
        if len(bonds1[cp1]) >  len(bondsM[m]) and len(bonds2[cp2]) == len(bondsM[m]):
            subset2.add(cp2)
            
        if len(bonds2[cp2]) >  len(bondsM[m]) and len(bonds1[cp1]) == len(bondsM[m]):
            subset1.add(cp1)
    
    return(subset1, subset2)
            

def get_mcs_atom_numbers(comp1, bMcs):
    MCSatoms = list(range(0, bMcs.GetNumAtoms()))
    MCSatomDict=[]
    
    for substruct in comp1.GetSubstructMatches(bMcs):
        d={} 
        for i, a in enumerate(MCSatoms):
            #if substruct[i] in available1:
            d[a] = substruct[i]
        MCSatomDict.append(d)
    return(MCSatomDict)

def make_bonds_dict(bondsComp2):
    d={}
    d2={}
    for x in bondsComp2:
        a1, a2, bondtype = x
        if a1 not in d.keys():
            d[a1]={'SINGLE':set(), 'DOUBLE':set(), 'AROMATIC':set(), 'TRIPLE':set()}
            d2[a1] = set()
        if a2 not in d.keys():
            d[a2]={'SINGLE':set(), 'DOUBLE':set(), 'AROMATIC':set(), 'TRIPLE':set()}
            d2[a2] = set()
            
        d[a1][bondtype].add(a2)
        d[a2][bondtype].add(a1)
        d2[a1].add(a2)
        d2[a2].add(a1)
    return(d, d2)               


        
def best_msc_RAs(matches, available, bonds, bondsM, bondTypeChanged):
# def best_msc_hits(matches, available, bonds, bondsM):
    ### find missing bonds between the compound and the MCS - use whole compounds, available atoms are only used for scoring
    ### if there are multiple matches choose the best one
    
    # for each matching MCS substructure find the following
    match_scores = [] #find the best of multiple matches
    match_r_bonds =[]
    match_r_atoms = []
    match_c_atoms = []
    match_match = []
    
    # find in the missing bonds for each substruct match
    # print('\nmatch')
    for match in matches:
        # Generate the atom match score
        av1 = set(match.values()).intersection(available)
        if len(av1)==0: continue
        
        # for each atom in the substructure
        r_bonds=[]
        r_atoms=[]
        for M_atom, M_bonds in bondsM.items():
            ## compare the mcs bonds and the full compound bonds using full compound numbering
            
            # look up the atom in the full compound and get the adjacent atoms
            comp1_adj = bonds[match[M_atom]]
            # get the adjacent atoms in the substructure
            mcs_adj = set([match[x] for x in M_bonds])
            
            if comp1_adj != mcs_adj: 
                # get the lost bonds
                lost_bonds = comp1_adj - mcs_adj
                
                # turn missing adjacent atoms into full bonds and store
                lost_bonds2 = [[match[M_atom], x] for x in lost_bonds ]
                r_bonds= r_bonds + lost_bonds2
                
                # save reacting atoms
                r_atoms.append(set([ y for x in lost_bonds2 for y in x])) 
        
        # deal with the bonds that changed type
        for diffBond in bondTypeChanged:
            if len(set(diffBond).intersection(av1)) ==2:
                r_bonds = r_bonds + [list(diffBond)]
                r_atoms = r_atoms + [set(diffBond)]

        match_scores.append(len(av1)/(len(r_bonds) +1))
        match_r_bonds.append(r_bonds)
        match_r_atoms.append(r_atoms)
        match_c_atoms.append(av1)

    if len(match_scores)>1:
        ind_best = match_scores.index(max(match_scores))
        match_scores =  [match_scores [ind_best]]
        match_r_bonds = [match_r_bonds[ind_best]]
        match_r_atoms = [match_r_atoms[ind_best]]
        match_c_atoms = [match_c_atoms[ind_best]]
        match_match = [matches[ind_best]]
    else:
        match_match = [matches[0]]
        
    if len(match_scores)>0:
        return(match_scores[0], match_r_bonds[0], match_r_atoms[0], match_c_atoms[0], match_match[0])
    else:
        return(None, None, None, None)

def scoreReaction(match_c_atoms, match_r_bonds):
    # # if there was multiple hits then find the best one
    match_scores = []
    for i, av1 in enumerate(match_c_atoms):
        r_bonds = match_r_bonds[i]
        match_scores.append(len(av1)/(len(r_bonds)+1))
        
    # if len(match_scores)>1:
    #     ind_best = match_scores.index(max(match_scores))
    #     match_scores =  [match_scores [ind_best]]
    #     match_r_bonds = [match_r_bonds[ind_best]]
    #     match_r_atoms = [match_r_atoms[ind_best]]
    #     match_c_atoms = [match_c_atoms[ind_best]]
    #     match_match = [matches[ind_best]]
        
    # if len(match_scores)>0:
    #     return(match_scores[0], match_r_bonds[0], match_r_atoms[0], match_c_atoms[0], matches[0])
    # else:
    #     return(None, None, None, None, None)

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

        fragDist = mergeDicts(fragDist, fd)
        bi1_filtered = mergeDicts(bi1_filtered, bi)
                
    #bi1_filtered={k:tuple(v) for k, v in bi1_filtered.items()}
     
    # sort fragDist by size
    #fragDist = {k: v for k, v in sorted(fragDist.items(), key=lambda item: item[1])}

    return(bi1_filtered, fragDist)


def mapSubstructs(comp1, bMcs):
    mapAt= comp1.GetSubstructMatches(bMcs) 
    mapAtoms = {}
    for mapping in mapAt:
        for i, x in enumerate(mapping):
            if x not in mapAtoms.keys():
                mapAtoms[x]=set()
            mapAtoms[x].add(i) 
    return(mapAtoms)

def getBondType(comp1, mapAtom1, bMcs):
    b = {x:list() for x in list(range(0,bMcs.GetNumAtoms()))}
    
    if mapAtom1:
        for bondComp in comp1.GetBonds(): 
            startAtom = bondComp.GetBeginAtomIdx()
            endAtom   = bondComp.GetEndAtomIdx()
            if startAtom in mapAtom1 and endAtom in mapAtom1:
                for y in mapAtom1[startAtom]:
                    b[y].append(str(bondComp.GetBondType()))
                for y in mapAtom1[endAtom]:
                    b[y].append(str(bondComp.GetBondType()))
            
    else:
        for x in comp1.GetBonds(): 
            b[x.GetBeginAtomIdx()].append(str(x.GetBondType()))
            b[x.GetEndAtomIdx()].append(str(x.GetBondType())) 
    return(b)

def changedBondsToRA(comp1, comp2, bMcs, atomMap1, atomMap2):
    mapAtom1 =mapSubstructs(comp1, bMcs)
    mapAtom2 =mapSubstructs(comp2, bMcs)  

    # get the bond type and allocate them into the best match 
    b1 = getBondType(comp1, mapAtom1, bMcs)
    b2 = getBondType(comp2, mapAtom2, bMcs)
    
    diff = {k1:[v1, b2[k1]] for k1, v1 in b1.items() if sorted(v1) != sorted(b2[k1]) }
    
    # map the compounds back to comp1 & comp2
    rAtoms1 = []
    rAtoms2 = []
    for k, v in diff.items():
        for cpd1, mp1 in  mapAtom1.items():
            if k in mp1:
                rAtoms1.append(cpd1)  #= v[0]
                continue
        for cpd2, mp2 in  mapAtom2.items():
            if k in mp2:
                rAtoms2.append(cpd2)  #= v[0]
                continue
    
    rBonds1 = set()
    for x in rAtoms1:
        rbonds = [tuple(sorted([x, y])) for y in atomMap1[x] if y in rAtoms1]
        # MNXR100083 due to a missing bond - linear to ring formation
        if not rbonds: print('problem atom', x)
        rBonds1 = rBonds1 | set(rbonds)

    rBonds2 = set()
    for x in rAtoms2:
        rbonds = [tuple(sorted([x, y])) for y in atomMap2[x] if y in rAtoms2]
        if not rbonds: print('problem atom', x)
        rBonds2 = rBonds2 | set(rbonds)        
    
                    
    return diff, rBonds1, rBonds2, rAtoms1, rAtoms2




#comps1 = ['C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC=C(C=C3)O']
#comps2 =  ['C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=C(C=C3)O)O']

#rfrags, adj, fps, sucesses = get_MCS2(comps1, comps2)
