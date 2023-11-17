#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 10:06:52 2022

@author: ruth
"""
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import numpy as np

def neutralize_atoms(mol):
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

def rot_ar_y(radi):
    return  np.array([[np.cos(radi), 0, np.sin(radi), 0],
                      [0, 1, 0, 0],
                      [-np.sin(radi), 0, np.cos(radi), 0],
                     [0, 0, 0, 1]], dtype=np.double)

#tyr = Chem.MolFromSmiles('C1=CC(=CC=C1CC(C(=O)O)N)O')
#ldopa = Chem.MolFromSmiles('C1=CC(=C(C=C1CC(C(=O)O)N)O)O')

for mol_names, smiles in comps.items():
    mol1 = Chem.MolFromSmiles(smiles[0])
    mol2 = Chem.MolFromSmiles(smiles[1]) 
    
    # 
    bi={}
    neutralize_atoms(smiles[0])
    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 3, bitInfo=bi )
    # prints = [(mol1, x, bi) for x in fp1.GetOnBits()]
    # Draw.DrawMorganBits(prints, molsPerRow=4, legends=[str(x) for x in fp1.GetOnBits()])
    # hits=[1607, 589, 310, 1754,  585, 1427]
    # prints = [(mol1, x, bi) for x in hits]
    # Draw.DrawMorganBits(prints, molsPerRow=3, legends=[str(x) for x in hits])
    
    
    bi={}
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 3, bitInfo=bi )
    # prints = [(mol2, x, bi) for x in fp2.GetOnBits()]
    # Draw.DrawMorganBits(prints, molsPerRow=4, legends=[str(x) for x in fp2.GetOnBits()])
    # hits=[ 1778, 1313, 310, 1754, 585,  1427]
    # prints = [(tyr, x, bi) for x in hits]
    # Draw.DrawMorganBits(prints, molsPerRow=3, legends=[str(x) for x in hits])
    
    
    # reacting bonds
    hit_ats = list(ldopa.GetSubstructMatch(tyr))
    hit_ats2 = []
    hit_bonds = []
    for bond in ldopa.GetBonds():
        aid1 = bond.GetBeginAtomIdx()
        aid2 = bond.GetEndAtomIdx()
        #if aid1==12 or aid2==12:crash
        if len(set([aid1, aid2]).intersection(set(hit_ats)))==1:
               hit_bonds.append(bond.GetIdx())
               hit_ats2.append(aid1)
               hit_ats2.append(aid2)
    
    colours = [(0,0,0.8)]*4
    bond_cols = {}
    for i, bd in enumerate(hit_bonds):
        bond_cols[bd] = colours[3 - i%4]
        
    atom_cols = {}
    for i, at in enumerate(hit_ats2): 
        atom_cols[at] = colours[i%4]
    
    d = rdMolDraw2D.MolDraw2DCairo(500, 500) # or MolDraw2DCairo to get PNGs
    rdMolDraw2D.PrepareAndDrawMolecule(d, ldopa, highlightAtoms=hit_ats2, highlightAtomColors=atom_cols,
                                       highlightBonds=hit_bonds, highlightBondColors=bond_cols)
    
    d.WriteDrawingText("highlight_changing_bonds/Ldopa_RB.png")
    
    # tyrosine
    
    
    d = rdMolDraw2D.MolDraw2DCairo(500, 500) # or MolDraw2DCairo to get PNGs
    rdMolDraw2D.PrepareAndDrawMolecule(d, tyr, highlightAtoms=[3], highlightAtomColors=atom_cols,
                                       highlightBonds=[], highlightBondColors=bond_cols)
    
    d.WriteDrawingText("highlight_changing_bonds/tyrosine_RB.png")