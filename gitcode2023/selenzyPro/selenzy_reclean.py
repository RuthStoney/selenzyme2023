#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 10:51:50 2023

@author: ruth
"""


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 10:53:03 2017

@author: jerrywzy, Pablo Carbonell
"""


import re
import os, subprocess
# import json
import csv
import argparse
#import quickRsim
#import quickRsim_mod3 as quickRsim
#import quickRsim_mod7 as quickRsim
import quickRsim_local2 as quickRsim
# import quickRsim_local3 as quickRsim
import numpy as np
import pandas as pd
from rdkit.Chem import AllChem, Draw
from datetime import datetime
from rxnmapper import RXNMapper


class preLoad(object):
    """ Container of precomputed data """
    def __init__(self):
        pass

    def fasta(self, datadir, ffile="seqs.fasta"):
        (sequence, names, descriptions, fulldescriptions, osource, pexistence) = readFasta(datadir, ffile)
        self.sequence = sequence
        self.names = names
        self.descriptions = descriptions
        self.fulldescriptions = fulldescriptions
        self.osource = osource
        self.pexistence = pexistence

    def fpData(self, datadir):
        with open(datadir+'ruthsprint.txt', 'a+') as handler:
            handler.write("in fpData" + str(quickRsim.fingerprint()) + "\n")
        with open(datadir+'ruthsprint.txt', 'a+') as handler:
            handler.write("in fpData" + datadir + "\n")

        self.fp = {}    
        for fpid in quickRsim.fingerprint():
            self.fp[fpid] = quickRsim.loadFingerprint(datadir, fpid)

        with open(datadir+'ruthsprint.txt', 'a+') as handler:
            handler.write("in fpData" + datadir + 'FP_MorgRF' + "\n")

        self.fpr = {}    
        for fpid in quickRsim.fingerprint():
            self.fpr[fpid] = quickRsim.loadFingerprint3(datadir, 'FP_MorgRF')
            
        # fp, fpn, fpp, fpfun = loadFingerprint(arg.datadir, arg.fp)
        # fpRF, fpnRF, fprRF, rfDict, rfdist = loadFingerprint3(arg.datadir, 'FP_MorgRF')

        with open(datadir+'ruthsprint.txt', 'a+') as handler:
            handler.write("in fpData done\n")

    def seqData(self, datadir, fl):        

        with open(datadir+'ruthsprint.txt', 'a+') as handler:
            handler.write("making MnxToUprot using " + str(fl[0]) + "\n")
        with open(os.path.join(datadir, fl[0])) as f:
            self.MnxToUprot = {}
            self.UprotToMnx = {}
            for row in f:
                mnxr, db, seqid, source, ec = row.split('\t')
                if db == 'uniprot':
                    if mnxr not in self.MnxToUprot:
                        self.MnxToUprot[mnxr] = set()
                    self.MnxToUprot[mnxr].add(seqid)
                    if seqid not in self.UprotToMnx:
                        self.UprotToMnx[seqid] = set()
                    self.UprotToMnx[seqid].add(mnxr)
        
        # with open(os.path.join(datadir, fl[1])) as f2:
        #     self.upclst = json.load(f2)
        
        # with open(os.path.join(datadir, fl[2])) as f3:
        #     self.clstrep= json.load(f3)

        self.seqorg = seqOrganism(datadir, fl[3])
        self.tax = readTaxonomy(datadir, fl[4])

        # self.seqorg = seqOrganism(datadir, fl[1])
        # self.tax = readTaxonomy(datadir, fl[2])

    def reacData(self, datadir, smf):
    #     """ Transitional downgrade mapping for mnx v3 to v2.
    #     In the future everything will be based on v3.
    #     """
        smiFile = os.path.join(datadir, smf[0])
        rxnRefFile = os.path.join(datadir, smf[1])
    #     rxnConsensus = os.path.join(datadir, smf[2])
        rxnProp = os.path.join(datadir, smf[3])
    #     rxnBrenda = os.path.join(datadir, smf[4])
    #     rxnSabiork = os.path.join(datadir, smf[5])
    #     rxnv3 = None
    #     if len(smf) > 6:
    #         rxnv3 = os.path.join(datadir, smf[6])
        self.smir = {}
        if os.path.exists(smiFile):
            self.smir = reactionSmiles(smiFile)
        self.rxnref = {}
        if os.path.exists(rxnRefFile):
            self.rxnref = reactionXref(rxnRefFile) #, rxnBrenda, rxnSabiork,rxnv3)
    #     self.rxndir = {}
    #     if os.path.exists(rxnConsensus):
    #         self.rxndir = readRxnCons(rxnConsensus)
        self.ecrxn = {}
        self.rxnec = {}
        if os.path.exists(rxnProp):
            self.ecrxn, self.rxnec = readRxnProp(rxnProp)
        self.ecsmi = ecSmiles(self.ecrxn, self.smir, self.rxnref)


def readData(datadir):
    with open(datadir+'ruthsprint.txt', 'a+') as handler:
        handler.write("in readdata\n")

    """ Read all data into memory """
    pc = preLoad()

    pc.fasta(datadir, 'seqs.fasta')
    pc.fpData(datadir)
    pc.seqData(datadir, ['reac_seqs.tsv', '_', '_', "seq_org.tsv", "org_lineage.csv"])
    # pc.reacData(datadir, ['reac_smi.csv', '', '', 'reac_prop.tsv'])
    pc.reacData(datadir, ['reac_smi.csv', 'reac_xref.tsv', '', 'reac_prop.tsv'])
    # pc.reacData(datadir, ['reac_smi.csv','reac_xref.tsv',"rxn_consensus_20160612.txt",
    #                       'reac_prop.tsv']) #, 'brenda-mnxref2.tsv', 'sabiork-mnxref2.tsv', 'reac_xref_v3.tsv'])
    return pc

def availableFingerprints():
    fp = quickRsim.fingerprint()

    return sorted(fp, key=lambda x: fp[x][4])

def sanitizeRxn(rxninfo, outrxn):
    """ It works both with the smiles string or a rxn file """
    try:
        if os.path.exists(rxninfo):
            if open(rxninfo).readline().startswith('$RXN'):
                rxn = AllChem.ReactionFromRxnFile(rxninfo)
            else:
                smarts =  open(rxninfo).readline().rstrip()
                try:
                    rxn = AllChem.ReactionFromSmarts(smarts, useSmiles=True)
                except:
                    rxn = AllChem.ReactionFromSmarts(smarts)
        smi = AllChem.ReactionToSmiles(rxn)
        with open(outrxn+'.smi', 'w') as handler:
            handler.write(smi)
        mdl = AllChem.ReactionToRxnBlock(rxn)
        with open(outrxn+'.rxn', 'w') as handler:
            handler.write(mdl)
        return smi
    except:
        return ''

def sanitizeSmarts(smarts, outrxn):
    """ It works both with the smiles string or a rxn file """
    try:
        try:
            rxn = AllChem.ReactionFromSmarts(smarts, useSmiles=True)
        except:
            rxn = AllChem.ReactionFromSmarts(smarts)
        smi = AllChem.ReactionToSmiles(rxn)
        with open(outrxn+'.smi', 'w') as handler:
            handler.write(smi)
        mdl = AllChem.ReactionToRxnBlock(rxn)
        with open(outrxn+'.rxn', 'w') as handler:
            handler.write(mdl)
        return smi
    except:
        return ''



def display_reaction(rxninfo, outfolder, outname, marvin=False):
    """ It works both with the smiles string or a rxn file """
    if marvin:
        try:
            cmd = ['molconvert', 'svg:w500', rxninfo]
            job = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = job.communicate()
            outimage = os.path.join(outfolder, outname+'.svg')
            return outimage, (600, 400)
        except:
            return '', (0,0)
    else:
        try:
            if os.path.exists(rxninfo):
                if open(rxninfo).readline().startswith('$RXN'):
                    rxn = AllChem.ReactionFromRxnFile(rxninfo)
                else:
                    smarts =  open(rxninfo).readline()
                    rxn = AllChem.ReactionFromSmarts(smarts)
                outimage = os.path.join(outfolder, outname+'.png')
                im = Draw.ReactionToImage(rxn)
                size = im.size
                im.save(outimage)
                return outimage, size
            else:
                rxn = AllChem.ReactionFromSmarts(rxninfo)
                outimage = os.path.join(outfolder, outname+'.png')
                im = Draw.ReactionToImage(rxn)
                size = im.size
                im.save(outimage)
                return outimage, size
        except:
            return '', (0,0)

def seqOrganism(datadir, fileSeqOrg):
    with open(datadir+'ruthsprint.txt', 'w+') as handler:
        handler.write("in seqOrganism in " + datadir + " reading" + fileSeqOrg +"\n") 

    seqorg = {}
    for line in open(os.path.join(datadir, fileSeqOrg)):
        row = line.rstrip().split('\t')
        seqorg[ row[0] ] = (row[1], row[2])
    with open(datadir+'ruthsprint.txt', 'a+') as handler:
        handler.write("out seqOrganism\n") 
    return seqorg

def readTaxonomy(datadir, fileLineage):
    tax = {}
    with open(os.path.join(datadir, fileLineage)) as handler:
        for row in csv.reader(handler):
            tax[row[0]] = row
    return tax

def taxDistance(tax, host, target):
    # tax, host, target = tax, host, seqorg[y][0]
    if host in tax and target in tax:
        hostLineage = set(tax[host])
        targetLineage = set(tax[target])
        distance = 1 + len(hostLineage ^ targetLineage)
        return distance
    else:
        return '-'

def seqScore(newscore=None):
    print('in seqScore')
    import string
    # Initial score
    # vdict = {
    #     string.ascii_uppercase[9]: ('Reaction similarity:', 50.0, True),
    #     string.ascii_uppercase[12]: ('Combined score', 0.0, True),
    #     string.ascii_uppercase[10]: ('RF similarity', 50.0, False),
    #     string.ascii_uppercase[8]: ('Sequence <a href="https://en.wikipedia.org/wiki/Conserved_sequence" target="_blank">conservation</a>:', 1.0, True),
    #     string.ascii_uppercase[4]: ('Sequence <a href="https://www.ncbi.nlm.nih.gov/taxonomy" target="_blank">taxonomic</a> distance:', -1.0, True),
    #     string.ascii_uppercase[7]: ('Uniprot <a href="http://www.uniprot.org/help/protein_existence" target="_blank">protein evidence</a>:', -0.1, True),      
    #     string.ascii_uppercase[14]: ('Percentage sheets:', 0.0, False),
    #     string.ascii_uppercase[15]: ('Percentage turns:', 0.0, False),
    #     string.ascii_uppercase[16]: ('Molecular weight:', 0.0, False),
    #     string.ascii_uppercase[17]: ('Isoelectric point:', 0.0, False)
    #     #string.ascii_uppercase[18]: ('Percentage polar amino acids:', 0.0, False)
    # }
    vdict = {
        string.ascii_uppercase[10]: ('Reaction similarity:', 0.0, False),
        string.ascii_uppercase[13]: ('Combined score', 50.0, True),
        string.ascii_uppercase[11]: ('RF similarity', 0.0, False),
        string.ascii_uppercase[9]: ('Sequence <a href="https://en.wikipedia.org/wiki/Conserved_sequence" target="_blank">conservation</a>:', 0.0, False),
        string.ascii_uppercase[4]: ('Sequence <a href="https://www.ncbi.nlm.nih.gov/taxonomy" target="_blank">taxonomic</a> distance:', -0.05, True),
        string.ascii_uppercase[8]: ('Uniprot <a href="http://www.uniprot.org/help/protein_existence" target="_blank">protein evidence</a>:', 0.0, False),      
        string.ascii_uppercase[15]: ('Percentage sheets:', 0.0, False),
        string.ascii_uppercase[16]: ('Percentage turns:', 0.0, False),
        string.ascii_uppercase[18]: ('Molecular weight:', 0.0, False),
        string.ascii_uppercase[19]: ('Isoelectric point:', 0.0, False)
        #string.ascii_uppercase[18]: ('Percentage polar amino acids:', 0.0, False)
    }
    nvdict = vdict.copy()
    # Reference order (alternatively, perhaps easier to keep table order)
    clist = [string.ascii_uppercase[x] for x in [10, 13, 11, 9, 4, 8, 15, 16, 18, 19]]
    # Update score if given and well-formed
    update = False
    if newscore is not None:
        nscore = [x[0] for x in newscore]
        if len( set(nscore) & set(clist) ) > 0:
            for val in vdict:
                nvdict[val] = (vdict[val][0], vdict[val][1], False)
            for val in newscore:
                if val[0] in nvdict:
                    try:
                        newval = float(val[1])
                        nvdict[val[0]] = (nvdict[val[0]][0], newval, True)
                        update = True
                    except:
                        continue
    if not update:
        nvdict = vdict
    score = []
    for x in clist:
        score.append( (x,) + nvdict[x] )
    return score

def updateScore(csvfile, score):
    """ Add or update score column and reorder """
    print('\nin update score')
    import string
    head, rows = read_csv(csvfile)
    data = pd.read_csv(csvfile)
    data.index = data.index + 1
    cols = data.columns.tolist()
    sco = pd.Series(np.zeros(len(data[cols[0]])), index=data.index)
    if 'Score' not in cols:
        data['Score'] = sco
        cols = ['Score'] + cols
        data = data[cols]
    colk = list(string.ascii_uppercase)
    
    for sc in score:
        try:
            coln = colk.index(sc[0])            
            val = sc[2]
            checked = sc[3]
            if checked:
                print('coln', coln, head[coln], '::::', sc[1])
                sco += val * data.iloc[:,coln]
        except:
            continue
        
        
    data['Score'] = sco
    data = data.sort_values('Score', ascending=False)
    updateMSA(os.path.dirname(csvfile), [[v] for v in data['Seq. ID']])
    data = data.reset_index(drop=True)
    data.index = data.index + 1
    data.rename_axis('Select', axis="columns")
    #data.to_csv(csvfile.replace('.csv', '_UDscore.csv'), quoting=csv.QUOTE_ALL, index=False)
    data.to_csv(csvfile, quoting=csv.QUOTE_ALL, index=False)
    print('score calculated\n')
    return data


def readFasta(datadir, fileFasta, limit=None):
    
    from Bio import SeqIO
    
    sequence = {}
    descriptions = {}
    fulldescriptions = {}
    osource = {}
    pexistence = {}
    names = []
    seen = set()
    seen_add = seen.add
    
    for seq_record in SeqIO.parse(os.path.join(datadir, fileFasta), "fasta"):
        ming = seq_record.id
        try:
            idonly = re.search(r'\|(.*?)\|',ming)
            x = idonly.group(1)
        except:
            x = ming
        fulldesc = seq_record.description
        desc = fulldesc.rsplit('OS=')[0]
        try:
            shortdesc = " ".join(desc.split()[1:])
        except:
            shortdesc = desc
        try:
            orgsource = fulldesc.rsplit('OS=')[1]
        except:
            orgsource = '-'
        shortos = orgsource.rsplit('GN=')[0]
        if ',' in shortdesc:
            y = shortdesc.replace(",", ";")
        else:
            y = shortdesc
        pe = orgsource.rsplit('PE=')[1].split(' ')[0]
        if x not in seen:
            names.append(x)
            seen_add(x)
        myseq = seq_record.seq
        sequence[x]=(myseq)
        descriptions[x]= y
        osource[x]=shortos
        fulldescriptions[x] = fulldesc
        pexistence[x] = pe
        # A practical limit hardcoded
        if limit is not None:
            try:
                if len(sequence) > limit:
                    break
            except:
                continue
    
    return (sequence, names, descriptions, fulldescriptions, osource, pexistence)
 
def readRxnCons(consensus):
    
    f = open(consensus, 'r')
    
    MnxDir = {}
    
    for line in f:
        splitdata = line.split()
        Mnx = splitdata[1]
        dirxn = splitdata[2]
        MnxDir[Mnx] = dirxn
        
    return (MnxDir)   
    
def getMnxSim(rxnInput, datadir, outdir, rxnm, drxn=0, fp='Morgan', pc=None):
    """ Commmand line arguments of quickRsim """


    args = [datadir, fp] + rxnInput + ['-out', os.path.join(outdir,'results_quickRsim.txt')]

    quickRsim.run( quickRsim.arguments(args), pc, rxnm )
    MnxSim = {}
    print('testing consensus file')
    
    try:
        MnxDirPref = readRxnCons(os.path.join(datadir, "rxn_consensus_20160612.txt"))
        MnxDirUsed = {}
    except:
        MnxDirPref ={}
        MnxDirUsed = {}
        with open(outdir+'_ruthsprint.txt', 'a+') as handler:
            handler.write("no consensus file\n")    
        print('no consensus file')
    EcNumber = {}
    MnxSimRF ={}
    SumScore ={}

    
    if drxn==1:
        print('in drxn == 1')
        with open(outdir+'_ruthsprint.txt', 'a+') as handler:
            handler.write("in getMnxSim drxn == 1\n")

        fileout = open(os.path.join(outdir, "/results_quickRsim.txt"), 'r')
        for line in fileout:
            splitdata = line.split()
            Mnx = splitdata[1]
            S1 = splitdata[2]
            S2 = splitdata[3]
            SMILES = splitdata[4]
            MnxSimRF[Mnx] = splitdata[5]
            EC = ''
            if len(splitdata) > 6:
                EC = splitdata[6]
            EcNumber[Mnx] = EC
            try:
                direction = MnxDirPref[Mnx] 
                if direction == '1':
                    MnxSim[Mnx]=S1
                    MnxDirUsed[Mnx]= '1'
                elif direction == '-1':
                    MnxSim[Mnx]=S2
                    MnxDirUsed[Mnx]= '-1'
                else:
                    if float(S1) >= float(S2):
                        MnxSim[Mnx]=S1
                        MnxDirUsed[Mnx]= '1'
                    elif  float(S2) > float(S1):
                        MnxSim[Mnx]=S2
                        MnxDirUsed[Mnx]= '-1'

            except KeyError:
                MnxDirPref[Mnx] = "N/A"   # for missing Keys?
                if S1 >= S2:
                    MnxSim[Mnx]=S1
                else:
                    MnxSim[Mnx]=S2
            try:
                if  1>=float(splitdata[5])>=0:
                    SumScore[Mnx] = (max([float(S1), float(S2)]) + float(splitdata[5]))/2
                else:
                    SumScore[Mnx] =  max([float(S1), float(S2)])
            except:
                SumScore[Mnx] = max([float(S1), float(S2)])       
        fileout.close()
            
        with open(outdir+'_ruthsprint.txt', 'a+') as handler:
            handler.write("exit getMnxSim drxn == 1\n")

        return (MnxSim, MnxDirPref, MnxDirUsed, SMILES, EcNumber, MnxSimRF, SumScore)
        
    else:
        print('in drxn == 0 ')
        with open(outdir+'_ruthsprint.txt', 'a+') as handler:
            handler.write("in getMnxSim drxn == 0\n\n")

        print(os.path.join(outdir, "results_quickRsim.txt"))
        fileout = open(os.path.join(outdir, "results_quickRsim.txt"), 'r')
        for line in fileout:
            splitdata = line.split()
            Mnx = splitdata[1]
            S1 = splitdata[2]
            S2 = splitdata[3]
            SMILES = splitdata[4]
            MnxSimRF[Mnx] = splitdata[5]
            EC = ''
            
            if len(splitdata) > 6:
                EC = splitdata[6]
            EcNumber[Mnx] = EC
            
            if float(S1) > float(S2):
                MnxDirUsed[Mnx]='1'
                MnxSim[Mnx]=S1  
            
            elif float(S2) > float(S1):
                MnxDirUsed[Mnx]='-1'
                MnxSim[Mnx]=S2
            elif float(S1) == float(S2):
                MnxDirUsed[Mnx]='0'
                MnxSim[Mnx]=S2   
                
            try:
                if  1>=float(splitdata[5])>=0:
                    SumScore[Mnx] = (max([float(S1), float(S2)]) + float(splitdata[5]))/2
                else:
                    SumScore[Mnx] =  max([float(S1), float(S2)])
            except:
                SumScore[Mnx] = max([float(S1), float(S2)])

        fileout.close()              

        with open(outdir+'_ruthsprint.txt', 'a+') as handler:
            handler.write("exit getMnxSim drxn == 0\n\n")
        
        return (MnxSim, MnxDirPref, MnxDirUsed, SMILES, EcNumber, MnxSimRF, SumScore)

def readRxnProp(rxnprop):
    ecrxn = {}
    rxnec = {}
    with open(rxnprop) as handler:
        for line in handler:
            if line.startswith('#'):
                continue
            row = line.rstrip().split('\t')
            rxnid = row[0]
            # ec = row[4].split(';')
            ec = row[3].split(';')
            if len(ec) > 0:
                for e in ec:
                    if e == '':
                        continue
                    if e not in ecrxn:
                        ecrxn[e] = set()
                    ecrxn[e].add(rxnid)
                    if rxnid not in rxnec:
                        rxnec[rxnid] = set()
                    rxnec[rxnid].add(e)
    return ecrxn, rxnec

def reactionXref(rxnRefFile, rxnBrenda=None, rxnSabiork=None, rxnv3=None):
    rxnref = {}
    for xref in (rxnRefFile, rxnBrenda, rxnSabiork):
        if xref is not None:
            with open(xref) as handler:
                for line in handler:
                    if line.startswith('#'):
                        continue
                    row = line.rstrip().split('\t')
                    rxnref[row[0]] = row[1]
                    rxnref['mnx:'+row[1]] = row[1]
    # Downgrade v3 to v2
    if rxnv3 is not None:
        for line in open(rxnv3):
            if line.startswith('#'):
                continue
            row = line.rstrip().split('\t')
            try:
                db, rid = row[0].split(':')
                if db == 'deprecated':
                    rxnref['mnx:'+row[1]] = rid
            except:
                continue
    return rxnref

def ecSmiles(ecrxn, rsmi, rxnref):
    smiec = {}
    for ec in ecrxn:
        for r in set(ecrxn[ec]) & set(rsmi):
            s = rsmi[r][0]
            if s not in smiec:
                smiec[s] = set()
            smiec[s].add(ec)
    for ec in ecrxn:
        """ Take the shortest SMARTS being unique for the EC class if possible """
        hits = sorted(set(ecrxn[ec]) & set(rsmi), key=lambda x: len(rsmi[x][0]))
        for h in hits:
            if len(rsmi[h]) > 0:
                rxnid = 'ec:'+ec
                rxnref[rxnid] = h
                if len(smiec[rsmi[h][0]]) == 1:
                       break            


def reactionSmiles(rxnSmilesFile):
    rsmi = {}
    with open(rxnSmilesFile) as handler:
        header = handler.readline()
        for line in handler:
            row = line.rstrip().split(',')
            if len(row) > 1:
                rid = row[0]
                rs1 = row[1]
                lr = rs1.split('>>')
                rs2 = lr[1]+'>>'+lr[0]
                rsmi[rid] = (rs1, rs2)
                
    return rsmi
   
def pepstats(file, outdir):
    print('pepstats in')
    outfile = os.path.join(outdir, "results.pepstats")
    args = ("pepstats -sequence {0} -outfile ".format(file) + outfile)
    os.system(args)
    
    f = open(outfile, "r")
    
    hydrop = {}
    weight = {}
    isoelec = {}
    polar = {}
    
    for line in f:
        if "PEPSTATS of" in line:
            splitdata = line.split()
            seq = splitdata[2].split('_')[0]
        elif "Molecular weight = " in line:
            splitdata = line.split()
            w = splitdata[3]
            weight[seq] = w
        elif "Isoelectric Point = " in line:
            splitdata = line.split()
            i = splitdata[3]
            isoelec[seq] = i
        elif "Polar	" in line:
            splitdata = line.split()
            percent = splitdata[3]
            polar[seq] = percent
            seq
    print('pepstats out')
    return (hydrop, weight, isoelec, polar)

def noAmbiguousSeqs(infile, outfile):
    print('preload Alphabet - deleted')
    """ Remove ambigous amino acid codes """
    from Bio.Data.IUPACData import protein_letters_1to3, extended_protein_values
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    #from Bio.Alphabet import generic_protein
    print('loaded Alphabet ')

    newrecords = []
    for record in SeqIO.parse(infile, "fasta"):
        newseq = ''
        for aa in record.seq:
            if aa not in protein_letters_1to3:
                newseq += extended_protein_values[aa][0]
            else:
                newseq += aa
        newrecords.append( SeqRecord(Seq(newseq), id = record.id, description = record.description) )
    SeqIO.write(newrecords, outfile, "fasta")
    print('done noAmbiguousSeqs')
    
def garnier(file, outdir):
    print('garnier in')
    fixfile = file+'.fix.fasta'
    noAmbiguousSeqs(file, fixfile)
    outfile = os.path.join(outdir, "garnier.txt")
    
    args = ("garnier -sequence {0} -outfile ".format(fixfile) + outfile)
    os.system(args)
    
    f = open(outfile, "r")

    helices = {}
    sheets = {}
    turns = {}
    coils = {}
    
    print('garnier loop')
    for line in f:
        if "Sequence:" in line:
            splitdata = line.split()
            seq = splitdata[2].split('_')[0]
        elif "percent:" in line:
            percents = line.split()
            h = percents[3]
            e = percents[5]
            t = percents[7]
            c = percents[9]
            helices[seq] = h
            sheets[seq] = e
            turns[seq] = t
            coils[seq] = c

    print('garnier out')
    return (helices, sheets, turns, coils)

def updateMSA(outdir, sortrows):
    with open(outdir+'_ruthsprint.txt', 'a+') as handler:
        handler.write("\nin updateMSA\n") 

    """ Update the fasta MSA file in the order of the scores """
    align_fasta = os.path.join(outdir, "sequences_aln.fasta")
    if not os.path.exists(align_fasta):
        with open(outdir+'_ruthsprint.txt', 'a+') as handler:
            handler.write("align_fasta did not exist\n")
        return
    fasta = {}

    with open(outdir+'_ruthsprint.txt', 'a+') as handler:
        handler.write("going through "  + align_fasta + "\n")

    with open(align_fasta) as handler:
        for line in handler:
            if line.startswith('>'):
                splitdata = line.split()
                upid = splitdata[0][1:]
                if len(upid.split('|')) > 1:
                    upid = upid.split('|')[1]
                elif len(upid.split('_')) > 1:
                    upid = upid.split('_')[0]
            if upid not in fasta:
                fasta[upid] = []

            fasta[upid].append(line)

    if len(fasta) == 0:
        with open(outdir+'_ruthsprint.txt', 'a+') as handler:
            handler.write("len fasta is 0\n")
        return
    outfile = align_fasta

    with open(outdir+'_ruthsprint.txt', 'a+') as handler:
        handler.write("writing " + outfile + "\n")

    with open(outfile, 'w') as hw:
        for row in sortrows:
            try:
                pid = row[0]
                if pid in fasta:
                    for seq in fasta[pid]:
                        hw.write(seq)
                else:
                    with open(outdir+'_ruthsprint.txt', 'a+') as handler:
                        handler.write(pid+" not in "  + str(align_fasta) + "\n")
            except:
                continue

def doMSA(finallistfile, outdir):
    print('doing doMSA')
    outfile = os.path.join(outdir, "sequences.score_ascii")
    outfile_html = os.path.join(outdir, "sequences.score_ascii.score_html")
    outfile_aln = os.path.join(outdir, "sequences.score_ascii.fasta_aln")
    align_html = os.path.join(outdir, "sequences_score.html")
    align_fasta = os.path.join(outdir, "sequences_aln.fasta")
    treefile = os.path.join(outdir, "sequences.dnd")
    args = ("t_coffee -in {0} -mode quickaln -output=score_ascii,fasta_aln,score_html -outfile ".format(finallistfile) +outfile+ " -newtree "+treefile)
    os.system(args)
    if os.path.exists(outfile_html):
        os.rename(outfile_html, align_html)
    if os.path.exists(outfile_aln):
        os.rename(outfile_aln, align_fasta)
    
    f = open(outfile, "r")

    cons = {}
    for line in f:
        if "   :  " in line:
            splitdata = line.split()
            upid = splitdata[0]
            if len(upid.split('|')) > 1:
                   upid = upid.split('|')[1]
            elif len(upid.split('_')) > 1:
                   upid = upid.split('_')[0]
            score = splitdata[2]
            cons[upid] = score

    return cons

def read_csv(csvfile):
    rows = []
    if os.path.exists(csvfile):
        with open(csvfile) as handler:
            cv = csv.reader(handler)
            head = next(cv)
            for row in cv:
                rows.append(row)
    return head, rows

def write_csv(csvfilepath, head, rows):
    with open (csvfilepath, 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL)
        writer.writerow(head)
        for r in rows:
            writer.writerow(r)



def sort_rows(rows, columns, head = ['Seq. ID', 'Description', 'Organism Source', 'Tax. distance', 'Rxn. ID','Rxn. ID - additional', 'EC Number', 'Uniprot protein evidence','Consv. Score','Rxn Sim.','Rxn Sim RF.','Direction Used','Combined score','% helices','% sheets','% turns','% coils', 'Mol. Weight', 'Isoelec. Point',
 'Polar %', 'Query','Hit']):
    
    for i in range(len(columns), 0, -1):
        key = columns[i-1]
        print('sorting by', head[abs(key)-1], key)
        if key != 0:
            if key < 0:
                try:
                    rows.sort(key = lambda x: -float(x[abs(key)-1]))
                except:
                    rows.sort(key = lambda x: x[abs(key)-1], reverse=True)
            else:
                try:
                    rows.sort(key = lambda x: float(x[key-1]))
                except:
                    rows.sort(key = lambda x: x[key-1])
    return rows


def sort_rows_df(rows, filt, head = ''):
    head = ['Seq. ID', 'Description', 'Organism Source', 'Tax. distance', 'Rxn. ID','Rxn. ID - additional', 
            'EC Number', 'Uniprot protein evidence','Consv. Score','Rxn Sim.','Rxn Sim RF.','Direction Used',
            'Combined score','% helices','% sheets','% turns','% coils', 'Mol. Weight', 'Isoelec. Point'
            ,'Polar %', 'Query','Hit']
    
    exp_numeric = ['Tax. distance', 'Uniprot protein evidence', 'Uniprot protein evidence',
       'Consv. Score', 'Rxn Sim.', 'Rxn Sim RF.', 'Direction Used',
       'Combined score', '% helices', '% sheets', '% turns', '% coils',
       'Mol. Weight', 'Isoelec. Point', 'Polar %']
    
    asc = [False if x <0 else True for x in filt]
    filt = [head[abs(x)-1] for x in list(filt)]
    print('sorting by', filt)
    
    
    check_num = []
    swap_nan = {}
    for i, x in enumerate(filt):
        if x in exp_numeric:
            check_num.append(x)
            if asc[i] == True:
                swap_nan[x] = np.inf
            else:
                swap_nan[x] = -np.inf
    
    sorted_data_test = pd.DataFrame(rows, columns = head)
    sorted_data_test[check_num] = sorted_data_test[check_num].apply(pd.to_numeric, errors = 'coerce')
    
    for c, sn in swap_nan.items():
        sorted_data_test[c] = sorted_data_test[c].replace(np.nan, sn)

    
    sorted_data_test2 = sorted_data_test.sort_values(filt, ascending = asc)
    #x = sorted_data_test2[['Seq. ID', 'Tax. distance', 'Combined score', 'Rxn Sim RF.']]
    
    return sorted_data_test2.values.tolist()

# def sort_rows(rows, columns, outdir):
#     with open(outdir+'_ruthsprint.txt', 'a+') as handler:
#         handler.write("\n in sort_rows number rows "+ str(len(rows)) + "\n")

#     # # remove any rows with blank data values
#     for i in range(len(columns), 0, -1):
#         key = columns[i-1]
#         # with open(outdir+'_ruthsprint.txt', 'a+') as handler:
#         #     for x in rows:
#         #         if x[abs(key)-1] == '-':
#         #             handler.write("sort_rows filtered " + str(key) + "  " + str(x) +"\n")


#         rows= [ x for x in rows if x[abs(key)-1] != '-' ]
#         # with open(outdir+'_ruthsprint.txt', 'a+') as handler:
#         #     handler.write("len rows after " +str(key)+ "   "+ str(len(rows)) + "\n")

#     for i in range(len(columns), 0, -1):
#         key = columns[i-1]
#         # print('\nkey', key)
#         with open(outdir+'_ruthsprint.txt', 'a+') as handler:
#             handler.write("in sort_rows key: "+ str(key) + "\n")
#         with open(outdir+'_ruthsprint.txt', 'a+') as handler:
#             handler.write("rows: "+ str(len(rows)) + "\n")
#         with open(outdir+'_ruthsprint.txt', 'a+') as handler:
#             handler.write("rows type: "+ str(type(rows)) + "\n")
#         with open(outdir+'_ruthsprint.txt', 'a+') as handler:
#             handler.write("end rows "+ str(rows[-4:]) + "\n")
#         with open(outdir+'_ruthsprint.txt', 'a+') as handler:
#             handler.write("columns "+ str(columns) + "\n")



#         if key != 0:
#             if key < 0:
#                 try:
#                     rows.sort(key = lambda x: -float(x[abs(key)-1]))
#                 except:
#                     rows.sort(key = lambda x: x[abs(key)-1], reverse=True)
#                     # print('e1')
#             else:
#                 try:
#                     rows.sort(key = lambda x: float(x[key-1]))
#                 except:
#                     rows.sort(key = lambda x: x[key-1])
#                     # print('e2')
#             # for x in rows[0:8]: print( x[11], x[4])
#             # print('end')
#             # for x in rows[4000:4005]: print( x[11], x[4])
#     return rows

def write_fasta(fastaFile, targets, pc, short=False, info=False, maxlength=25):
    def shorten(s, maxlength=maxlength):
        ns = re.sub('\s+', '_', s)
        ns = re.sub('\|','_', ns)
        if len(ns) > maxlength:
            ns = ''.join(list(ns)[0:maxlength])+'...'
        return ns
    with open(fastaFile, "w") as f:
        for t in targets:
            try: 
                seq = pc.sequence[t]
                if info:
                    desc = ''
                    try:
                        desc = pc.descriptions[t]
                    except:
                        pass
                    org = ''
                    try:
                        org = pc.seqorg[t][1]
                    except:
                        pass
                    ecl = set()
                    if t in pc.UprotToMnx:
                        for r in pc.UprotToMnx[t]:
                            if r in pc.rxnec:
                                ecl |= pc.rxnec[r]
                    ec = ';'.join( sorted(ecl) )
                    st = '_'.join( [t, shorten(desc,30), shorten(org), ec] )
                    print ('>{0} \n{1}'.format(st, seq), file=f)
                elif short:
                    print ('>{0} \n{1}'.format(t, seq), file=f)
                else:
                    fdesc = pc.fulldescriptions[t]
                    print ('>{0} \n{1}'.format(fdesc, seq), file=f)
            except KeyError:
                pass

def short_fasta(fastafile):
    dirname = os.path.dirname(fastafile)
    basename = os.path.basename(fastafile)
    shortfile = 'short_'+basename
    shortname = os.path.join(dirname, shortfile)
    with open(fastafile) as fasta, open(shortname, 'w') as fastaw:
        for line in fasta:
            if line.startswith('>'):
                head = line.split('OS=')
                if len(head[0].split('|')) >= 2:
                    seqid = head[0].split('|')[1]
                    line = '>'+seqid+'\n'
            fastaw.write(line)
    return shortfile


def extend_sequences(initialfastafile, fastafile, workfolder, noMSA):
    """ Extend the fasta file """
    csvfile = os.path.join(workfolder, 'selenzy_results.csv')
    shortfile = short_fasta(os.path.join(workfolder, fastafile))
    try:
        # TO do: Check that does not assume Uniprot format...
        (sequence, names, descriptions,
         fulldescriptions, osource, pexistence) = readFasta(workfolder, fastafile, limit=1000)
        
        if len(sequence) == 0:
            return csvfile

        (hydrop, weight, isoelec, polar, helices,
         sheets, turns, coils) = sequence_properties(os.path.join(workfolder, shortfile))
        fasta = open(os.path.join(workfolder, initialfastafile)).readlines()
        fasta.extend(open(os.path.join(workfolder, fastafile)).readlines())
        with open(os.path.join(workfolder, initialfastafile), 'w') as handler:
            for line in fasta:
                handler.write(line)
        shortfile = short_fasta(os.path.join(workfolder, initialfastafile))
        if not noMSA:
            cons = conservation_properties(os.path.join(workfolder, shortfile))
        else:
            cons = ({},{})
        
    except:
        return csvfile
    # Extend the csvfile
    head, rows = read_csv(csvfile)
    # Update the conservation scores
    if not noMSA:
        for r in rows:
            seqid = r[head.index('Seq. ID')]
            if seqid in cons:
                r[head.index('Consv. Score')] = cons[seqid]
    for k in range(0, len(fulldescriptions)):
        n = names[k]
        try:
            conservation = cons[n]
        except:
            conservation = 0
        try:
            h = helices[n]
            e = sheets[n]
            t = turns[n]
            c = coils[n]
        except:
            h = '-'
            e = '-'
            t = '-'
            c = '-'
        try:
            w = weight[n]
            i = isoelec[n]
            pol = polar[n]
        except:
            w = '-'
            i = '-'
            pol = '-'
        try:
            description = descriptions[n]
            source = osource[n]
            ext = pexistence[n]
        except:
            description = '-'
            source = '-'
        row = [0.0, n, description, source, -1,
               '-', '-', ext, conservation, 1, 1, 1,
               h, e, t, c, w, i, pol, '-', '-']
        rows.append( row )
    write_csv(csvfile, head, rows)
    return csvfile

def sequence_properties(fastaShortNameFile):
    #analysis of final list of sequences
    (hydrop, weight, isoelec, polar) = pepstats(fastaShortNameFile, os.path.dirname(fastaShortNameFile))
    (helices, sheets, turns, coils) = garnier(fastaShortNameFile,  os.path.dirname(fastaShortNameFile))
    return hydrop, weight, isoelec, polar, helices, sheets, turns, coils

def conservation_properties(fastaFile):
    cons = doMSA(fastaFile,  os.path.dirname(fastaFile))
    return cons

def analyse(rxnInput, targ, datadir, outdir, csvfilename, pdir=0, host='83333', NoMSA=True, pc=None, fp='Morgan'):
    # targ, datadir, outdir, csvfilename, pdir, host, NoMSA, pc, fp = arg.tar, arg.datadir, arg.outdir, arg.outfile, 'dc.csv', '83333',  True,  None, 'Morgan'
    
    with open(outdir+'_ruthsprint.txt', 'w+') as handler:
        handler.write("in Selenzy_update\n starting analyse to make file " + csvfilename +"\n")

    # ##### my code
    rxnm = RXNMapper()    

    now=datetime.now()
    ct=now.strftime("%H:%M:%S")

    with open(outdir+'_ruthsprint.txt', 'a+') as handler:
        handler.write("running selenzyme at "+ ct +"\n" +
            'input: ' + 'rxnInput ' + str(rxnInput) +"\n" +
            'targ: '+str(targ) +"\n" +
            'datadir: ' + str(datadir)+"\n" +
            'outdir: ' + str(outdir) +"\n" +
            'csvfilename: ' +str(csvfilename) +"\n" +
            'pdir: '+ str(pdir)+"\n" +
            'NoMSA: '+ str(NoMSA) +"\n" +
            'pc: ' + str(pc)+"\n" +
            'fp: ' + str(fp))

    #####

    datadir = os.path.join(datadir)
    outdir = os.path.join(outdir)
    
#    rxnname = os.path.splitext(rxn)[0]
#    csvname = rxnname.rsplit('/', 1)[-1]
    
    if csvfilename: 
        csvfilename = csvfilename
    else:
        csvfilename = "results_selenzy.csv"
    print ("Acquiring databases...")
    if pc is None:
        print ("Reading Data...")
        pc = readData(datadir) 
        print('data in')

    with open(outdir+'_ruthsprint.txt', 'a+') as handler:
        handler.write("databases acquired\n") 
    
    print ("Running quickRsim...pdir =", pdir)
    try:
        (MnxSim, MnxDirPref, MnxDirUsed, Smiles, EcNumber, MnxSimRF, SumScore) = getMnxSim(rxnInput, datadir, outdir, rxnm, pdir, fp, pc)
        #print('SumScore2', SumScore)
    except:
        with open(outdir+'_ruthsprint.txt', 'a+') as handler:
            handler.write("LOST MnxSim\n" + '   '.join([str(x) for x in [rxnInput, datadir, outdir, rxnm, pdir, fp, pc]])  )

        return False, pc

    with open(outdir+'_ruthsprint.txt', 'a+') as handler:
        handler.write("back in selenzy.analyse \n") #"+rxnInput+"

    # crash

    # sequence = pc.sequence
    # names = pc.names
    descriptions = pc.descriptions
    fulldescriptions = pc.fulldescriptions
    osource = pc.osource
    pexistence = pc.pexistence
    seqorg = pc.seqorg
    tax = pc.tax
    MnxToUprot = pc.MnxToUprot
    # upclst = pc.upclst
    # clstrep = pc.clstrep
    smir = pc.smir


    print('\ntop reaction RF')
    top_reactions1 = []
    list_mnx1 = sorted(MnxSimRF, key=MnxSimRF.__getitem__, reverse=True)
    for x in list_mnx1:
        if MnxSimRF[x] != 'nan':
            if x in MnxToUprot:
                top_reactions1.append(x)
                if len(top_reactions1)==10: break
    for x in top_reactions1: print(x, MnxSimRF[x], str(round(float(MnxSimRF[x]), 2)), MnxSim[x])
 
    print('\ntop reactions')    
    top_reactions2 = []
    list_mnx2 = sorted(MnxSim, key=MnxSim.__getitem__, reverse=True)
    for x in list_mnx2:
        if MnxSim[x] != 'nan':
            if x in MnxToUprot:
                top_reactions2.append(x)
                if len(top_reactions2)==10: break
    for x in top_reactions2: print(x, MnxSim[x], str(round(float(MnxSim[x]), 2)), MnxSimRF[x])    
    


    with open(outdir+'_ruthsprint.txt', 'a+') as handler:
        handler.write("seqorg "+str(seqorg)[0:10]+"\ndone variables!\n")

    print('sorting reactions before selecting the top ', targ)
    #list_mnx = sorted(MnxSim, key=MnxSim.__getitem__, reverse=True)  #allow user to manipulate window of initial rxn id list
    list_mnx = sorted(SumScore, key=SumScore.__getitem__, reverse=True)  #allow user to manipulate window of initial rxn id list
    print ("Creating initial MNX list...")
    targets = set()
    UprotToMnx = {}
    UprotToMnxOther = {}
    
    with open(outdir+'_ruthsprint.txt', 'a+') as handler:
        handler.write("about to loop through all mnxr e.g."+str(list_mnx[0:5])+"\n")
        for x in list_mnx[0:5]:
            handler.write(x + '  ' + str(MnxSim[x])  +"\n")
            
    # first creating fasta file, f, for further data extraction
    for i, x in enumerate(list_mnx):
        up = MnxToUprot.get(x)  
        # if x == 'MNXR111576':crash
        #if i<30: print(MnxSim[x], MnxSimRF[x], SumScore[x], x, EcNumber[x])
        if up is not None:
            for y in up:
                if len(targets) >= int(targ):    # allow user to manipulate desired number of entries for resulting table
                    break
                else:
                    # added if statement
                    if y not in UprotToMnx.keys():
                        UprotToMnx[y] = x
                    else:
                        if y not in UprotToMnxOther: 
                            UprotToMnxOther[y] = []
                        UprotToMnxOther[y].append(x+'_'+ MnxSim[x] + '_' + MnxSimRF[x])
                    targets.add(y)
 
        else:
            if len(UprotToMnx) == 0:
                print('missing enzymes for top reactions', x, MnxSim[x], SumScore[x])
            with open(outdir+'_ruthsprint.txt', 'a+') as handler:
                handler.write("missing from MNxToUprot (reac_.tsv)"+ str(x) +" \n")            
            continue
        
    with open(outdir+'_ruthsprint.txt', 'a+') as handler:
        handler.write("got "+str(len(targets))+" targets\n" + str(UprotToMnx) + "\n\n")



    print('processing fasta....')
    fastaFile = os.path.join(outdir, "sequences.fasta")
    write_fasta(fastaFile, targets, pc)
    # Avoid issues with sequence ids
    fastaShortNameFile = os.path.join(outdir, "seqids.fasta")
    fastaInfoNameFile = os.path.join(outdir, "seqinfo.fasta")
    write_fasta(fastaShortNameFile, targets, pc, short=True)
    write_fasta(fastaInfoNameFile, targets, pc, info=True)
    

    print('sequence properties....', NoMSA)
    (hydrop, weight, isoelec, polar, helices, sheets, turns, coils) = sequence_properties(fastaShortNameFile)
    if not NoMSA:
        cons = conservation_properties(fastaInfoNameFile)
    else:
        cons = {}
    
    print ("Acquiring sequence properties...")
    # final table, do all data and value storing before this!
    tdist = {}
    rows = []

    with open(outdir+'_ruthsprint.txt', 'a+') as handler:
        handler.write("\n\n\nGetting properties from fasta\n"+str(len(targets)))


    for y in targets:
        with open(outdir+'_ruthsprint.txt', 'a+') as handler:
            handler.write(y + "\n")
        try:
            # Essential sequence information
            if y not in descriptions.keys():
                with open(outdir+'_ruthsprint.txt', 'a+') as handler:
                    handler.write(str(y)+" not in descriptions\n") 
                desc = 'NA'
            else: desc = descriptions[y]

            if y not in fulldescriptions.keys():
                with open(outdir+'_ruthsprint.txt', 'a+') as handler:
                    handler.write(str(y)+" not in fulldescriptions\n")
                fdesc = 'NA'
            else: fdesc = fulldescriptions[y]

            if y not in osource.keys():
                if y in pc.seqorg.keys():
                    org = pc.seqorg[y][1]
                else:
                    org = 'NA'
                with open(outdir+'_ruthsprint.txt', 'a+') as handler:
                    handler.write(str(y)+" not in osource\n")                    
            else: org = osource[y]

            if y not in pexistence.keys():
                with open(outdir+'_ruthsprint.txt', 'a+') as handler:
                    handler.write(str(y)+" not in pexistence\n")
                ext = 'NA'
            else: ext = pexistence[y]

            if y not in UprotToMnx.keys():
                with open(outdir+'_ruthsprint.txt', 'a+') as handler:
                    handler.write(str(y)+" not in UprotToMnx\n")
                print('how why', y)
                mnx = 'NA'
            else: mnx = UprotToMnx[y]

            if y not in UprotToMnxOther.keys():
                mnxOther = '-'
            else: mnxOther = '|'.join(UprotToMnxOther[y])
            
            if mnx not in MnxSim.keys():
                with open(outdir+'_ruthsprint.txt', 'a+') as handler:
                    handler.write(str(y)+" not in MnxSim\n")


            # desc = descriptions[y]
            # fdesc = fulldescriptions[y]
            # org = osource[y]
            # ext = pexistence[y]
            # mnx = UprotToMnx[y]
            rxnsimpre = float(MnxSim[mnx])
            rxnsim = float("{0:.5f}".format(rxnsimpre))

            with open(outdir+'_ruthsprint.txt', 'a+') as handler:
                handler.write(str(y)+" got data step 1\n")

            # Non-essential sequence information
            # cn = 0
            # repid = 0
            # try:
            #     cn = upclst.get(y)
            #     repid = clstrep[cn]
            # except:
            #     cn = 0
            #     repid = y
            try:
                conservation = float(cons[y])
            except:
                conservation = 0.0
            try:
                ecid = EcNumber[mnx]
            except:
                ecid = ''
            try:
                h = helices[y]
                e = sheets[y]
                t = turns[y]
                c = coils[y]
            except:
                h = '-'
                e = '-'
                t = '-'
                c = '-'
            try:
                w = weight[y]
                i = isoelec[y]
                pol = polar[y]
            except:
                w = '-'
                i = '-'
                pol = '-'
            try:
                rxndirused = MnxDirUsed[mnx]
                # rxndirpref = MnxDirPref[mnx]
            except:
                rxndirused = 1
                # rxndirpred = 1
            try:
                mnxSimRF = MnxSimRF[mnx]
            except:
                mnxSimRF = float("nan")
            try:
                sumScore = SumScore[mnx]
            except:
                sumScore = rxnsim

            mnxSmiles = ''
            if mnx in smir:
                if rxndirused == 1:
                    mnxSmiles = smir[mnx][0]
                else:
                    mnxSmiles = smir[mnx][1]
            
            if org not in tdist:
                if y in seqorg:
                    tdist[org] = taxDistance(tax, host, seqorg[y][0])
                    # print('in seqorg', tdist[org], seqorg[y][0])
                    with open(outdir+'_ruthsprint.txt', 'a+') as handler:
                        handler.write(str(tdist[org])+" tdist[org]\n")                    
                else:
                    #tdist[org] = '-'
                    tdist[org] = float('Inf')
                    # print('NOPE')
                    # with open(outdir+'_ruthsprint.txt', 'a+') as handler:
                    #     handler.write(str(org)+" INVENTING tdist[org]\n")

            rows.append( (y, desc, org, tdist[org], mnx, mnxOther, ecid, ext, conservation, rxnsim, mnxSimRF, rxndirused, sumScore, h, e, t, c, w, i, pol, Smiles, mnxSmiles) )

            with open(outdir+'_ruthsprint.txt', 'a+') as handler:
                handler.write(str(y)+" end of iteration - appended to rows\n")
       
        except KeyError:
            print(y)
            with open(outdir+'_ruthsprint.txt', 'a+') as handler:
                handler.write("KeyError outer "+str(y)+"\n")
            pass

    with open(outdir+'_ruthsprint.txt', 'a+') as handler:
        handler.write("\nend of massive loop to add data to rows. no rows"+ str(len(rows)) +"\n\nreplacing empty vals with NA\n")

    # this will filter out files that dont have values in columns 4 (tdist[org]), -9 (rxndirused) and -10(rxndirpref)
    #sortrows = sort_rows(rows, (-10, -9, 4), outdir)
    #sortrows = sort_rows(rows, (-11, -8, 3), outdir)

    head = ('Seq. ID','Description', 'Organism Source', 'Tax. distance', 'Rxn. ID', 'Rxn. ID - additional', 'EC Number', 'Uniprot protein evidence', 'Consv. Score',
            'Rxn Sim.', 'Rxn Sim RF.', "Direction Used", "Combined score",
            '% helices', '% sheets', '% turns', '% coils', 'Mol. Weight', 'Isoelec. Point', 'Polar %','Query', 'Hit')       

    #sortrows = sort_rows(rows, (-12, 5), outdir)
    sortrows = sort_rows_df(rows, (-13, 4), head)
    # new filter rs
    # sortrows = [['NA' if (y == '-' or y == '' ) else y for y in x] for x in sortrows] 

    with open(outdir+'_ruthsprint.txt', 'a+') as handler:
        handler.write("sorted rows\n")

    updateMSA(outdir, sortrows)
    with open(outdir+'_ruthsprint.txt', 'a+') as handler:
        handler.write("updatedMSA\n")



    #csvfilename = csvfilename.replace('.csv', '_MNXR189840.csv')
    write_csv(os.path.join(outdir, csvfilename), head, sortrows)
    write_csv(os.path.join(outdir, csvfilename.replace('.csv', '2.csv')), head, sortrows)

    with open(outdir+'_ruthsprint.txt', 'a+') as handler:
        handler.write("CSV file created " + outdir+ csvfilename + "\n")

    print ("CSV file created.", outdir+ csvfilename)

    return True, pc
    
    
def arguments(args=None):
    parser = argparse.ArgumentParser(description='SeqFind script for Selenzy')
    parser.add_argument('rxn', 
                        help='Input reaction [default = rxn file]')
    parser.add_argument('datadir',
                        help='specify data directory for required databases files, please end with slash')
    parser.add_argument('outdir',
                        help='specify output directory for all output files, including final CSV file, please end with slash')
    
    parser.add_argument('-tar', type=float, default=200,
                        help='Number of targets to display in results [default = 20]')
    parser.add_argument('-d', type=float, default=0,
                        help='Use similiarity values for preferred reaction direction only [default=0 (OFF)]')

    parser.add_argument('-outfile',
                        help='specify non-default name for CSV file output')
    parser.add_argument('-NoMSA', action='store_true',
                        help='Do not compute MSA/conservation scores')
    parser.add_argument('-smarts', action='store_true',
                        help='Input is a reaction SMARTS string')
    parser.add_argument('-smartsfile', action='store_true',
                        help='Input is a reaction SMARTS file')
    parser.add_argument('-host', type=str, default='83333',
                        help='Host organism taxon id [default: E. coli]')
    arg = parser.parse_args(args=args)
    return arg

if __name__ == '__main__':
    arg = arguments()
    
    # arg = arguments(['O=C([O-])CCC(=O)C(=O)[O-].NC(CC(=O)[O-])C(=O)O>>O=C([O-])CC(=O)C(=O)[O-].NC(CCC(=O)[O-])C(=O)O', 
    #     '/home/ruth/code/update_selenzyme/run_folder/selenzyme/selenzyPro/data/',
    #     '/home/ruth/code/update_selenzyme/run_folder/selenzyme/selenzyPro/uploads/',
    #     '-outfile', 'selenzy_results.csv' ] )

    # arg = arguments(['C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@@]43C)[C@@H]1CC[C@]2(O)C(=O)CO>>C[C@]12C[C@H](O)[C@H]3[C@@H](CCC4=CC(=O)CC[C@@]43C)[C@@H]1CC[C@]2(O)C(=O)CO', 
    #     '/home/ruth/code/update_selenzyme/run_folder_min_dec/data/',
    #     '/home/ruth/code/update_selenzyme/run_folder/selenzyme/selenzyPro/uploads/',
    #     '-outfile', 'selenzy_results.csv' ] )

    
    newpath = os.path.join(arg.outdir)
    if not os.path.exists(newpath):
        os.makedirs(newpath)
        with open(os.path.join(arg.datadir)+'ruthsprint.txt', 'w+') as f: 
            f.write('new file\n')
    
    if arg.smarts is not None:
        rxnInput = ['-smarts', arg.rxn]
    elif arg.smartsfile:
        rxnInput = ['-smartsfile', arg.rxn]
    else:
        rxnInput = ['-rxn', arg.rxn]

    analyse(rxnInput, arg.tar, arg.datadir, arg.outdir, arg.outfile, arg.d, arg.host, NoMSA=arg.NoMSA)

    # head, rows = read_csv(arg.outdir +'selenzy_results.csv')       
    #score = [('J', 'Reaction similarity:', 50.0, True), ('M', 'Combined score', 0.0, False), ('K', 'RF similarity', 50.0, False), ('I', 'Sequence <a href="https://en.wikipedia.org/wiki/Conserved_sequence" target="_blank">conservation</a>:', 1.0, True), ('E', 'Sequence <a href="https://www.ncbi.nlm.nih.gov/taxonomy" target="_blank">taxonomic</a> distance:', -1.0, True), ('H', 'Uniprot <a href="http://www.uniprot.org/help/protein_existence" target="_blank">protein evidence</a>:', -0.1, True), ('O', 'Percentage sheets:', 0.0, False), ('P', 'Percentage turns:', 0.0, False), ('Q', 'Molecular weight:', 0.0, False), ('R', 'Isoelectric point:', 0.0, False)]
    #updateScore(csvfile, score)
    # # remove rows with empty values
    # filt = [-13, 4]
    # filt = (-11, 4)
    # sortrows = sort_rows_df(rows, filt) # , head
    # sorted_data_test = pd.DataFrame(sortrows, columns = head)
    # sorted_data_test2 = sorted_data_test.sort_values(['Rxn Sim RF.', 'Tax. distance'], ascending = [False, True])
    # sorted_data_test2 = sorted_data_test.sort_values(['Rxn Sim RF.'], ascending = [False])
    # x = sorted_data_test2[['Seq. ID', 'Tax. distance', 'Combined score', 'Rxn Sim RF.']]
    # y = x.sort_index()
                        

    
    
    
       




                 
# python quickRsim_local4_debug.py \
#     /home/ruth/code/update_selenzyme/run_folder/selenzyme/selenzyPro/data/ \
#     Morgan \
#     -smarts '[H]Oc1c(OC([H])([H])[H])c([H])c(C([H])=C([H])C(=O)OC2([H])OC([H])(C([H])([H])O[H])C([H])(O[H])C([H])(O[H])C2([H])O[H])c([H])c1OC([H])([H])[H]>>[H]OC(=O)C([H])([H])C([H])(OC(=O)C([H])=C([H])c1c([H])c(OC([H])([H])[H])c(O[H])c(OC([H])([H])[H])c1[H])C(=O)O[H].[H]OC([H])([H])C([H])(O[H])C([H])(O[H])C([H])(O[H])C([H])(O[H])C([H])=O' \
#     -chem /home/ruth/code/update_selenzyme/selenzyme_update/data_google_cloud_edited/data_test/chem_prop.tsv \
#     -out /home/ruth/code/update_selenzyme/run_folder/RSquickRsim_ring_MTrue.txt 


# python gitcode/selenzyPro/selenzy_reclean.py \
#     '[H]Oc1c(OC([H])([H])[H])c([H])c(C([H])=C([H])C(=O)OC2([H])OC([H])(C([H])([H])O[H])C([H])(O[H])C([H])(O[H])C2([H])O[H])c([H])c1OC([H])([H])[H]>>[H]OC(=O)C([H])([H])C([H])(OC(=O)C([H])=C([H])c1c([H])c(OC([H])([H])[H])c(O[H])c(OC([H])([H])[H])c1[H])C(=O)O[H].[H]OC([H])([H])C([H])(O[H])C([H])(O[H])C([H])(O[H])C([H])(O[H])C([H])=O' \
#     data/ \
#     selenzyme2/selenzyPro/uploads/ \
#     -outfile selenzy_results.csv







