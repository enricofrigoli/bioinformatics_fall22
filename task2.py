# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 13:47:54 2022

@author: lotte
"""
import Bio.SeqIO as SeqIO
from Bio import pairwise2
from datetime import datetime

humDir = r"C:\Users\lotte\Documents\FAU\BioInf\human.fa"
mouseDir = r"C:\Users\lotte\Documents\FAU\BioInf\mouse.fa"

muromonab = SeqIO.to_dict(SeqIO.parse(r"C:\Users\lotte\Documents\FAU\BioInf\muromonab.fasta", "fasta"))
muromonabSeq = (muromonab['muromonab']).seq

# bevacizumab, caplacizumab

def humanization_model(seq,humDir,mouseDir):
    print(datetime.now())
    humDict = SeqIO.to_dict(SeqIO.parse(humDir, "fasta"))
    mouseDict = SeqIO.to_dict(SeqIO.parse(mouseDir, "fasta")) 
    
    humSequences = [s.seq for s in humDict.values()]
    mouseSequences = [s.seq for s in mouseDict.values()]
    
    print('starting pairwise alignment human dataset\n')    
    humAlignments = [pairwise2.align.globalxx(seq,humS) for humS in humSequences]
    print("done pairwise alignment human dataset\n")
   
    print(datetime.now())
    
    print('starting pairwise alignment mouse dataset\n')
    mouseAlignments = [pairwise2.align.globalxx(seq,mouseS) for mouseS in mouseSequences]
    print("done pairwise alignment mouse dataset\n")
    
    print(datetime.now())
         
    return humAlignments, mouseAlignments

def get_highest_score(humAlignments,mouseAlignments):

    hMax = 0
    mMax = 0
        
    for hA,mA in zip(humAlignments,mouseAlignments):
        for align in hA:
            if align.score > hMax: 
                hMax = align.score
                hSeq = align.seqB
        for align in mA:
            if align.score > mMax:
                mMax = align.score
                mSeq = align.seqB

    if hMax > mMax:
        return {'dataset': 'human', 'seq': hSeq, 'score': hMax}
    return {'dataset': 'mouse', 'seq': mSeq, 'score': mMax}        

# h,m = humanization_model(testSeq, humDir, mouseDir)
# high = get_highest_score(h, m)

hMuromonabSeq, mMuromonabSeq = humanization_model(muromonabSeq, humDir, mouseDir)
highMuromonabSeq = get_highest_score(hMuromonabSeq, mMuromonabSeq)
print(highMuromonabSeq)
