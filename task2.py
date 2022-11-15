# -*- coding: utf-8 -*-
import Bio.SeqIO as SeqIO
from Bio import pairwise2
import os
from datetime import datetime

# get the directories
# used files have to be stored in same directory as the code
cwd = os.getcwd()
humDir = os.path.join(cwd,"human.fa")
mouseDir = os.path.join(cwd, "mouse.fa")
muromonab = os.path.join(cwd, "muromonab.fasta")
bevacizumab = os.path.join(cwd, "bevacizumab.fasta")
caplacizumab = os.path.join(cwd, "caplacizumab.fasta")

# get muromonab sequence
muromonab = SeqIO.to_dict(SeqIO.parse(muromonab, "fasta"))
muromonab = list(muromonab.values())
muromonab = muromonab[0].seq

# get bevacizumab sequence
bevacizumab = SeqIO.to_dict(SeqIO.parse(bevacizumab, "fasta"))
bevacizumab = list(bevacizumab.values())
bevacizumab = bevacizumab[0].seq

# get caplacizumab sequence
caplacizumab = SeqIO.to_dict(SeqIO.parse(caplacizumab, "fasta"))
caplacizumab = list(caplacizumab.values())
caplacizumab = caplacizumab[0].seq


# local alignment of given sequence with all sequences given sequences of humans and mice
# use global alignment if sequences have roughly same length else use local
def alignment(seq,humDir,mouseDir):
    print(datetime.now())
    
    # get human and mouse dicts
    humDict = SeqIO.to_dict(SeqIO.parse(humDir, "fasta"))
    mouseDict = SeqIO.to_dict(SeqIO.parse(mouseDir, "fasta")) 
    
    # get sequences
    humSequences = [s.seq for s in humDict.values()]
    mouseSequences = [s.seq for s in mouseDict.values()]
    
    # start alignment for all given human sequences
    print('starting pairwise alignment human dataset\n')    
    humAlignments = [pairwise2.align.localxx(seq,humS) for humS in humSequences]
    print("done pairwise alignment human dataset\n")  
    print(datetime.now())
    
    # start alignment for all given mouse sequences
    print('starting pairwise alignment mouse dataset\n')
    mouseAlignments = [pairwise2.align.localxx(seq,mouseS) for mouseS in mouseSequences]
    print("done pairwise alignment mouse dataset\n")    
    print(datetime.now())
         
    return humAlignments, mouseAlignments

# find the highest score of all possible alignments
def get_highest_score(humAlignments,mouseAlignments):

    hMax = 0
    mMax = 0
      
    # find highest score of human alignments and highest score of mouse alignments
    for hA,mA in zip(humAlignments,mouseAlignments):
        for align in hA:
            if align.score > hMax: 
                hMax = align.score
                hSeq = align.seqB
        for align in mA:
            if align.score > mMax:
                mMax = align.score
                mSeq = align.seqB
    
    # compare human and mouse highscore and return highest
    if hMax > mMax:
        return {'dataset': 'human', 'seq': hSeq, 'score': hMax}
    return {'dataset': 'mouse', 'seq': mSeq, 'score': mMax}        


hMuromonabSeq, mMuromonabSeq = alignment(muromonab, humDir, mouseDir)
highMuromonabSeq = get_highest_score(hMuromonabSeq, mMuromonabSeq)
print("highest score for alignment with muromonab: "+ str(highMuromonabSeq))

hBevacizumabSeq, mBevacizumabSeq = alignment(bevacizumab, humDir, mouseDir)
highBevacizumabSeq = get_highest_score(hBevacizumabSeq, mBevacizumabSeq)
print("highest score for alignment with bevacizumab: " + str(highBevacizumabSeq))

hCaplacizumabSeq, mCaplacizumabSeq = alignment(caplacizumab, humDir, mouseDir)
highCaplacizumabSeq = get_highest_score(hCaplacizumabSeq, mCaplacizumabSeq)
print("highest score for alignment with caplacizumab: " + str(highCaplacizumabSeq))
