from collections import Counter
from settings import *
import re, time, sys, os, platform, subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

settings = getSettings()


def levenshteinDistance(s1, s2):
    matrix = np.zeros ((len(s1)+1, len(s2)+1)) 
    matrix [0:len(s1)+1, 0] = [x for x in range (0, len(s1)+1)] 
    matrix [0,0: len(s2)+1] = [y for y in range (0, len(s2)+1)] 

    for x in range(1, len(s1)+1): # iterating rows, we start from 1 because 0 was preset before
        for y in range(1, len(s2)+1): 			# if characters are equal, we do nothing
            if s1[x-1] == s2[y-1]: 
                matrix [x,y] = matrix[x-1,y-1]
            else: 
                matrix [x,y] = min(
                    matrix[x-1,y] + 1,
                    matrix[x-1,y-1] + 1,
                    matrix[x,y-1] + 1
                )
    return (matrix[len(s1), len(s2)])
    

def countMutantSeqOccurence(fastaFile, startPatternToDetect = settings["seqOccurenceCountingMutagenesis"]["countMutation_seq_start"], endPatternToDetect = settings["seqOccurenceCountingMutagenesis"]["countMutation_seq_end"], sourceFolder = "fasta/pre-treated/MutagenesisData/"):
    MutantSeqDict= {}
    MutantSeqList = []
    
    file_number = int(fastaFile[fastaFile.find("T0"):fastaFile.find("T0")+3][-1])
    print(file_number)
    print("getting data")
    
    for record in SeqIO.parse(sourceFolder+ fastaFile, "fasta-2line"):
        seq = str(record.seq)
        mainSeq=""

        if re.search(startPatternToDetect, seq) != None and re.search(endPatternToDetect, seq) != None :
            startIndex = re.search(startPatternToDetect, seq).start() # The real start should be -9 from the start
            endIndex = re.search(endPatternToDetect, seq).end()
            mainSeq = seq[startIndex:endIndex+ len(endPatternToDetect)]
        # MutantSeqList.append(SeqRecord(Seq(product), id = record.id, name = record.name, description=""))
        
        if not mainSeq in MutantSeqDict.keys(): 
            MutantSeqDict[mainSeq] = 1
        else: 
            MutantSeqDict[mainSeq] += 1

    print("getting dict data finish")

    #  sorting the dictionnary
    
    sortedMutantSeqList = []
    sortedMutantSeqDict = {}    
     
    for key, val in MutantSeqDict.items(): sortedMutantSeqList.append((val, key))
    sortedMutantSeqList.sort(reverse=True, key=lambda a: a[0])    
    for val, key in sortedMutantSeqList: sortedMutantSeqDict[key] = val    
    
    return (sortedMutantSeqList, MutantSeqDict, sortedMutantSeqDict, fastaFile)


def writeSeqMutagenesisFasta(SequenceDict, fastaFile, destinationFolder = "fasta/treated/MutagenesisData/"):
    seqId = 0
    recordList = []
    file_number = fastaFile[fastaFile.find("T0"):fastaFile.find("T0")+3]
    
    for seq, nbFound in SequenceDict.items():
        if nbFound >= 1 :
            seqName = "MutantSequence" + str(seqId) + "_"+ file_number  + "_N="+ str(nbFound)
            recordList.append(SeqRecord(Seq(seq), id = seqName, name= "", description=""))
            seqId += 1
    
    if not os.path.isdir(destinationFolder) : os.mkdir(destinationFolder)
    SeqIO.write(recordList, destinationFolder + fastaFile.replace(".fasta", "_MutagenesisSequences.fasta"), "fasta-2line")



def countMutation(fastaFile, sourceFolder = "fasta/alignment/MutagenesisData/"):
    # nbMut = lambda ref_seq, mut_seq: sum(ei != ej for ei, ej in zip(ref_seq, mut_seq))
    
    mutationPosCountDict = {}
    mutationCountList = []
    mutationTypeCountDict = {"A":0,"T":0,"C":0,"G":0, "insertion":0, "deletion":0}
    
    recordList = list(SeqIO.parse(sourceFolder + fastaFile, "fasta"))
    
    refSeq = str(recordList[0].seq).upper()
    
    for i in range(1, len(refSeq)+1): 
        mutationPosCountDict[i] = 0

    for record in recordList:
        # nbMutation = nbMut(refSeq, record.seq)
        seq = str(record.seq).upper()
        nbMutation = 0
        for i in mutationPosCountDict.keys():
            if refSeq[i-1] != seq[i-1]:
                nbMutation+=1
                mutationPosCountDict[i] +=1
                if refSeq[i-1] =="-": mutationTypeCountDict["insertion"]+=1
                elif seq[i-1] =="-": mutationTypeCountDict["deletion"]+=1
                else: mutationTypeCountDict[seq[i-1]] +=1
        mutationCountList.append(nbMutation)

    return (mutationCountList, mutationPosCountDict, mutationTypeCountDict)


# function plotSeqLengthDistribution
# This fonction plot as an histogram the sequence length distribution of a fasta file
# @param
# @seqLengthList, the list of sequence length to plot
def plotSeqLengthDistribution(seqLengthList, fastaFile, destinationFolder = "figures/MutagenesisData/sequence_length_distribution/"):
    print("Plotting sequence length distribution for ", fastaFile)
    
    if not os.path.isdir(destinationFolder) : os.makedirs(destinationFolder)
    
    plt.figure(figsize=(15,8), facecolor='white')
    plt.title("Sequence length distribution")
    plt.suptitle("for " + fastaFile.replace(".", " ").replace("fasta", ""))
    plt.xlim([min(seqLengthList)-1, max(seqLengthList)+1])
    plt.xlabel('Read Length')
    plt.ylabel('Number of Reads')
    plt.hist(seqLengthList, alpha=0.5, facecolor='purple')
    # plt.show()
    plt.savefig(destinationFolder + fastaFile.replace(".fasta", "_lengthDistribution.png"))
    plt.close()
    print("Treatment finish for ", fastaFile)

# function plotMutationDistribution
# This fonction plot as an histogram mutation count distribution of a fasta file
# @param
# @mutationCountList, the list of mutation count to plot
def plotMutationDistribution(mutationCountList, fastaFile, destinationFolder = "figures/MutagenesisData/mutation_count_distribution/"):
    print("Plotting mutation count distribution for ", fastaFile)
    if not os.path.isdir(destinationFolder) : os.makedirs(destinationFolder)
    plt.figure(figsize=(15,8), facecolor='white')
    plt.title("Mutation count distribution")
    plt.suptitle("for " + fastaFile.replace(".", " ").replace("fasta", ""))
    plt.xlim([min(mutationCountList)-1, max(mutationCountList)+1])
    plt.xlabel('number of mutation')
    plt.ylabel('Number of mutation')
    plt.hist(mutationCountList, alpha=0.5, facecolor='purple')
    # plt.show()
    plt.savefig(destinationFolder + "mutation_count_distribution_" + fastaFile.replace(".fasta", ".png"))
    plt.close()