from collections import Counter
from fileinput import filelineno
from statistics import mean
from settings import *
import re, time, sys, os, platform, subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from matplotlib import pyplot as plt
from matplotlib.patches import Patch
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

        if len(mainSeq) < 140 : continue
        
        if not mainSeq in MutantSeqDict.keys() and mainSeq != "": 
            MutantSeqDict[mainSeq] = 1
        elif mainSeq in MutantSeqDict.keys()  and mainSeq != "": 
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



def countMutationForMutagenesisData(fastaFile, sourceFolder = "fasta/alignment/MutagenesisData/", destinationFolder = "results/MutagenesisData/countMutation/"):
    if not os.path.isdir(destinationFolder) : os.mkdir(destinationFolder)
    # nbMut = lambda ref_seq, mut_seq: sum(ei != ej for ei, ej in zip(ref_seq, mut_seq))
    file_number = fastaFile[fastaFile.find("T0"):fastaFile.find("T0")+3]
    print("counting mutation for : ", file_number )
    mutationPosCountDict = {}
    mutationCountList = []
    mutationTypeCountDict= {
        "insertion":0, 
        "deletion":0, "A->C":0, 
"A->G":0, 
"A->T":0,
"C->A":0,
"C->G":0,
"C->T":0,
"G->A":0,
"G->C":0,
"G->T":0,
"T->A":0, 
"T->C":0, 
"T->G":0, 
}
    recordList = list(SeqIO.parse(sourceFolder + fastaFile, "fasta"))
    
    refSeq = str(recordList[0].seq).upper()
    startGapePos = {}
    mutationTypeByPosDict = {}
    for i in range(1, len(refSeq)+1): mutationPosCountDict[i] = 0
    for i in range(1, len(refSeq)+1): 
        mutationTypeByPosDict[i] = {"insertion":0, "deletion":0,"A->C":0, "A->G":0, "A->T":0,"C->A":0,"C->G":0,"C->T":0,"G->A":0,"G->C":0,"G->T":0,"T->A":0, "T->C":0, "T->G":0, }
    for record in SeqIO.parse(sourceFolder + fastaFile, "fasta"):
        if record.id == "SunY_sequence": continue
        seq = str(record.seq).upper()

        if not record.id in startGapePos .keys() : startGapePos [record.id] =0
        for i  in range(len(seq)): 
            if seq[i] == "-" : startGapePos[record.id] +=1
            else : break
        nbOccurence = int(record.id[record.id.rfind("="):].replace("=",""))
        
        nbMutation = 0
        for i in mutationPosCountDict.keys():
            if refSeq[i-1] != seq[i-1] and refSeq[i-1] !="N" and  seq[i-1] !="N":
                nbMutation+= 1
                mutationPosCountDict[i] +=1 * nbOccurence
                if refSeq[i-1] =="-": 
                    mutationTypeCountDict["insertion"]+=1 * nbOccurence
                    mutationTypeByPosDict[i]["insertion"]+=1 * nbOccurence
                elif seq[i-1] == "-": 
                    mutationTypeCountDict["deletion"]+=1 * nbOccurence
                    mutationTypeByPosDict[i]["deletion"]+=1 * nbOccurence
                else: 
                    mutationType = refSeq[i-1] + "->"+ seq[i-1]
                    if not mutationType in mutationTypeCountDict.keys(): 
                        mutationTypeCountDict[mutationType] =1*nbOccurence
                        mutationTypeByPosDict[i][mutationType] =1*nbOccurence
                    else: 
                        mutationTypeCountDict[mutationType] +=1*nbOccurence
                        mutationTypeByPosDict[i][mutationType] +=1*nbOccurence
        mutationCountList.extend([nbMutation] *nbOccurence)

    print("writting log file")
    log = "mutation count for : " + fastaFile + "\n"
    log += "the avrage number of mutation pear sequence  in the file is : " + str(round(sum(mutationCountList)/len(mutationCountList), 2)) + " mutation pear sequence with " + str(len(mutationCountList)) + "\n"
    log += "the avrage number of mutation pear nucleotyde in the file is : " + str(sum(mutationCountList)/len(mutationCountList)/len(refSeq)) + "\n"
    log += "which represent : " + str(sum(mutationCountList)/len(mutationCountList)/len(refSeq)*100) + " %" + "\n"
    log+= "the type of position and their number are as followed : \n"
    for k,v in mutationTypeCountDict.items(): log+= k + " : " + str(v) + " which represent " +  str(round(v / sum(mutationTypeCountDict.values()) *100, 2)) + " % of the mutation \n"
    log+= "please look at the graph for more infos on the count of mutation by position"
    # i = 0
    # for k,v in startGapePos.items():
        # print(k + " : " + str(v))
        # i+=1
        # if i>10: break
    # print(mean(startGapePos.values()))
    with open(destinationFolder+ file_number + "count_mutation_log.txt", "w") as logFile:
        logFile.write(log)

    return (mutationCountList, mutationPosCountDict, mutationTypeCountDict, mutationTypeByPosDict)


# function plotSeqLengthDistribution
# This fonction plot as an histogram the sequence length distribution of a fasta file
# @param
# @seqLengthList, the list of sequence length to plot
def plotSeqLengthDistributionForMutagenesisData(seqLengthList, fastaFile, destinationFolder = "figures/MutagenesisData/sequence_length_distribution/"):
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
def plotMutationDistributionForMutagenesisData(mutationCountList, fastaFile, destinationFolder = "figures/MutagenesisData/mutation_count_distribution/"):
    if ".aln" in fastaFile: fastaFile= fastaFile.replace(".aln", ".fasta")
    print("Plotting mutation count distribution for ", fastaFile)
    if not os.path.isdir(destinationFolder) : os.makedirs(destinationFolder)
    mutationCountList.sort()
    tmp = Counter(mutationCountList)
    mutationCounter = {}
    # print(tmp)
    for k,v in tmp.items() : mutationCounter[str(k)] = round(v/sum(tmp.values())*100, 2)
    plt.figure(figsize=(15,8), facecolor='white')
    plt.title("Mutation count distribution")
    plt.suptitle("for " + fastaFile.replace(".", " ").replace("fasta", ""))
    # plt.xlim([min(tmp.keys())-1, max(tmp.keys())+1])
    plt.xlim([0, 15])
    plt.xlabel('number of mutation')
    plt.ylabel('percentage of sequence with this number of mutation')
    # plt.hist(mutationCountList, alpha=0.5, facecolor='purple')
    plt.bar(mutationCounter.keys(), mutationCounter.values())
    # plt.show()
    plt.savefig(destinationFolder + "mutation_count_distribution_" + fastaFile.replace(".fasta", ".png"))
    plt.close()



def plotMutationPosDistributionForMutagenesisData(mutationPosCountDict, fastaFile, destinationFolder = "figures/MutagenesisData/mutation_count_distribution_by_position/"):
    if ".aln" in fastaFile: fastaFile= fastaFile.replace(".aln", ".fasta")
    print("Plotting mutation pos count distribution for ", fastaFile)
    if not os.path.isdir(destinationFolder) : os.makedirs(destinationFolder)
    tmp={}
    for k, v in mutationPosCountDict.items(): tmp[str(k)] = round(v/sum(mutationPosCountDict.values())*100, 2)
    # for k, v in tmp.items(): print(k, str(v))
    plt.figure(figsize=(15,8), facecolor='white')
    plt.title("Mutation count distribution by position")
    plt.suptitle("for " + fastaFile.replace(".", " ").replace("fasta", ""))
    plt.xlim([min(mutationPosCountDict.keys())-1, 170])
    plt.xticks(np.arange(-1, 171, step=5))
    plt.xlabel('Position dans la sequence de sunwise')
    plt.ylabel('Percentage of mutation in the file')
    plt.bar(tmp.keys(), tmp.values())
    # plt.show()
    plt.savefig(destinationFolder + "mutation_count_by_position_distribution_" + fastaFile.replace(".fasta", ".png"))
    plt.close()


def plotMutationTypeByPosDistributionForMutagenesisData(mutationTypeByPosDict, fastaFile, destinationFolder = "figures/MutagenesisData/mutation_count_distribution_type_by_position/"):
    if ".aln" in fastaFile: fastaFile= fastaFile.replace(".aln", ".fasta")
    print("Plotting mutation pos count distribution for ", fastaFile)
    if not os.path.isdir(destinationFolder) : os.makedirs(destinationFolder)
    tmp={}
    typePosDict = {}
    for k, v in mutationTypeByPosDict.items():
        sigmaV = sum(v.values())
        for k2,v2 in v.items():
            if not k2 in typePosDict.keys(): 
                if sigmaV >  0 :  typePosDict[k2] = [round(v2/sum(v.values())*100, 2)]
                else : typePosDict[k2]  = [0]
            else: 
                if sigmaV > 0 : typePosDict[k2].append(round(v2/sum(v.values())*100, 2))
                else: typePosDict[k2].append(0)
            # tmp[str(k)] = {k2: round(v2/sum(v.values())*100, 2)}
    # for k, v in tmp.items(): print(k, str(v))
    mutationTypeColorDict= {
        "insertion": "black", 
        "deletion":"grey",
        "A->C":"#e62e00", 
"A->G":"#ff5c33", 
"A->T":"#ff8566",
"C->A":"#3366ff",
"C->G":"#809fff",
"C->T":"#b3c6ff",
"G->A":"#339966",
"G->C":"#66cc99",
"G->T":"#9fdfbf",
"T->A":"#993366", 
"T->C":"#d279a6", 
"T->G":"#e6b3cc", 
}
    # plt.figure(figsize=(15,8), facecolor='white')
    fig, ax = plt.subplots()
    ax.set_title("Mutation count distribution by position for each type of mutation")
    # ax.set("for " + fastaFile.replace(".", " ").replace("fasta", ""))
    ax.set_xlim([0, 170])
    ax.set_xticks(np.arange(-1, 171, step=5))
    ax.set_xlabel('Position dans la sequence de sunY')
    ax.set_ylabel('Percentage of mutation in the file')
    patches = [Patch(color=v, label=k) for k, v in mutationTypeColorDict.items()]
    ax.legend(title='Type of mutation', labels=mutationTypeColorDict.keys(), handles=patches, bbox_to_anchor=(1.04, 0.5), loc='center left', borderaxespad=0, fontsize=15, frameon=False)
    for k, v in typePosDict.items():
        ax.bar(mutationTypeByPosDict.keys(), v, color = mutationTypeColorDict[k])
    # plt.show()
    fig.savefig(destinationFolder + "mutation_type_by_position_distribution_" + fastaFile.replace(".fasta", ".png"),bbox_inches='tight')
    plt.close()

def checkSeqLength(fastaFile):

    lengthList = []
    for record in SeqIO.parse(fastaFile, "fasta-2line"): lengthList.append(len(record.seq))

    print(mean(lengthList))
    print(min(lengthList))
    print(max(lengthList))

def generateSunyMutant():
    base= ["A", "T", "C", "G"]
    baseSeq = settings["refSeqSequence"]
    mutantList = []
    for i in range(len(baseSeq)):
        for nuc in base:
            if nuc != baseSeq[i]:
                seqName = "sunY_mutant_pos_" + str(i+1) + "_" + baseSeq[i] + "_to_" + nuc
                mutedSeq = baseSeq[0:i] + nuc + baseSeq[i+1:]
                mutantList.append(SeqRecord(Seq(mutedSeq), id= seqName, name="", description=""))

    SeqIO.write(mutantList, "fasta/sunyMutan.fasta", "fasta-2line")

def checkMutanPresence(fileList, mutanFasta, fileSourcePath = "fasta/treated/MutagenesisData/"):
    mutanList = {}
    for record in SeqIO.parse(mutanFasta, "fasta-2line"):
        mutanList [record.id] = str(record.seq)

    fileInfos = {}
    for file in fileList:
        print("starting check for ", file)
        # fileInfos[file] = {"nbMutanFound":0}
        fileInfos[file] = {"nbMutanFound":0, "nbSeqInFile":0}
        listSeq = []
        for record in SeqIO.parse(fileSourcePath  + file, "fasta-2line"):
            listSeq.append(str(record.seq))

        fileInfos[file]["nbSeqInFile"] = len(listSeq)
        
        for seq in mutanList.values():
            fileInfos[file]["nbMutanFound"] += sum(seq in s for s in listSeq) 
                    # fileInfos[file]["nbMutanFound"]+=1
                

            log = "Suny mutan count \n"
            for fileName, d in fileInfos.items():
                log += str(d["nbMutanFound"]) + " sunY mutan have been found in the file " + fileName + " for a total of " + str(d["nbSeqInFile"]) + " sequence \n"
                log += "it represent " + str(d["nbMutanFound"]/ d["nbSeqInFile"] *100) + " % of the sequence \n"

            with open("results/sunYMutanSearch.txt", "w") as logFile: logFile.write(log)

# checkSeqLength("fasta/treated/MutagenesisData/L447T06.R1_pre-treated_MutagenesisSequences.fasta")

generateSunyMutant()