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



    

def treatFasta(fastaFile, sourceFolder = "fasta/raw_fasta/", destinationFolder = "fasta/pre-treated/"):
    refSeq= SeqRecord(Seq(settings["refSeqSequence"]), id = settings["refSeqId"], name = "", description= "")
    seqLengthList = [len(refSeq.seq)]
    seqList = [refSeq]
    sub2Only = 0
    sunyBiginningOnly = 0
    turn = fastaFile[fastaFile.find("T0"):fastaFile.find("T0")+3]
    print(turn)
    haveBoth= 0
    for record in SeqIO.parse(sourceFolder+ fastaFile, "fasta"):
        haveSuny = False
        haveSub2= False
        
        seq = str(record.seq)
        if settings["Suny_beginning"] in seq and int(turn[-1]) <=5: 
            haveSuny=True
            seq = seq[seq.find(settings["Suny_beginning"]) -5:]

        if settings["Suny_beginning2"] in seq and int(turn[-1]) >5: 
            haveSuny=True
            seq = seq[seq.find(settings["Suny_beginning2"])  -5:]
        if settings["sub2"] in seq and int(turn[-1]) <=5: 
            haveSub2=True
            seq = seq[:seq.rfind(settings["sub2"]) + len(settings["sub2"])]

        if settings["exon"] in seq and int(turn[-1]) > 5: 
            haveSub2=True
            seq = seq[:seq.rfind(settings["exon"]) + len(settings["exon"])]
        
        record.seq = Seq(seq)
        
        if len(seq) > 100 and int(turn[-1]) <= 5:
            if haveSuny and not haveSub2: sunyBiginningOnly+=1
            if not haveSuny and haveSub2: sub2Only+=1
            if haveSub2 and haveSuny: 
                haveBoth +=1
                seqList.append(record)
                seqLengthList.append(len(seq))

        elif len(seq) > 100 and int(turn[-1]) > 5:
            if haveSuny and not haveSub2: sunyBiginningOnly+=1
            if not haveSuny and haveSub2: sub2Only+=1 
            if haveSuny: 
                haveBoth +=1
                seqList.append(record)
                seqLengthList.append(len(seq))

    log = "the file " + fastaFile + " after treatement have " + str(len(seqList)) + " sequences \n"
    # log+ = "the file " + fastaFile + " after treatement have " + str(len(seqList)+ len(seqList2)) + " sequences"
    log+="including : \n" + str(sunyBiginningOnly) + " sequences have the suny beginning motif only \n"
    log+=str(sub2Only) + " sequences have the sub2 motif only \n"
    log+=str(haveBoth) + " sequences have the sub2 and suny motif \n"
    if not os.path.isdir(destinationFolder) : os.makedirs(destinationFolder)
    if not os.path.isfile(destinationFolder+"log.txt"): logFile = open(destinationFolder+ "log.txt", "w").close()
    with open(destinationFolder+"log.txt", "a") as logFile:
        logFile.write(log)

    print("write fasta")
    SeqIO.write(seqList, destinationFolder + fastaFile.replace(".fasta", "_pre-treated.fasta"), "fasta-2line")
    return seqLengthList


def plotSeqLengthDistribution(seqLengthList, fastaFile, destinationFolder = "figures/sequence_length_distribution/"):
    print("plotting sequence length distribution for ", fastaFile)
    if not os.path.isdir(destinationFolder) : os.makedirs(destinationFolder)
    plt.figure(figsize=(15,8), facecolor='white')
    plt.title("Sequence length distribution")
    plt.suptitle("lenth distribution for " + fastaFile.replace(".", " ").replace("fasta", ""))
    plt.xlim([min(seqLengthList)-1, max(seqLengthList)+1])
    plt.xlabel('Read Length')
    plt.ylabel('Number of Reads')
    plt.hist(seqLengthList, alpha=0.5, facecolor='purple')
    # plt.show()
    plt.savefig(destinationFolder + fastaFile.replace(".fasta", "_lengthDistribution.png"))
    plt.close()
    print("treatement finish for ", fastaFile)

def plotSeqCountDistribution(seqCountList, fastaFile, destinationFolder = "figures/variant_count_distribution/"):
    seqCountList.sort()
    tmp2 = Counter(seqCountList)
    print("plotting sequence count distribution for ", fastaFile)
    if not os.path.isdir(destinationFolder ) : os.makedirs(destinationFolder )

    plt.title("Sequence count distribution")
    plt.suptitle("sequence count distribution for " + fastaFile.replace(".", " ").replace("fasta", ""))
    plt.xlim([min(tmp2.keys())-1, 300])
    plt.xlabel('Copy number by variants')
    plt.ylabel('Count')
    plt.yscale('log')
    plt.plot(tmp2.keys(), tmp2.values(), "b")
    # plt.show()
    plt.savefig(destinationFolder  + fastaFile.replace(".fasta", "_sequenceCountDistribution.png"))
    plt.close()


def plotSeqSubstratCountDistribution(seqCountList, fastaFile, destinationFolder = "figures/variant_count_distribution/"):
    seqCountList.sort()
    tmp2 = Counter(seqCountList)
    print("plotting sequence count distribution for ", fastaFile)
    if not os.path.isdir(destinationFolder) : os.makedirs(destinationFolder)

    plt.title("Sequence substrat count distribution")
    plt.suptitle("sequence substrat count distribution for " + fastaFile.replace(".", " ").replace("fasta", ""))
    plt.xlim([min(tmp2.keys())-1, 300])
    plt.xlabel('Copy number by variants')
    plt.ylabel('Count')
    plt.yscale('log')
    plt.plot(tmp2.keys(), tmp2.values(), "b")
    # plt.show()
    plt.savefig(destinationFolder + fastaFile.replace(".fasta", "_sequenceSubstratCountDistribution.png"))
    plt.close()

def plotSeqProduitCountDistribution(seqCountList, fastaFile, destinationFolder = "figures/variant_count_distribution/"):
    seqCountList.sort()
    tmp2 = Counter(seqCountList)
    print("plotting sequence count distribution for ", fastaFile)
    if not os.path.isdir(destinationFolder) : os.makedirs(destinationFolder)

    plt.title("Sequence produit  count distribution")
    plt.suptitle("sequence produit  count distribution for " + fastaFile.replace(".", " ").replace("fasta", ""))
    plt.xlim([min(tmp2.keys())-1, 300])
    plt.xlabel('Copy number by variants')
    plt.ylabel('Count')
    plt.yscale('log')
    plt.plot(tmp2.keys(), tmp2.values(), "b")
    # plt.show()
    plt.savefig(destinationFolder + fastaFile.replace(".fasta", "_sequenceSeqProduitCountDistribution.png"))
    plt.close()


def countMutation(fastaFile, startPaternToDetect1_5 = settings["Suny_beginning"], endPaternToDetect1_5 = settings["sub2Strict"], startPaternToDetect6_9 = settings["Suny_beginning"], endPaternToDetect6_9 = settings["exon"], sourceFolder = "fasta/pre-treated/"):
    mutationDict = {}
    substratDict = {}
    seqProduitDict= {}
    substratDictCount = {}
    substratList = []
    seqRecordDict= {}
    nbMut = lambda suny_seq, mutSeq: sum(ei != ej for ei, ej in zip(suny_seq, mutSeq))
    hamingDistanceList= []
    levenshteinDistanceList= []
    turn = int(fastaFile[fastaFile.find("T0"):fastaFile.find("T0")+3][-1])
    print("getting data")
    
    for record in SeqIO.parse(sourceFolder+ fastaFile, "fasta-2line"):
        seq = str(record.seq)
        mainSeq=""
        substrat=""
        seqProduit=""
        if turn <= 5:
            if re.search(startPaternToDetect1_5, seq) != None:
                startIndex = re.search(startPaternToDetect1_5, seq).start()
                mainSeq= seq[startIndex:startIndex +142]

            if re.search(endPaternToDetect1_5, seq) != None :
                substrat = seq[seq.find(mainSeq) + len(mainSeq):]


            if re.search(startPaternToDetect1_5, seq) != None and re.search(endPaternToDetect1_5, seq) != None :
                startIndex = re.search(startPaternToDetect1_5, seq).start()
                endIndex = re.search(endPaternToDetect1_5, seq).end()
                # endIndex = re.search(settings["sunY_end"], seq).end()
                seqProduit= seq[startIndex:endIndex]
            substratList.append(SeqRecord(Seq(substrat), id = record.id, name=record.name, description=""))
        

            if not mainSeq in mutationDict.keys():
                mutationDict[mainSeq] = 1
                substratDict[mainSeq]= [substrat] 
                # seqRecordDict[mainSeq] = [record] 
                hamingDistanceList.append(nbMut(mainSeq, settings["refSeqSequence"][4:-2])) 
                # levenshteinDistanceList.append(levenshteinDistance(mainSeq, settings["refSeqSequence"][4:-2]))
            else: 
                mutationDict[mainSeq] +=1
                # seqRecordDict[mainSeq].append(record)
                if not substrat in substratDict[mainSeq] : substratDict[mainSeq].append(substrat)

            if not seqProduit in seqProduitDict.keys(): seqProduitDict[seqProduit] = 1
            else: seqProduitDict[seqProduit] += 1
            if not substrat in substratDictCount.keys(): substratDictCount[substrat] = 1
            else: substratDictCount[substrat] +=1


        # elif turn > 5: mainSeq= seq[seq.find(startPaternToDetect)-5:seq.find(startPaternToDetect) + 137]

    print("geting dict data finish")

    #  sorting the dictionnary
    sorted = []
    sortedSubstrat = []
    sortedSeqProduit= []
    for key, val in mutationDict.items(): sorted.append((val, key))
    for key, val in substratDictCount.items(): sortedSubstrat.append((val, key))
    for key, val in seqProduitDict.items(): sortedSeqProduit.append((val, key))
    sorted.sort(reverse=True, key=lambda a: a[0])
    sortedSubstrat.sort(reverse=True, key=lambda a: a[0])
    sortedSeqProduit.sort(reverse=True, key=lambda a: a[0])
    newMutationDict = {}
    newSubstratDict = {}
    newSeqProduitDict= {}
    for val, key in sorted : newMutationDict[key] = val
    for val, key in sortedSubstrat : newSubstratDict[key] = val
    for val, key in sortedSeqProduit: newSeqProduitDict[key] = val

    print("writing csv")
    csvData = {"sequence":newMutationDict.keys(), "seqCount":newMutationDict.values(), "hamingDistance": hamingDistanceList}
    csvData2 = {"sequence":newSubstratDict.keys(), "seqCount":newSubstratDict.values()}
    csvData3 = {"sequence":newSeqProduitDict.keys(), "seqCount":newSeqProduitDict.values()}
    df = pd.DataFrame.from_dict(csvData)
    df2 = pd.DataFrame.from_dict(csvData2)
    df3 = pd.DataFrame.from_dict(csvData3)
    df= df.sort_values(by=["seqCount"], ascending=False)
    df2= df2.sort_values(by=["seqCount"], ascending=False)
    df3= df3.sort_values(by=["seqCount"], ascending=False)
    if not os.path.isdir("results_csv") : os.makedirs("results_csv")
    df.to_csv("results_csv/"+ fastaFile[0:fastaFile.rfind(".")+1]  + "_readCount.csv", index=False,   sep=",")
    df2.to_csv("results_csv/"+ fastaFile[0:fastaFile.rfind(".")+1]  + "_substratCount.csv", index=False,   sep=",")
    df3.to_csv("results_csv/"+ fastaFile[0:fastaFile.rfind(".")+1]  + "_seqProduitCount.csv", index=False,   sep=",")
    print(len(df))
    plotSeqCountDistribution(list(newMutationDict.values()), fastaFile)
    plotSeqSubstratCountDistribution(list(newSubstratDict.values()), fastaFile)
    plotSeqProduitCountDistribution(list(newSeqProduitDict.values()), fastaFile)

    return (newMutationDict, substratDict, newSubstratDict, newSeqProduitDict, substratList, fastaFile)

def writeSubstratOnly(substratList, fastaFile, destinationFolder = "fasta/treated/product_By_Read/"):
    if not os.path.isdir(destinationFolder): os.mkdir(destinationFolder)
    SeqIO.write(substratList,destinationFolder+ fastaFile.replace(".fasta", "_substratOnly.fasta"), "fasta-2line")

def writeSeqMoreThan10TimeFasta(newMutationDict, substratDict, substratDict2, fastaFile, destinationFolder = "fasta/treated/Sequence_Occurence/", destinationFolder2 = "fasta/treated/Product_Occurence_By_Sequence/"):
    seqId = 0
    recordList = []
    substratList = []
    for seq, nbFound in newMutationDict.items():
        if nbFound >= 10 :
            seqName = "seq_" + str(seqId)+ "_N="+ str(nbFound)
            recordList.append(SeqRecord(Seq(seq), id = seqName, name= "", description=""))
            substratList.append(SeqRecord(Seq(seq), id = seqName, name= "", description=""))
            substratID=0
            for record in substratDict[seq]:
                substratName = seqName + "_substrat_" + str(substratID)
                substratList.append(SeqRecord(Seq(record), id = substratName, name ="", description=""))
                substratID+=1

            seqId+=1
    
    if not os.path.isdir(destinationFolder) : os.mkdir(destinationFolder)
    if not os.path.isdir(destinationFolder2) : os.mkdir(destinationFolder2)
    SeqIO.write(recordList,destinationFolder + fastaFile.replace(".fasta", "_more10SeqPresence.fasta"), "fasta-2line")
    SeqIO.write(substratList,destinationFolder2 + fastaFile.replace(".fasta", "_more10SeqPresence_substratPossibilities.fasta"), "fasta-2line")



def writeSeqSubstratMoreThan10TimeFasta(substratDict,fastaFile, destinationFolder = "fasta/treated/product_Occurence/"):
    seqId = 0
    recordList = []
    for seq, nbFound in substratDict.items():
        if nbFound >= 10 :
            seqName = "seq_substrat" + str(seqId)+ "_N="+ str(nbFound)
            recordList.append(SeqRecord(Seq(seq), id = seqName, name= "", description=""))

            seqId+=1
    
    if not os.path.isdir(destinationFolder) : os.mkdir(destinationFolder)
    SeqIO.write(recordList,destinationFolder + fastaFile.replace(".fasta", "_more10SeqSubstratPresence.fasta"), "fasta-2line")

def writeSeqProduitMoreThan10TimeFasta(produitDict, fastaFile, destinationFolder = "fasta/treated/Sequence_With_Sub2_Occurence/"):
    seqId = 0
    recordList = []
    for seq, nbFound in produitDict.items():
        if nbFound >= 10 :
            seqName = "seq_produit" + str(seqId)+ "_N="+ str(nbFound)
            recordList.append(SeqRecord(Seq(seq), id = seqName, name= "", description=""))

            seqId+=1
    
    if not os.path.isdir(destinationFolder) : os.mkdir(destinationFolder)
    SeqIO.write(recordList,destinationFolder + fastaFile.replace(".fasta", "_more10SeqProduitPresence.fasta"), "fasta-2line")



def main():
    if os.path.isfile("fasta/pre-treated/log.txt"): os.remove("fasta/pre-treated/log.txt")

    # fileList = [file for file in os.listdir("fasta/raw_fasta/") if "trimmed" in file and "R1" in file]
    # for file in fileList:
        # seqLengthList = treatFasta(file)
        # plotSeqLengthDistribution(seqLengthList, file)

    print("switching to mutation count")
    if not os.path.isdir("fasta/treated/") : os.mkdir("fasta/treated/")
    fileList = [file for file in os.listdir("fasta/pre-treated/") if "R1" in file]
    for file in fileList:
        # newMutationDict, substratDict, newSubstratDict, newSeqProduitDict, substratList, fastaFile= countMutation(file,endPaternToDetect1_5=settings["exon_N_sub2"])
        newMutationDict, substratDict, newSubstratDict, newSeqProduitDict, substratList, fastaFile= countMutation(file)
        # writeSeqMoreThan10TimeFasta(newMutationDict,substratDict,  newSubstratDict, file)
        # writeSeqSubstratMoreThan10TimeFasta(newSubstratDict, file)
        writeSeqProduitMoreThan10TimeFasta(newSeqProduitDict, file)
        # writeSubstratOnly(substratList, fastaFile)

if __name__ == '__main__':
    main()