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
    
# function preTreatmentFasta
# this function pretreat the raw fasta file in order 
# to filter sequence of a minimal length
# and that contain both of  the start and end patern
# this function also provide a log file that indicate the amount of sequence filtered and the amount of sequence with one or the other patern
# @param,
# @fastaFile, The raw fasta file to treat 
# @sourceFolder, the source folder where to find fasta file, 
# @destinationFolder, the destination folder where to write the pre-treated fasta
def preTreatmentFasta(fastaFile, sourceFolder = settings["pretreatment"]["sourceFolder"], destinationFolder = settings["pretreatment"]["destinationFolder"]):
    
    refSeq= SeqRecord(Seq(settings["refSeqSequence"]), id = settings["refSeqId"], name = "", description= "")
    seqLengthList = [len(refSeq.seq)]
    seqList = [refSeq]
    
    EndMotifOnly = 0
    BeginningMotifOnly = 0
    BothMotif= 0
    
    file_number = fastaFile[fastaFile.find("T0"):fastaFile.find("T0")+3]
    print(file_number)

    for record in SeqIO.parse(sourceFolder+ fastaFile, "fasta"):
        haveBeginningMotif = False
        haveEndMotif= False
        
        seq = str(record.seq)
        if settings["pretreatment"]["pre_treatment_start_1"] in seq and int(file_number[-1]) <=5: 
            haveBeginningMotif=True
            seq = seq[seq.find(settings["pretreatment"]["pre_treatment_start_1"]) -5:]

        if settings["pretreatment"]["pre_treatment_start_2"] in seq and int(file_number[-1]) >5: 
            haveBeginningMotif=True
            seq = seq[seq.find(settings["pretreatment"]["pre_treatment_start_2"])  -5:]

        if settings["pretreatment"]["pre_treatment_end_1"] in seq and int(file_number[-1]) <=5: 
            haveEndMotif=True
            seq = seq[:seq.rfind(settings["pretreatment"]["pre_treatment_end_1"]) + len(settings["pretreatment"]["pre_treatment_end_1"])]

        if settings["pretreatment"]["pre_treatment_end_2"] in seq and int(file_number[-1]) > 5: 
            haveEndMotif=True
            seq = seq[:seq.rfind(settings["pretreatment"]["pre_treatment_end_2"]) + len(settings["pretreatment"]["pre_treatment_end_2"])]
        
        record.seq = Seq(seq)
        
        if len(seq) > settings["pretreatment"]["length_threshold"] and int(file_number[-1]) <= 5:
            if haveBeginningMotif and not haveEndMotif: BeginningMotifOnly+=1
            if not haveBeginningMotif and haveEndMotif: EndMotifOnly+=1
            if haveEndMotif and haveBeginningMotif: 
                BothMotif +=1
                seqList.append(record)
                seqLengthList.append(len(seq))

        elif len(seq) > settings["pretreatment"]["length_threshold"] and int(file_number[-1]) > 5:
            if haveBeginningMotif and not haveEndMotif: BeginningMotifOnly+=1
            if not haveBeginningMotif and haveEndMotif: EndMotifOnly+=1 
            if haveBeginningMotif: 
                BothMotif +=1
                seqList.append(record)
                seqLengthList.append(len(seq))

    log = "The file " + fastaFile + " after pre-treatment has " + str(len(seqList)) + " sequences \n"
    log+= "including : \n" + str(BeginningMotifOnly) + " sequences have the Beginning motif only \n"
    log+= str(EndMotifOnly) + " sequences have the End motif only \n"
    log+= str(BothMotif) + " sequences have both the Beginning and End motif. \n"
    log+= "The Beginning motif is %s and the End motif is %s." %(settings["pretreatment"]["pre_treatment_start_1"], settings["pretreatment"]["pre_treatment_end_1"])
    log+= "-------------\n"

    if not os.path.isdir(destinationFolder) : os.makedirs(destinationFolder)
    if not os.path.isfile(destinationFolder+"log.txt"): logFile = open(destinationFolder+ "log.txt", "w").close()
    with open(destinationFolder+"log.txt", "a") as logFile:
        logFile.write(log)

    print("write fasta")
    SeqIO.write(seqList, destinationFolder + fastaFile.replace(".fasta", "_pre-treated.fasta"), "fasta-2line")
    return seqLengthList

# function plotSeqLengthDistribution
# This fonction plot as an histogram the sequence length distribution of a fasta file
# @param
# @seqLengthList, the list of sequence length to plot
def plotSeqLengthDistribution(seqLengthList, fastaFile, destinationFolder = "figures/sequence_length_distribution/"):
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

def plotSeqCountDistribution(seqCountList, fastaFile, destinationFolder = "figures/variant_count_distribution/"):
    seqCountList.sort()
    tmp2 = Counter(seqCountList)
    print("Plotting sequence count distribution for ", fastaFile)
    
    if not os.path.isdir(destinationFolder ) : os.makedirs(destinationFolder )

    plt.title("Count distribution")
    plt.suptitle("for " + fastaFile.replace(".", " ").replace("fasta", ""))
    plt.xlim([min(tmp2.keys())-1, 300])
    plt.xlabel('Number of Reads by variants')
    plt.ylabel('Count')
    plt.yscale('log')
    plt.plot(tmp2.keys(), tmp2.values(), "b")
    # plt.show()
    plt.savefig(destinationFolder  + fastaFile.replace(".fasta", "_sequenceCountDistribution.png"))
    plt.close()


def plotSubstrateCountDistribution(seqCountList, fastaFile, destinationFolder = "figures/product_count_distribution/"):
    seqCountList.sort()
    tmp2 = Counter(seqCountList)
    
    print("Plotting Substrate count distribution for ", fastaFile)
    
    if not os.path.isdir(destinationFolder) : os.makedirs(destinationFolder)

    plt.title("Substrate count distribution")
    plt.suptitle(" for " + fastaFile.replace(".", " ").replace("fasta", ""))
    plt.xlim([min(tmp2.keys())-1, 300])
    plt.xlabel('Number of Reads of products')
    plt.ylabel('Count')
    plt.yscale('log')
    plt.plot(tmp2.keys(), tmp2.values(), "b")
    # plt.show()
    plt.savefig(destinationFolder + fastaFile.replace(".fasta", "_SubstrateCountDistribution.png"))
    plt.close()

def plotSeqProductCountDistribution(seqCountList, fastaFile, destinationFolder = "figures/variant_with_product_count_distribution/"):
    seqCountList.sort()
    tmp2 = Counter(seqCountList)
    
    print("Plotting sequence count distribution for ", fastaFile)
    
    if not os.path.isdir(destinationFolder) : os.makedirs(destinationFolder)

    plt.title("Sequence with product count distribution")
    plt.suptitle("for " + fastaFile.replace(".", " ").replace("fasta", ""))
    plt.xlim([min(tmp2.keys())-1, 300])
    plt.xlabel('Number of Reads by variants/products')
    plt.ylabel('Count')
    plt.yscale('log')
    plt.plot(tmp2.keys(), tmp2.values(), "b")
    # plt.show()
    plt.savefig(destinationFolder + fastaFile.replace(".fasta", "_SequenceProductCountDistribution.png"))
    plt.close()


def countSeqOccurences(fastaFile, startPatternToDetect = settings["seqOccuranceCounting"]["countMutation_seq_start"], endPatternToDetect = settings["seqOccuranceCounting"]["countMutation_seq_Substrate"], sourceFolder = "fasta/pre-treated/"):
    seqDict = {}
    productDict = {}
    seqProductDict= {}
    productDictCount = {}
    productList = []
    seqRecordDict= {}
    
    nbMut = lambda ref_seq, mut_seq: sum(ei != ej for ei, ej in zip(ref_seq, mut_seq))
    hamingDistanceList= []
    levenshteinDistanceList= []
    
    file_number = int(fastaFile[fastaFile.find("T0"):fastaFile.find("T0")+3][-1])
    print(file_number)
    print("getting data")
    
    for record in SeqIO.parse(sourceFolder+ fastaFile, "fasta-2line"):
        seq = str(record.seq)
        mainSeq=""
        product=""
        seqProduct=""

        if file_number > 5:
            if re.search(startPatternToDetect, seq) != None:
                startIndex = re.search(startPatternToDetect, seq).start()
                mainSeq = seq[startIndex : startIndex + settings["seqOccuranceCounting"]["countMutation_length"]]

            if re.search(endPatternToDetect, seq) != None :
                product = seq[seq.find(mainSeq) + len(mainSeq):]


            if re.search(startPatternToDetect, seq) != None and re.search(endPatternToDetect, seq) != None :
                startIndex = re.search(startPatternToDetect, seq).start()
                endIndex = re.search(endPatternToDetect, seq).end()
                # endIndex = re.search(settings["sunY_end"], seq).end()
                seqProduct= seq[startIndex:endIndex]
            productList.append(SeqRecord(Seq(product), id = record.id, name=record.name, description=""))
        

            if not mainSeq in seqDict.keys():
                seqDict[mainSeq] = 1
                productDict[mainSeq]= [product] 
                # seqRecordDict[mainSeq] = [record] 
                hamingDistanceList.append(nbMut(mainSeq, settings["refSeqSequence"][4:-2])) 
                # levenshteinDistanceList.append(levenshteinDistance(mainSeq, settings["refSeqSequence"][4:-2]))
            else: 
                seqDict[mainSeq] +=1
                # seqRecordDict[mainSeq].append(record)
                if not product in productDict[mainSeq] : productDict[mainSeq].append(product)

            if not seqProduct in seqProductDict.keys(): seqProductDict[seqProduct] = 1
            else: seqProductDict[seqProduct] += 1
            if not product in productDictCount.keys(): productDictCount[product] = 1
            else: productDictCount[product] +=1


        # elif file_number > 5: mainSeq= seq[seq.find(startPatternToDetect)-5:seq.find(startPatternToDetect) + 137]

    print("geting dict data finish")

    #  sorting the dictionnary
    
    sorted = []
    sortedProduct= []
    sortedSeqProduct= []
    
    for key, val in seqDict.items(): sorted.append((val, key))
    for key, val in productDictCount.items(): sortedProduct.append((val, key))
    for key, val in seqProductDict.items(): sortedSeqProduct.append((val, key))
    
    sorted.sort(reverse=True, key=lambda a: a[0])
    sortedProduct.sort(reverse=True, key=lambda a: a[0])
    sortedSeqProduct.sort(reverse=True, key=lambda a: a[0])
    
    newSeqDict= {}
    newProductDict= {}
    newSeqProductDict= {}
    
    for val, key in sorted : newSeqDict[key] = val
    for val, key in sortedProduct: newProductDict[key] = val
    for val, key in sortedSeqProduct: newSeqProductDict[key] = val

    print("writing csv")
    
    csvData = {"sequence":newSeqDict.keys(), "seqCount":newSeqDict.values(), "hamingDistance": hamingDistanceList}
    csvData2 = {"sequence":newProductDict.keys(), "seqCount":newProductDict.values()}
    csvData3 = {"sequence":newSeqProductDict.keys(), "seqCount":newSeqProductDict.values()}
    
    df = pd.DataFrame.from_dict(csvData)
    df2 = pd.DataFrame.from_dict(csvData2)
    df3 = pd.DataFrame.from_dict(csvData3)
    
    df= df.sort_values(by=["seqCount"], ascending=False)
    df2= df2.sort_values(by=["seqCount"], ascending=False)
    df3= df3.sort_values(by=["seqCount"], ascending=False)
    
    if not os.path.isdir("results_csv") : os.makedirs("results_csv")
    
    df.to_csv("results_csv/"+ fastaFile[0:fastaFile.rfind(".")+1]  + "_readCount.csv", index=False,   sep=",")
    df2.to_csv("results_csv/"+ fastaFile[0:fastaFile.rfind(".")+1]  + "_SubstrateCount.csv", index=False,   sep=",")
    df3.to_csv("results_csv/"+ fastaFile[0:fastaFile.rfind(".")+1]  + "_SeqProductCount.csv", index=False,   sep=",")
    print(len(df))
    
    # plotSeqCountDistribution(list(newSeqDict.values()), fastaFile)
    # plotSubstrateCountDistribution(list(newProductDict.values()), fastaFile)
    # plotSeqProductCountDistribution(list(newSeqProductDict.values()), fastaFile)

    return (newSeqDict, productDict, newProductDict, newSeqProductDict, productList, fastaFile)

def writeSubstrateReads(SubstrateList, fastaFile, destinationFolder = "fasta/treated/product_By_Read/"):
    if not os.path.isdir(destinationFolder): os.mkdir(destinationFolder)
    SeqIO.write(SubstrateList,destinationFolder+ fastaFile.replace(".fasta", "_SubstrateReadsList.fasta"), "fasta-2line")

def writeSeqFasta(newMutationDict, SubstrateDict, SubstrateDict2, fastaFile, destinationFolder = "fasta/treated/Sequence_Occurence/", destinationFolder2 = "fasta/treated/Product_Occurence_By_Sequence/"):
    seqId = 0
    recordList = []
    SubstrateList = []
    file_number = fastaFile[fastaFile.find("T0"):fastaFile.find("T0")+3][-1]
    
    for seq, nbFound in newMutationDict.items():
        if nbFound >= settings["seqOccuranceCounting"]["minimumOccurence"]:
            seqName = "seq_" + str(seqId)+ "_"+ file_number + "_N="+ str(nbFound)
            recordList.append(SeqRecord(Seq(seq), id = seqName, name= "", description=""))
            SubstrateList.append(SeqRecord(Seq(seq), id = seqName, name= "", description=""))
            SubstrateID=0
            for record in SubstrateDict[seq]:
                SubstrateName = seqName + "_Substrate_" + str(SubstrateID)
                SubstrateList.append(SeqRecord(Seq(record), id = SubstrateName, name ="", description=""))
                SubstrateID+=1

            seqId+=1
    
    if not os.path.isdir(destinationFolder) : os.mkdir(destinationFolder)
    if not os.path.isdir(destinationFolder2) : os.mkdir(destinationFolder2)
    SeqIO.write(recordList,destinationFolder + fastaFile.replace(".fasta", "SeqOccurence.fasta"), "fasta-2line")
    SeqIO.write(SubstrateList,destinationFolder2 + fastaFile.replace(".fasta", "SubstratePossibilitiesbySequence.fasta"), "fasta-2line")



def writeSubstrateFasta(SubstrateDict,fastaFile, destinationFolder = "fasta/treated/product_Occurence/"):
    seqId = 0
    recordList = []
    file_number = fastaFile[fastaFile.find("T0"):fastaFile.find("T0")+3]
    for seq, nbFound in SubstrateDict.items():
        if nbFound >= settings["seqOccuranceCounting"]["minimumOccurence"] :
            seqName = "seq_Substrate" + str(seqId) + "_"+ file_number + "_N="+ str(nbFound)
            recordList.append(SeqRecord(Seq(seq), id = seqName, name= "", description=""))

            seqId+=1
    
    if not os.path.isdir(destinationFolder) : os.mkdir(destinationFolder)
    SeqIO.write(recordList,destinationFolder + fastaFile.replace(".fasta", "_ProductOccurance.fasta"), "fasta-2line")

def writeSeqProductFasta(productDict, fastaFile, destinationFolder = "fasta/treated/Sequence_With_Sub2_Occurence/"):
    seqId = 0
    recordList = []
    file_number = fastaFile[fastaFile.find("T0"):fastaFile.find("T0")+3]
    for seq, nbFound in productDict.items():
        if nbFound >= settings["seqOccuranceCounting"]["minimumOccurence"] :
            seqName = "seq_product" + str(seqId) + "_"+ file_number  + "_N="+ str(nbFound)
            recordList.append(SeqRecord(Seq(seq), id = seqName, name= "", description=""))

            seqId+=1
    
    if not os.path.isdir(destinationFolder) : os.mkdir(destinationFolder)
    SeqIO.write(recordList,destinationFolder + fastaFile.replace(".fasta", "_Sequence_with_product_Occurence.fasta"), "fasta-2line")



def countMutation(fastaFile, sourceFolder = "alignment/"):
    # nbMut = lambda ref_seq, mut_seq: sum(ei != ej for ei, ej in zip(ref_seq, mut_seq))
    mutationPosCountDict= {}
    mutationCountList =[]
    mutationTypeCountDict = {"A":0,"T":0,"C":0,"G":0, "insertion":0, "deletion":0}
    recordList = list(SeqIO.parse(sourceFolder + fastaFile, "fasta"))
    refSeq = str(recordList[0].seq).upper()
    for i in range(1, len(refSeq)+1): mutationPosCountDict[i] =0
    for record in recordList:
        # nbMutation = nbMut(refSeq, record.seq)
        seq= str(record.seq).upper()
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

# function plotMutationDistribution
# This fonction plot as an histogram mutation count distribution of a fasta file
# @param
# @mutationCountList, the list of mutation count to plot
def plotMutationDistribution(mutationCountList, fastaFile, destinationFolder = "figures/mutation_count_distribution/"):
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