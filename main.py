from settings import *
from treatment import *
from treatment_mutagenesis import *
from preTreatment import *
from alignment import *
import re, time, sys, os, platform, subprocess

def main():
    settings = getSettings()
    
    # Pre-treatment
    if settings["pretreatment"]["preTreatFile"]:
        if os.path.isfile("fasta/pre-treated/log.txt"): os.remove("fasta/pre-treated/log.txt")
        if not os.path.isdir("figures/") : os.mkdir("figures/")
        if not os.path.isdir("results_csv/") : os.mkdir("results_csv/")
        if not os.path.isdir("results/") : os.mkdir("results/")
 
        print("starting pre-treatment")
        fileList = [file for file in os.listdir("fasta/raw_fasta/") if not "trimmed" in file and "R1" in file and file.endswith(".fasta")]
        for file in fileList:
            seqLengthList = preTreatmentFasta(file)
        
    print("switching to sequence occurence count function")

    # Directed Evolution

    if not os.path.isdir("fasta/treated/") : os.mkdir("fasta/treated/")
    if not os.path.isdir("fasta/treated/DirectedEvolutionData/") : os.mkdir("fasta/treated/DirectedEvolutionData/")

    # fileList = [file for file in os.listdir("fasta/pre-treated/DirectedEvolutionData/") if "R1" in file]
    
    # for file in fileList:
        # newMutationDict, SubstrateDict, newSubstrateDict, newSeqProductDict, SubstrateList, fastaFile= countSeqOccurences(file)
        # writeSeqFasta(newMutationDict,SubstrateDict, file)
        # writeSubstrateFasta(newSubstrateDict, file)
        # writeSeqProductFasta(newSeqProductDict, file)
        # writeSubstrateReads(SubstrateList, fastaFile)

    # Mutagenesis

    if not os.path.isdir("fasta/treated/MutagenesisData/") : os.mkdir("fasta/treated/MutagenesisData/")

    # fileList = [file for file in os.listdir("fasta/pre-treated/MutagenesisData/") if "R1" in file and file.endswith(".fasta")]
    
    # for file in fileList:
        # sortedMutantSeqList, MutantSeqDict, sortedMutantSeqDict, fastaFile = countMutantSeqOccurence(file)
        # writeSeqMutagenesisFasta(sortedMutantSeqDict, file)

    # alignedFile = align("L447T06.R1_pre-treated_MutagenesisSequences.fasta", 10000000, sourcePath =  "fasta/treated/MutagenesisData/", destinationPath = "fasta/alignment/MutagenesisData/")
    if not os.path.isdir("fasta/alignment/") : os.mkdir("fasta/alignment/")
    # fileList = [file for file in os.listdir("fasta/treated/MutagenesisData/") if file.endswith(".fasta")]
    
    # for file in fileList: 
        # alignedFile = align(file, 10000, sourcePath =  "fasta/treated/MutagenesisData/", destinationPath = "fasta/alignment/MutagenesisData/")

    if not os.path.isdir("figures/MutagenesisData/") : os.mkdir("figures/MutagenesisData/")
    fileList = [file for file in os.listdir("fasta/alignment/MutagenesisData") ]
    for file in fileList:
        if file.endswith("log.txt") : os.remove("fasta/alignment/MutagenesisData/"+ file)

    # mutationCountList, mutationPosCountDict, mutationTypeCountDict= countMutationForMutagenesisData("L447T09.R1_pre-treated_MutagenesisSequences_align.fasta", sourceFolder="fasta/alignment/MutagenesisData/")  
    fileList = [file for file in os.listdir("fasta/alignment/MutagenesisData/")  if file.endswith(".fasta") or file.endswith(".aln")]
    
    for file in fileList:
        # mutationCountList, mutationPosCountDict, mutationTypeCountDict= countMutationForMutagenesisData(file, sourceFolder="fasta/alignment/MutagenesisData/")
        mutationCountList, mutationPosCountDict, mutationTypeCountDict, mutationTypeByPosDict = countMutationForMutagenesisData(file)
        plotMutationDistributionForMutagenesisData(mutationCountList, file)
        plotMutationPosDistributionForMutagenesisData(mutationPosCountDict, file)
        plotMutationTypeByPosDistributionForMutagenesisData(mutationTypeByPosDict, file)

if __name__ == '__main__':
    main()