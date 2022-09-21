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
 
        print("starting pre-treatment")

        fileList = [file for file in os.listdir("fasta/raw_fasta/") if "trimmed" in file and "R1" in file]
        for file in fileList:
            seqLengthList = preTreatmentFasta(file)
        
    print("switching to sequence occurence count function")

    # Directed Evolution

    if not os.path.isdir("fasta/treated/") : os.mkdir("fasta/treated/")
    # if not os.path.isdir("fasta/treated/DirectedEvolutionData/") : os.mkdir("fasta/treated/DirectedEvolutionData/")

    fileList = [file for file in os.listdir("fasta/pre-treated/DirectedEvolutionData/") if "R1" in file]
    
    for file in fileList:
        plotSeqLengthDistribution(seqLengthList, file)
        plotSeqCountDistribution(list(newSeqDict.values()), fastaFile)
        plotSubstrateCountDistribution(list(newProductDict.values()), fastaFile)
        plotSeqProductCountDistribution(list(newSeqProductDict.values()), fastaFile)

        # newMutationDict, SubstrateDict, newSubstrateDict, newSeqProductDict, SubstrateList, fastaFile= countSeqOccurences(file, startPatternToDetect=settings["pretreatment"]["pre_treatment_start_2"] , endPatternToDetect=settings["pretreatment"]["pre_treatment_end_2"] )
        newMutationDict, SubstrateDict, newSubstrateDict, newSeqProductDict, SubstrateList, fastaFile= countSeqOccurences(file)
        writeSeqFasta(newMutationDict,SubstrateDict, newSubstrateDict, file)
        writeSubstrateFasta(newSubstrateDict, file)
        writeSeqProductFasta(newSeqProductDict, file)
        writeSubstrateReads(SubstrateList, fastaFile)

    # Mutagenesis

    if not os.path.isdir("fasta/treated/MutagenesisData/") : os.mkdir("fasta/treated/MutagenesisData/")

    fileList = [file for file in os.listdir("fasta/pre-treated/MutagenesisData/") if "R1" in file]
    
    for file in fileList:
        # plotSeqLengthDistribution(seqLengthList, file)
        sortedMutantSeqList, MutantSeqDict, sortedMutantSeqDict, fastaFile = countMutantSeqOccurence(file)
        writeSeqMutagenesisFasta(sortedMutantSeqDict, file)


    if not os.path.isdir("fasta/alignment/") : os.mkdir("fasta/alignment/")
    fileList = [file for file in os.listdir("fasta/treated/MutagenesisData/") if file.endswith(".fasta")]
    
    for file in fileList: 
        alignedFile = align(file,1000000000, sourcePath =  "fasta/treated/MutagenesisData/", destinationPath = "fasta/alignment/MutagenesisData/")

    if not os.path.isdir("figures/MutagenesisData/") : os.mkdir("figures/MutagenesisData/")
    fileList = [file for file in os.listdir("fasta/alignment/MutagenesisData") ]
    
    for file in fileList: 
        mutationCountList, mutationPosCountDict, mutationTypeCountDict= countMutation(file)
        plotMutationDistribution(mutationCountList, file)

if __name__ == '__main__':
    main()