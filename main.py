from fileinput import filelineno
from settings import *
from treatment import *
from alignment import *
import re, time, sys, os, platform, subprocess

def main():
    settings = getSettings()
    if os.path.isfile("fasta/pre-treated/log.txt"): os.remove("fasta/pre-treated/log.txt")

    # fileList = [file for file in os.listdir("fasta/raw_fasta/") if "trimmed" in file and "R1" in file]
    # for file in fileList:
        # seqLengthList = preTreatmentFasta(file)
        # plotSeqLengthDistribution(seqLengthList, file)

    print("switching to sequence occurence count function")

    if not os.path.isdir("fasta/treated/") : os.mkdir("fasta/treated/")

    # fileList = [file for file in os.listdir("fasta/pre-treated/") if "R1" in file]
    
    # for file in fileList:
        # newMutationDict, SubstrateDict, newSubstrateDict, newSeqProductDict, SubstrateList, fastaFile= countSeqOccurences(file, startPatternToDetect=settings["pretreatment"]["pre_treatment_start_2"] , endPatternToDetect=settings["pretreatment"]["pre_treatment_end_2"] )
        # newMutationDict, SubstrateDict, newSubstrateDict, newSeqProductDict, SubstrateList, fastaFile= countSeqOccurences(file)
        # writeSeqFasta(newMutationDict,SubstrateDict,  newSubstrateDict, file)
        # writeSubstrateFasta(newSubstrateDict, file)
        # writeSeqProductFasta(newSeqProductDict, file)
        # writeSubstrateReads(SubstrateList, fastaFile)

    fileList = [file for file in os.listdir("fasta/treated/Sequence_With_Sub2_Occurence/") if file.endswith(".fasta")]
    for file in fileList: 
        alignedFile =align(file,1000000000)

    fileList = [file for file in os.listdir("alignment/") ]
    for file in fileList: 
        mutationCountList, mutationPosCountDict, mutationTypeCountDict= countMutation(file)
        plotMutationDistribution(mutationCountList, file)

if __name__ == '__main__':
    main()