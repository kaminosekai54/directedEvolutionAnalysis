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

def align(fastafile, limit, sourcePath = "fasta/treated/Sequence_With_Sub2_Occurence/", destinationPath="fasta/alignment/"):
    refSeq= SeqRecord(Seq(settings["refSeqSequence"]), id = settings["refSeqId"], name = "", description= "")
    recordList = [refSeq]
    for record in SeqIO.parse(sourcePath + fastafile, "fasta"):
        if len(recordList)+1 < limit : recordList.append(record)
        else : break

    algo = "-align"
    if len(recordList) >= 200 : algo = "-super5"
    if not os.path.isdir(destinationPath): os.mkdir(destinationPath)
    alignedFile= destinationPath + fastafile.replace(".fasta", "_align.fasta")
    tmpFile = destinationPath + fastafile
    SeqIO.write(recordList, tmpFile, "fasta")
    print("tmp file written")
    osName = platform.system()
    mafftExe = ""
    if osName == "Linux": mafftExe = settings["linuxMafftPath"]
    # elif osName == "Windows": mafftExe =os.path.abspath(settings["winMafftPath"])
    elif osName == "Windows": mafftExe =os.path.abspath(settings["winMafftPath"])

    mafft_cline = mafftExe + " --auto --out "+ os.path.abspath(alignedFile) + " " + os.path.abspath(tmpFile)
    print("Sometimes the subproccess get locked for no reason. If it happens, please quit and run the following command manually :")
    print(mafft_cline)
    print("Processing alignment")
    # child= subprocess.Popen(str(mafft_cline), stdout = subprocess.PIPE, stderr=subprocess.PIPE,          shell = (sys.platform!="win32"))
    child= subprocess.check_output(str(mafft_cline),text=True)
    # child.wait()
    
    if os.path.isfile(alignedFile):print("Alignment finished : the aligned file is " + alignedFile)
    else: print("It seems that an error occured, the alignement could not be done.")
    
    if os.path.isfile(tmpFile): os.remove(tmpFile)
    return alignedFile

# def main():
    # align("L447T01.R1.fasta", 2000)


# if __name__ == '__main__':
    # main()