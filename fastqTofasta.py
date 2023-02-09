import gzip
import re, os, time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO

def fastqToFasta(fastqFile, sourceFolder = "raw_data/", destinationFolder = "fasta/raw_fasta/"):
    if not os.path.isdir(destinationFolder) : os.makedirs(destinationFolder)

    fastaFile = destinationFolder+ fastqFile.replace(".fastq.gz", ".fasta")
    if os.pathth.isfile(fastaFile) : return 
    with gzip.open(sourceFolder + fastqFile, "rt") as fastq:
        SeqIO.convert(fastq, "fastq", fastaFile, "fasta")
    # os.remove(fastqFile)
    print("convertion finish for ", fastqFile)

def convertAll(path = "raw_data/"):
    sourceFolder = str(input("please enter the path to the source folder containing the fastQ files"))
    destinationFolder = str(input("please enter the path where to store the converted fasta"))
    fileList = [file for file in os.listdir(sourceFolder) if file.endswith(".fastq.gz")]
    for file in fileList:
        fastqToFasta(file, sourceFolder=sourceFolder, destinationFolder=destinationFolder)

convertAll()