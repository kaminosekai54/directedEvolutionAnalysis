from settings import *
import  os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from treatment import *



def preTreatmentFasta(fastaFile, sourceFolder = settings["pretreatment"]["sourceFolder"], destinationFolder = settings["pretreatment"]["destinationFolder"]):
    if not os.path.isdir(destinationFolder) : os.mkdir(destinationFolder)
    file_number = fastaFile[fastaFile.find("T0"):fastaFile.find("T0")+3]
    print(file_number)
    if int(file_number[-1]) <=5:
        destinationFolder += '/DirectedEvolutionData/'
        seqLengthList = preTreatDirectedEvolutionData(fastaFile, destinationFolder=destinationFolder)
        if settings["pretreatment"]["plotSeqLengthDistribution"] : plotSeqLengthDistribution(seqLengthList, fastaFile, destinationFolder = "figures/DirectedEvolutionData/sequence_length_distribution/")
    elif int(file_number[-1]) > 5:
        destinationFolder += '/MutagenesisData/'
        seqLengthList=preTreatMutagenesisData(fastaFile, destinationFolder=destinationFolder)
        if settings["pretreatment"]["plotSeqLengthDistribution"] : plotSeqLengthDistribution(seqLengthList, fastaFile, destinationFolder = "figures/MutagenesisData/sequence_length_distribution/")


# function preTreatDirectedEvolutionData
# this function pretreat the raw fasta file in order 
# to filter sequence of a minimal length
# and that contain both of  the start and end patern
# this function also provide a log file that indicate the amount of sequence filtered and the amount of sequence with one or the other patern
# @param,
# @fastaFile, The raw fasta file to treat 
# @sourceFolder, the source folder where to find fasta file, 
# @destinationFolder, the destination folder where to write the pre-treated fasta
def preTreatDirectedEvolutionData(fastaFile, sourceFolder = settings["pretreatment"]["sourceFolder"], destinationFolder = settings["pretreatment"]["destinationFolder"]):
    
    refSeq= SeqRecord(Seq(settings["refSeqSequence"]), id = settings["refSeqId"], name = "", description= "")
    seqLengthList = [len(refSeq.seq)]
    seqList = [refSeq]
    
    EndMotifOnly = 0
    BeginningMotifOnly = 0
    BothMotif = 0
    
    file_number = fastaFile[fastaFile.find("T0"):fastaFile.find("T0")+3]
    print(file_number)
    baseLengthFile = len(list(SeqIO.parse(sourceFolder+ fastaFile, "fasta")))

    for record in SeqIO.parse(sourceFolder+ fastaFile, "fasta"):
        haveBeginningMotif = False
        haveEndMotif= False
        
        seq = str(record.seq)
        if settings["pretreatment"]["pre_treatment_start_1"] in seq: 
            haveBeginningMotif=True
            seq = seq[seq.find(settings["pretreatment"]["pre_treatment_start_1"]) -5:]

        if settings["pretreatment"]["pre_treatment_end_1"] in seq : 
            haveEndMotif=True
            seq = seq[:seq.rfind(settings["pretreatment"]["pre_treatment_end_1"]) + len(settings["pretreatment"]["pre_treatment_end_1"])]

        record.seq = Seq(seq)
        
        if len(seq) > settings["pretreatment"]["length_threshold"] and int(file_number[-1]) <= 5:
            if haveBeginningMotif and not haveEndMotif: BeginningMotifOnly+=1
            if not haveBeginningMotif and haveEndMotif: EndMotifOnly+=1
            if haveEndMotif and haveBeginningMotif: 
                BothMotif +=1
                seqList.append(record)
                seqLengthList.append(len(seq))


    log = "The file " + fastaFile + " after pre-treatment has " + str(len(seqList)) + " sequences \n"
    log+= "including : \n" + str(BeginningMotifOnly) + " sequences have the Beginning motif only \n Which represent " + str(round(BeginningMotifOnly/baseLengthFile*100, 2)) + " % of the raw fasta file \n"
    log+= str(EndMotifOnly) + " sequences have the End motif only \n Which represent " + str(round(EndMotifOnly/baseLengthFile*100, 2)) + " % of the raw fasta file \n"
    log+= str(BothMotif) + " sequences have both the Beginning and End motif. \nWhich represent " + str(round(BothMotif/baseLengthFile*100, 2)) + " % of the raw fasta file \n"
    log+= "The Beginning motif is %s and the End motif is %s." %(settings["pretreatment"]["pre_treatment_start_1"], settings["pretreatment"]["pre_treatment_end_1"])
    log+= "-------------\n"

    if not os.path.isdir(destinationFolder) : os.makedirs(destinationFolder)
    if not os.path.isfile(destinationFolder+"log.txt"): logFile = open(destinationFolder+ "log.txt", "w").close()
    with open(destinationFolder+"log.txt", "a") as logFile:
        logFile.write(log)

    print("write fasta")
    SeqIO.write(seqList, destinationFolder + fastaFile.replace(".fasta", "_pre-treated.fasta"), "fasta-2line")
    return seqLengthList



# function preTreatMutagenesisData
# this function pretreat the raw fasta file in order 
# to filter sequence of a minimal length
# and that contain both of  the start and end patern
# this function also provide a log file that indicate the amount of sequence filtered and the amount of sequence with one or the other patern
# @param,
# @fastaFile, The raw fasta file to treat 
# @sourceFolder, the source folder where to find fasta file, 
# @destinationFolder, the destination folder where to write the pre-treated fasta
def preTreatMutagenesisData(fastaFile, sourceFolder = settings["pretreatment"]["sourceFolder"], destinationFolder = settings["pretreatment"]["destinationFolder"]):
    
    refSeq= SeqRecord(Seq(settings["refSeqSequence"]), id = settings["refSeqId"], name = "", description= "")
    seqLengthList = [len(refSeq.seq)]
    seqList = [refSeq]
    
    EndMotifOnly = 0
    BeginningMotifOnly = 0
    BothMotif = 0
    
    file_number = fastaFile[fastaFile.find("T0"):fastaFile.find("T0")+3]
    print(file_number)
    baseLengthFile = len(list(SeqIO.parse(sourceFolder+ fastaFile, "fasta")))
    for record in SeqIO.parse(sourceFolder+ fastaFile, "fasta"):
        haveBeginningMotif = False
        haveEndMotif= False
        
        seq = str(record.seq)

        if settings["pretreatment"]["pre_treatment_start_2"] in seq: 
            haveBeginningMotif=True
            seq = seq[seq.find(settings["pretreatment"]["pre_treatment_start_2"])  :]
            # seq = seq[seq.find(settings["pretreatment"]["pre_treatment_start_2"])  -9:]


        if settings["pretreatment"]["pre_treatment_end_2"] in seq: 
            haveEndMotif=True
            seq = seq[:seq.rfind(settings["pretreatment"]["pre_treatment_end_2"]) + len(settings["pretreatment"]["pre_treatment_end_2"])]
        
        record.seq = Seq(seq)
        

        if len(seq) > settings["pretreatment"]["length_threshold"]:
            if haveBeginningMotif and not haveEndMotif: BeginningMotifOnly+=1
            if not haveBeginningMotif and haveEndMotif: EndMotifOnly+=1 
            if haveBeginningMotif and haveEndMotif: 
                BothMotif +=1
                seqList.append(record)
                seqLengthList.append(len(seq))

    log = "The file " + fastaFile + " after pre-treatment has " + str(len(seqList)) + " sequences \n"
    log+= "including : \n" + str(BeginningMotifOnly) + " sequences have the Beginning motif only \n Which represent " + str(round(BeginningMotifOnly/baseLengthFile*100, 2)) + " % of the raw fasta file \n"
    log+= str(EndMotifOnly) + " sequences have the End motif only \n Which represent " + str(round(EndMotifOnly/baseLengthFile*100, 2)) + " % of the raw fasta file \n"
    log+= str(BothMotif) + " sequences have both the Beginning and End motif. \nWhich represent " + str(round(BothMotif/baseLengthFile*100, 2)) + " % of the raw fasta file \n"
    log+= "The Beginning motif is %s and the End motif is %s." %(settings["pretreatment"]["pre_treatment_start_2"], settings["pretreatment"]["pre_treatment_end_2"])
    log+= "-------------\n"

    if not os.path.isdir(destinationFolder) : os.makedirs(destinationFolder)
    if not os.path.isfile(destinationFolder+"log.txt"): logFile = open(destinationFolder+ "log.txt", "w").close()
    with open(destinationFolder+"log.txt", "a") as logFile:
        logFile.write(log)

    print("write fasta")
    SeqIO.write(seqList, destinationFolder + fastaFile.replace(".fasta", "_pre-treated.fasta"), "fasta-2line")
    return seqLengthList