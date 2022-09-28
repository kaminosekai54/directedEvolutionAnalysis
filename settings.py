settings = {
    "refSeqId": "SunY_sequence",
"refSeqSequence": "GTGGTAAAATCTGCCTAAACGGGGAAACTCTCACTGAGACAATCCCGTGCTAAATCAGCAATAGCTGTAAATGCCTAACGACTACACGGTAGACAACTCTAAGAGTTGAAGGTATAGTCTAAACTGCAAGGTGACTTGCAGATATCGG",

    "pretreatment" :{ # this variable are concerning only the pretreatment part
    "sourceFolder":"fasta/raw_fasta/", 
    "destinationFolder":"fasta/pre-treated/", 
"pre_treatment_start_1":"TAAAATCTGCCTAAACGG", # SunY beginning from nucleotide 5
"pre_treatment_start_2":"GTGGTAAAAT", # SunY beginning from nucleotide 9
# "pre_treatment_start_2":"TCTGCCTAAACGG", # SunY beginning from nucleotide 9
"pre_treatment_end_1":"AACTTCAAATATCTTCGGAACTCA", # Sub2
# "pre_treatment_end_2":"TCCCTATCTCTCTAT", # exon_OT
"pre_treatment_end_2":"CAGATATCGG", # exon_OT
"length_threshold" : 100,
"preTreatFile" : False,
"plotSeqLengthDistribution" : True,
    },

"seqOccurenceCounting" : { # this part concerned only the counting and the writting of the output file for Directed Evolution data
"countMutation_seq_start": "TAAAATCTGCCTAAACGG", # SunY beginning from nucleotide 5
"countMutation_seq_end": "CAGATATCGG", # SunY end
"countMutation_seq_Substrate": "CAGATATCGGAACTTCAAATATCTTCGGAACTCA", # SunY end with Sub2
"countMutation_length" : 142,
"minimumOccurence" : 10,
},
"seqOccurenceCountingMutagenesis" : { # this part concerned only the counting and the writting of the output file for Mutagenesis data
# "countMutation_seq_start": "TCTGCCTAAACGG", # SunY beginning from nucleotide 9
"countMutation_seq_start": "GTGGTAAAAT", # SunY beginning from nucleotide 9
# "countMutation_seq_end": "TCCCTATCTCTCTAT", # exon_OT
"countMutation_seq_end": "CAGATATCGG", # exon_OT
"countMutation_length" : 142,
"minimumOccurence" : 1,
},

"winMafftPath" : "mafft-win/mafft.bat",
"linuxMafftPath" : "/usr/bin/mafft",

#"Suny_beginning" : "TAAAATCTGCCTAAACGG",
#"Suny_beginning2" : "TCTGCCTAAACGG",
#"sub2": "AACTTCAAATATCTTCGGAACTCA",
#"sub2Strict": "CAGATATCGGAACTTCAAATATCTTCGGAACTCA",
#"exon_OT": "TCCCTATCTCTCTAT",
#"exon_CL": "TGGTGACATGCAGATATCGG",
#"exon_N_sub2" : (r"CAGATATCGGTCCCTATCTCTCTAT[A,C,C,G]{1,20}AACTTCAAATATCTTCGGAACTCA"),
}


def getSettings():
    return settings 