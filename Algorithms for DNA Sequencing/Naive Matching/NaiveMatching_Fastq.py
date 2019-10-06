__author__  = "Filippo Bergeretti"
__date__    = "06 October 2019"


import wget
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def getFile(url, filename):
    '''Genome downloader'''
    wget.download(url, filename)
    print("\n\nGenome downloaded as", filename)
    return filename


#====================================== Naive Matching Algorithm Functions ===========================================

def readGenome(fasta):
    '''Genome data cleaner'''
    genome = ''
    with open(fasta, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
            else:
                print(line)
    return genome

def reverseComplement(seq):
    '''Give the reverse complement (t) of a DNA sequence (seq) using a dictionary of complementation
    A->T; C->G; G->C; T->A'''
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in seq:
        t = complement[base] + t
    return t

def naive(p, t):
    '''naive exact matching algorithm'''
    occurrences = []
    for i in range(len(t) - len(p) + 1):    # loop over alignments
        match = True
        for j in range(len(p)):             # loop over characters
            if t[i+j] != p[j]:              # compare characters on P
                match = False
                break
        if match:
            occurrences.append(i)           # all chars matched; record
    return occurrences

def naive_mismatch(p, t):
    '''naive matching algorithm with 2 mismatches allowed.
    looking for approximate matches for P itself, not its reverse complement.'''
    occurrences = []
    for i in range(len(t) - len(p) + 1):    # loop over alignments
        match = True
        mis = 0                             # added mismatch countdown
        for j in range(len(p)):             # loop over characters
            if t[i + j] != p[j]:            # compare characters on P
                mis += 1
                if mis > 2:
                    match = False
                    break
        if match:
            occurrences.append(i)           # all chars matched; record
    return occurrences

def Naive_Revcomp(P, t):
    '''naive exact matching algorithm that is strand-aware.
    If P is ACT, your function should find occurrences of both ACT and its reverse complement AGT in T.'''
    P_rc = reverseComplement(P)
    if P != P_rc:
        rocc = naive(P_rc, t)
        occ = naive(P, t)
        totalocc = occ + rocc
        return totalocc
    else:
        occ = naive(P,t)
        return occ


#========================================== FASTQ Analysis Functions ===============================================


def readFastq(filename):
    '''Function called in order to obtain a list of sequences and a list of qualities from a single FASTQ file'''
    sequences = []
    qualities = []
    names     = []
    with open(filename) as fh:
        while True:
            nam = fh.readline().rstrip()                # read name line
            nam = nam[:14]
            seq = fh.readline().rstrip()                # read base sequence
            fh.readline()                               # skip placeholder line
            qual = fh.readline().rstrip()               # base quality line
            if len(seq) == 0:
                break
            names.append(nam)
            sequences.append(seq)
            qualities.append(qual)
    return names, sequences, qualities                  # Output lists


def Phred33toQ(qualities):
    '''Function that turns ASCII-encoded quality values into Q. returns an integer'''
    return ord(qualities)-33

def QString(element):
    '''Function that returns Q values for elements of a list. Returns a list of integers'''
    qstring = []
    for character in element:
        q = Phred33toQ(character)
        qstring.append(q)
    return qstring

def Average(lst):
    '''Function to get average of a list'''
    return sum(lst) / len(lst)

def MinQ(read):
    '''Function to get the minimum q value and its offset in order to understand
    which sequencing cycle had problems.'''
    min_q       = min(int(s) for s in read)
    min_index   = read.index(min_q)
    return min_index

def FastqAnalyze(filename):
    '''Fastq analyzer, not using multidimensional lists for debugging and accessibility purposes.
    Pandas dataframe used to help visualize data and improve decision making during coding time.
    Remember that a sequencing cycle corresponds to a particular offset in all the reads'''
    output          = readFastq(filename)
    name_list       = output[0]
    qual_list       = output[2]
    qlist           = []
    qavg            = []
    seq_len         = []
    minidx          = []
    # Dataframe list creation
    for element in qual_list:
        qlist.append(QString(element))                  # Converting to Q values(int)
    for i in range(len(qlist)):
        qavg.append(Average(qlist[i]))                  # Q average value list
        minidx.append(MinQ(qlist[i]))                   # Q minimum value index list
        seq_len.append(len(qlist[i]))                   # seq length
    # Dictionary for Dataframe
    dict = {'Sequence id': name_list,'Q average': qavg,'Sequence length': seq_len,'Min Q index': minidx}
    df = pd.DataFrame(dict)
    print("Fastq Q values analysis for ", filename, " sequences:\n====================================================")
    print(df.sort_values(by=['Min Q index']))           # Sorting for min Q index value, first hint
    # Plotting
    plot_array = np.array(minidx)                       # Setting up a NumPy array
    counts = np.bincount(plot_array)
    highest = np.argmax(counts)                         # Offset problem result
    plt.hist(plot_array, bins = 100)
    plt.title('Offsets')
    plt.xlabel("Sequencing Cycle")
    plt.ylabel("Frequency")
    plt.savefig('OffsetHistogram.png')
    print("\nThe sequencing cycle number ", highest, " has the problem.")
    return qlist



#===================================================== ANSWERS ========================================================

# Set URL for genome download and P,T values. No I/O interface on this one
url                 = "https://d28rh4a8wq0iu5.cloudfront.net/ads1/data/lambda_virus.fa"
Pattern             = "AGGT"
# Function callings
genome_fasta        = getFile(url, "lambda.fasta")
genome              = readGenome(genome_fasta)
occurrences         = Naive_Revcomp(Pattern, genome)
occurrences_mism    = naive_mismatch(Pattern, genome)
print('offset of leftmost occurrence: %d' % min(occurrences))
print('# occurrences: %d' % len(occurrences))
print('\nAllowing up to 2 mismatches:\noffset of leftmost occurrence: %d' % min(occurrences_mism))
print('# occurrences: %d' % len(occurrences_mism))

# Set URL for FASTQ file download
url_2               = "https://d28rh4a8wq0iu5.cloudfront.net/ads1/data/ERR037900_1.first1000.fastq"
# Function calling
human               = getFile(url_2, "human.fastq")
analysis            = FastqAnalyze(human)

# Cleaning
if os.path.exists("lambda.fasta") and ("human.fastq"):
  os.remove("lambda.fasta")
  os.remove("human.fastq")
else:
  print("Error")
