''' Dataset: *.fasta
     A repeat is a substring of a DNA sequence that occurs in multiple copies (more than one) somewhere in the sequence.
     Although repeats can occur on both the forward and reverse strands of the DNA sequence, we will only consider
     repeats on the forward strand here. Also we will allow repeats to overlap themselves. For example, the sequence
     ACACA contains two copies of the sequence ACA - once at position 1 (index 0 in Python), and once at position 3.
     Given a length n, your program should be able to identify all repeats of length n in all sequences in the FASTA
     file. Your program should also determine how many times each repeat occurs in the file, and which is the most
     frequent repeat of a given length.

    Exam questions:
    1)  Find the most frequently occurring repeat of length n in all sequences. How many times does it occur in all?
    2)  Find all repeats of length 12 in the input file. Let's use Max to specify the number of copies
        of the most frequent repeat of length 12. How many different 12-base sequences occur Max times?
    3)  Which one of the following repeats of length 7 has a maximum number of occurrences? (a,b,c or d)

        pandas
        '''

__author__ = "Filippo Bergeretti"
__date__ = "07 June 2019"

# Imports
import pandas as pd

# Lists
seq_headers_list = []
seq_list = []
seq_idx_list = []
seq_lengths = []
seq_joined = []
row_list = []
total_kmers = []
max_count_list = []
patternlist = []
datalist = []
kmer_freq = []
n_kmers = []
global_patternlist = []
global_max_count = []

# Opening file
with open("dna2.fasta", "r+") as file:
    for row in file:
        row_list.append(row.strip())

# Cleaning data and listing strings from text rows.
    i = 0
    while i < len(row_list):
        if ">" in row_list[i]:
            seq_headers_list.append(row_list[i])
        if ">" in row_list[i - 1]:
            seq_list.append(row_list[i])
        if ">" not in row_list[i - 1] and ">" not in row_list[i]:
            seq_list.append(row_list[i])
            space = ""
            joined = space.join(seq_list[-2:])
            del seq_list[-2:]
            seq_list.append(joined)

        i += 1

# Sequences obtained as list
print("\nSequence list from *.fasta:\n", seq_list, "\n")
for seq in seq_list:
    seq_idx_list.append(seq_list.index(seq))
print("Indexes: ", seq_idx_list, "\n")

# Sequence id list
x = len(seq_headers_list)
print("\nNumber of records in file: \t", x)
# Sequence length list
for position in range(len(seq_list)):
    l = len(seq_list[position])
    seq_lengths.append(l)




# Looping sequence list to find k-mers

k = 6
i = 0
newdict = {}

for text in seq_list:
    print("\n\nScanning k-mers in sequence with index", seq_idx_list[i], ": \n", seq_list[i], )

    def Count(text, pattern):
        count = 0
        text_length = len(text)
        pattern_length = len(pattern)
        for i in range(text_length - pattern_length + 1):  # +1 due to indexing starting from 0
            if text[i: i + pattern_length] == pattern:
                count += 1
        return count


    def FrequentWords(text, k):
        frequency_patterns = []  # How many times a pattern is encountered
        pattern_counts = []  # number of kmers found
        text_length = len(text)
        print("DNA lenght = " + str(text_length) + " bp.")

        # For cycle, finding the patterns
        for i in range(text_length - k + 1):  # +1 due to indexing starting from 0
            pattern = text[i: i + k]  # Saving pattern
            number = Count(text, pattern)  # Calling the Count function again to count +1
            pattern_counts.append(number)  # Adding the +1 to the total number of kmers

        max_count = max(pattern_counts)  # Selecting the most frequent kmers
        print("\t\t\t===================================================================================\n\t\t\tThe most frequent value is " + str(max_count))

        print("\n\nMost frequent pattern(s):")
        for i in range(text_length - k + 1):  # Array appearance count for most frequent kmers
            if pattern_counts[i] == max_count:
                frequency_patterns.append(text[i: i + k])
        frequency_patterns = set(frequency_patterns)  # Generating a set of the characters in
        patternlist = list(frequency_patterns)
        return list(patternlist), max_count  # Result of the function is a list

    kmer_data = FrequentWords(text, k)
    print(kmer_data[0])
    print(kmer_data[1])
    global_patternlist.append(kmer_data[0])
    global_max_count.append(kmer_data[1])
    i+=1





print(global_patternlist)
print(global_max_count)



# intialize dataset
datalist.append(seq_headers_list)
datalist.append(seq_list)
datalist.append(seq_lengths)
df = pd.DataFrame({'Length': seq_lengths,'N k-mer': global_max_count, 'k-mer sequences' : global_patternlist }, index=range(len(seq_list)))#,'global_max_count': kmer_freq, 'N k-mer': n_kmers}, index=range(len(seq_list)))
print(df)















