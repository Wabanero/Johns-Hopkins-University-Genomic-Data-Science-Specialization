''' Dataset: *.fasta
    Exam questions:
    1)  How many records are in the file?
    2)  What are the lengths of the sequences in the file? What is the longest sequence and
        what is the shortest sequence? Is there more than one longest or shortest sequence? What are their identifiers?
    3)  What is the length of the longest ORF in the file?
        What is the identifier of the sequence containing the longest ORF?
        For a given sequence identifier, what is the longest ORF contained
        in the sequence represented by that identifier?
        What is the starting position of the longest ORF in the sequence that contains it?

        Regex Module
        BioPython
        Pandas
        '''

__author__ = "Filippo Bergeretti"
__date__ = "01 June 2019"

# Imports
import pandas as pd
import regex as re
from Bio import SeqIO

# Lists
seq_headers_list = []
seq_list = []
seq_lengths = []
seq_joined = []
row_list = []
datalist = []
datalist_orf = []
seq_idx_list = []
max_count_list = []

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

    # 1)    How many records are in the file?
    x = len(seq_headers_list)
    print("\nNumber of records in file: \t", x)

    # 2a)    What are the lengths of the sequences in the file?
    for position in range(len(seq_list)):
        y = seq_list[position]
        print("\nSequence ", (position + 1), "\n", seq_list[position], "Length : \t", len(y))
    print("\n\n=================================================================================\n\n")

    # 2b)    What is the longest sequence and what is the shortest sequence?
    for position in range(len(seq_list)):
        l = len(seq_list[position])
        seq_lengths.append(l)

    # intialize dataset
    datalist.append(seq_headers_list)
    datalist.append(seq_list)
    datalist.append(seq_lengths)
    df = pd.DataFrame({'Sequence id': seq_headers_list, 'Length': seq_lengths}, index=range(len(seq_list)))
    df_descending_length = df.sort_values('Length', ascending=False)
    print(df_descending_length)

    longest_seq_list = []
    shortest_seq_list = []

    for column in df_descending_length[['Length']]:
        columnSeriesObj = df_descending_length[column]
        print('Colunm Name : \n\t', column)
        print('Column Contents : \n\t', columnSeriesObj.values, "\n\n")
        maxl = int(df_descending_length.iat[0, 1])
        minl = int(df_descending_length.iat[-1, 1])
        print("single max length: ",df_descending_length.iat[0, 1])  # single max length @ df column "Length" (column 1)
        print("single min length: ",df_descending_length.iat[-1, 1])  # single min length @ df column "Length" (column 1)
        print("===============================================================\n\n")

    # Returning index of longest and shortest sequence(s). To have number of lists, add 1 to indexes
    longest_seq_idx = df_descending_length.index[df_descending_length['Length'] == maxl].tolist()
    shortest_seq_idx = df_descending_length.index[df_descending_length['Length'] == minl].tolist()
    print("\nThe longest sequence is sequence at index: \t\t", longest_seq_idx, " ", (df_descending_length.iat[0, 1]),
          "\t\tbp long.")
    print("\nThe shortest sequence is sequence at index: \t", shortest_seq_idx, " ", (df_descending_length.iat[-1, 1]),
          "\t\tbp long.\n\n\n")

    # 3)  What is the length of the longest ORF in the file?
    #     What is the identifier of the sequence containing the longest ORF?
    orf_length = []
    orf_record = []
    orf_seq = []
    orf_pos = []
    orf_frame = []
    records = SeqIO.parse('dna2.fasta', 'fasta')
    print("======\n ORFs\n======")
    for record in records:
        for strand, seq in (1, record.seq), (-1, record.seq.reverse_complement()):
            for frame in range(3):
                index = frame
                while index < len(record) - 6:
                    match = re.match('(ATG(?:\S{3})*?T(?:AG|AA|GA))', str(seq[index:]))
                    if match:
                        orf = match.group()
                        index += len(orf)
                        if len(orf) > 1000:
                            orf_length.append(len(orf))
                            orf_record.append(record.id)
                            orf_seq.append(orf)
                            orf_frame.append(frame)

                            pos = str(record.seq).find(orf) + 1
                            orf_pos.append(pos)
                            print("{}...{} - length {}, strand {}, frame {}, pos {}, name {}".format \
                                      (orf[:6], orf[-3:], len(orf), strand, frame, pos, record.id))
                    else:
                        index += 3
    print(orf_length)
    max_orf_value = max(orf_length)
    max_orf_index = orf_length.index(max_orf_value)
    x = orf_length.index(max(orf_length))  # Longest ORF index in ORF lengths list
    print("\nLongest ORF: \n", "Length : \t", max(orf_length), "\t\tId: \t", orf_record[x], "\nSeq :",
          str(orf_seq[x])[:10], "...", str(orf_seq[x])[-10:])

    # 3b)   What is the starting position of the longest ORF in the sequence that contains it?
    #       Creating dataframe for rapid access
    print("Position: ", orf_pos[x],"\n\n")
    datalist_orf.append(orf_length)
    datalist_orf.append(orf_record)
    datalist_orf.append(orf_pos)
    datalist_orf.append(orf_frame)

    df2 = pd.DataFrame({'Sequence id': orf_record, 'ORF Length': orf_length,'ORF Position': orf_pos,'Frame': orf_frame, }, index=range(len(orf_record)))
    df2_descending_length = df2.sort_values(['ORF Length','Frame'], ascending=False)
    print(df2_descending_length)


