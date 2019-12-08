__author__  = "Filippo Bergeretti"
__date__    = "08 December 2019"

import itertools
import random


def randomString(stringLength,n):
    """Generate a random string of fixed length """
    seqlist = []
    i = 0
    while i < n:
        bases = "ATCG"
        seq = ''.join([random.choice(bases)for n in range(stringLength)])
        seqlist.append(seq)
        i+=1
    return seqlist


def Overlap(a, b, min_length):
    '''function to find the overlap between two DNA sequences (strings a,b) with a specific minimum overlap'''
    start = 0                                   # Index start
    while True:
        start = a.find( b[:min_length], start)  # Look in a for the prefix of min length of b
        if start == -1:
            return 0                            # Returns 0 if no overlap detected
        if b.startswith(a[start:]):             # Verify that the prefix of b = suffix of a starting at position start.
            return len(a) -start                # The len(a)- start is then the length of the overlap
        start += 1                              # Increment start so that it doesn't assess the same location over again


# ========================================= MAIN FUNCTION ==============================================================

def Scs(sequences):
    '''implementing the shortest common super string using just a brute force method to check all the possible
    combinations of reads. Using the overlap function and the permutations function from previous practicals.
    n! iterations with n = number of sequences'''

    shortest_superstring = None

    for ss_permutation in itertools.permutations(sequences):        # Every possible permutation of strings(sequences)
        superstring = ss_permutation[0]                             # First superstring = first seq of permutation
        for i in range(len(sequences)-1):                           # Looping through seq list, find overlaps
            a = ss_permutation[i]                                   # Seq a
            b = ss_permutation[i+1]                                 # Seq b
            min_length = 1
            overlap_length = Overlap(a, b, min_length)              # Overlap function called between a, b
            superstring += ss_permutation[i+1][overlap_length:]     # append the part of the next string which doesn't overlap (suffix)
        if shortest_superstring is None or len(superstring) < len(shortest_superstring):  # Testing shortest
            shortest_superstring = superstring

    return shortest_superstring


dna_len = 8
n = 3
seq = randomString(dna_len, n)
print(dna_len, "bp sequences:")
print(' '.join(seq))
shortest = Scs(seq)
print("Shortest common superstring is:\t",shortest, len(shortest), "bp")
