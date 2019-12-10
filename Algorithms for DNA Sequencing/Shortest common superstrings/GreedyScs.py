__author__  = "Filippo Bergeretti"
__date__    = "09 December 2019"

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

def PickMaxOverlap(reads, k):
    '''Return a pair of reads from the list with a maximal suffix/prefix overlap >= k.  Returns overlap length 0 if
    there are no such overlaps.'''
    reada, readb = None, None                           # Reads definition
    best_overlap_length = 0                             # Best overlap length initialised to zero
    for a, b in itertools.permutations(reads, 2):       # For each pair of reads a-b, calculate the overlap length
        overlap_length = Overlap(a, b, min_length=k)
        if overlap_length > best_overlap_length:        # If a better value than previous one is found, save it.
            reada = a                                   # Saving pair of reads
            readb = b
            best_overlap_length = overlap_length
    return reada, readb, best_overlap_length



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

def GreedyScs(sequences, k):
    '''Greedy shortest-common-superstring merge. Repeat until no edges (overlaps of length >= k) remain. The algorithm
    will make a series of decisions and at each decision point, it will choose the option that reduces the length of the
    eventual superstring the most. Picking the longest overlap seems to make sense because the longer the overlaps
    between the strings, the shorter the final string will be. The greedy algorithm can sometimes report a superstring,
    a string that does contain all of the input strings, but it's not necessarily the shortest common superstring.'''


    read_a, read_b, olen = PickMaxOverlap(sequences, k)     # Calculate reads A and B, with the best overlap from sequences
    while olen > 0:                                         # As long as overlap length > 0, replacing A and B, in the set of sequences, with their overlap
        sequences.remove(read_a)
        sequences.remove(read_b)
        sequences.append(read_a + read_b[olen:])            # Appending A plus the suffix of read B
        read_a, read_b, olen = PickMaxOverlap(sequences, k) # Recalculate with the greatest overlap obtained.
    return ''.join(sequences)


dna_len     = 4
n           = 10
min_overlap = 1
seq         = randomString(dna_len, n)
print(dna_len, "bp sequences:")
print(' '.join(seq))
shortest = GreedyScs(seq, min_overlap)
print("Shortest common superstring by greedy implementation is:\t",shortest, len(shortest), "bp")
