__author__  = "Filippo Bergeretti"
__date__    = "10 December 2019"

import os
os.environ["PATH"] += os.pathsep + 'D:/Programmi/Graphviz/bin'
from graphviz import Digraph
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


def DeBruijn(st, k):
    """ Return a list holding, for each k-mer, its left
        k-1-mer and its right k-1-mer in a pair """
    edges = []
    nodes = set()
    for i in range(len(st) - k + 1):
        edges.append((st[i:i+k-1], st[i+1:i+k]))
        nodes.add(st[i:i+k-1])
        nodes.add(st[i+1:i+k])
    return nodes, edges


# ========================================= MAIN FUNCTION ==============================================================

def Visualizer(st, k):
    """ Visualize a directed multigraph using graphviz """
    nodes, edges = DeBruijn(st, k)
    lbl = (str(st)+ '\n'+ str(k)+ ' k-mers')
    dot = Digraph(name='test', encoding='utf-8')
    for node in nodes:
        dot.node(node)

    for src, dst in edges:
        dot.edge(src, dst)
    dot.attr(label=lbl)
    dot.attr(fontsize='10')
    dot.render('test-output/round-table.gv', view=True)

    dot_str = 'digraph "DeBruijn graph" {\n'
    for node in nodes:
        dot_str += '  %s [label="%s"] ;\n' % (node, node)
    for src, dst in edges:
        dot_str += '  %s -> %s ;\n' % (src, dst)
    return dot_str + '}\n'


dna_len     = 12
n           = 1
kmer_length = 3
seq         = randomString(dna_len, n)
sequence    = seq[n-1]
print(dna_len, "bp sequence:")
print(' '.join(seq))
result = Visualizer(sequence, kmer_length)
print(result)
