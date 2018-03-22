import mummerTools as mt
from subprocess import call
import nuclseqTools

with open ("../canu1.7/tig938.fasta") as s1:
    for line in s1:
        if line.startswith(">"):
            continue
        else:
            seq1 = line
with open ("../canu1.7/tig996.fasta") as s2:
    for line in s2:
        if line.startswith(">"):
            continue
        else:
            seq2 = line
#seq1 = nuclseqTools.reverse_complement(seq1)

delta, reflen, querylen = mt.align(seq1,seq2)

alignment,dir = mt.findDirectionRef(delta, reflen, querylen)
print(dir)
'''
coords = alignment[0].split(" ")

ref_coords = (int(coords[0])-1, int(coords[1]))

query_coords = (0, int(coords[3]))

seq1 = seq1[ref_coords[0]:ref_coords[1]]
seq2 = seq2[query_coords[0]:query_coords[1]]

merged_seq = nuclseqTools.createConsensus(alignment[1:],seq1,seq2)
print(merged_seq)
'''
