import mummerTools as mt
from subprocess import call
import nuclseqTools as nt

with open ("../canu1.7/ref.fasta") as s1:
    for line in s1:
        if line.startswith(">"):
            continue
        else:
            seq1 = line
with open ("../canu1.7/query.fasta") as s2:
    for line in s2:
        if line.startswith(">"):
            continue
        else:
            seq2 = line
#seq1 = nuclseqTools.reverse_complement(seq1)

delta, reflen, querylen = mt.align(seq1,seq2)

#alignment,dir = mt.findDirectionRef(delta, reflen, querylen)
#print(dir)

alignment = mt.findAlignment(delta, reflen)
newseq = []

if alignment != None:

    # Collect the coordinates
    coords = alignment[0].split(" ")
    ref_coords = (int(coords[0])-1, int(coords[1]))
    query_coords = (0, int(coords[3]))

    # Add beginning of ref sequence to the new contig (newseq),
    # up until overlap begins.
    newseq.append(seq1[:ref_coords[0]])
    # Then collect the subsequences that align
    alrefseq = seq1[ref_coords[0]:ref_coords[1]]
    alqueryseq = seq2[query_coords[0]:query_coords[1]]

    # Add "consensus" of aligned sequence to new contig
    alignment_ints = alignment[1:] # First entry in alignment is the header containing coordinates
    newseq.append(nt.createConsensus(alignment_ints,alrefseq,alqueryseq))
    newseq.append(seq2[query_coords[1]:])
    newseq = "".join(newseq)
    with open("out.fasta", "w") as out:
        out.write(">seq\n"+newseq)

    '''
coords = alignment[0].split(" ")

ref_coords = (int(coords[0])-1, int(coords[1]))

query_coords = (0, int(coords[3]))

seq1 = seq1[ref_coords[0]:ref_coords[1]]
seq2 = seq2[query_coords[0]:query_coords[1]]

merged_seq = nuclseqTools.createConsensus(alignment[1:],seq1,seq2)
print(merged_seq)
'''
