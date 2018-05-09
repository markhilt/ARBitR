#!/usr/bin/env python
# -*- coding: utf-8 -*-

import mappy as mp
import argparse
import re

parser = argparse.ArgumentParser(description="Extract split mapping hits.")
parser.add_argument("fasta", help="Reference fasta file. Required.", type = str)
parser.add_argument("fastq", help="Input fastq file. Required.", type = str)
args = parser.parse_args()

"""
def breakScaffold(aligned_read_seq, aligned_read_coords):
    '''
    Given an aligned read sequence and
    Might be good to make a consensus of all reads aligning here men orkar fan inte nu.
    '''
    read_seq = aligned_read_seq[aligned_read_coords[0]:aligned_read_coords[1]]
    newseq = "".join([ scaffold_seq[:scaffold_coords[0]], read_seq, scaffold_seq[scaffold_coords[1]:] ])
    return newseq
"""



def chopCigar(cigar):
    '''
    Given a cigar string, chops it up into a list.
    '''
    cig_list = []
    pos = 0
    for idx, c in enumerate(cigar):
        if c.isalpha():
            cig_list.append(cigar[pos:idx+1])
            pos  = idx + 1
    return cig_list

def findReadCoords(cig_list,ref_aln_start,ref_scf_coords,read_aln_start):
    '''
    Given a cigar string, reference scaffold coordinates and read sequence, calculates the
    coordinates in the read sequence that correspond to the scaffold N's.
    '''
    pos = ref_aln_start # track position in reference
    readpos = read_aln_start # track position in read
    readcoords = [0,0]
    for cg in cig_list:

        if pos < ref_scf_coords[0]:
            # If no indels, increment both positions by number of steps
            if cg[-1] == "M":
                pos += int(cg[:-1])
                readpos += int(cg[:-1])
            # If insertion in read, increase readpos by number of steps
            elif cg[-1] == "I":
                readpos += int(cg[:-1])
            # If deletion in read, increase reference pos by number of steps
            elif cg[-1] == "D":
                pos += int(cg[:-1])

        # Check if first scaffold coordinate was passed
        elif pos >= ref_scf_coords[0]:
            if readcoords[0] == 0:
                # If scaffold starting coordinate was passed and
                # readcoords hasn't been uptaed, check by how many steps and
                # update
                overstep_start = pos - ref_scf_coords[0]
                readcoords[0] = readpos - overstep_start

            # Also check if scaffold ending coordinate was passed
            if pos >= ref_scf_coords[1]:
                # If so, check by how much. Rest will align, because of N's
                overstep_end = pos - ref_scf_coords[1]
                readcoords[1] = readpos - overstep_end
                break

            # If it wasn't passed, keep going through cigar string until it is passed
            else:
                # If no indels, increment both positions by number of steps
                if cg[-1] == "M":
                    pos += int(cg[:-1])
                    readpos += int(cg[:-1])
                # If insertion in read, increase readpos by number of steps
                elif cg[-1] == "I":
                    readpos += int(cg[:-1])
                # If deletion in read, increase reference pos by number of steps
                elif cg[-1] == "D":
                    pos += int(cg[:-1])

    return readcoords


def createFixedSequence(old_seq, fills):
    '''
    Given a list of fill sequences and coords where they should go, fills the old sequence
    with the fills and returns the fixed sequence. Fills is a list of (start, stop, seq) tuples,
    where start and stop are coordinates and seq is the sequence that goes between them
    '''
    fills = sorted(fills)
    # First break old_seq into segments
    segments = []
    prev = 0
    for f in fills:
        segments.append( old_seq[prev:f[0]] )
        segments.append(f[2])
        prev = f[1]

        # If at last entry, also append remaining old sequence
        if f == fills[-1]:
            segments.append( old_seq[f[1]:] )

    newseq = "".join(segments)
    return newseq

def find_scaffold(seq):
    '''
    This shit will not be needed later
    '''
    Ns = [n.span() for n in re.finditer(r"((N)\2{9,})", seq)] # Convert span of hits to list
    #Ns = [n.span() for n in re.finditer(r"[N]+", seq)] # Convert span of hits to list
    return Ns

def getnames():
    refnames = []
    with open(args.fasta, "r") as fasta:
        for line in fasta:
            line = line.strip()
            if line.startswith(">"):
                refnames.append(line[1:])
    return refnames

def main():
    a = mp.Aligner(args.fasta, preset="map-pb")

    #tig164 = a.seq("tig00000164", 247307, 253307)

    #if not a: raise Exception("ERROR: failed to load/build index file '{}'".format(args.fasta))
    refnames = getnames()
    scaffolds = {} # {tig: <iterator of RE hits>}
    """
    testref = "ATGAGAATGCGAGATTTNNNNNNNNNNNATGCCTAGAGCTGAGGATGAGCT"
    testseq = "AGAGCGTAGACTTAGTGATGAtTGCACGTGAATGCGCAGGCTGAGCGATGAGCT"
    facit = [14,28]

    #testcigar = "3M1D6M1I4M4I14M1D21M"
    testcigar = "3M2D3M1I3M1D1I2M7I15M1I1M1D2M1D6M1I8M"
    scf_coords = (17,28)
    read_start = 0
    ref_start = 3
    testcoords = findReadCoords(testcigar,ref_start,scf_coords,read_start)
    fill = testseq[testcoords[0]:testcoords[1]]
    f = [(scf_coords[0],scf_coords[1], fill)]
    print(createFixedSequence(testref,f))
    """
    for tig in refnames:
        seq = str( a.seq(tig) )
        scaffolds[tig] = find_scaffold(seq)
    '''
    for k,v in scaffolds.items():
        for s in v:
            print(k, s)
            #    tig = "linked_contig_25"
            #    start = 350708
            #    stop = 350718
    '''
    #scaffolds["linked_contig_25_pilon"] = [(44267, 44405)]

    # Map reads, and shout if mapping is found across scaffold coords
    mapping = {} # To fill with reads mapping across scaffolded regions
    fillings = {} # To fill with sequences
    for name, seq, qual in mp.fastx_read(args.fastq):

        # Skip reads shorter than 5kb; won't be useful
        if len(seq) > 19000:
            for h in a.map(seq):

                # Look for alignments across specific region
                if h.ctg in scaffolds.keys():
                    # if mapped read is on a contig in our dict, check coords
                    for s in scaffolds[h.ctg]:
                        #start = s.span()[0]
                        #stop = s.span()[1]
                        start = int(s[0])
                        stop = int(s[1])

                        # For now, just take first read that aligns across the scaffold
                        if int(h.r_st) < start and int(h.r_en) > stop:
                            cig_list = chopCigar(h.cigar_str)

                            # If on the - strand, revcomp and reverse cigar
                            if h.strand == -1:
                                seq = mp.revcomp(seq)
                                cig_list = cig_list[::-1]

                            readcoords = findReadCoords(cig_list,h.r_st,s, h.q_st)
                            #print(readcoords)
                            fill = seq[readcoords[0]:readcoords[1]]

                            if h.ctg not in fillings:
                                fillings[h.ctg] = [(start, stop, fill)]
                            else:
                                ls = fillings[h.ctg]
                                ls.append( (start, stop, fill) )
                                fillings[h.ctg] = ls

                            # Remove from scaffolds dict
                            ls = scaffolds[h.ctg]
                            ls.remove(s)
                            scaffolds[h.ctg] = ls
                            if scaffolds[h.ctg] == []:
                                del scaffolds[h.ctg]

    print("Done mapping. Filling scaffolds.")

    # Create new sequences and write report
    with open("test.fasta", "w") as fastaout:
        with open("report_test.txt", "w") as reportout:
            for r,f in fillings.items():
                fasta = "".join( [">", r, "\n", createFixedSequence(a.seq(r),f)] )
                fastaout.write(fasta+"\n")
                for scf in f:
                    report_line = "Fixed scaffold {}: {}-{}\n".format(r,str(scf[0]),str(scf[1]))
                    reportout.write(report_line)

            for scf_tig,scf_coords in scaffolds.items():
                for entry in scf_coords:
                    report_line = "Remaining scaffold {}: {}-{}\n".format(scf_tig,str(entry[0]),str(entry[1]))
                    reportout.write(report_line)





if __name__ == "__main__":
    main()

# python miniasm_bridge.py ../assemblies/comb_asm1-7_merged_scaffolded_pil.fasta ../reads/canu-pb-nano1.7.trimmedReads.fasta
