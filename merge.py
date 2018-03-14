#!/home/markus/miniconda3/bin/python3.6
# -*- coding: utf-8 -*-

import argparse
import numpy as np
import time
import pysam
import timeit
import os
from subprocess import call
from scipy.stats import t
import formatting_tools as ft
import calcESD as esd
import graph_building as gb

parser = argparse.ArgumentParser(description="Reads a bam file, creates links between contigs based on linked read information, and outputs a .gfa.")
parser.add_argument("input_bam", help="Input bam file. Required.", type = str)
parser.add_argument("-i", "--input_fasta", help="Input fasta file for contig merging. Optional. If not specified, will only output linkage .gfa and .gv for manual merging.", type = str)
parser.add_argument("-s","--region_size", help="Size of region of contig start and end to collect barcodes from. [20000]", default = 20000, type = int)
parser.add_argument("-n","--barcode_number", help="Minimum number of shared barcodes to create link. [1]", default = 1, type = int)
parser.add_argument("-f","--barcode_fraction", help="Minimum fraction of shared barcodes to create link. [0.01]", default = 0.01, type = float)
parser.add_argument("-q","--mapq", help="Mapping quality cutoff value. [60]", default = 60, type = int)
parser.add_argument("-o","--output", help="Prefix for output files.", type = str)
args = parser.parse_args()

#### Functions

# Given reads from samfile, return GEM barcodes aligned to this region in a set
def getGEMs(reads):
    BC_list = []
    for read in reads:
        if read.has_tag('BX') == True and read.mapping_quality >= args.mapq:
            BC = read.get_tag("BX")
            BC_list.append(BC)
    BC_list = set(BC_list)
    return BC_list

# To split up the contig into windows
def getWindows(contig):
    length = contig_dict[contig]
    start = (0,args.region_size)

    # If contig is shorter than region_size, just take all of it to "end" as well
    if (length - args.region_size) < 0:
        endstart = 0
    else:
        endstart = length - args.region_size
    end = (endstart,length)
    return [start,end]

# To count the number of reads aligned to a region
def countReads(contig,coords_to_check):
    cov_arr = samfile.count_coverage(contig, coords_to_check[0],coords_to_check[1])
    cov = sum([sum(cov_arr[x]) for x in range(0,4,1)]) / 50
    return cov


def trimFasta(contig_side_to_trim):
    '''
    Given a fasta entry with side to trim, returns new start and end coordinates
    '''
    # Input format: "tig-side"
    tig = contig_side_to_trim[:-1]
    side = contig_side_to_trim[-1]
    coords_at_a_time = 50
    mincov = 50

    if tig in fastafile.references:
        cov = 0
        coords_to_trim = 0

        if side == "s":
            coords_to_check = [0,coords_at_a_time]
            while cov < mincov:
                cov = countReads(tig,coords_to_check)

                if cov <= mincov:
                    coords_to_trim += coords_at_a_time # Update with new end coordinate
                    coords_to_check = [x+coords_at_a_time for x in coords_to_check] # And move to next 10 coordinates

        elif side == "e":
            stop = trimmed_fasta_coords[tig][1]
            coords_to_check = [stop-coords_at_a_time,stop]
            while cov <= mincov:
                cov = countReads(tig,coords_to_check)

                if cov <= mincov:
                    coords_to_trim = coords_to_trim - coords_at_a_time    # Update with new end coordinate
                    coords_to_check = [x-coords_at_a_time for x in coords_to_check] # And move to next 10 coordinates

    return (coords_to_trim)

def reverse_complement(nuclstring):
    rev_comped = ""
    for l in reversed(nuclstring):
        if l == "A" or l == "a":
            rev_comped += "T"
        elif l == "T" or l == "t":
            rev_comped += "A"
        elif l == "C" or l == "c":
            rev_comped += "G"
        elif l == "G" or l == "g":
            rev_comped += "C"
        elif l == "N" or l == "n":
            rev_comped += "N"
    return rev_comped

def createConsensus(delta,string1,string2):
    '''
    Given two overlapping nucleotide strings (excluding overhangs) and alignment information,
    returns the merged sequence
    delta = [int,int,int, ...]
    strings = "ATCG"
    '''
    new_string1 = string1
    new_string2 = string2
    start = 0
    for i in delta:
        if i > 0:
            start += i
            new_string2 = new_string2[:start-1] + "." + new_string2[start-1:]

        else:
            i = -i
            start += i
            new_string1 = new_string1[:start-1] + "." + new_string1[start-1:]

    # Write the new string. Insertions are favoured over deletions.
    # N's are put at mismatches
    cons = ""
    for i in range(0,len(new_string1),1):
        if new_string1[i] == new_string2[i]:
            cons += new_string1[i]
        elif new_string1[i] == ".":
            cons += new_string2[i]
        elif new_string2[i] == ".":
            cons += new_string1[i]
        else:
            cons += "N"

    return cons

samfile = pysam.AlignmentFile(args.input_bam, "rb")
contig_list = samfile.references
contig_dict = {}
i = 0

# Create output file name if nothing specified
if args.output:
    outfilename = args.output
elif "/" in args.input_bam:
    outfilename = args.input_bam.split("/")[-1].split(".bam")[0]
else:
    outfilename = args.input_bam.split(".bam")[0]

# Create dict of contigs and their lengths, save also to a gfa file
gfa_header = []
while i < len(contig_list):
    contig_dict[contig_list[i]] = samfile.lengths[i]
    gfa_header.append("S\t{0}\t*\tLN:i:{1}".format(contig_list[i],samfile.lengths[i]))
    i += 1

# Main.
localtime =  localtime = time.asctime( time.localtime(time.time()) )
starttime = timeit.default_timer()
print("\nStarting merge pipeline. {0}\n".format(localtime))
GEMlist = {}
GEMcomparison = {}
bclib_list = []

# First step is to collect all barcodes (passing -q cutoff) that are aligned to each contigs first
# and last regions (-l)
print("Starting barcode collection. Found {0} contigs.\n".format(len(contig_list)))
donecontigs = 0
for contig in contig_list:

    # Report progress every 20 windows
    if donecontigs in range(0,100000000,20):
        print("Completed: {0}% ({1} out of {2})".format( str(round(donecontigs / len(contig_list) * 100, 2) ), donecontigs, len(contig_list)) )

    # Create windows from first and last X kb from each contig
    windowlist = getWindows(contig)

    # Collect barcodes from every window
    for window in windowlist:
        reads = samfile.fetch(contig, window[0], window[1])
        GEMs = getGEMs(reads)

        # For the merge pipeline, windows equal start and end 20kb of contig
        if window == windowlist[0]:
            region = contig + "s"
        elif window == windowlist[1]:
            region = contig + "e"

        if GEMs:
            GEMlist[region] = GEMs

    donecontigs += 1

samfile.close()

print("\nBarcode listing complete. Starting pairwise comparisons...\n")

# Second step is to compare the barcodes in every region to all other regions
donewindows = 0
prev_contig = ""
GEMlistlen = len(GEMlist)
for region in GEMlist:
    contig = region[:-1]

    # Report progress every 20 windows
    if donewindows in range(0,100000000,20):
        print("Completed: {0}% ({1} out of {2})".format( str(round(donewindows / GEMlistlen * 100, 2) ),donewindows, GEMlistlen))

    lib1 = GEMlist[region] # Reference region's barcodes
    nested_dict = {} # Each entry in the original dictionary is another dictionary

    for region2 in GEMlist:
        lib2 = GEMlist[region2] # Comparison region's barcodes

        shared = set(lib1).intersection(lib2) # Find shared ones
        lib1length = int(len(lib1)) # Get number of barcodes in this region
        lib2length = int(len(lib2))
        sharedlen = int(len(shared))
        totallength = lib1length + lib2length - sharedlen # The total number of barcodes to consider is this

        # To find the fraction of shared barcodes, NEEDS SOME FORM OF NORMALIZATION?
        if totallength != 0:
            fraction = sharedlen / totallength
        else:
            fraction = 0

        # Fill the nested dictionary
        nested_dict[region2] = (fraction,sharedlen)

    prev_contig = contig
    donewindows += 1
    GEMcomparison[region] = nested_dict

# Save gfa header
gfa_header = "H\tVN:Z:avid/links\n" + "\n".join(gfa_header)

# Write the output comparison matrix
print("\nWriting comparison matrix to {0}.txt\n".format(outfilename) )
with open(outfilename + ".txt", "w", encoding = "utf-8") as outfile:
    outfile.write(ft.formatTable(GEMcomparison))

print("chromQC started at {0}.\n".format(localtime))
completed_windows = 0

# Third step is to create edges based on outlying fraction values.
# If a certain window has an outlier to another window, an edge is
# created
windows = sorted(GEMcomparison.keys())
totwin = len(windows)
print("Number of windows: "+str(totwin))
presentwindows = {}
edges = {}
for region in GEMcomparison:
    contig = region[:-1]
    window = region[-1]

    fractions = {}
    counts = {}

    # Report progress every 100 windows
    if completed_windows in range(0,10000000,100):
        progress = round( (completed_windows/totwin)*100, 2)
        print("Progress: {0} % of windows completed ({1} out of {2})".format(str(progress),completed_windows,totwin))

    # Calculate outliers from the comparisons of window k to all other windows
    outliers = esd.getOutliers(GEMcomparison[region])
    outliers_short = {}

    # Remove fractions less than -f, including lower outliers
    # and outliers with too few shared barcodes (-n)
    for k,v in outliers.items():
        if v[0] > args.barcode_fraction and v[1] > args.barcode_number:
            outliers_short[k] = v

    # Sort outliers_short dict by fraction value into the list sorted_outliers
    sorted_outliers = [(k, outliers_short[k]) for k in sorted(outliers_short, key=outliers_short.get, reverse=True)]

    # Skip to next line if there are no good links
    if len(sorted_outliers) == 0:
        completed_windows += 1
        continue

    # Several hits will be given edges in the gfa, if present
    for l in sorted_outliers:
        link = l[0]
        ltig = link[:-1]
        lwin = link[-1]

        # Skip hits to self
        if ltig == contig:
            continue

        # Check if current edge was added in the other direction,
        # don't add in again
    #    if region in edges[link]:
    #        continue

        # If not, and the region is not yet in the dict,
        # add the region as key and link as value in a list
        # If the region has been added to the dict previously,
        # append new link to the list
    #    else:
        if region in edges:
            linkedlist = edges[region]

            linkedlist.append(link)
            edges[region] = linkedlist

        else:
            edges[region] = [link]

    completed_windows += 1

# Write gfa
print("\nWriting gfa to {0}.gfa.\n".format(outfilename))

with open(outfilename + ".gfa", "w", encoding = "utf-8") as gfa:
    gfa.write(gfa_header+"\n")
    gfalist = ft.formatGFA(edges)
    for i in gfalist:
        gfa.write(i)
'''
# Write tsv
with open (outfilename+".tsv", "w", encoding = "utf-8") as tsvout:
    tsvlist = ft.formatTSV(edges)
    for i in tsvlist:
        tsvout.write(i)
'''

if not args.input_fasta:
    print("No fasta file specified for merging. Pipeline finished.")

# Third step is to align ends of contigs, following unique edges,
# merge them into a "consensus" sequence, and write a new fasta with
# merged contigs.
else:
    fastafile = pysam.FastaFile(args.input_fasta)
    samfile = pysam.AlignmentFile(args.input_bam, "rb")
    print("Found fasta file for merging: {0}".format(args.input_fasta))

    # Testing reformatting of edges from dict to list of tuples
    edges_list = []
    for i in edges:
        for a in edges[i]:
            # If edge was not already added in the opposite direction, add it
            if ( (a,i) ) not in edges_list:
                edges_list.append( (i, a) )

    # Reformat unique edges to a dict {"linked_contig_X":[ ("old_contig_1f/r"), ("old_contig_2f/r"), ... ]}
    global linked_contigs
    global done_edges
    linked_contigs = {}
    done_edges = []
    counter = 1

    # Time to start linking contigs together based on the linkage graph
    # First find nodes that will start or end linked contigs
    linked_contig_ends = gb.findStartingNode(edges)

    # Then go through these and build one contig between two entries
    for c in linked_contig_ends:
        if c[:-1] not in done_edges:
            linked_tig_ID = "linked_contig_"+str(counter)
            linked_contig = gb.create_linked_contig(c)
            done_edges.append(linked_contig[0][:-1])
            done_edges.append(linked_contig[-1][:-1])
            linked_contigs[linked_tig_ID] = linked_contig
            counter += 1

    for i in linked_contigs:
        print(i,linked_contigs[i])
        '''
    # Next trim fasta sequences to be merged in low coverage regions

    # trimmed_fasta_coords is a dict with coords to keep from original fasta
    # Format: {contig: [start_coord, end_coord]}
    # Start by filling with old coords, which will then be changed
    trimmed_fasta_coords = {}
    for i in range(0,len(fastafile.references),1):
        trimmed_fasta_coords[fastafile.references[i]] = [0, fastafile.lengths[i]]

    # Then find new coordinates for all sequences to merge
    print("Trimming contig edges.: {0}".format(args.input_fasta))
    for edge in linked_contigs.values():
        for i in edge:
            tig = i[:-1]
            # Trim only sides where links are formed
            # If first connected node, trim last coordinate
            if i == edge[0]:
                trimmed_fasta_coords[tig] = [trimmed_fasta_coords[tig][0], trimmed_fasta_coords[tig][1] + trimFasta(tig+"e")]
            # If last connected node, trim first coordinate
            elif i == edge[-1]:
                trimmed_fasta_coords[tig] = [trimmed_fasta_coords[tig][0] + trimFasta(tig+"s"), trimmed_fasta_coords[tig][1]]
            # If node is connected at both sides, trim both
            else:
                trimmed_fasta_coords[tig] = [trimmed_fasta_coords[tig][0] + trimFasta(tig+"s"), trimmed_fasta_coords[tig][1] + trimFasta(tig+"e")]

    # Then walk through the linked_contigs dict, gradually building each contig from left to right.
    # Sequences to be merged are saved in merged_seqs in the format
    # {new_linked_contig_ID: "ATCG..." }
    merged_seqs = {}
    unaligned_linked_contigs = {} # Collect linked contigs where no alignment was found
    contigs_to_exclude = [] # Collect contigs that were merged into new linked contigs
    print("Building linked contigs.")
    for linked_tig in linked_contigs:
        # Go through every linked contig in the list up until the next to last;
        # at this step the last contig will be merged into the sequence
        tig_counter = 0 # To keep track of where in the linked contig list we are
        tig_started = False # To check if merged contig was initiated

        # List of nucleotide strings to merge.
        # Format: [contig1_left_overhang, contig1-contig2_aligned_consensus, contig2_right_overhang ...]
        newseq = []
        for t in linked_contigs[linked_tig]:
            # If starting_tig was not created before, i.e. if starting from the first linked contig,
            # assign these variables
            if tig_started == False:
                tig_started = True
                starting_tig = linked_contigs[linked_tig][0]
                refname, refdirection = starting_tig[:-1], starting_tig[-1]
                ref_fasta = fastafile.fetch(reference=refname, start=trimmed_fasta_coords[refname][0], end=trimmed_fasta_coords[refname][1])
                # Reverse complement if necessary
                if refdirection == "r":
                    ref_fasta = reverse_complement(ref_fasta)

            # If tig_started == True, i.e. this is the second or more run through the loop,
            # there is already a starting point, i.e. newseq.
            # Keep aligning and merging linked contigs to this one.
            else:
                ref_fasta = newref # last entry of previous iteration
                refdirection = ""

            reflen = len(ref_fasta) # Collect length of reference sequence

            # If at the last contig in the list, just append it to newseq and then break
            if tig_counter == len(linked_contigs[linked_tig]) - 1:
                newseq.append(ref_fasta)
                break

            next_tig = linked_contigs[linked_tig][tig_counter+1]
            queryname, querydirection = next_tig[:-1], next_tig[-1]
            query_fasta = fastafile.fetch(reference=queryname, start=trimmed_fasta_coords[queryname][0], end=trimmed_fasta_coords[queryname][1])

            # Reverse complement if necessary
            if querydirection == "r":
                query_fasta = reverse_complement(query_fasta)

            # Write trimmed contigs to temp files for mummer alignment in working directory
            with open("ref.fa","w",encoding = "utf-8") as fastaout:
                fastaout.write(">"+refname+refdirection+"\n")
                fastaout.write(ref_fasta)
            with open("query.fa","w",encoding = "utf-8") as fastaout:
                fastaout.write(">"+queryname+querydirection+"\n")
                fastaout.write(query_fasta)

            # Call mummer. Select overlap based on nucmer
            call(["nucmer", "ref.fa", "query.fa"])

            # Select alignment with correct start and end coordinates for both contigs
            with open("out.delta","r",encoding = "utf-8") as delta:
                passed = False # To check if the expected alignment was found or not
                d = [] # Reconstruct delta in this list

                for line in delta:
                    line = line.strip()
                    fields = line.split(" ")

                    # Filter out correct alignment in the delta
                    # New alignments have 7 fields
                    if len(fields) == 7:
                        alignment_ref = (int(fields[0]),int(fields[1]))
                        alignment_query = (int(fields[2]),int(fields[3]))

                        # Because we revcomped earlier, all alignments will start somewhere in ref and end in last coord of ref,
                        # and start at beginning of query and end somewhere in query
                        if alignment_ref[1] == reflen and alignment_query[0] == 1:
                            # Correct alignment found. Start going through it.
                            passed = True

                    # Reconstruct delta
                    if len(fields) == 1 and passed == True:
                        deltaint = int(line)
                        if deltaint != 0:
                            d.append(deltaint)

                        else:
                            break

                # If the correct alignment was found, append sequences to newseq
                if passed == True:
                    # If this is the first time through the loop,
                    # add beginning of ref sequence to the new contig (newseq),
                    # up until overlap begins. Else it has been added before as
                    # query.
                    newseq.append(ref_fasta[:alignment_ref[0]])
                    # Then collect the subsequences that align
                    alrefseq = ref_fasta[alignment_ref[0]:alignment_ref[1]]
                    alqueryseq = query_fasta[alignment_query[0]:alignment_query[1]]

                    # Add "consensus" of aligned sequence to new contig
                    newseq.append(createConsensus(d,alrefseq,alqueryseq))
                    # And then assign the remainder of query contig as new reference
                    newref = query_fasta[alignment_query[1]:]

                # If the correct alignment was not found, instead scaffold by inserting 10*N:
                else:
                    newseq.append(ref_fasta)
                    newseq.append("NNNNNNNNNN")
                    newref = query_fasta

                # Exclude the tigs that were just added to newseq
                if refname not in contigs_to_exclude:
                    contigs_to_exclude.append(refname)
                contigs_to_exclude.append(queryname)

                # Go to next contig in the list
                tig_counter += 1

        # After iterating through every element in linked_contigs[linked_tig][:-1],
        # merge newseq and add to the merged_seqs dict
        merged_seqs[linked_tig] = "".join(newseq)

    print("Writing merged fasta to {0}.fasta".format(outfilename))
    with open(outfilename+".fasta","w",encoding = "utf-8") as fastaout:
        for i in merged_seqs:
            fastaout.write(">"+i+"\n")
            fastaout.write(merged_seqs[i]+"\n")

        for i in fastafile.references:
            if i not in contigs_to_exclude:
                fastaout.write(">"+i+"\n")
                fastaout.write(fastafile.fetch(reference=i)+"\n")

    print("Writing report to report.txt")
    with open("report.txt","w",encoding = "utf-8") as report:
        for i in linked_contigs:
            report.write(i)
            for j in linked_contigs[i]:
                report.write("\t"+j)
            report.write("\n")

print("Done.")

'''
