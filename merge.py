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

# Included modules
import formatting_tools as ft
import calcESD as esd
import graph_building as gb
import nuclseqTools as nt
import mummerTools as mt

parser = argparse.ArgumentParser(description="Reads a bam file, creates links between contigs based on linked read information, and outputs a .gfa.")
parser.add_argument("input_bam", help="Input bam file. Required.", type = str)
parser.add_argument("-i", "--input_fasta", help="Input fasta file for contig merging. Optional. If not specified, will only output linkage .gfa and .gv for manual merging.", type = str)
parser.add_argument("-s","--region_size", help="Size of region of contig start and end to collect barcodes from. [20000]", default = 20000, type = int)
parser.add_argument("-n","--barcode_number", help="Minimum number of shared barcodes to create link. [1]", default = 1, type = int)
parser.add_argument("-f","--barcode_fraction", help="Minimum fraction of shared barcodes to create link. [0.01]", default = 0.01, type = float)
parser.add_argument("-q","--mapq", help="Mapping quality cutoff value. [60]", default = 60, type = int)
parser.add_argument("-o","--output", help="Prefix for output files.", type = str)
args = parser.parse_args()

def getGEMs(reads):
    '''
    Given reads from samfile, return GEM barcodes aligned to this region in a set
    '''
    BC_list = []
    for read in reads:
        if read.has_tag('BX') == True and read.mapping_quality >= args.mapq:
            BC = read.get_tag("BX")
            BC_list.append(BC)
    BC_list = set(BC_list)
    return BC_list

def getWindows(contig):
    '''
    To split up the contig into windows
    '''
    length = contig_dict[contig]
    start = (0,args.region_size)

    # If contig is shorter than region_size, just take all of it to "end" as well
    if (length - args.region_size) < 0:
        endstart = 0
    else:
        endstart = length - args.region_size
    end = (endstart,length)
    return [start,end]

def countReads(contig,coords_to_check):
    '''
    To count the number of reads aligned to a region
    '''
    cov_arr = samfile.count_coverage(contig, coords_to_check[0],coords_to_check[1])
    cov = sum([sum(cov_arr[x]) for x in range(0,4,1)]) / 50
    return cov

def trimFasta(contig_side_to_trim):
    '''
    Given a fasta entry with side to trim, returns new start and end coordinates
    Input format: "tigs" or "tige"
    '''

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

    # Reformat edges from dict to list of tuples
    gb.edges_as_list(edges)

    # Build linked contigs in a dict {"linked_contig_X":[ ("old_contig_1f/r"), ("old_contig_2f/r"), ... ]}
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
            # If direction of node is unknown (last character in string == "u"),
            # trim both sides
            # If first connected node, trim last coordinate
            if i == edge[0] and i[-1] != "u":
                trimmed_fasta_coords[tig] = [trimmed_fasta_coords[tig][0], trimmed_fasta_coords[tig][1] + trimFasta(tig+"e")]
            # If last connected node, trim first coordinate
            elif i == edge[-1] and i[-1] != "u":
                trimmed_fasta_coords[tig] = [trimmed_fasta_coords[tig][0] + trimFasta(tig+"s"), trimmed_fasta_coords[tig][1]]
            # If node is connected at both sides, trim both
            else:
                trimmed_fasta_coords[tig] = [trimmed_fasta_coords[tig][0] + trimFasta(tig+"s"), trimmed_fasta_coords[tig][1] + trimFasta(tig+"e")]

    # Before merging can begin, the direction of any contigs of unknown
    # orientation needs to be determined. If the direction cannot be determined,
    # the linked contig is broken.
    print("\nFinding orientation of unknown contigs...")
    linked_contigs_oriented = {} # This dict will be filled with the oriented contigs
    for linked_tig,tig_list in linked_contigs.items():
        tig_counter = 0 # To keep track of position in tig_list
        c = 0 # Count new linked contigs formed by breaking earlier ones
        new_tig_list = [] # This list will be filled with the oriented contigs
        new_linked_tig_ID = ".".join([linked_tig,str(c)])
        for t in tig_list:
            tig = t[:-1]
            orientation = t[-1]
            # If orientation is not "u", add entries as before
            if orientation != "u":
                new_tig_list.append(t)
            else:
                # If orientation is unknown, use mummerTools module to find it
                # The contig to compare to depends on where in the linked contig we are
                # If at the first contig to be linked, compare to the second
                if t == tig_list[0]:
                    ref_fasta = fastafile.fetch(reference=tig, start=trimmed_fasta_coords[tig][0], end=trimmed_fasta_coords[tig][1])
                    ctig = tig_list[1][:-1]
                    corientation = tig_list[1][-1]
                    query_fasta = fastafile.fetch(reference=ctig, start=trimmed_fasta_coords[ctig][0], end=trimmed_fasta_coords[ctig][1])
                    if corientation == "r":
                        query_fasta = nt.reverse_complement(query_fasta)
                    delta, reflen, querylen = mt.align(ref_fasta,query_fasta)
                    alignment,dir = mt.findDirectionRef(delta, reflen, querylen)
                    if dir != None:
                        tig_ori = tig+dir
                        new_tig_list.append(tig_ori)

                    # If direction still cannot be determined, continue iterating through tig_list
                    else:
                        continue

                # If anywhere else, compare to the surrounding contigs
                else:
                    # First the previous contig
                    ctig = tig_list[tig_counter-1][:-1]
                    corientation = tig_list[tig_counter-1][-1]
                    ref_fasta = fastafile.fetch(reference=ctig, start=trimmed_fasta_coords[ctig][0], end=trimmed_fasta_coords[ctig][1])
                    if corientation == "r":
                        ref_fasta = nt.reverse_complement(ref_fasta)
                    query_fasta = fastafile.fetch(reference=tig, start=trimmed_fasta_coords[tig][0], end=trimmed_fasta_coords[tig][1])
                    delta, reflen, querylen = mt.align(ref_fasta,query_fasta)
                    alignment,dir = mt.findDirectionQuery(delta, reflen, querylen)
                    if dir != None:
                        tig_ori = tig+dir
                        new_tig_list.append(tig_ori)

                    # If direction still cannot be determined, and at the last element in the tig_list,
                    # just continue without adding the last element
                    elif t == tig_list[-1]:
                        continue

                    # If direction cannot be determined from the previous contig, instead try the next one,
                    # if not already at the end of the tig_list
                    else:
                        ref_fasta = fastafile.fetch(reference=tig, start=trimmed_fasta_coords[tig][0], end=trimmed_fasta_coords[tig][1])
                        ctig = tig_list[tig_counter+1][:-1]
                        corientation = tig_list[tig_counter+1][-1]
                        query_fasta = fastafile.fetch(reference=ctig, start=trimmed_fasta_coords[ctig][0], end=trimmed_fasta_coords[ctig][1])
                        if corientation == "r":
                            query_fasta = nt.reverse_complement(query_fasta)
                        delta, reflen, querylen = mt.align(ref_fasta,query_fasta)
                        alignment,dir = mt.findDirectionRef(delta, reflen, querylen)
                        if dir != None:
                            tig_ori = tig+dir
                            new_tig_list.append(tig_ori)

                        # If dir == None, break the linked contig
                        # If there is anything in new_tig_list, add it to the new dict
                        else:
                            if len(new_tig_list) > 1:
                                linked_contigs_oriented[new_linked_tig_ID] = new_tig_list
                                c += 1
                                new_linked_tig_ID = ".".join([linked_tig,str(c)])
                                new_tig_list = []
            tig_counter += 1
        # After finishing iteration through each tig_list, and
        # there is anything in new_tig_list, add it to the new dict
        if len(new_tig_list) > 1:
            linked_contigs_oriented[new_linked_tig_ID] = new_tig_list
            c += 1

    print("Writing pre-oriented report to {0}.preor.report".format(outfilename))
    with open(outfilename+".preor.report","w",encoding = "utf-8") as report:
        for i in linked_contigs:
            report.write(i)
            for j in linked_contigs[i]:
                report.write("\t"+j)
            report.write("\n")

    print("Writing report to {0}.report".format(outfilename))
    with open(outfilename+".report","w",encoding = "utf-8") as report:
        for i in linked_contigs_oriented:
            report.write(i)
            for j in linked_contigs_oriented[i]:
                report.write("\t"+j)
            report.write("\n")

    # For the contigs that were now directed, find new trim coordinates if they are at the beginning
    # or the end of the linked contig
    # TODO

    # Then walk through the linked_contigs dict, gradually building each contig from left to right.
    # Sequences to be merged are saved in merged_seqs in the format
    # {new_linked_contig_ID: "ATCG..." }
    merged_seqs = {}
    contigs_to_exclude = [] # Collect contigs that were merged into new linked contigs
    print("Building linked contigs.")

    for linked_tig_ID,tig_list in linked_contigs_oriented.items():
        # Go through every linked contig in the dict to align and merge the sequences
        # contained in the values
        tig_counter = 0 # To keep track of where in the linked contig list we are

        # newseq is a List of nucleotide strings to merge.
        # This is the most important variable for the rest of this script
        # Format: [contig1_left_overhang, contig1-contig2_aligned_consensus, contig2_right_overhang ...]
        newseq = []
        starting_tig = None
        for next_tig in tig_list:
            # Exclude the tig that is being added, so it's not included in the final fasta
            contigs_to_exclude.append(next_tig[:-1])

            # If starting_tig was not created before, i.e. if starting from the first contig,
            # assign these variables
            if starting_tig == None:
                # Collect first contig
                starting_tig = next_tig
                # Collect name and direction of first contig
                refname, refdirection = starting_tig[:-1], starting_tig[-1]
                # And sequence
                ref_fasta = fastafile.fetch(reference=refname, start=trimmed_fasta_coords[refname][0], end=trimmed_fasta_coords[refname][1])
                # Reverse complement if necessary
                if refdirection == "r":
                    ref_fasta = nt.reverse_complement(ref_fasta)

            # If tig_started == True, i.e. this is the second or more run through the loop,
            # there is already a starting point, i.e. newref.
            # Keep aligning and merging linked contigs to this one.
            # Reference is then the last contig visited in the iteration and query the current
            else:
                # This path is taken when two sequences are to be aligned

                # Collect next tig to be merged into the linked contig
                queryname, querydirection = next_tig[:-1], next_tig[-1]
                query_fasta = fastafile.fetch(reference=queryname, start=trimmed_fasta_coords[queryname][0], end=trimmed_fasta_coords[queryname][1])

                # Reverse complement if necessary
                if querydirection == "r":
                    query_fasta = nt.reverse_complement(query_fasta)

                # Give trimmed contig sequences to mummer to align. Mummer is currently wrapped
                # in the mummerTools module
                delta, reflen, querylen = mt.align(ref_fasta,query_fasta)
                alignment = mt.findAlignment(delta, reflen)

                # If the expected alignment was found, append sequences to newseq
                if alignment != None:

                    # Collect the coordinates
                    coords = alignment[0].split(" ")
                    ref_coords = (int(coords[0])-1, int(coords[1]))
                    query_coords = (0, int(coords[3]))

                    # Add beginning of ref sequence to the new contig (newseq),
                    # up until overlap begins.
                    newseq.append(ref_fasta[:ref_coords[0]])
                    # Then collect the subsequences that align
                    alrefseq = ref_fasta[ref_coords[0]:ref_coords[1]]
                    alqueryseq = query_fasta[query_coords[0]:query_coords[1]]

                    # Add "consensus" of aligned sequence to new contig
                    alignment_ints = alignment[1:] # First entry in alignment is the header containing coordinates
                    newseq.append(nt.createConsensus(alignment_ints,alrefseq,alqueryseq))
                    # And then assign the remainder of query contig as new reference
                    ref_fasta = query_fasta[query_coords[1]:]

                # If the correct alignment was not found, instead scaffold by inserting 10*N:
                else:
                    newseq.append(ref_fasta)
                    newseq.append("NNNNNNNNNN")
                    ref_fasta = query_fasta

                if next_tig == tig_list[-1]:
                    # When at the last element in tig_list, there is nothing left to align.
                    # Merge newseq and add to the merged_seqs dict
                    newseq.append(ref_fasta)
                    merged_seqs[linked_tig_ID] = "".join(newseq)

            # Go to next contig in the list
            tig_counter += 1

    print("\nWriting merged fasta to {0}.fasta".format(outfilename))
    with open(outfilename+".fasta","w",encoding = "utf-8") as fastaout:
        for i in merged_seqs:
            fastaout.write(">"+i+"\n")
            fastaout.write(merged_seqs[i]+"\n")

        for i in fastafile.references:
            if i not in contigs_to_exclude:
                fastaout.write(">"+i+"\n")
                fastaout.write(fastafile.fetch(reference=i)+"\n")

print("Done.")
