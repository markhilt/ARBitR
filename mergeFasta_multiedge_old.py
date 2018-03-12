#!/usr/bin/python3.6
# -*- coding: utf-8 -*-

import argparse
import numpy as np
import time
import pysam
import timeit
import os
from subprocess import call
from scipy.stats import t

parser = argparse.ArgumentParser(description="Reads and trims a fasta file for regions without coverage.")
parser.add_argument("input_bam", help="Input bam file. Required.", type = str)
parser.add_argument("-i", "--input_fasta", help="Input fasta file for contig merging. Optional. If not specified, will only output linkage .gfa and .gv for manual merging.", type = str)
parser.add_argument("-o","--output", help="Prefix for output files.", type = str)
args = parser.parse_args()


# To count the number of reads aligned to a region
def countReads(contig,coords_to_check):
    cov_arr = samfile.count_coverage(contig=contig, start=coords_to_check[0],stop=coords_to_check[1])
    cov = sum([sum(cov_arr[x]) for x in range(0,4,1)]) / 50
    return cov

# Given a fasta entry with side to trim, returns new start and end coordinates
def trimFasta(contig_side_to_trim):
    # Input format: "tig-side"
    tig = contig_side_to_trim.split("-")[0]
    side = contig_side_to_trim.split("-")[1]
    coords_at_a_time = 50
    mincov = 10

    if tig in fastafile.references:
        cov = 0
        coords_to_trim = 0

        if side == "start":
            coords_to_check = [0,coords_at_a_time]
            while cov < mincov:
                cov = countReads(tig,coords_to_check)

                if cov <= mincov:
                    coords_to_trim += coords_at_a_time # Update with new end coordinate
                    coords_to_check = [x+coords_at_a_time for x in coords_to_check] # And move to next 10 coordinates

        elif side == "end":
            stop = trimmed_fasta_coords[tig][1]
            coords_to_check = [stop-coords_at_a_time,stop]
            while cov <= mincov:
                cov = countReads(tig,coords_to_check)

                if cov <= mincov:
                    coords_to_trim = coords_to_trim - coords_at_a_time    # Update with new end coordinate
                    coords_to_check = [x-coords_at_a_time for x in coords_to_check] # And move to next 10 coordinates

    return (coords_to_trim)

# Given a region in the format tigXXX-side,
# returns all edges from input region
def findLink(region):
    # Collect all edges matching input region
    current_edges = []
    for (a,b) in edges_list:
        if a == region:
            current_edges.append(b)
        elif b == region:
            current_edges.append(a)
    return current_edges

# starting_node is a string of the format tig-start/end, list_of_edges a list of len 3
# If successful, return two strings: same and diff
# same is the middle node which has two edges and
# diff the second node
# If unsuccessful return None,None
def find3way(starting_node,list_of_edges):
    if len(list_of_edges) != 3:
        return None
    else:
        node1 = list_of_edges[0].split("-")[0]
        node2 = list_of_edges[1].split("-")[0]
        node3 = list_of_edges[2].split("-")[0]
        if node1 == node2:
            same = node1
            diff = list_of_edges[2]
        elif node1 == node3:
            same = node1
            diff = list_of_edges[1]
        elif node2 == node3:
            same = node2
            diff = list_of_edges[0]
        else:
            same = None

        # If two of the three edges are in fact to the same node,
        # they would be expected to have two edges each:
        # 1 to starting_node and 1 to the next node (diff)
        if same != None:
            # Collect edges from both ends of the middle node
            links_same_start = findLink("-".join( [same,"start"] ))
            links_same_end = findLink("-".join( [same,"end"] ))

            # First check that both have len == 2,
            # i.e. there are no other nodes involved
            if not len(links_same_start) == 2 or not len(links_same_end) == 2:
                return None

            # Then check that i and diff are in both
            if starting_node in links_same_start and diff in links_same_start \
            and starting_node in links_same_end and diff in links_same_end:
                # If so, the middle node seems ok. Last thing
                # to do is to check diff for any other edges
                links_diff = findLink(diff)

                # Expected len == 3: both ends of middle node,
                # and starting_node
                if len(links_diff) != 3:
                    return None

                if starting_node in links_diff \
                and "-".join( [same,"start"] ) in links_diff \
                and "-".join( [same,"end"] ) in links_diff:
                    return (same,diff)

# Given a node with only one edge, returns all linked nodes to this one
def create_linked_contig(starting_region):
    tig = starting_region.split("-")[0]
    side = starting_region.split("-")[1]
    if side == "end":
        orientation = "f"
    elif side == "start":
        orientation = "r"
    else:
        orientation = "u"

    links = [ "".join([tig,orientation]) ]
    i = findLink(starting_region)
    if len(i) == 1:
        i = i[0]
    elif len(i) == 3:
        s = find3way(starting_region, i)
        same, diff = s[0], s[1]
        links.append( "".join([same,"u"]) ) # Direction of middle node is unknown
        i = diff

    while i != None:
        tig = i.split("-")[0]
        side = i.split("-")[1]
        if side == "start":
            orientation = "f"
            opposite = "end"
        elif side == "end":
            orientation = "r"
            opposite = "start"
        else:
            orientation = "u"
            opposite = "u"

        links.append( "".join([tig,orientation]) )
        # Look at opposite side of current node for unqiue edges
        next_edge = "-".join([tig,opposite])

        l = findLink(next_edge)
        # If only one matching edge, check that connected node
        # has no other edges
        if len(l) == 1 and len(findLink(l[0])) == 1:
            i = l[0]

        # If 3 edges, there is a chance two of them are to the same
        # node, on opposite sides, and the last to another node
        # This means that the one where both ends are connected
        # belongs in between.
        elif len(l) == 3:
            s = find3way(i,l)
            if s != None:
                same, diff = s[0], s[1]
                # Congrats, you have found a hit!
                links.append( "".join([same,"u"]) ) # Direction of middle node is unknown
                i = diff
            else:
                i = None

        # If there are two edges, and the connect to the same node at opposite ends,
        # add this node in unknown orientation and break
        elif len(l) == 2 and l[0].split("-")[0] == l[1].split("-")[0]:
            links.append( "".join( [l[0].split("-")[0],"u"] ) )
            i = None

        else:
            i = None

    return links

def reverse_complement(nuclstring):
    rev_comped = []
    for l in reversed(nuclstring):
        if l == "A" or l == "a":
            rev_comped.append("T")
        elif l == "T" or l == "t":
            rev_comped.append("A")
        elif l == "C" or l == "c":
            rev_comped.append("G")
        elif l == "G" or l == "g":
            rev_comped.append("C")
        elif l == "N" or l == "n":
            rev_comped.append("N")
    return "".join(rev_comped)

# Reads a delta file opened before and sent to this.
# Returns
#def readDelta(delta):


# Given two overlapping nucleotide strings (excluding overhangs) and alignment information,
# returns the merged sequence
def createConsensus(delta,string1,string2):
    # delta = [int,int,int, ...]
    # strings = "ATCG"
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
    cons = []
    for i in range(0,len(new_string1),1):
        if new_string1[i] == new_string2[i]:
            cons.append(new_string1[i])
        elif new_string1[i] == ".":
            cons.append(new_string2[i])
        elif new_string2[i] == ".":
            cons.append(new_string1[i])
        else:
            cons.append("N")

    return "".join(cons)

# If running into a contig where the direction is unknown
# eg tig1u while going through the linked_contigs dict,
# instead rely on alignment information as to which direction
# the unknown one should be merged into the surrounding contigs.
# Give a delta file of alignment between two sequences, they query of
# unknown orientation
# Return ordered contigs based on alignment
# If only one alignment is found, the direction can still be inferred.
# If no alignment is found, return None.
def findDirection(delta):
    ori = ""
    score = 10000
    passed = False
    d = []
    for line in delta:
        line = line.strip()
        print(line)
        fields = line.split(" ")

        # Collect sequence lengths
        if line.startswith(">"):
            reflen = fields[2]
            querylen = fields[3]

        # Filter out correct alignment in the delta
        # New alignments have 7 fields
        if len(fields) == 7:
            alignment_ref = (int(fields[0]),int(fields[1]))
            alignment_query = (int(fields[2]),int(fields[3]))

            # Now we will look for alignments where the ref ends on
            # its' last coordinate aka reflen
            if alignment_ref[1] == reflen:

                # One possible alignment found. Check in which direction query aligns
                # If score is lower than previous best hit, use new hit
                if alignment_query[0] == 1 and fields[4] < score:
                    # Orientation is ref_f query_f
                    score = fields[4]
                    ori = "f"
                    passed = True

                elif alignment_query[0] == querylen and fields[4] < score:
                    # Orientation is ref_f query_r
                    score = fields[4]
                    ori = "r"
                    passed = True

        # Reconstruct delta
        if len(fields) == 1 and passed == True:
            deltaint = int(line)
            if deltaint != 0:
                d.append(deltaint)

            else:
                passed = False

    if ori != "":
        return (d, ori, alignment_ref, alignment_query)
    else:
        return None


def findDirectionStarting(delta):
    # First get trimming coordinates
    tig1 = contig1_side.split("-")[0]
    if contig1_side.split("-")[1] == "start":
        trim1 = [trimmed_fasta_coords[tig1][0] + trimFasta( "-".join([contig1,"start"])), trimmed_fasta_coords[tig1][1]]
    else:
        trim1 = [0, trimmed_fasta_coords[tig1][1] + trimFasta( "-".join([contig1,"end"]))]

    tig2 = contig2_side.split("-")[0] # Side doesn't matter for contig2 - trim both sides anyway
    trim2 = [trimmed_fasta_coords[tig2][0] + trimFasta( "-".join([contig2,"start"])), trimmed_fasta_coords[tig2][1] + trimFasta( "-".join([contig2,"end"]) )]

    tig3 = contig3_side.split("-")[0]
    if contig3_side.split("-")[1] == "start":
        trim3 = [trimmed_fasta_coords[tig3][0] + trimFasta( "-".join([contig3,"start"])), trimmed_fasta_coords[tig3][1]]
    else:
        trim3 = [0, trimmed_fasta_coords[tig3][1] + trimFasta( "-".join([contig3,"end"]))]

    # Then align with mummer
    ref_fasta = fastafile.fetch(reference=refname, start=trimmed_fasta_coords[refname][0], end=trimmed_fasta_coords[refname][1])
    # First get trimming coordinates
    # Trim both sides of starting tig
    trim1 = [trimmed_fasta_coords[starting_tig][0] + trimFasta( "-".join([starting_tig,"start"])), trimmed_fasta_coords[starting_tig][1] + trimFasta( "-".join([starting_tig,"end"]) )]

    tig2 = linked_tig.split("-")[0]
    if linked_tig.split("-")[1] == "start":
        RC = False
        trim2 = [trimmed_fasta_coords[tig2][0] + trimFasta( "-".join([tig2,"start"])), trimmed_fasta_coords[tig2][1]]
    else:
        RC = True
        trim2 = [0, trimmed_fasta_coords[tig2][1] + trimFasta( "-".join([tig2,"end"]))]

    # Then align with mummer. Use starting as ref and linked_tig as query
    # As we will go left to right starting from starting_tig,
    # reverse complement ref_fasta if edge is from the end if the linked_tig
    ref_fasta = fastafile.fetch(reference=starting_tig, start=trim1[0], end=trim1[1])
    reflen = len(ref_fasta)
    query_fasta = fastafile.fetch(reference=tig2, start=trim2[0], end=trim2[1])
    if RC == True:
        query_fasta = reverse_complement(query_fasta)

    # Write trimmed contigs to temp files for mummer alignment in working directory
    with open("ref.fa","w",encoding = "utf-8") as fastaout:
        fastaout.write(">"+tig2+"\n")
        fastaout.write(ref_fasta)
    with open("query.fa","w",encoding = "utf-8") as fastaout:
        fastaout.write(">"+starting_tig+"\n")
        fastaout.write(query_fasta)

    # Call mummer. Select overlap based on nucmer
    call(["nucmer", "ref.fa", "query.fa"])

    with open("out.delta","r",encoding = "utf-8") as delta:
        passed = False # To check if the expected alignment was found or not
        ori = []
        for line in delta:
            line = line.strip()
            fields = line.split(" ")

            # Filter out correct alignment in the delta
            # New alignments have 7 fields
            if len(fields) == 7:
                alignment_ref = (int(fields[0]),int(fields[1]))
                alignment_query = (int(fields[2]),int(fields[3]))

                # Because we revcomped earlier, we will now look for
                # alignments starting on coord 1 of query and either first or last coord
                # of ref
                if alignment_query[0] == 1:
                    score = fields[4]
                    # One possible alignment found. Check in which direction ref aligns
                    if alignment_ref[1] == reflen:
                        if RC == False:
                            # Orientation is ref_f query_f
                            ori.append( ("f",score) )
                        elif RC == True:
                            ori.append( ("r",score) )

                    elif alignment_ref[1] == 1:
                        if RC == False:
                            # Orientation is ref_r query_f
                            ori.append( ("r",score) )
                        elif RC == True:
                            ori.append( ("f",score) )

        # If aligned to both sides, use the one with best score
        # i.e. minimal amount of indels/mismatches
        if len(ori) > 1:
            if ori[0][1] > ori[1][1]:
                return ori[1][0]
            elif ori[0][1] < ori[1][1]:
                return ori[0][0]
            else:
                return None


# Create output file name if nothing specified
if args.output:
    outfilename = args.output
elif "/" in args.input_bam:
    outfilename = args.input_bam.split("/")[-1].split(".bam")[0]
else:
    outfilename = args.input_bam.split(".bam")[0]

if not args.input_fasta:
    print("No fasta file specified for merging. Pipeline finished.")

else:
    fastafile = pysam.FastaFile(args.input_fasta)
    samfile = pysam.AlignmentFile(args.input_bam, "rb")
    print("Found fasta file for merging: {0}".format(args.input_fasta))

    edges = { \
    "tig00000180-end": ["tig00000204-start","tig00000204-end", "tig00000519-start"], \
    "tig00000519-start": ["tig00000204-start","tig00000204-end"], \
    "tig00000319-end": ["tig00000094-start"], \
    "tig00000094-end": ["tig00000519-end","tig00000074-end"], \
    "tig00000430-start": ["tig00000039-end","tig00000039-start"], \
    "tig00000430-end": ["tig00000460-end"] \
    }

    # Testing reformatting of edges from dict to list of tuples
    edges_list = []
    for i in edges:
        for a in edges[i]:
            # If edge was not already added in the opposite direction, add it
            if ( (a,i) ) not in edges_list:
                edges_list.append( (i, a) )

    # Reformat unique edges to a dict {"linked_contig_X":[ ("old_contig_1f/r"), ("old_contig_2f/r"), ... ]}
    linked_contigs = {}
    done_edges = []
    all_connected_nodes = set()
    duplicates = []
    counter = 1 # To keep track of names of new contigs
    for i in edges_list:
        # Collect all nodes in a single set
        all_connected_nodes.update(i)

    # Create dict of coordinates to keep from all fasta entries
    # Format: {contig: [start_coord, end_coord]}
    # Start by filling with old coords, which will then be changed
    trimmed_fasta_coords = {}
    for i in range(0,len(fastafile.references),1):
        trimmed_fasta_coords[fastafile.references[i]] = [0, fastafile.lengths[i]]

    ##### Find starting node #####
    # A starting node needs to have the following:
    # 1. A unique edge to another node, which also only
    # connects back to the starting node.
    # 2. The opposite end of the starting node either has no edges,
    # or 2, 3 different, or more than 3 edges (making it impossible to determine the path)
    # If it has 3 edges of which 2 are to the same node (in different ends), it
    # might be possible to determine a path (depending on edges on the 3 connected nodes)
    for i in all_connected_nodes:
        # First check that there is only one edge, and connecting node has
        # only one (the same) edge
        # Special case if there are three edges: then check that find3way does
        # not return None
        if ( len(findLink(i)) == 1 \
        and len(findLink(findLink(i)[0])) == 1 ) \
        or ( len(findLink(i)) == 3 \
        and find3way(i,findLink(i)) != None ):
            # If we got to here, we have a candidate starting node with good edges
            # Now check if opposite end of current node has 1 edge, whose node has 1 edge,
            # or the special case above again, we are in the middle of a region to merge
            # If so, continue
            itig = i.split("-")[0]
            iside = i.split("-")[1]
            if iside == "start":
                iopposite = "end"
            else:
                iopposite = "start"

            opposite_end = "-".join([itig,iopposite])

            if ( len(findLink(opposite_end)) == 1 \
            and  len(findLink(findLink(opposite_end))) == 1) \
            or ( len(findLink(opposite_end)) == 3 \
            and find3way(opposite_end,findLink(opposite_end)) != None) \
            or ( len(findLink(opposite_end)) == 2 \
            and findLink(opposite_end)[0].split("-")[0] == findLink(opposite_end)[1].split("-")[0]):
                continue

            # Starting point found.
            # Check if it was added in opposite orientation
            if i.split("-")[0] in done_edges:
                continue

            # Now check for the special case where both ends of current node
            # connect to the same node end
            if ( len(findLink(i)) == 1 \
            and len(findLink(opposite_end)) == 1 \
            and findLink(i)[0].split("-")[0] == findLink(opposite_end)[0].split("-")[0]):
                #i = "-".join([itig,"u"])
                idir = findDirectionStarting(itig, findLink(i))
                if idir != None:
                    i = "-".join( [itig,idir] )
                else:
                    continue

            # Starting node identified
            linked_tig_ID = "linked_contig_"+str(counter)
            linked_contig = create_linked_contig(i)
            print("Starting node: "+i)
            done_edges.append(linked_contig[-1][:-1])
            linked_contigs[linked_tig_ID] = linked_contig
            counter += 1

    for i in linked_contigs:
        print(i,linked_contigs[i])

    # Next find new coordinates for all sequences to merge
    print("Trimming contig edges.: {0}".format(args.input_fasta))
    for edge in linked_contigs.values():
        for i in edge:
            tig = i[:-1]
            # Trim only sides where links are formed
            # If first connected node in f or last in r, trim last coordinate
            if (i == edge[0] and i[-1] == "f") \
            or (i == edge[-1] and i[-1] == "r"):
                trimmed_fasta_coords[tig] = [trimmed_fasta_coords[tig][0], trimmed_fasta_coords[tig][1] + trimFasta(tig+"-end")]
            # If last connected node in f or first in r, trim first coordinate
            elif (i == edge[-1] and i[-1] == "f") \
            or (i == edge[0] and i[-1] == "r"):
                trimmed_fasta_coords[tig] = [trimmed_fasta_coords[tig][0] + trimFasta(tig+"-start"), trimmed_fasta_coords[tig][1]]
            # If node is connected at both sides or side is unknown, trim both
            else:
                trimmed_fasta_coords[tig] = [trimmed_fasta_coords[tig][0] + trimFasta(tig+"-start"), trimmed_fasta_coords[tig][1] + trimFasta(tig+"-end")]

    # Then walk through the linked_contigs dict, gradually building each contig from left to right.
    # Sequences to be merged are saved in merged_seqs in the format
    # {new_linked_contig_ID: "ATCG..." }
    merged_seqs = {}
    leftovers = {} # Collect leftovers of broken contigs
    contigs_to_exclude = [] # Collect contigs that were merged into new linked contigs
    print("Building linked contigs.")
    for linked_tig in linked_contigs:
        # Go through every linked contig in the list up until the next to last;
        # at this step the last contig will be merged into the sequence
        tig_counter = 0 # To keep track of where in the linked contig list we are
        tig_started = False # To check if merged contig was initiated
        newseq = [] # List of nucleotide strings to merge.
        # Format: [contig1_left_overhang, contig1-contig2_aligned_consensus, contig2_right_overhang ...]

        for t in linked_contigs[linked_tig]:
            # If at the last contig in the list, just append it to newseq and then break
            if tig_counter == len(linked_contigs[linked_tig]) - 1:
                newseq.append(ref_fasta)
                break

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
                # For now assume that refdirection != u

            # If tig_started == True, i.e. this is the second or more run through the loop,
            # there is already a starting point, i.e. newseq.
            # Keep aligning and merging linked contigs to this one.
            else:
                ref_fasta = newref # remainder of query from previous iteration
                refdirection = ""

            reflen = len(ref_fasta) # Collect length of reference sequence

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

                # The following assumes that two sequences in forward-forward were aligned
                # If direction is unknown, determine direction from delta instead
                if querydirection == "u":
                    alignment = findDirection(delta)
                    if alignment != None:
                        d, ori, alignment_ref, alignment_query = alignment[0], alignment[1], alignment[2], alignment[3]
                        passed = True

                        # If orientation is reverse,
                        # reverse complement sequence and
                        if ori == "r":
                            unaligned_query_fasta = query_fasta[:alignment_query[1]]
                            unaligned_query_fasta = reverse_complement(unaligned_query_fasta)
                            alqueryseq = query_fasta[alignment_query[1]:alignment_query[0]]
                            alqueryseq = reverse_complement(aligned_query_fasta)

                            # See below for explanation of this
                            newseq.append(ref_fasta[:alignment_ref[0]])
                            alrefseq = ref_fasta[alignment_ref[0]:alignment_ref[1]]
                            newseq.append(createConsensus(d,alrefseq,alqueryseq))
                            newref = unaligned_query_fasta
                            tig_counter += 1
                            continue

                    # If direction cannot be determined,
                    # break the current contig
                    else:
                        leftovers[linked_tig] = linked_contigs[linked_tig][tig_counter:]
                        break

                else:
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
                    # Add beginning of ref sequence to the new contig (newseq),
                    # up until overlap begins.
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
