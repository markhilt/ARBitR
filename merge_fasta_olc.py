#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
.. module:: merge_fasta
    :synopsis: merge_fasta implements functions related to fasta file operations
     during the anvil pipeline.

Copyright (c) 2019, Johannesson lab
Licensed under the GPL3 license. See LICENSE file.
"""

import numpy as np
import pysam
from scipy.stats import t
import mappy as mp

# Included modules
import nuclseqTools as nt
import misc

# In dev:
class Overlapgraph:
    """Simple overlap graph
    Nodes are contigs, edges are overlaps between contigs in the form
    (node1, node1_ori, node2, node2_ori, mappy.Alignment).
    """
    def __init__(self, nodes = [], edges = []):
        self.nodes = nodes # [node1, node2, ...]
        self.edges = edges # [(node1, node1_ori, node2, node2_ori, mappy.Alignment), ... ]
        self.paths = [] # Paths through the graph. Called externally.
        self.partial_paths = [] # Also keep partial paths
        self.incomplete_paths = []  # Combinations of partial paths with None at
                                    # the gap position

    def __str__(self):
        return "Overlapgraph([{0} nodes, {1} edges, {2} complete paths, {3} incomplete_paths])".format( len(self.nodes), \
                                                                    len(self.edges), \
                                                                    len(self.paths), \
                                                                    len(self.incomplete_paths))

    def addNode(self, node):
        '''
        Adds node to graph
        '''
        self.nodes.append(node)

    def addEdges(self, edges):
        '''
        Adds edges (list) to graph.
        '''
        self.edges = self.edges + edges

    def nodes():
        '''
        Returns all nodes in graph
        '''
        return self.nodes

    def edges():
        '''
        Returns all edges in graph
        '''
        return self.edges

    def reverse_direction(self,dir):
        '''
        Reverse direction of dir.
        '''
        if dir == "+":
            return "-"
        else:
            return "+"

    def outgoing(self,node, direction):
        '''
        Returns all outgoing edges from node
        '''
        return  [ edge for edge in self.edges if edge[0] == node and edge[1] == direction ]

    def traverse(self, start, start_direction, target, target_direction, path = [], edge = None):
        '''
        Appends all possible paths between the nodes start and target to self.paths
        '''
        if edge != None:
            path = path + [edge]
        if start == target and start_direction == target_direction:
            self.paths.append(path)
        else:
            outgoing_edges = self.outgoing(start,start_direction)
            if outgoing_edges == [] and path != [] and start != target:
                self.partial_paths.append(path)
            for edge in outgoing_edges:
                if edge not in path:
                    self.traverse(edge[2], edge[3], target, target_direction, path, edge)

    def reverse_path(self, path):
        # Reverse the given path
        reversed = []
        for step in path[::-1]:
            for out in self.outgoing(step[2], self.reverse_direction(step[3])):
                if out[2] == step[0] and out[3] == self.reverse_direction(step[1]):
                    reversed.append(out)
        return reversed

    def resolve_incomplete_paths(self, start, start_ori, target, target_ori):
        '''Find possible incomplete paths.

        Description:
            Sometimes a complete path between start and target cannot be found,
            but there are still extensions from start and/or target. In these
            cases we also want to include these contigs in the final scaffold,
            and insert the gap where approproate.

        Args:
            start (str): starting contig name.
            start_ori (str): starting contig orientation.
            target (str): target contig name.
            target_ori (str): target contig orientation.
            incomplete_paths (merge_fasta.Overlapgraph.incomplete_paths):

        Returns:
            list: possible paths.
        '''
        # The simplest path is just start and target contigs with a gap in between,
        # use this as a starting point.
        self.incomplete_paths.append( [ (start, start_ori, target, target_ori, None) ] )

        # I think we first need to sort out paths starting from start and ones
        # starting from target
        forward_paths, reverse_paths = [],[]
        for path in self.partial_paths:
            if path[0][0] == start:
                forward_paths.append(path)
            elif path[0][0] == target:
                reverse_paths.append(path)

        # Then we can try to combine them in different ways
        # Only use each contig once.
        for path in forward_paths:
            used_contigs = [ p for p in path for p in (p[0], p[2])]

            if len(reverse_paths) > 0:
                # Try every member of reverse_paths
                for rev_path in reverse_paths:
                    # Find if there is a reverse path that doesn't use any of the same
                    # contigs as the forward one.
                    rev_path_contigs = [ p for p in rev_path for p in (p[0], p[2])]
                    if set(used_contigs).isdisjoint(rev_path_contigs) == False: # Check for overlaps
                        # If not, it's a possible path
                        #ipdb.set_trace()
                        rev_path = self.reverse_path(rev_path)

                        # TODO: find why sometimes only the alignment is found
                        # in one direction
                        if len(rev_path) > 0:
                            self.incomplete_paths.append(   path + \
                                                            [(path[-1][2], \
                                                            path[-1][3], \
                                                            rev_path[0][0], \
                                                            rev_path[0][1], \
                                                            None)] + rev_path )
            else:
                self.incomplete_paths.append( path + \
                                            [(path[-1][2], \
                                            path[-1][3], \
                                            target, \
                                            target_ori, \
                                            None)] )

        # We also need to do the same for the extensions from target,
        # in case there are no extensions from start
        for path in reverse_paths:
            used_contigs = [ p for p in path for p in (p[0], p[2])]

            if len(forward_paths) > 0:
                # Try every member of forward_paths
                for for_path in forward_paths:
                    # Find if there is a reverse path that doesn't use any of the same
                    # contigs as the forward one.
                    for_path_contigs = [ p for p in rev_path for p in (p[0], p[2])]
                    if not set(used_contigs).isdisjoint(for_path_contigs): # Check for overlaps
                        # If not, it's a possible path
                        # Sort out directions and append
                        rev_path = self.reverse_path(path)
                        if len(rev_path) > 0:
                            self.incomplete_paths.append(   for_path + \
                                                        [(for_path[-1][2], \
                                                        for_path[-1][3], \
                                                        rev_path[0][0], \
                                                        rev_path[0][1], \
                                                        None)] + rev_path )
            else:
                rev_path = self.reverse_path(path)
                if len(rev_path) > 0:
                    self.incomplete_paths.append( [(start, \
                                                    start_ori, \
                                                    rev_path[0][0], \
                                                    rev_path[0][1], None)] \
                                                    + rev_path )

def findOverlap(seq1,seq2):
    '''Finds and returns all prefix/suffix overlaps between the given two sequences.

    '''

    #### First find reference suffix overlaps. We will treat seq1 as reference
    # Write seq1 to fasta file
    with open("tmp.fasta", "w") as tmpfasta:
        tmpfasta.write(">tmp\n")
        tmpfasta.write(seq1)

    # Build index from seq1
    idx = mp.Aligner(fn_idx_in="tmp.fasta", preset="map-pb")#k=19, w=10, scoring=[1,4,6,26])

    # Align
    alignments = idx.map(seq2)

    # Iterate over alignment and search for overlapping ends
    suffix_overlaps, prefix_overlaps = [], []
    if alignments:
        for aln in alignments:

            # Filter for alignments starting within the first 1kb or ending
            # within the last 1 kb of seq1, and starting within the first or
            # last 1 kb of seq2
            # We also need to control that the alignment is in the right direction

            # Overlap is at suffix of reference
            if aln.r_en > aln.ctg_len - 1000:

                # Overlap is at prefix of query, in this case it must be in
                # forward orientation
                if aln.q_st < 1000 and aln.strand == 1:
                    suffix_overlaps.append(aln)

                # Overlap is at suffix of query, in this case it must be in
                # reverse orientation
                if aln.q_en > len(seq2)-1000 and aln.strand == -1:
                    suffix_overlaps.append(aln)

    ### Then do the same for prefix overlaps, by first reverse complementing
    # seq1
    seq1 = nt.reverse_complement(seq1)
    with open("tmp.fasta", "w") as tmpfasta:
        tmpfasta.write(">tmp\n")
        tmpfasta.write(seq1)

    idx = mp.Aligner(fn_idx_in="tmp.fasta", preset="map-pb")#k=19, w=10, scoring=[1,4,6,26])
    alignments = idx.map(seq2)
    prefix_overlaps = []
    if alignments:
        for aln in alignments:

            # Overlap is at suffix of (revcomped) reference
            if aln.r_en > aln.ctg_len - 1000:

                # Overlap is at prefix of query, in this case it must be in
                # reverse orientation
                if aln.q_st < 1000 and aln.strand == 1:
                    prefix_overlaps.append(aln)

                # Overlap is at suffix of query, in this case it must be in
                # forward orientation
                if aln.q_en > len(seq2)-1000 and aln.strand == -1:
                    prefix_overlaps.append(aln)

    return suffix_overlaps, prefix_overlaps

def findOverlaps(fasta):
    '''Compute all prefix/suffix overlaps of given sequences.

    Args:
        fasta: fasta in dict format where keys are headers and values sequences
    Returns:
        list: list of tuples describing the overlaps.
    '''
    overlaps = []
    for tig1, seq1 in fasta.items():
        for tig2, seq2 in fasta.items():

            # Avoid self-alignments
            if tig1 != tig2:

                # Find overlaps between suffix of seq1 and prefix or suffix of seq2
                suffix_ovls, prefix_overlaps = findOverlap(seq1, seq2)
                for ovl in suffix_ovls:
                    if ovl.strand > 0:
                        strand = "+"
                    else:
                        strand = "-"
                    overlaps.append( (tig1, "+", tig2, strand, ovl) )

                for ovl in prefix_overlaps:
                    if ovl.strand > 0:
                        strand = "+"
                    else:
                        strand = "-"
                    overlaps.append( (tig1, "-", tig2, strand, ovl) )

    return overlaps

def countReads(contig,coords_to_check, region_length):
    '''
    To count the number of reads aligned to a region
    '''
    cov_arr = samfile.count_coverage(contig, coords_to_check[0],coords_to_check[1])
    cov = sum([sum(cov_arr[x]) for x in range(0,4,1)]) / region_length
    return cov

def trimFasta(trimmed_fasta_coords, contig_side_to_trim):
    '''
    Given a fasta entry with side to trim, returns new start and end coordinates
    Input format: "tigs" or "tige"

    Returns:
        int: number of coordinates to trim off either start or end side of
            input contig. Negative if end side.
    '''

    tig = contig_side_to_trim[:-1]
    side = contig_side_to_trim[-1]
    coords_at_a_time = 50
    mincov = 50


    if tig in fastafile.references:
        ctg_len = trimmed_fasta_coords[tig][1]
        cov = 0
        coords_to_trim = 0

        if side == "s":
            coords_to_check = [0,coords_at_a_time]
            cov = countReads(tig,coords_to_check, coords_at_a_time)

            while cov <= mincov:
                coords_to_trim += coords_at_a_time # Update with new end coordinate
                coords_to_check = [x+coords_at_a_time for x in coords_to_check] # And move to next coordinates

                # If contig boundary was passed without coverage being enough,
                # use all of it
                if coords_to_trim >= ctg_len - coords_at_a_time:
                    coords_to_trim = 0
                    break

                cov = countReads(tig,coords_to_check, coords_at_a_time)

        elif side == "e":
            coords_to_check = [ctg_len-coords_at_a_time,ctg_len]
            cov = countReads(tig,coords_to_check, coords_at_a_time)

            while cov <= mincov:
                coords_to_trim -= coords_at_a_time    # Update with new end coordinate
                coords_to_check = [x-coords_at_a_time for x in coords_to_check] # And move to next coordinates

                # If contig boundary was passed without coverage being enough,
                # use all of it
                if abs(coords_to_trim) >= ctg_len - coords_at_a_time:
                    coords_to_trim = 0
                    break

                cov = countReads(tig,coords_to_check, coords_at_a_time)

    return coords_to_trim

def trimSequences(paths):
    '''Trim away low quality regions of input sequences

    Description:
        Because de novo assembled contigs often end in low quality regions
        that are of too poor sequence to find good overlaps between, we want to
        trim input contigs of regions where reads don't map. Only trim regions
        where there is a potential overlap, i.e. NOT at the start of the first
        contig and end of the last contig in a path.

    Args:
        paths (list): list of lists. Each nested list contains ordered
        graph_building.Junction objects.
    Returns:
        dict: trimmed_fasta_coords. Keys: input contig headers, values:
            start and end coordinates to keep.
    '''
    # trimmed_fasta_coords is a dict with coords to keep from original fasta
    # Format: {contig: [start_coord, end_coord]}
    # Start by filling with old coords, which will then be changed
    trimmed_fasta_coords = {}
    for idx, ctg in enumerate(fastafile.references):
        trimmed_fasta_coords[ctg] = [0, fastafile.lengths[idx]]

    # Then find new coordinates for all sequences to merge
    for path in paths:
        for junction in path:
            if junction.start != None:
                start_tig, start_side = junction.start[:-1], junction.start[-1]
            else:
                start_tig, start_side = None, None
            if junction.target != None:
                target_tig, target_side = junction.target[:-1], junction.target[:1]
            else:
                target_tig, target_side = None, None
            connections = junction.connections

            # Trim the sides of contigs where a junction is formed,
            # and don't trim places where there are no junctions.
            if start_side == "s":
                trimmed_fasta_coords[start_tig] =   [trimmed_fasta_coords[start_tig][0] + \
                                                    trimFasta(trimmed_fasta_coords, junction.start), \
                                                    trimmed_fasta_coords[start_tig][1]]
            elif start_side == "e":
                trimmed_fasta_coords[start_tig] =   [trimmed_fasta_coords[start_tig][0], \
                                                    trimmed_fasta_coords[start_tig][1] + \
                                                    trimFasta(trimmed_fasta_coords, junction.start)]
            if target_side == "s":
                trimmed_fasta_coords[target_tig] =  [trimmed_fasta_coords[target_tig][0] + \
                                                    trimFasta(trimmed_fasta_coords, junction.target), \
                                                    trimmed_fasta_coords[target_tig][1]]
            elif target_side == "e":
                trimmed_fasta_coords[target_tig] =  [trimmed_fasta_coords[target_tig][0], \
                                                    trimmed_fasta_coords[target_tig][1] + \
                                                    trimFasta(trimmed_fasta_coords, junction.target)]

            # Also trim everything in connections
            for conn in connections:
                trimmed_fasta_coords[conn] = [trimmed_fasta_coords[conn][0] + \
                                            trimFasta(trimmed_fasta_coords, conn+"s"), \
                                            trimmed_fasta_coords[conn][1] + \
                                            trimFasta(trimmed_fasta_coords, conn+"e")]
    return trimmed_fasta_coords

def chop_cigar(cigar):
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

def reverse_cigar(cigar_list):
    '''
    Reverses the input list of cigar operations.
    '''
    reversed = []
    for cig in cigar_list[::-1]:
        if cig[-1] == "M":
            reversed.append(cig)

        elif cig[-1] == "D":
            reversed.append(cig[:-1]+"I")
        elif cig[-1] == "I":
            reversed.append(cig[:-1]+"D")
    return reversed

def shortestPath(paths):
    '''Find the shortest path.

    Description: Each member of paths is a list, containing starting contig,
        starting contig direction, aligned contig, aligned contig direction,
        and alignment information in a mappy.Alignment object. The path with
        the one with the longest overlap distance is returned.
    Args:
        paths: graph_building.Overlapgraph.paths
    Returns:
        List: The path that contains the longest overlap distance.
    '''
    best_path = paths[0]

    # Allow for incomplete paths
    best_ovl = 0
    best_ovl += sum( [step[4].blen for step in best_path if step[4] != None ] )

    for path in paths:
        ovl_dist = 0
        ovl_dist += sum( [step[4].blen for step in path if step[4] != None] )
        if ovl_dist > best_ovl:
            best_path = path
            best_ovl = ovl_dist
    return best_path

def createConsensus(step,string1,string2):
    '''
    Given two overlapping nucleotide strings and mappy alignment information,
    merges and returns the nucleotide strings.
    '''
    ref_ori, target_ori = step[1], step[3]
    aln = step[4]

    # Reverse complement sequences if in - direction
    # If reference is in -, the overlap was already calculated
    # from the revcomp. I.e. no need to recalculate position and cigar string.
    if ref_ori == "-":
        string1 = nt.reverse_complement(string1)
    if target_ori == "-":
        string2 = nt.reverse_complement(string2)

    # If there is no alignment, merge by gap insertion
    if aln == None:
        output_string = [string1,"NNNNNNNNNN", string2]
        query_pos = 0
        gap = True

    else:
        cig_list = chop_cigar(aln.cigar_str)
        # Track position in both strings
        ref_pos, query_pos = aln.r_st, aln.q_st

        if target_ori == "-":
            cigar_list = reverse_cigar(cig_list)
            query_pos = len(string2) - aln.q_en

        # First add sequence of string1 up until alignment start
        output_string = [string1[:ref_pos]]

        # Then walk through the aligned region in the cigar string,
        # gradually building the sequence. Because of the tendency of PacBio
        # to miss some bases, we will use the extra base at every indel position
        for cig in cig_list:

            # If strings match, there is no problem
            # Use whatever sequence, they are the same
            # Then move both position trackers
            if cig[-1] == "M":
                output_string.append( string1[ref_pos:ref_pos+int(cig[:-1])] )
                ref_pos += int(cig[:-1])
                query_pos += int(cig[:-1])

            # If insertion in query, add extra bases from query and increase position
            elif cig[-1] == "I":
                output_string.append( string2[query_pos:query_pos+int(cig[:-1])] )
                query_pos += int(cig[:-1])

            # If deletion in query, add extra bases from reference and increase position
            elif cig[-1] == "D":
                output_string.append( string1[ref_pos:ref_pos+int(cig[:-1])] )
                ref_pos += int(cig[:-1])

        # After iterating over the whole alignment, add remaining bases from
        # query sequence
        output_string.append( string2[query_pos:] )
        gap = False

    return ''.join(output_string), query_pos, gap

def mergeSeq(path):
    '''Traverse the path to merge visited contigs.

    Description:
        Create a combined sequence of all contigs included in path. A merge
        will be performed at overlaps, and if there is no overlap, i.e. the
        <mappy.Alignment object> is None, a gap is introduced.
    Args:
        path (list): List of steps in the path. A step looks like:
            ('tig1_name', 'tig1_dir', 'tig2_name', 'tig2_dir', <mappy.Alignment object>),
            where dirs are either "+" or "-"
    Returns:
        str: merged_sequence, the merged sequence
        list: included_contigs, names of included contigs
    '''
    merged_sequence = ""
    included_contigs = []
    coords_to_trim1 = 0
    coords_to_trim2 = 0
    included_contigs.append(path[0][0])
    n_gaps, n_merges = 0,0

    for ovl in path:

        #print("\nNow merging: ", ovl, str(ovl[4]))
        consensus, query_pos, gap_bool = createConsensus(ovl, trimmed_fasta[ovl[0]], trimmed_fasta[ovl[2]])

        # Add bases from previous step's q_st
        merged_sequence = merged_sequence[:-coords_to_trim1] + consensus[coords_to_trim2:]

        coords_to_trim1 = len(trimmed_fasta[ovl[2]]) - query_pos
        coords_to_trim2 = query_pos
        included_contigs.append(ovl[2])

        # Count
        if gap_bool == True:
            n_gaps += 1
        else:
            n_merges += 1

    return merged_sequence, included_contigs, n_gaps, n_merges

def combine_paths(linkpath, contigs):
    '''Find if each junction in linkpath can be merged by alignment or gap
    introduction. Return a path where each step is a tuple describing the
    potential overlap.

    Args:
        linkpath (list): list of graph_building.Junction objects to be combined
            into a scaffold.
        contigs (dict): fasta in dict format. Keys: contig names, values: trimmed
            sequences.

    Returns:
        path (list): List of steps in the path. Each step is a tuple looking like:
            ('tig1_name', 'tig1_dir', 'tig2_name', 'tig2_dir', <mappy.Alignment object>),
            where dirs are either "+" or "-". If <mappy.Alignment object> is None,
            scaffold by gap introduction.
        all_edges (list): List of all edges that were formed.
    '''
    filled_path = []
    all_edges = []
    for junction in linkpath:

        to_overlap = {  junction.start[:-1]: contigs[junction.start[:-1]], \
                        junction.target[:-1]: contigs[junction.target[:-1]]}
        for conn in junction.connections:
            to_overlap[conn] = contigs[conn]

        # Build an overlap graph
        edges = findOverlaps(to_overlap)
        all_edges = all_edges + edges
        graph = Overlapgraph(list(to_overlap.keys()),edges)

        # Find which orientation to look for
        if junction.start[-1] == "s":
            to_overlap[junction.start[:-1]] = nt.reverse_complement(to_overlap[junction.start[:-1]])
            ref_dir = "-"
        else:
            ref_dir = "+"

        if junction.target[-1] == "e":
            target_dir = "-"
        else:
            target_dir = "+"

        # Find a path between start and target
        graph.traverse(junction.start[:-1],ref_dir, junction.target[:-1], target_dir)

        if graph.paths != []:
            bestpath = shortestPath(graph.paths)
        else:
            # Find partial paths starting from target
            def reverse_ori(ori):
                return "+" if ori == "-" else "-"

            graph.traverse(junction.target[:-1],reverse_ori(target_dir), junction.start[:-1], reverse_ori(ref_dir))
            graph.resolve_incomplete_paths(junction.start[:-1], ref_dir, junction.target[:-1], target_dir)
            bestpath = shortestPath(graph.incomplete_paths)
        filled_path.append(bestpath)
    filled_path = [ step for partial_path in filled_path for step in partial_path]

    return filled_path, all_edges

def build_scaffolds(paths):
    scaffold_sequences, scaffold_correspondences = {}, {}
    all_edges = []
    n_gaps, n_merges = 0,0
    misc.printstatus("Number of paths: "+str(len(paths)))

    for idx, path in enumerate(paths):
        misc.printstatusFlush("[ SCAFFOLDING ]\t" + misc.reportProgress(idx+1, len(paths)))

        # Skip unknown ending contigs, for now
        # Also avoid two listcomps...
        if None in [ junction.start for junction in path] \
        or None in [ junction.target for junction in path]:
            continue

        # Collect all relevant sequences from fasta
        linked_contigs = [ [junction.start[:-1], junction.target[:-1]] + \
                            junction.connections for junction in path]
        linked_contigs = [ step for partial_path in linked_contigs for step in partial_path]
        to_overlap = { k:v for k,v in trimmed_fasta.items() if k in linked_contigs}

        # Start overlapping
        filled_path, edges = combine_paths(path, to_overlap)
        all_edges.append(edges)

        # Create scaffold
        scaffold_sequence, included, ng, nm = mergeSeq(filled_path)
        scaffold_sequences["scaffold_"+str(idx)] = scaffold_sequence
        scaffold_correspondences["scaffold_"+str(idx)] = included
        n_gaps += ng
        n_merges += nm

    misc.printstatus("[ SCAFFOLDING ]\t" + misc.reportProgress(idx+1, len(paths)))

    return scaffold_sequences, scaffold_correspondences, all_edges, n_gaps, n_merges

def main(input_fasta, input_bam, paths):
    '''Controller for merge_fasta.

    Args:
        input_fasta (str): Path to fasta file to create scaffolds from.
        input_bam (str): Path to bam file of reads mapped to input_fasta.
        paths (list): list of lists containing graph_building.Junction objects
            describing the paths inferred previously during the pipeline.
    Returns:
        dict: scaffolded fasta to output. Keys: fasta headers. Values: the
            resulting sequence.
        dict: correspondence, which contigs went into which scaffold.
    '''

    global samfile
    global fastafile
    fastafile = pysam.FastaFile(input_fasta)
    samfile = pysam.AlignmentFile(input_bam, "rb")

    # Get trim coordinates based on read mappings in samfile
    misc.printstatus("Trimming contig ends...")
    trimmed_fasta_coords = trimSequences(paths)

    # Trim fasta
    global trimmed_fasta
    trimmed_fasta = {}
    for tig in samfile.references:
        trimmed_fasta[tig] = fastafile.fetch(   reference=tig, \
                                                start=trimmed_fasta_coords[tig][0], \
                                                end=trimmed_fasta_coords[tig][1])
    samfile.close()
    fastafile.close()

    # Start finding overlaps
    misc.printstatus("Creating scaffolds...")
    scaffold_sequences, scaffold_correspondences, all_edges, n_gaps, n_merges = build_scaffolds(paths)
    all_edges = [edge for ls in all_edges for edge in ls]

    misc.printstatus("Scaffolding completed.")
    misc.printstatus("Number of aligned merges: {}".format(str(n_merges)))
    misc.printstatus("Number of gaps introduced: {}".format(str(n_gaps)))
    #complete_overlap_graph = Overlapgraph(list(trimmed_fasta.keys()), all_edges)
    #writeGfa(complete_overlap_graph)

    # Collect contigs that were not put into a scaffold
    misc.printstatus("Collecting leftover sequences.")
    used_contigs = [ tig for value in scaffold_correspondences.values() for tig in value]
    leftover_contigs = [ ctg for ctg in trimmed_fasta.keys() if ctg not in used_contigs ]

    for idx,tig in enumerate(leftover_contigs):
        scaffold_sequences["unplaced_contig_"+str(idx)] = trimmed_fasta[tig]
        scaffold_correspondences["unplaced_contig_"+str(idx)] = [tig]

    return scaffold_sequences, scaffold_correspondences
