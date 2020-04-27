#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""merge_fasta

merge_fasta implements functions related to fasta file and nucleotide string
operations during the anvil pipeline.

Copyright (c) 2020, Markus Hiltunen
Licensed under the GPL3 license. See LICENSE file.
"""

import numpy as np
import pysam
from scipy.stats import t
import mappy as mp

import nuclseqTools as nt
import misc

class Overlapgraph:
    """Simple overlap graph
    Nodes are contigs, edges are overlaps between contigs in the form
    (node1, node1_ori, node2, node2_ori, mappy.Alignment).
    """
    def __init__(   self, nodes = set(), \
                    edges = set()):
        self.nodes = nodes # {node1, node2, ...}
        self.edges = edges # {(node1, node1_ori, node2, node2_ori, mappy.Alignment), ... }
        self.paths = [] # Paths through the graph. Called externally.
        self.partial_paths = [] # Also keep partial paths
        self.incomplete_paths = []  # Combinations of partial paths with None at
                                    # the gap position

    def __str__(self):
        return "Overlapgraph([\
                {0} nodes, \
                {1} edges, \
                {2} complete paths, \
                {3} incomplete_paths])".format( \
                len(self.nodes), \
                len(self.edges), \
                len(self.paths), \
                len(self.incomplete_paths))

    def addNode(self, node):
        '''Adds node to graph
        '''
        self.nodes.add(node)

    def addEdges(self, edges):
        '''Adds edges (list) to graph.
        '''
        self.edges.add(edges)

    def nodes():
        '''Returns all nodes in graph
        '''
        return self.nodes

    def edges():
        '''Returns all edges in graph
        '''
        return self.edges

    def outgoing(self, node, direction):
        '''Returns all outgoing edges from node
        '''
        return [ edge for edge in self.edges if edge[0] == node and edge[1] == direction ]

    def traverse(self, start, start_direction, target, target_direction):
        '''Appends all possible paths between the nodes start and target to self.paths

        Description:
            Look for paths connecting start and target nodes. Every combination
            of edges is considered, and partial paths, i.e. extensions from
            start that don't reach target, are also kept.

        Args:
            start (str): starting contig name.
            start_direction (str): starting contig orientation.
            target (str): target contig name.
            target_direction (str): target contig orientation.
        '''
        def __visit(n, n_dir, path):
            visited = [v[0] for v in path]

            # Check if target is reached yet.
            if n == target:
                # Also check that direction is correct. If wrong, there might be
                # a partial_path to keep
                if n_dir == target_direction:
                    self.paths.append(path)
                elif len(path) > 1:
                    self.partial_paths.append(path[:-1])

            # Check that n has not been visited before in the path
            elif n not in visited:
                # Look for outgoing edges from the current node
                outgoing = self.outgoing(n, n_dir)
                # If there aren't any, save partial path and return
                if len(outgoing) == 0 and path != []:
                        self.partial_paths.append(path)
                else:
                    # If there is at least one, continue traversal, creating
                    # branches if more than one outgoing edge
                    for edge in outgoing:
                        if edge[2] not in visited:
                            __visit(edge[2], edge[3], path+[edge])
                        # If edge[2] was visited before, there might be a partial
                        # path to save
                        elif path != []:
                            self.partial_paths.append(path)

        __visit(start, start_direction, [])

    def reverse_path(self, path):
        '''Reverse the given path
        '''
        def reverse_direction(dir):
            return "-" if dir == "+" else "+"
        reversed = []
        for step in path[::-1]: # Iterate over path in reverse
            # Iterate over outgoing edges
            for out in self.outgoing(step[2], reverse_direction(step[3])):
                # Look for the edge we're looking for
                if out[2] == step[0] and out[3] == reverse_direction(step[1]):
                    reversed.append(out)
                    break
                # If the reverse is not in the graph, we need to return None
                if step == path[0]:
                    return None
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
        '''
        # The simplest path is just start and target contigs with a gap in between.
        # Add as incomplete path
        self.incomplete_paths.append( [(start, start_ori, target, target_ori, None)] )

        # I think we first need to sort out paths starting from start and ones
        # starting from target
        forward_paths, reverse_paths = [],[]
        for path in self.partial_paths:
            if path[0][0] == start:
                forward_paths.append(path)
            elif path[0][0] == target:
                reverse_paths.append(path)

        # Then we can try to combine them in different ways
        for path in forward_paths:
            # Append the incomplete_path with the gap right before target
            self.incomplete_paths.append(path + \
                                        [(path[-1][2], \
                                        path[-1][3], \
                                        target, \
                                        target_ori, \
                                        None)] )

            # Then look through the reverse paths to see if there is a possible
            # combination, i.e. one that doesn't use any of the same contigs.
            if len(reverse_paths) > 0:
                # Try every member of reverse_paths
                for rev_path in reverse_paths:
                    # Find if there is a reverse path that doesn't use any of the same
                    # contigs as the forward one.
                    used_contigs = set([ p for p in path for p in (p[0], p[2])])
                    rev_path_contigs = set([ p for p in rev_path for p in (p[0], p[2])])
                    if not set(used_contigs).isdisjoint(rev_path_contigs): # Make sure there are no overlaps
                        # If no overlaps between the two sets, it's a possible path
                        rev_path = self.reverse_path(rev_path)

                        # TODO: find why sometimes only the alignment is found
                        # in one direction
                        if rev_path:
                            self.incomplete_paths.append(   path + \
                                                            [(path[-1][2], \
                                                            path[-1][3], \
                                                            rev_path[0][0], \
                                                            rev_path[0][1], \
                                                            None)] + rev_path )

        # We also need to do the same for the extensions from target,
        # in case there are no extensions from start
        for path in reverse_paths:

            if len(forward_paths) > 0:
                # Try every member of forward_paths
                for for_path in forward_paths:
                    # Find if there is a reverse path that doesn't use any of the same
                    # contigs as the forward one.
                    used_contigs = set([ p for p in path for p in (p[0], p[2])])
                    for_path_contigs = set([ p for p in for_path for p in (p[0], p[2])])
                    if not set(used_contigs).isdisjoint(for_path_contigs): # Make sure there are no overlaps
                        # If no overlaps between the two sets, it's a possible path
                        # Sort out directions and append
                        rev_path = self.reverse_path(path)
                        if rev_path:
                            self.incomplete_paths.append(for_path + \
                                                        [(for_path[-1][2], \
                                                        for_path[-1][3], \
                                                        rev_path[0][0], \
                                                        rev_path[0][1], \
                                                        None)] + rev_path )
            else:
                rev_path = self.reverse_path(path)
                if rev_path:
                    self.incomplete_paths.append(   [(start, \
                                                    start_ori, \
                                                    rev_path[0][0], \
                                                    rev_path[0][1], None)] \
                                                    + rev_path )

def findOverlap(seq1,seq2):
    '''Find overlaps between two sequences.

    Description:
        Use mappy to find overlaps between all prefix-suffix pairs of the
        two given nucleotide sequences.

    Args:
        seq1 (str): First nucleotide sequence, serving as reference.
        seq2 (str): Second nucleotide sequence, serving as query.

    Returns:
        list: suffix_overlaps, list of overlaps from the suffix of seq1
        list: prefix_overlaps, list of overlaps from the suffix of revcomped seq1
    '''
    alignment_preset = "map-pb"

    #### First find reference suffix overlaps. We will treat seq1 as reference
    # Write seq1 to temporary fasta file
    with open("tmp.fasta", "w") as tmpfasta:
        tmpfasta.write(">tmp\n")
        tmpfasta.write(seq1)

    # Build index from seq1
    #idx = mp.Aligner(fn_idx_in="tmp.fasta", preset="map-pb")#k=19, w=10, scoring=[1,4,6,26])
    idx = mp.Aligner(fn_idx_in="tmp.fasta", preset=alignment_preset)

    # Align
    alignments = idx.map(seq2)

    # Iterate over alignments and search for overlapping ends
    suffix_overlaps, prefix_overlaps = [], []
    if alignments:
        for aln in alignments:

            # Filter for alignments starting within the first 1kb or ending
            # within the last 1 kb of seq1, and starting within the first or
            # last 1 kb of seq2
            # We also need to control that the alignment is in the right direction

            # Overlap is at suffix of reference. Ideally it ends at the last
            # coordinate of reference, but because contig ends usually have
            # poor quality, the alignment might not reach that far. We can still
            # use the overlap, as long as there is an overhang from the query
            # sequence.
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

    #idx = mp.Aligner(fn_idx_in="tmp.fasta", preset="map-pb")#k=19, w=10, scoring=[1,4,6,26])
    idx = mp.Aligner(fn_idx_in="tmp.fasta", preset=alignment_preset)
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
    '''Compute all prefix/suffix overlaps of given sequences. Return only the
    longest overlap for each prefix-suffix combination.

    Args:
        fasta: fasta in dict format where keys are headers and values sequences
    Returns:
        list: list of tuples describing the overlaps.
    '''
    overlaps = {}
    for tig1, seq1 in fasta.items():
        for tig2, seq2 in fasta.items():

            # Avoid self-alignments
            if tig1 != tig2:

                # Find overlaps between suffix of seq1 and prefix or suffix of seq2
                suffix_ovls, prefix_overlaps = findOverlap(seq1, seq2)
                for ovl in suffix_ovls:
                    strand = "+" if ovl.strand > 0 else "-"

                    # If current prefix-suffix pair not in overlaps, append it
                    if (tig1, "+", tig2, strand) not in overlaps:
                        overlaps[(tig1, "+", tig2, strand)] = ovl
                    else:
                        # Else check if the new overlap has more matching
                        # bases, in that case use the new overlap
                        old_mlen = overlaps[(tig1, "+", tig2, strand)].mlen
                        if ovl.mlen > old_mlen:
                            overlaps[(tig1, "+", tig2, strand)] = ovl

                # prefix_overlaps are actually suffix overlaps from the revcomp
                # of seq1
                for ovl in prefix_overlaps:
                    strand = "+" if ovl.strand > 0 else "-"

                    # If current prefix-suffix pair not in overlaps, append it
                    if (tig1, "-", tig2, strand) not in overlaps.keys():
                        overlaps[(tig1, "-", tig2, strand)] = ovl
                    else:
                        # Else check if the new overlap has more matching
                        # bases, in that case use the new overlap
                        old_mlen = overlaps[(tig1, "-", tig2, strand)].mlen
                        if ovl.mlen > old_mlen:
                            overlaps[(tig1, "-", tig2, strand)] = ovl

    # Reformat to list of tuples
    return [tuple(list(k)+[ovl]) for k, ovl in overlaps.items()]

def countReads(contig,coords_to_check, region_length):
    '''
    To count the number of reads aligned to a region
    '''
    cov_arr = samfile.count_coverage(contig, coords_to_check[0],coords_to_check[1])
    cov = sum([sum(cov_arr[x]) for x in range(0,4,1)]) / region_length
    return cov

def trimFasta(trimmed_fasta_coords, contig_side_to_trim, mincov):
    '''
    Given a fasta entry with side to trim, returns new start and end coordinates
    Input format: "tigs" or "tige"

    Returns:
        int: number of coordinates to trim off either start or end side of
            input contig. Negative if end side.
    '''

    tig = contig_side_to_trim[:-1]
    side = contig_side_to_trim[-1]
    step_len = 50

    if tig in fastafile.references:
        ctg_len = trimmed_fasta_coords[tig][1]
        cov = 0
        coords_to_trim = 0

        if side == "s":
            coords_to_check = [0,step_len]
            cov = countReads(tig,coords_to_check, step_len)

            while cov <= mincov:
                coords_to_trim += step_len # Update with new end coordinate
                coords_to_check = [x+step_len for x in coords_to_check] # And move to next coordinates

                # If contig boundary was passed without coverage being enough,
                # use all of it
                if coords_to_trim >= ctg_len - step_len:
                    coords_to_trim = 0
                    break

                cov = countReads(tig,coords_to_check, step_len)

        elif side == "e":
            coords_to_check = [ctg_len-step_len,ctg_len]
            cov = countReads(tig,coords_to_check, step_len)

            while cov <= mincov:
                coords_to_trim -= step_len    # Update with new end coordinate
                coords_to_check = [x-step_len for x in coords_to_check] # And move to next coordinates

                # If contig boundary was passed without coverage being enough,
                # use all of it
                if abs(coords_to_trim) >= ctg_len - step_len:
                    coords_to_trim = 0
                    break

                cov = countReads(tig,coords_to_check, step_len)

    return coords_to_trim

def trimSequences(paths, mincov):
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
        mincov (int): Trim contig ends with lower average coverage than this
            value
    Returns:
        dict: trimmed_fasta_coords. Keys: input contig headers, values:
            start and end coordinates to keep, in addition to True or False
            for start and end if they were trimmed or not.
    '''
    # trimmed_fasta_coords is a dict with coords to keep from original fasta
    # Format: {contig: [start_coord, end_coord, bool, bool]}
    # Start by filling with old coords, which will then be changed
    trimmed_fasta_coords = {}
    for idx, ctg in enumerate(fastafile.references):
        trimmed_fasta_coords[ctg] = [0, fastafile.lengths[idx], False, False]

    # Then find new coordinates for all sequences to merge
    for idx, path in enumerate(paths):
        if idx in range(0,100000000,5):
            misc.printstatusFlush("[ TRIMMING ]\t" + misc.reportProgress(idx+1, len(paths)))

        for junction in path:
            if junction.start != None:
                start_tig, start_side = junction.start[:-1], junction.start[-1]
            else:
                start_tig, start_side = None, None
            if junction.target != None:
                target_tig, target_side = junction.target[:-1], junction.target[-1]
            else:
                target_tig, target_side = None, None
            connections = junction.connections

            # Trim the sides of contigs where a junction is formed,
            # and don't trim places where there are no junctions.
            if start_side == "s" \
            and trimmed_fasta_coords[start_tig][2] == False:
                trimmed_fasta_coords[start_tig] =   [trimmed_fasta_coords[start_tig][0] + \
                                                    trimFasta(trimmed_fasta_coords, junction.start, mincov), \
                                                    trimmed_fasta_coords[start_tig][1], \
                                                    True, \
                                                    trimmed_fasta_coords[start_tig][3]]
            elif start_side == "e" \
            and trimmed_fasta_coords[start_tig][3] == False:
                trimmed_fasta_coords[start_tig] =   [trimmed_fasta_coords[start_tig][0], \
                                                    trimmed_fasta_coords[start_tig][1] + \
                                                    trimFasta(trimmed_fasta_coords, junction.start, mincov), \
                                                    trimmed_fasta_coords[start_tig][2], \
                                                    True]
            if target_side == "s" \
            and trimmed_fasta_coords[target_tig][2] == False:
                trimmed_fasta_coords[target_tig] =  [trimmed_fasta_coords[target_tig][0] + \
                                                    trimFasta(trimmed_fasta_coords, junction.target, mincov), \
                                                    trimmed_fasta_coords[target_tig][1], \
                                                    True, \
                                                    trimmed_fasta_coords[target_tig][3]]
            elif target_side == "e" \
            and trimmed_fasta_coords[target_tig][3] == False:
                trimmed_fasta_coords[target_tig] =  [trimmed_fasta_coords[target_tig][0], \
                                                    trimmed_fasta_coords[target_tig][1] + \
                                                    trimFasta(trimmed_fasta_coords, junction.target, mincov), \
                                                    trimmed_fasta_coords[target_tig][2], \
                                                    True]

            # Also trim everything in connections
            for conn in connections:
                if not trimmed_fasta_coords[conn][2] == True \
                or not trimmed_fasta_coords[conn][3] == True:
                    trimmed_fasta_coords[conn] = [trimmed_fasta_coords[conn][0] + \
                                                trimFasta(trimmed_fasta_coords, conn+"s", mincov), \
                                                trimmed_fasta_coords[conn][1] + \
                                                trimFasta(trimmed_fasta_coords, conn+"e", mincov), \
                                                True, True]
    misc.printstatus("[ TRIMMING ]\t" + misc.reportProgress(idx+1, len(paths)))

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
    '''Find the path with the longest overlap distance.

    Description: Each member of paths is a list, containing starting contig,
        starting contig direction, aligned contig, aligned contig direction,
        and alignment information in a mappy.Alignment object. The path with
        the longest overlap distance is returned.
    Args:
        paths: graph_building.Overlapgraph.paths
    Returns:
        List: The path that contains the longest overlap distance.
    '''
    # Start at the first path
    best_path = paths[0]

    # Allow for incomplete paths
    best_ovl = sum( [step[4].blen for step in best_path if step[4] != None ] )
    for path in paths[1:]:
        ovl_dist = 0
        ovl_dist += sum( [step[4].blen for step in path if step[4] != None] )
        if ovl_dist > best_ovl:
            best_path = path
            best_ovl = ovl_dist
    return best_path

def createConsensus(step,string1,string2, gapsize):
    '''Given two overlapping nucleotide strings and mappy alignment information,
    merges and returns the nucleotide strings.

    Args:
        step: step in the path. Looks like:
        ("ctg1", "ctg1_orientation", "ctg2", "ctg2_orientation", <mappy.Alignment object>)
        string1: first nucleotide string
        string2: second nucleotide string
        gapsize (int): Gap size for scaffolding.
    Returns:
        Str: The nucleotide sequence starting at the alignment r_st, if any,
        otherwise at the start of the gap.
        Int: query_pos, query position.
        Bool: gap, if there is a gap in the sequence or not.
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
    if not aln:
        output_string = ["N"*gapsize, string2]
        gap = True
        ovl_len = 0

    else:
        cig_list = chop_cigar(aln.cigar_str)
        # Track position in both strings
        ref_pos, query_pos = aln.r_st, aln.q_st

        if target_ori == "-":
            cigar_list = reverse_cigar(cig_list)
            query_pos = len(string2) - aln.q_en

        output_string = []

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
        ovl_len = len(''.join(output_string))
        output_string.append( string2[query_pos:] )
        gap = False

    return ''.join(output_string), gap, ovl_len

def mergeSeq(path, gapsize):
    '''Traverse the path to merge visited contig sequences.

    Description:
        Create a combined sequence of all contigs included in path. A merge
        will be performed at overlaps, and if there is no overlap, i.e. the
        <mappy.Alignment object> is None, a gap is introduced. The sequence is
        built by iteration over the path. During each iteration, if there is
        an overlap, the merged sequence is first trimmed for the number of bases
        corresponding to the reference overlap length, and then extended,
        starting at the overlap beginning and including all of the query sequence.
        If there is no overlap, the merged sequence is not trimmed and the
        extension begins at the gap.
    Args:
        path (list): List of steps in the path. Each step is a tuple looking like:
            ('tig1_name', 'tig1_dir', 'tig2_name', 'tig2_dir', <mappy.Alignment object>),
            where dirs are either "+" or "-".
        gapsize (int): Gap size for scaffolding.
    Returns:
        str: merged_sequence, the merged/gapped nucleotide sequence
        list: included_contigs, names of included contigs
        int: n_gaps, number of gaps
        int: n_merges, number of aligned merges
        list: bed_coords, coordinates of the different merged feaures
    '''

    # Initiate the merged sequence by adding the first reference sequence
    # Control for reverse orientation
    first_step = path[0]
    merged_sequence = trimmed_fasta[first_step[0]] if first_step[1] == "+" \
    else nt.reverse_complement(trimmed_fasta[first_step[0]])
    bed_coords = [( "0", str(len(merged_sequence)), first_step[0], \
                    "0", first_step[1], "0", str(len(merged_sequence)), \
                    "0,0,255")]


    # Keep track of which sequences that we've used
    included_contigs = [first_step[0]]
    reference_position = 0
    query_position = 0
    n_gaps, n_merges = 0,0 # Count how many gaps and merges we make

    for ovl in path:
        # Next, iterate over the steps in the path. At each path, add sequence
        # corresponding to the alignment and remaining query sequence, or if there
        # is no alignment, add sequence corresponding to the gap and the whole
        # query sequence.
        r_name, q_name = ovl[0], ovl[2]
        r_ori, q_ori = ovl[1], ovl[3]
        aln = ovl[4]
        seq1, seq2 = trimmed_fasta[r_name], trimmed_fasta[q_name]

        # If there is an alignment at this step in the path, remove bases from
        # the merged sequence from where the alignment starts and onwards
        if aln:
            bases_to_trim = len(trimmed_fasta[r_name]) - aln.r_st
            if bases_to_trim > 0:
                merged_sequence = merged_sequence[:-bases_to_trim]
            else:
                # If we got here it means that the whole reference sequence is
                # in the alignment. In this case, we cannot simply remove the
                # number of bases in the merged sequence that correspond to the
                # query length, because the query may already have been aligned
                # into the merged sequence. Instead, recalculate the alignment from
                # the merged sequence.
                suffix_overlaps, prefix_overlaps = findOverlap(merged_sequence, trimmed_fasta[q_name])
                # We are only interested in suffix_overlaps in the expected
                # orientation of query
                exp_ori = 1 if q_ori == "+" else -1
                suffix_overlaps = [sovl for sovl in suffix_overlaps if sovl.strand == exp_ori]

                # If more than one overlap, use the longest one
                if len(suffix_overlaps) > 0:
                    best_ovl = suffix_overlaps[0]
                    for sovl in suffix_overlaps[1:]:
                        ovl_dist = sovl.blen
                        if ovl_dist > best_ovl.blen:
                            best_ovl = sovl

                    r_name, r_ori = "tmp", "+"
                    seq1 = merged_sequence

                else:
                    aln = None

        consensus, gap_bool, ovl_len = createConsensus( (r_name, r_ori, \
                                                        q_name, q_ori, aln), \
                                                        seq1, seq2, gapsize)

        # Count & update bed coords
        if gap_bool == True:
            n_gaps += 1
        else:
            n_merges += 1
            bed_coords.append( (str(len(merged_sequence)+1), \
                                str(len(merged_sequence)+ovl_len+1), \
                                "overlap", "0", q_ori, \
                                str(len(merged_sequence)+1), \
                                str(len(merged_sequence)+ovl_len+1), "255,0,0"))

        # Extend merged_sequence by the new consensus
        merged_sequence += consensus
        included_contigs.append(q_name)
        # Append seq2 to bed_coords
        bed_coords.append( (str(len(merged_sequence)-len(seq2)+1), \
                            str(len(merged_sequence)+1), \
                            q_name, "0", q_ori, \
                            str(len(merged_sequence)-len(seq2)+1), \
                            str(len(merged_sequence)+1), "0,0,255" ) )

    return merged_sequence, included_contigs, n_gaps, n_merges, bed_coords

def combine_paths(linkpath):
    '''Find if each junction in linkpath can be merged by alignment or gap
    introduction. Return a path where each step is a tuple describing the
    potential overlap.

    Args:
        linkpath (list): list of graph_building.Junction objects to be combined
            into a scaffold.

    Returns:
        filled_path (list): List of steps in the path. Each step is a tuple looking like:
            ('tig1_name', 'tig1_dir', 'tig2_name', 'tig2_dir', <mappy.Alignment object>),
            where dirs are either "+" or "-". If <mappy.Alignment object> is None,
            scaffold by gap introduction.
        all_edges (list): List of all edges that were formed.
    '''
    filled_path = []
    all_edges = []
    for junction in linkpath:
        # Take a different route if there is a None in the junction
        if junction.start and junction.target:
            # Extract all sequences in this junction to a separate dict
            to_overlap = {  junction.start[:-1]: trimmed_fasta[junction.start[:-1]], \
                            junction.target[:-1]: trimmed_fasta[junction.target[:-1]]}
            for conn in junction.connections:
                to_overlap[conn] = trimmed_fasta[conn]

            ref_dir = "+" if junction.start[-1] == "e" else "-"
            target_dir = "-" if junction.target[-1] == "e" else "+"
            # Find which direction to traverse the graph in.
            if ref_dir == "-":
                to_overlap[junction.start[:-1]] = nt.reverse_complement(to_overlap[junction.start[:-1]])

            # Find all overlaps between sequences in question and if any are found,
            # create an overlap graph to describe them.
            edges = findOverlaps(to_overlap)
            if edges:
                all_edges = all_edges + edges
                graph = Overlapgraph(   set(to_overlap.keys()), \
                                        set(edges))

                # Try to find a path between start and target
                # Partial paths are also saved in the graph
                graph.traverse( junction.start[:-1], ref_dir, \
                                junction.target[:-1], target_dir)

                if len(graph.paths) == 1:
                    # If one path, use it right away
                    bestpath = list(graph.paths)[0]

                elif len(graph.paths) > 1:
                    # If there is more than one path, find the shortest path
                    bestpath = shortestPath(graph.paths)

                else:
                    # If no path, try to find incomplete paths
                    # Find partial paths also starting from target
                    def reverse_ori(ori):
                        return "+" if ori == "-" else "-"

                    graph.traverse( junction.target[:-1], \
                                    reverse_ori(target_dir), \
                                    junction.start[:-1], \
                                    reverse_ori(ref_dir))
                    graph.resolve_incomplete_paths( junction.start[:-1], \
                                                    ref_dir, \
                                                    junction.target[:-1], \
                                                    target_dir)
                    bestpath = shortestPath(graph.incomplete_paths)

                filled_path.extend(bestpath)
            else:
                # If no edges, put None at overlap and drop any connections
                filled_path.append( (junction.start[:-1], ref_dir, \
                                    junction.target[:-1], target_dir, None))

        # If there is a None in the junction
        else:
            if junction.start:
                to_overlap = {junction.start[:-1]: trimmed_fasta[junction.start[:-1]]}
                for conn in junction.connections:
                    to_overlap[conn] = trimmed_fasta[conn]
                ref_dir = "+" if junction.start[-1] == "e" else "-"
                # Find which direction to traverse the graph in.
                if ref_dir == "-":
                    to_overlap[junction.start[:-1]] = nt.reverse_complement(to_overlap[junction.start[:-1]])
                edges = findOverlaps(to_overlap)
                if edges:
                    all_edges = all_edges + edges
                    graph = Overlapgraph(   set(to_overlap.keys()), \
                                            set(edges))

                    # Try to find partial paths extending from start
                    # We can use graph.traverse for this by defining a target
                    # that's not in the graph, i.e. complete paths will be
                    # impossible to find
                    graph.traverse( junction.start[:-1], ref_dir, \
                                    "not_a_target", "+")

                    if graph.partial_paths:
                        bestpath = shortestPath(graph.partial_paths)
                        filled_path.extend(bestpath)

            elif junction.target:
                #if junction.target == "8e":
                #    import ipdb; ipdb.set_trace()
                # Same as above, except in the opposite direction
                to_overlap = {junction.target[:-1]: trimmed_fasta[junction.target[:-1]]}
                for conn in junction.connections:
                    to_overlap[conn] = trimmed_fasta[conn]
                #target_dir = "-" if junction.target[-1] == "e" else "+"
                target_dir = "+" if junction.target[-1] == "e" else "-"
                if target_dir == "-":
                    to_overlap[junction.target[:-1]] = nt.reverse_complement(to_overlap[junction.target[:-1]])
                edges = findOverlaps(to_overlap)
                if edges:
                    all_edges = all_edges + edges
                    graph = Overlapgraph(   set(to_overlap.keys()), \
                                            set(edges))
                    graph.traverse( junction.target[:-1], target_dir, \
                                "not_a_target", "+")
                    if graph.partial_paths:
                        # If there is a partial_path, we need to reverse it
                        # and add it to the beginning of filled_path
                        reversed_paths = [graph.reverse_path(path) for path in graph.partial_paths]
                        # Remove None's, in separate step to avoid calculating
                        # non-existent paths twice
                        reversed_paths = [p for p in reversed_paths if p]

                        if reversed_paths:
                            bestpath = shortestPath(reversed_paths)
                            filled_path = bestpath + filled_path

    return filled_path, all_edges

def build_scaffolds(paths, gapsize):
    scaffold_sequences, scaffold_correspondences = {}, {}
    all_edges = []
    n_gaps, n_merges = 0,0
    misc.printstatus("Number of paths: "+str(len(paths)))
    bed = {} # To collect bed coordinates

    for idx, path in enumerate(paths):
        misc.printstatusFlush("[ SCAFFOLDING ]\t" + misc.reportProgress(idx+1, len(paths)))

        # Collect all relevant sequences from fasta
        linked_contigs = [ [junction.start[:-1], junction.target[:-1]] + \
                            junction.connections for junction in path \
                            if junction.start and junction.target]
        linked_contigs = [ step for partial_path in linked_contigs for step in partial_path]

        # Start overlapping
        filled_path, edges = combine_paths(path)
        # It is possible that there is no filled_path, in the case that the
        # path had a single junction which had a None at junction.start or
        # junction.target and no overlaps were found. In this case, continue.
        if filled_path:
            all_edges.extend(edges)

            # Create scaffold
            scaffold_sequence, included, ng, nm, bed_coords = mergeSeq(filled_path, gapsize)
            scaffold_sequences["scaffold_"+str(idx)] = scaffold_sequence
            scaffold_correspondences["scaffold_"+str(idx)] = included
            bed["scaffold_"+str(idx)] = bed_coords
            n_gaps += ng
            n_merges += nm

    misc.printstatus("[ SCAFFOLDING ]\t" + misc.reportProgress(idx+1, len(paths)))

    return scaffold_sequences, scaffold_correspondences, all_edges, n_gaps, n_merges, bed

def main(input_fasta, input_bam, paths, mincov, gapsize):
    '''Controller for merge_fasta.

    Args:
        input_fasta (str): Path to fasta file to create scaffolds from.
        input_bam (str): Path to bam file of reads mapped to input_fasta.
        paths (list): list of lists containing graph_building.Junction objects
            describing the paths inferred previously during the pipeline.
        mincov (int): Minimum average coverage for trimming.
        gapsize (int): Gap size when scaffolding by gap introduction.
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
    trimmed_fasta_coords = trimSequences(paths, mincov)

    # Trim fasta
    global trimmed_fasta
    trimmed_fasta = {}
    for tig in samfile.references:
        trimmed_fasta[tig] = fastafile.fetch(reference=tig, \
                                             start=trimmed_fasta_coords[tig][0], \
                                             end=trimmed_fasta_coords[tig][1])
    samfile.close()
    fastafile.close()

    # Start finding overlaps
    misc.printstatus("Creating scaffolds...")

    scaffold_sequences, scaffold_correspondences, all_edges, \
    n_gaps, n_merges, bed = build_scaffolds(paths, gapsize)

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

    return scaffold_sequences, scaffold_correspondences, bed
