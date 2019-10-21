#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
.. module:: barcode_collection
    :synopsis: barcode_collection implements methods to search the input
    bam file for reads and extract their barcodes.

Copyright (c) 2019, Johannesson lab
Licensed under the GPL3 license. See LICENSE file.
"""

import numpy as np
import pysam
from scipy.stats import t

import misc

def collectGEMs(window, mapq):
    """
    Collects the barcodes from the given window
    """
    tig, start, stop = window[0], window[1], window[2]
    reads = samfile.fetch(tig, start, stop)
    BC_set = set()
    occurrences = set()
    for read in reads:
        if read.has_tag('BX') == True and read.mapping_quality >= mapq:
            BC = read.get_tag("BX")

            # This part is to only add barcodes that are present in more than one
            # read, to avoid sporadic hits.
            if BC in occurrences:
                BC_set.add(BC)
            else:
                occurrences.add(BC)

    return BC_set

"""
def getWindows(contig, region_size, contig_length):
    '''
    To split up the contig into windows
    '''
    start = (0,region_size)

    # If contig is shorter than region_size, just take all of it to "end" as well
    if (contig_length - region_size) < 0:
        endstart = 0
    else:
        endstart = contig_length - region_size
    end = (endstart,contig_length)
    return [start,end]
"""

def getWindows(region_size, contig_dict):
    '''Generate windows from contig_dict
    Args:
        region_size (int): User specified window size.
        contig_dict (dict): Contigs to use and their lengths.
    Yields:
        iterator of windows
    '''

    for contig, length in contig_dict.items():
        if length <= region_size:
            # If contig is shorter than region_size, yield only one window of
            # the whole contig
            yield (contig+"a", 0, length)

        else:
            # Otherwise yield two windows of the contig: one for the start
            # and one for the end regions of the contig
            yield (contig+"s", 0, region_size)
            yield (contig+"e", length - region_size, length)

def main(input_bam, contig_dict, region_size=20000, mapq=60):
    global samfile
    samfile = pysam.AlignmentFile(input_bam, "rb")
    GEMlist = {} # Inappropriately named "list"

    # First step is to collect all barcodes (passing -q cutoff) that are aligned
    # to each contigs first and last regions (-l)
    misc.printstatus("Starting barcode collection. Found {0} contigs.".format(len(contig_dict.keys())))

    # Generate windows
    windows = getWindows(region_size, contig_dict)

    # Iterate over windows to collect barcodes sets
    for idx, window in enumerate(windows):
        # Unpack variables, for readability
        region, contig, start, end = window[0], window[0][:-1], window[1], window[2]

        # Print progress. Number of windows is dependent on if running on
        # backbone or on small contigs.
        if idx in range(0,100000000,20):
            if region[-1] == "a":
                misc.printstatusFlush("[ BARCODE COLLECTION ]\t"+misc.reportProgress(idx, len(contig_dict.keys())))
            else:
                misc.printstatusFlush("[ BARCODE COLLECTION ]\t"+misc.reportProgress(idx, len(contig_dict.keys())*2))

        # Collect barcodes from the window
        GEMs = collectGEMs( (contig, start, end), mapq)

        # If at least 100 barcodes in list, use it
        if len(GEMs) > 100:
            GEMlist[region] = GEMs

    '''
    for contig in contig_dict:

        # Report progress every 20 windows
        if donecontigs in range(0,100000000,20):
            misc.printstatusFlush("[ BARCODE COLLECTION ]\t"+misc.reportProgress(donecontigs, len(contig_dict.keys())))

        # Collect the barcodes from the start and end of contig and
        # put them into the GEMlist dict
        windowlist = getWindows(contig, region_size, contig_dict[contig])

        # Collect barcodes from every window
        for window in windowlist:
            GEMs = collectGEMs( (contig, window[0], window[1]), mapq)

            if window == windowlist[0]:
                region = contig + "s"
            elif window == windowlist[1]:
                region = contig + "e"

            # If at least 100 barcodes in list, use it
            if len(GEMs) > 100:
                GEMlist[region] = GEMs

        donecontigs += 1
    '''

    if region[-1] == "a":
        misc.printstatus("[ BARCODE COLLECTION ]\t"+misc.reportProgress(len(contig_dict.keys()), len(contig_dict.keys())))
    else:
        misc.printstatus("[ BARCODE COLLECTION ]\t"+misc.reportProgress(len(contig_dict.keys())*2, len(contig_dict.keys())*2))
    samfile.close()

    return GEMlist
