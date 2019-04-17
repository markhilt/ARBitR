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
    BC_list = []
    for read in reads:
        if read.has_tag('BX') == True and read.mapping_quality >= mapq:
            BC = read.get_tag("BX")
            BC_list.append(BC)
    BC_list = set(BC_list)

    return BC_list

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

def main(input_bam, region_size, mapq, contig_dict):
    global samfile
    samfile = pysam.AlignmentFile(input_bam, "rb")
    GEMlist = {} # Inappropriately named "list"

    # First step is to collect all barcodes (passing -q cutoff) that are aligned
    # to each contigs first and last regions (-l)
    misc.printstatus("Starting barcode collection. Found {0} contigs.".format(len(contig_dict.keys())))
    donecontigs = 0
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

            if GEMs:
                GEMlist[region] = GEMs

        donecontigs += 1

    misc.printstatus("[ BARCODE COLLECTION ]\t"+misc.reportProgress(donecontigs, len(contig_dict.keys())))

    samfile.close()

    return GEMlist
