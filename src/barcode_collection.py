#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""barcode_collection

barcode_collection implements methods to search the input
bam file for reads and extract their barcodes.

Copyright (c) 2020, Markus Hiltunen
Licensed under the GPL3 license. See LICENSE file.
"""

import numpy as np
import pysam
from scipy.stats import t
from collections import defaultdict

import misc

def collectGEMs(window, mapq, quant):
    """Collect the barcodes from the given window, if they fulfill some criteria.
    Args:
        window (str): Genomic region to collect barcodes from.
        mapq (int): Mapping quality cutoff to consider a read.
        quant (int): Cutoff for number of reads per barcode.
    Returns:
        set: BC_set, set of barcodes from this region.
    """
    tig, start, stop = window[0], window[1], window[2]
    reads = samfile.fetch(tig, start, stop)
    BC_set = set()
    occurrences = defaultdict(int)
    for read in reads:
        #if read.has_tag('BX') and read.mapping_quality >= mapq:
        if read.has_tag('BX') \
        and read.mapping_quality >= mapq \
        and read.is_proper_pair \
        and not read.is_qcfail:
        #and not read.is_duplicate \
            BC = read.get_tag("BX")
            occurrences[BC] += 1

            if occurrences[BC] >= quant:
                BC_set.add(BC)

    return BC_set

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

def main(input_bam, contig_dict, region_size=20000, mapq=60, bc_quant = 2):
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
                misc.printstatusFlush("[ BARCODE COLLECTION ]\t" + \
                misc.reportProgress(idx, len(contig_dict.keys())))
            else:
                misc.printstatusFlush("[ BARCODE COLLECTION ]\t" + \
                misc.reportProgress(idx, len(contig_dict.keys())*2))

        # Collect barcodes from the window
        GEMs = collectGEMs( (contig, start, end), mapq, bc_quant)

        # If at least 100 barcodes in list, use it
        if len(GEMs) > 100:
            GEMlist[region] = GEMs

    if region[-1] == "a":
        misc.printstatus("[ BARCODE COLLECTION ]\t" + \
        misc.reportProgress(len(contig_dict.keys()), len(contig_dict.keys())))
    else:
        misc.printstatus("[ BARCODE COLLECTION ]\t" + \
        misc.reportProgress(len(contig_dict.keys())*2, len(contig_dict.keys())*2))
    samfile.close()

    return GEMlist
