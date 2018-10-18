#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pysam
import calcESD
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Reads through a bam file produced by \
                                            longranger align to collect GEM barcodes \
                                            in sliding windows of a genome. \
                                            Then compares barcodes of every window \
                                            to all other windows to detect possible misassemblies.")

parser.add_argument("input_bam", help="Input bam file. Required.", type = str)
parser.add_argument("-w","--windowsize", help="Size of sliding windows. [10000]", default = 10000, type = int)
parser.add_argument("-s","--stepsize", help="Step size of where to start new window. [9000]", default = 9000, type = int)
parser.add_argument("-q","--mapq", help="Mapping quality cutoff value. [60]", default = 60, type = int)
parser.add_argument("-m","--expected_molecule_size", help="Expected mean molecule size that went into Chromium sequencing (bp). [45000]", default = 45000, type = int)
parser.add_argument("-o","--output", help="Prefix for output files.", type = str)
args = parser.parse_args()

def getOut():
    '''
    Creates a prefix for output files
    '''
    if args.output:
        outfilename = args.output
    elif "/" in args.input_bam:
        outfilename = args.input_bam.split("/")[-1].split(".bam")[0]
    else:
        outfilename = args.input_bam.split(".bam")[0]
    return outfilename

def collectRefLen(samfile):
    '''
    Returns contigs and their lengths as a dict
    '''
    i = 0
    contig_dict = {}
    contig_list = samfile.references
    while i < len(contig_list):
        contig_dict[contig_list[i]] = samfile.lengths[i]
        i += 1

    return contig_dict

def createWindows(contig_dict):
    '''
    To split up each contig in contig_dict into windows
    '''
    contig_windows = []

    for k,v in contig_dict.items():
        startcoord = 0
        length = v

        while startcoord < length:
            contig_windows.append( (k, startcoord, startcoord + args.windowsize) )
            startcoord = startcoord + args.stepsize

    return set(contig_windows)

def collectGEMs(window):
    """
    Collects the barcodes from the given window
    """
    tig, start, stop = window[0], window[1], window[2]
    reads = samfile.fetch(tig, start, stop)
    BC_list = []
    for read in reads:
        if read.has_tag('BX') == True and read.mapping_quality >= args.mapq:
            BC = read.get_tag("BX")
            BC_list.append(BC)
    BC_list = set(BC_list)

    return BC_list

def recodeGEMbase(gembase):
    """
    Given a dict of windows and their barcodes as strings, returns
    """

def reportProgress(current,total):
    return "Completed: {0}% ({1} out of {2})".format( str(round( (current / total) * 100, 2)), current, total)

def createBarcodeDatabase(windows):
    """
    Goes through every window in windows, calls collectGEMs
    to collect barcodes from each window. Returns a dict as { (tig, start-coord, stop-coord): ["barcode1", "barcode2", ... ] }
    """
    GEMbase = {}
    done = 0
    for w in windows:
        if done in range(0,100000000,300):
            print(reportProgress(done, len(windows)))

        GEMbase[w] = collectGEMs(w)
        done += 1

    return GEMbase

def compareGEMlibs_np(lib1,lib2):
    '''
    Compares two numpy arrays of barcodes indices, collects all shared ones and
    counts them. Returns fraction of shared barcodes.
    '''
    shared = np.intersect1d(lib1,lib2) # Find shared ones
    totallength = len(lib1) + len(lib2) - len(shared) # Total number of unshared barcodes

    # Find the fraction of shared barcodes
    if totallength != 0:
        fraction = len(shared) / totallength
    else:
        fraction = 0

    return fraction

def compareGEMlibs_str(lib1,lib2):
    '''
    Compares two sets of barcodes, collects all shared ones and
    counts them. Returns fraction of shared barcodes.
    '''
    shared = lib1.intersection(lib2) # Find shared ones
    totallength = len(lib1) + len(lib2) - len(shared) # Total number of unshared barcodes

    # Find the fraction of shared barcodes
    if totallength != 0:
        fraction = len(shared) / totallength
    else:
        fraction = 0

    return fraction

def pairwiseComps(GEMbase):
    """
    Does all pairwise comparisons of barcodes using numpy and pandas
    """
    # First create a 2-dim array of length len(GEMbase) x len(GEMbase),
    # which we will fill with values of shared barcode fractions
    arr = np.zeros( [len(GEMbase.keys()), len(GEMbase.keys())] )
    windows = sorted(GEMbase.keys())
    GEMcomparison = pd.DataFrame(arr, columns=windows,index=windows) # Collect in a pandas dataframe for easy indexing
    done = 0
    total_windows = len(GEMbase.keys())

    while len(GEMbase.keys()) > 1:

        # Report progress every 20 windows
        if done in range(0,100000000,20):
            print(reportProgress(done, total_windows))

        # Pop out an entry in the dict
        reg = GEMbase.popitem()
        lib1 = reg[1]

        # And compare this entry to all other entries, successively filling the dataframe GEMcomparison
        # at both data points
        for k,v in GEMbase.items():
            lib2 = v
            shared_frac = compareGEMlibs_str(lib1,lib2)
            GEMcomparison[reg[0]][k] = shared_frac
            GEMcomparison[k][reg[0]] = shared_frac

        done += 1

    return GEMcomparison

def evaluateOutliers(win,outl):
    '''
    Looks through the outliers in the dict "outl" after suspicious hits,
    compared to the window
    '''
    # For now, consider the 7 closest windows to either side to be expected
    # This will depend on window size and mean length of DNA molecules that
    # went into the Chromium sequencing
    tig = win[0]
    start = win[1]
    stop = win[2]

    # For now, skip the windows close to contig ends
    if start < args.expected_molecule_size \
    or stop > ref_lengths[tig] - args.expected_molecule_size:
        return []

    # Calculate coordinates of the windows that are expected to surround the given window
    expected_windows = []

    for n in range(1,7):
        expected_windows.append( (tig, start + args.stepsize * n, start + args.stepsize * n + args.windowsize) )
        expected_windows.append( (tig, start - args.stepsize * n, start - args.stepsize * n + args.windowsize) )

    suspicious = []
    for o in outl:
        # If outlier is from a different contig, a possible misassembly is detected
        if o[0] != tig:
            suspicious.append(o)

        # If contig is the same, but the window number is not close, something
        # strange is also going on
        elif o not in expected_windows:
            suspicious.append(o)

        # Otherwise this outlier passed

    return suspicious

def main():
    """
    Steps:
    1. Read input files, split genome into windows, etc
    2. Collect the barcodes from every window
    3. Compare every window to all other windows
    4. Collect outliers, that is, regions if significantly enriched barcode fractions
    5. Check for suspicious windows
    6. Write output
    """
    # First step
    print("Starting chromqc pipeline")
    global samfile, ref_lengths
    samfile = pysam.AlignmentFile(args.input_bam, "rb")
    ref_lengths = collectRefLen(samfile)
    windows = createWindows(ref_lengths)
    outfilename = getOut()

    # 2nd step
    print("\nStarting barcode collection")
    GEMbase = createBarcodeDatabase(windows)
    samfile.close()

    # 3rd step
    print("Done\n\nStarting pairwise comparisons")
    GEMcomparison = pairwiseComps(GEMbase) # GEMcomparison will be a pandas dataframe

    # 4th and 5th steps
    print("Done\n\nStarting outlier detection")
    colnames = list(GEMcomparison.columns)

    output = []
    for col in GEMcomparison:
        fracs = np.asarray(GEMcomparison[col])
        outliers = calcESD.getOutliers_QC(fracs,colnames,10)
        outliers_short = {}

        # Remove lower outliers and 0's
        for k,v in outliers.items():
            if v > np.mean(fracs):
                outliers_short[k] = v

        suspects = evaluateOutliers(col, outliers_short)

        if suspects != []:
            output.append("{}:{}-{} ({})".format(col[0],col[1],col[2], ref_lengths[col[0]]))
            for k,v in outliers_short.items():
                output.append("\t{}\t{}".format(str(k), str(v)))

    print("\nDone. Writing output to {0}.txt.\n".format(outfilename))
    with open(outfilename+".misasm.txt","w",encoding = "utf-8") as out:
        out.write("\n".join(output))


if __name__ == "__main__":
    main()
