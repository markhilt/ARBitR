#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
anvil.py
Version 0.1
Author: Markus Hiltunen
E-mail: markus.hiltunen@ebc.uu.se

Description: Controller of the AnVIL script collection. This script controls the
parameters and workflow of AnVIL.

LICENSING
"""

import argparse
import numpy as np
import pysam
from scipy.stats import t
import time
import os

# Included modules
import formatting_tools
#import validate_runner
import graph_building
import barcode_collection
import misc
import merge_fasta

parser = argparse.ArgumentParser(description="Reads a bam file, creates links \
                                            between contigs based on linked read \
                                            information, and outputs a .gfa.")
parser.add_argument("input_bam", \
                    help="Input bam file. Required.", \
                    type = str)
parser.add_argument("-i", "--input_fasta", \
                    help="Input fasta file for contig merging. Optional. \
                    If not specified, will only output linkage graph in \
                    .gfa and .tsv format.", \
                    type = str)
parser.add_argument("-s","--region_size", \
                    help="Size of region of contig start and end to collect \
                    barcodes from. [20000]", \
                    default = 20000, \
                    type = int)
parser.add_argument("-n","--barcode_number", \
                    help="Minimum number of shared barcodes to create link. [1]", \
                    default = 1, \
                    type = int)
parser.add_argument("-f","--barcode_fraction", \
                    help="Minimum fraction of shared barcodes to create link. [0.01]", \
                    default = 0.01, \
                    type = float)
parser.add_argument("-q","--mapq", \
                    help="Mapping quality cutoff value. [60]", \
                    default = 60, \
                    type = int)
parser.add_argument("-o","--output", \
                    help="Prefix for output files.", \
                    type = str)
args = parser.parse_args()

def getOut():
    '''
    Creates a prefix for output files.
    '''
    if args.output:
        outfilename = args.output
    elif "/" in args.input_bam:
        outfilename = args.input_bam.split("/")[-1].split(".bam")[0]
    else:
        outfilename = args.input_bam.split(".bam")[0]
    return outfilename

def writeGfa(outfilename, gfa_header, graph):
    '''
    Writes the graph in gfa format.
    '''
    gfa_header = "H\tVN:Z:AnVIL/link_graph\n" + "\n".join(gfa_header)
    with open(outfilename + ".gfa", "w", encoding = "utf-8") as gfa:
        gfa.write(gfa_header+"\n")
        gfa.write(formatting_tools.formatGFA(graph))

# Deprecated
def writeGfa_from_dict(gfa_header, graph):
    '''
    Writes the graph in gfa format.
    '''
    gfa_header = "H\tVN:Z:AnVIL/link_graph\n" + "\n".join(gfa_header)
    with open(outfilename + ".gfa", "w", encoding = "utf-8") as gfa:
        gfa.write(gfa_header+"\n")
        gfalist = formatting_tools.formatGFA(graph)
        for i in gfalist:
            gfa.write(i)

def writeFasta(outfilename, linked_scaffolds):
    with open(outfilename+".fasta","w",encoding = "utf-8") as fastaout:
        for k,v in linked_scaffolds.items():
            fastaout.write(">"+k+"\n")
            fastaout.write(v+"\n")

def writeTSV(gfa_header, graph):
    '''
    Writes the graph in tsv format, compatible with LINKS.
    TODO
    '''

def writePaths(outfilename, scaffolds):
    '''
    Writes the usable paths from the graph
    '''
    with open(outfilename+".paths.txt", "w") as pathsout:
        for k,v in scaffolds.items():
            pathsout.write("{}\t{}\n".format(k,v))

def main():
    misc.printstatus("Starting AnVIL.")
    outfilename = getOut() # Create a prefix for output files

    global input_contigs
    samfile = pysam.AlignmentFile(args.input_bam, "rb")
    input_contigs = samfile.references
    samfile.close()

    # First step is to do the initial check through of the bam file,
    # to find any region where there is a reduction in barcode continuity
    misc.printstatus("Searching input bam file for suspicious regions: {}".format(args.input_bam))
    #suspicious_regions = validate_runner.main(args.input_bam)

    # Second step is to collect the barcodes from the input bam file
    misc.printstatus("Collecting barcodes.")
    GEMlist, gfa_header = barcode_collection.main(  args.input_bam, \
                                                    args.region_size, \
                                                    args.mapq)
    # Also add the suspicious region's barcodes
    # TODO

    # Third step is to build the link graph based on the barcodes
    misc.printstatus("Creating link graph.")
    graph, paths = graph_building.main(input_contigs, \
                                    GEMlist, \
                                    args.barcode_number, \
                                    args.barcode_fraction)

    misc.printstatus("Writing graph to {}.gfa.".format(outfilename))
    writeGfa(outfilename, gfa_header, graph)
    misc.printstatus("Writing paths to {}.paths.txt.".format(outfilename))
    writePaths(outfilename, paths)

    if args.input_fasta:
        if os.path.isfile(args.input_fasta):
            # If user gave an assembly fasta file, use this for merging
            misc.printstatus("Found fasta file for merging: {}".format(args.input_fasta))
            new_scaffolds = merge_fasta.main(args.input_fasta, args.input_bam, paths)
            misc.printstatus("Writing merged fasta to {0}.fasta".format(outfilename))
            writeFasta(outfilename,new_scaffolds)
        else:
            raise Exception("Fasta file not found: ".format(args.input_fasta))

    else:
        # else finish
        misc.printstatus("No fasta file specified for merging. Pipeline finished.")

    if os.path.isfile("tmp.fasta"):
        os.remove("tmp.fasta")
    misc.printstatus("AnVIL successfully completed!\n")

if __name__ == "__main__":
    main()
