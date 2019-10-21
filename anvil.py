#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
anvil.py
Version 0.1
Author: Markus Hiltunen
E-mail: markus.hiltunen@ebc.uu.se

Description: Controller of the AnVIL script collection. This script controls the
parameters and workflow of AnVIL.

Copyright (c) 2019, Johannesson lab
Licensed under the GPL3 license. See LICENSE file.
"""

import time
import os

import argparse
import numpy as np
import pysam
from scipy.stats import t

import formatting_tools
import graph_building_olc
import barcode_collection_olc
import misc
import merge_fasta_olc
import fill_junctions

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
parser.add_argument("-m","--molecule_size", \
                    help="Estimated mean molecule size that went into Chromium \
                    sequencing. Linked reads spanning a distance larger than \
                    this size should be rare. [45000]", \
                    default = 45000, \
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
        outfilename = args.input_bam.split("/")[-1].split(".bam")[0]+".anvil"
    else:
        outfilename = args.input_bam.split(".bam")[0]+".anvil"
    return outfilename

# Deprecated
def formatContigs(samfile):
    '''
    Creates a dict where keys are contig names and values their lengths
    '''
    i = 0
    gfa_header = []
    contig_dict = {}
    contig_list = samfile.references
    while i < len(contig_list):
        contig_dict[contig_list[i]] = samfile.lengths[i]
        gfa_header.append("S\t{0}\t*\tLN:i:{1}".format(contig_list[i],samfile.lengths[i]))
        i += 1

    return contig_dict, gfa_header

def writeGfa(outfilename, input_contig_lengths, graph):
    '''
    Writes the graph in gfa format.
    '''

    with open(outfilename + ".gfa", "w", encoding = "utf-8") as gfa:
        gfa.write("H\tVN:Z:AnVIL/link_graph\n")
        for k,v in input_contig_lengths.items():
            gfa.write("S\t{0}\t*\tLN:i:{1}\n".format(k,str(v)))
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

# Not used
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

    # Unpack arguments
    region_size = args.region_size
    molecule_size = args.molecule_size
    mapq = args.mapq
    barcode_number = args.barcode_number
    barcode_fraction = args.barcode_fraction

    if region_size > molecule_size:
        misc.printstatus("Larger --region_size than --molecule_size detected. Using default values instead.")
        region_size, molecule_size = 20000, 45000

    outfilename = getOut() # Create a prefix for output files
    samfile = pysam.AlignmentFile(args.input_bam, "rb")
    input_contig_lengths = dict(zip(samfile.references, samfile.lengths))
    samfile.close()

    # Split dataset into backbone and small contigs
    misc.printstatus("Collecting contigs.")
    backbone_contig_lengths = { ctg:length for ctg, length in input_contig_lengths.items() if length > molecule_size}
    small_contig_lengths = {k:input_contig_lengths[k] for k in input_contig_lengths.keys() - backbone_contig_lengths.keys()}

    # First step is to collect the barcodes for the backbone graph
    misc.printstatus("Collecting barcodes for linkgraph.")
    GEMlist = barcode_collection_olc.main(  args.input_bam, \
                                            backbone_contig_lengths, \
                                            region_size, \
                                            mapq)

    # Second step is to build the link graph based on the barcodes
    misc.printstatus("Creating link graph.")
    backbone_graph = graph_building_olc.main(backbone_contig_lengths, \
                                            GEMlist, \
                                            barcode_number, \
                                            barcode_fraction)

    misc.printstatus("Writing link graph to {}.backbone.gfa.".format(outfilename))
    writeGfa(outfilename+".backbone", backbone_contig_lengths, backbone_graph)

    # Third step is to traverse the graph and build paths
    misc.printstatus("Finding paths.")
    backbone_graph.unambiguousPaths() # Fill graph.paths

    # Fourth step is to collect the barcodes from the input bam file,
    # this time for the small contigs
    misc.printstatus("Collecting barcodes from short contigs.")
    GEMlist = barcode_collection_olc.main(  args.input_bam, \
                                            small_contig_lengths, \
                                            molecule_size, \
                                            mapq)

    # Fifth step is to pull in the short contigs into the linkgraph junctions,
    # if they have
    # Sixth step is to fill the junctions in the backbone_graph
    paths = fill_junctions.fillJunctions(backbone_graph, GEMlist)

    '''
    # Fifth step is to build the full link graph based on the barcodes
    misc.printstatus("Creating full link graph.")
    graph = graph_building_olc.main(input_contig_lengths, \
                                    GEMlist, \
                                    args.barcode_number, \
                                    args.barcode_fraction)

    misc.printstatus("Writing complete link graph to {}.gfa.".format(outfilename))
    writeGfa(outfilename, input_contig_lengths, graph)
    '''

    writePaths(outfilename+".pre-merge", {str(idx):path for idx, path in enumerate(paths)})

    if os.path.isfile(args.input_fasta):
        # If user gave an assembly fasta file, use this for merging
        new_scaffolds, scaffold_correspondence = merge_fasta_olc.main(args.input_fasta, args.input_bam, paths)
        misc.printstatus("Found fasta file for merging: {}".format(args.input_fasta))
        misc.printstatus("Writing merged fasta to {0}.fasta".format(outfilename))
        writeFasta(outfilename,new_scaffolds)
        writePaths(outfilename+".correspondence", scaffold_correspondence)

    else:
        misc.printstatus("No fasta file found for merging. Pipeline finished.")

    if os.path.isfile("tmp.fasta"):
        os.remove("tmp.fasta")

    misc.printstatus("AnVIL successfully completed!\n")

if __name__ == "__main__":
    main()
