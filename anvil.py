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
import graph_building
import barcode_collection
import misc
import merge_fasta
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
parser.add_argument("-F","--barcode_fraction", \
                    help="Minimum fraction of shared barcodes to create a link. \
                    [0.01]", \
                    default = 0.01, \
                    type = float)
parser.add_argument("-f","--barcode_factor", \
                    help="Factor to determine outliers. [39]", \
                    default = 39, \
                    type = int)
parser.add_argument("-q","--mapq", \
                    help="Mapping quality cutoff value for linkgraph. [60]", \
                    default = 60, \
                    type = int)
parser.add_argument("-Q","--short_mapq", \
                    help="Mapping quality cutoff value for pulling in short contigs. [20]", \
                    default = 20, \
                    type = int)
parser.add_argument("-c","--coverage", \
                    help="Coverage cutoff for trimming contig ends. [20]", \
                    default = 20, \
                    type = int)
parser.add_argument("-g","--gapsize", \
                    help="Gapsize for building scaffolds. [100]", \
                    default = 100, \
                    type = int)
parser.add_argument("-b","--bc_quantity", \
                    help="Cutoff for number of reads per barcode. [3]", \
                    default = 3, \
                    type = int)
parser.add_argument("-o","--output", \
                    help="Prefix for output files.", \
                    type = str)
args = parser.parse_args()

def getOut():
    '''Creates a prefix for output files.
    '''
    if args.output:
        outfilename = args.output
    elif "/" in args.input_bam:
        outfilename = args.input_bam.split("/")[-1].split(".bam")[0]+".anvil"
    else:
        outfilename = args.input_bam.split(".bam")[0]+".anvil"
    return outfilename

def writeGfa(outfilename, input_contig_lengths, graph):
    '''Writes the link graph in gfa format.
    '''
    with open(outfilename + ".gfa", "w", encoding = "utf-8") as gfa:
        gfa.write("H\tVN:Z:AnVIL/link_graph\n")
        for k,v in input_contig_lengths.items():
            gfa.write("S\t{0}\t*\tLN:i:{1}\n".format(k,str(v)))
        gfa.write(formatting_tools.formatGFA(graph))

def writeFasta(outfilename, linked_scaffolds):
    '''Writes fasta output.
    '''
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

def writeBed(outfilename, bed_dict):
    ''' Write a bed file from dict.
    '''
    with open(outfilename+".bed", "w") as out:
        for k,v in bed_dict.items():
            for feature in v:
                out.write(k+"\t"+"\t".join(feature)+"\n")

def main():
    misc.printstatus("Starting AnVIL.")

    # Unpack arguments
    region_size = args.region_size
    molecule_size = args.molecule_size
    mapq = args.mapq
    short_mapq = args.short_mapq
    barcode_factor = args.barcode_factor
    barcode_fraction = args.barcode_fraction
    mincov = args.coverage
    bc_quantity = args.bc_quantity
    gapsize = args.gapsize

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
    GEMlist = barcode_collection.main(  args.input_bam, \
                                            backbone_contig_lengths, \
                                            region_size, \
                                            mapq, \
                                            bc_quantity)

    # Second step is to build the link graph based on the barcodes
    misc.printstatus("Creating link graph.")
    backbone_graph = graph_building.main(backbone_contig_lengths, \
                                            GEMlist, \
                                            barcode_factor, \
                                            barcode_fraction)

    misc.printstatus("Writing link graph to {}.backbone.gfa.".format(outfilename))
    writeGfa(outfilename+".backbone", backbone_contig_lengths, backbone_graph)

    # Third step is to traverse the graph and build paths
    misc.printstatus("Finding paths.")
    backbone_graph.unambiguousPaths() # Fill graph.paths
    misc.printstatus("Found {} paths.".format(len(backbone_graph.paths)))
    writePaths(outfilename+".pre-fill", {str(idx):path for idx, path in enumerate(backbone_graph.paths)})

    # Fourth step is to collect the barcodes from the input bam file,
    # this time for the small contigs
    misc.printstatus("Collecting barcodes from short contigs.")
    GEMlist = barcode_collection.main(  args.input_bam, \
                                            small_contig_lengths, \
                                            molecule_size, \
                                            short_mapq)

    # Fifth step is to pull in the short contigs into the linkgraph junctions,
    # if they have
    # Sixth step is to fill the junctions in the backbone_graph
    paths = fill_junctions.fillJunctions(backbone_graph, GEMlist)

    writePaths(outfilename+".pre-merge", {str(idx):path for idx, path in enumerate(paths)})

    if os.path.isfile(args.input_fasta):
        # If user gave an assembly fasta file, use this for merging
        misc.printstatus("Found fasta file for merging: {}".format(args.input_fasta))
        new_scaffolds, scaffold_correspondence, bed = merge_fasta.main(args.input_fasta, args.input_bam, paths, mincov, gapsize)
        misc.printstatus("Writing merged fasta to {0}.fasta".format(outfilename))
        writeFasta(outfilename,new_scaffolds)
        writePaths(outfilename+".correspondence", scaffold_correspondence)
        writeBed(outfilename, bed)

    else:
        misc.printstatus("No fasta file found for merging. Pipeline finished.")

    if os.path.isfile("tmp.fasta"):
        os.remove("tmp.fasta")

    misc.printstatus("AnVIL successfully completed!\n")

if __name__ == "__main__":
    main()
