#!/usr/bin/python3.4
# -*- coding: utf-8 -*-

''' For creating a bclib from reads aligned to a specified contig and coordinates in a bam file.'''

import argparse
import pysam
import os

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input", help="Input bam file", type = str)
parser.add_argument("-r","--region", help="Region of interest in the format contig_start-end", type = str)
parser.add_argument("-o","--output", help="Output prefix", type = str)
args = parser.parse_args()

# Return GEM barcodes aligned to this region in a set
def getGEMs(reads):
    BC_list = []
    for read in reads:
        if read.has_tag('BX') == True and read.mapping_quality >= 30:
            BC = read.get_tag("BX")
            BC_list.append(BC)
    BC_list = set(BC_list)
    return BC_list

samfile = pysam.AlignmentFile(args.input, "rb")

contig = args.region.split("_")[0]
start = int(args.region.split("_")[1].split("-")[0])
end = int(args.region.split("_")[1].split("-")[1])
reads = samfile.fetch(contig, start, end)

GEMs = getGEMs(reads)

if args.output:
    outfile = args.output + ".bclib"
else:
    outfile = args.region + ".bclib"

with open(outfile, mode = "w") as out:
    for i in GEMs:
        out.write(i+"\n")
