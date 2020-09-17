#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""convert_fastq.py

Author: Markus Hiltunen
E-mail: markus.hiltunen@ebc.uu.se

Description: Converts an stLFR fastq file, produced from stLFR_read_demux,
into a format compatible with ARBitR.

Copyright (c) 2020, Markus Hiltunen
Licensed under the GPL3 license. See LICENSE file.
"""

import argparse

parser = argparse.ArgumentParser(description="Converts an stLFR fastq file, \
                                produced from stLFR_read_demux, \
                                into a format compatible with ARBitR. The \
                                converted fastq is printed in standard out.")
parser.add_argument("input_fastq", \
                    help="Input fastq file. Required.", \
                    type = str)
args = parser.parse_args()

def main():
    with open(args.input_fastq, "r") as fq:
        for idx, line in enumerate(fq):
            line = line.strip()
            if line.startswith("@") and idx % 4 == 0: # Change every 4th line
                barcode = line.split("#")[1].split("/")[0] # Grab the right substring
                print(line.split("\t")[0] + " BX:Z:" + barcode)
            else:
                print(line)

if __name__ == "__main__":
    main()
