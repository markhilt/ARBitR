#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""misc

misc implements various functions used during the ARBitR pipeline.

Copyright (c) 2020, Markus Hiltunen
Licensed under the GPL3 license. See LICENSE file.
"""

import time

def printstatus(msg):
    '''
    Prints current step in the workflow.
    '''
    now = time.asctime( time.localtime(time.time()) )
    print("[{}]\t{}".format(now, msg))

def printstatusFlush(msg):
    '''
    Prints current step in the workflow.
    '''
    now = time.asctime( time.localtime(time.time()) )
    print("[{}]\t{}".format(now, msg), end="\r", flush=True)

def reportProgress(current,total):
    return "Completed: {0}% ({1} out of {2})".format( str(round( (current / total) * 100, 2)), current, total)


# Deprecated
def readPaths(paths):
    '''
    Reads externally provided paths file for scaffolding.

    Input file contains one path per line with contigs to create a scaffold
    from tab delimited and with the orientation (f/r) appended as the last character,
    e.g. contig1f
    '''

    with open(paths, "r") as pths:
        out_paths = {}
        counter = 1
        for line in pths:
            line = line.strip()
            path = line.split("\t")
            out_paths["scaffold_"+str(counter)] = path
            counter += 1
    return out_paths
