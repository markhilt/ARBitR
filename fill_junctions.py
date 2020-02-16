#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
.. module:: fill_junctions.py

Copyright (c) 2019, Johannesson lab
Licensed under the GPL3 license. See LICENSE file.
"""

import numpy as np
import pandas as pd

import graph_building
import misc

def fillJunctions(backbone_graph, GEMlist):
    '''
    Args:
        backbone_graph (Linkgraph)
        GEMlist (dict)
    Returns:
        list: list of paths with junctions filled.
    '''

    filled_junction_paths = []
    barcode_factor = 15 # hardcoded, for now

    # Iterate over paths and every junction in the path
    # Create a barcode comparison of the junction and all small contigs
    # Use ESD to determine which small contigs have a significant amount
    # of shared barcodes - put these into the junction
    for idx, path in enumerate(backbone_graph.paths):
        # Report progress every 100 windows
        if idx in range(0,10000000,1):
            misc.printstatusFlush("[ PATH FILLING ]\t" + misc.reportProgress(idx+1, len(backbone_graph.paths)))

        filled_path = []

        # Check outgoing edges from both start and target in full_graph.
        # If they are connected to both sides, add them to junction.
        for junction in path:
            tigs, fractions = zip(*[(k, graph_building.compareGEMlibs(junction.barcodes, v)) for k,v in GEMlist.items()])
            fracs = pd.Series(fractions, index = tigs)
            fracs = fracs[fracs > 0]
            outliers = graph_building.calcOutliers(fracs, barcode_factor)

            # Old outlier method:
            #outliers = esd.getOutliers_QC(np.array(fractions),tigs,10)

            # Add any outliers to junction.connections
            filled_path.append( graph_building.Junction(junction.start, junction.target, junction.connections + [ o[:-1] for o in list(outliers.index)] ))

        filled_junction_paths.append(filled_path)
    misc.printstatus("[ PATH FILLING ]\t" + misc.reportProgress(idx+1, len(backbone_graph.paths)))

    return filled_junction_paths
