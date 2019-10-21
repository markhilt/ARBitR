#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
.. module:: graph_building
    :synopsis: graph_building describes the Linkgraph class and implements
    functions related to graph building during the anvil pipeline.

Copyright (c) 2019, Johannesson lab
Licensed under the GPL3 license. See LICENSE file.
"""

import numpy as np

import graph_building_olc
import calcESD as esd
import misc

# Deprecated
def evaluateConnections(nodes, backbone_contigs):
    '''
    Check that given list of nodes fulfill some requirements.
    '''
    # Both nodes on all contigs need to be connected. Other nodes we don't
    # care about. They also need to be absent in the backbone graph.
    nodes = [n[:-1] for n in nodes if n[:-1] not in backbone_contigs]
    duplicates = set([x for x in nodes if nodes.count(x) > 1 ])
    return list(duplicates)


def fillJunctions(backbone_graph, GEMlist):
    '''
    Args:
        backbone_graph (Linkgraph)
        GEMlist (dict)
    Returns:
        list: list of paths with junctions filled.
    '''

    filled_junction_paths = []
    """
    backbone_contigs = list(set( [junction.start[:-1] for path in \
                            backbone_graph.paths for junction in \
                            path if junction.start != None] + \
                            [junction.target[:-1] for path in \
                            backbone_graph.paths for junction in \
                            path if junction.target != None] ))
    """

    # Iterate over paths and every junction in the path
    # Create a barcode comparison of the junction and all small contigs
    # Use ESD to determine which small contigs have a significant amount
    # of shared barcodes - put these into the junction
    for idx, path in enumerate(backbone_graph.paths):
        # Report progress every 100 windows
        if idx in range(0,10000000,1):
            misc.printstatusFlush("[ PATH FILLING ]\t" + misc.reportProgress(idx, len(backbone_graph.paths)))

        filled_path = []

        # Check outgoing edges from both start and target in full_graph.
        # If they are connected to both sides, add them to junction.
        for junction in path:
            tigs, fractions = zip(*[(k, graph_building_olc.compareGEMlibs(junction.barcodes, v)) for k,v in GEMlist.items()])

            """
            # Mask all 0's
            #fractions = np.ma.masked_equal(fractions,0)
            fractions = np.array(fractions)
            tigs = [i for j, i in enumerate(tigs) if j not in fractions[fractions != 0]]
            fractions = fractions[fractions != 0]
            """

            outliers = esd.getOutliers_QC(np.array(fractions),tigs,10)

            #import ipdb; ipdb.set_trace()

            # Add any outliers to junction.connections
            filled_path.append( graph_building_olc.Junction(junction.start, junction.target, junction.connections + [ o[:-1] for o in list(outliers.keys())] ))

            """
            connections =   list(set(junction.connections +
                            evaluateConnections(full_graph.getNodes(junction.start), backbone_contigs) + \
                            evaluateConnections(full_graph.getNodes(junction.target), backbone_contigs)))
            filled_path.append( graph_building_olc.Junction(junction.start, junction.target, connections) )
            """

        filled_junction_paths.append(filled_path)


    return filled_junction_paths
