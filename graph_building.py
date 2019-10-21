#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
.. module:: graph_building
    :synopsis: graph_building describes the Linkgraph class and implements
    functions related to graph building during the anvil pipeline.

Copyright (c) 2019, Johannesson lab
Licensed under the GPL3 license. See LICENSE file.
"""

import mappy as mp
import numpy as np
import pandas as pd

import nuclseqTools as nt
import calcESD as esd
import misc

class Junction:
    '''
    This class describes an AnVIL "junction". This is a part of a Linkgraph path
    where there is a known start and target node and their orientations are known.
    Then there can be one or more connected nodes of unknown orientation.
    '''
    def __init__(self, start_node, target_node, connected = [], barcodes = {}):
        self.start = start_node # Name of starting contig
        self.target = target_node # Name of target contig
        self.connections = connected # List of names of connected contigs
        self.barcodes = barcodes # Set of barcodes

    def __repr__(self):
        return "start: {0}, target: {1}, connections: {2}".format(self.start, self.target, self.connections)

    def __str__(self):
        return "Junction([{0},{1},{2}])".format(self.start, self.target, self.connections)

    def start():
        '''Returns starting node and orientation.
        '''
        return self.start

    def target():
        '''Returns target node and orientation.
        '''
        return self.target

    def connections():
        '''Return connections.
        '''
        return self.connections

    def barcodes():
        '''Returns barcodes.
        '''
        return self.barcodes

    def addConnection(node):
        '''Adds node to the list of connections.
        '''
        self.connections.append(node)

class Linkgraph:
    def __init__(self, nodes = [], edges = []):
        '''
        Constructor for the Linkgraph class. Nodes are ends of the input contigs,
        e.g. contig1s means the start of a contig named contig1. For now we assume
        that contigs are correctly assembled, which means that there is always
        a single edge between contig1s and contig1e. Drawing this edge in the graph
        is thus redundant. Instead we use the function opposite to traverse this
        edge.

        '''
        self.nodes = nodes  # [node1, node2, ...]

        # Draw edges in both directions
        # [(node1, node2, set(node1 barcodes), set(node1 barcodes)), ...] i.e. node1->node2 etc
        self.edges = list(set(edges + [(a[1], a[0], a[2]) for a in edges]))
        self.paths = [] # Paths as lists of junctions

    def __str__(self):
        return "Linkgraph([nodes: {0},edges: {1}])".format(len(self.nodes), len(self.edges))

    def addNode(node):
        '''
        Adds node to graph
        '''
        self.nodes.append(node)

    def addEdge(edge):
        '''
        Adds edge to graph. Edges have the format (node1, node2, weight),
        where nodes are directed, i.e. node1 is contig1s and node2 contig2s.
        This means that (node1, node2, weight) is the same as (node2, node1, weight).
        '''
        reverse_edge = (edge[1], edge[0], edge[2])
        self.edges.append(edge)
        self.edges.append(reverse_edge)

    def nodes():
        '''
        Returns all nodes in graph
        '''
        return self.nodes

    def edges():
        '''
        Returns all edges in graph
        '''
        return self.edges

    def opposite(self,node):
        '''
        Finds and returns the opposite end of a given node of format tigs, tige or tigu
        '''
        tig = node[:-1]
        side = node[-1]
        if side == "s":
            opposite = "e"
        elif side == "e":
            opposite = "s"
        else:
            opposite = "u"
        return "".join([tig,opposite])

    def getEdges(self,node):
        '''
        Returns all outgoing edges from the given node.
        '''
        outgoing_edges = []
        for e in self.edges:
            if e[0] == node:
                outgoing_edges.append(e)
        return outgoing_edges

    def getNodes(self,node):
        '''
        Returns all connected nodes to the given node.
        '''
        connected_nodes = []
        for e in self.edges:
            if e[0] == node:
                connected_nodes.append(e[1])
        return connected_nodes

    def extend(self,node):
        '''
        Returns a connected node to the given node, if there is
        an unambiguous connection. The connected node must have
        no other incoming edges.
        There are three possibilities:
        1. One outgoing edge.
        2. Two outgoing edges to s and e of the same contig.
        3. Three outgoing edges, two to s and e of the same contig and
        one to an additional node.

        Returns:
            Junction: If a single good connection is found, returns a junction
                describing this.
            None: if no good connection is found.
        '''

        if node not in self.nodes:
            raise Exception("Linkgraph error: node {} not in graph".format(node))

        else:
            # Edges are ordered, so we move from left to right along the edge
            connections = self.getNodes(node)

            # If only one outgoing edge, we continue to check the connected node
            # for a single edge back to the original node
            if len(connections) == 1 and len(self.getNodes(connections[0])) == 1:
                return Junction(node, connections[0], [], GEMs[node]|GEMs[connections[0]])

            elif len(connections) == 2 and connections[0][:-1] == connections[1][:-1]:
                # If there are two outgoing edges to the same contig,
                # check if there is anywhere to continue to. If not, return
                # outgoing_edges
                if len(self.getNodes(connections[0])) == 1 \
                and len(self.getNodes(connections[1])) == 1:
                    return Junction(node,None, [connections[0][:-1]], GEMs[node])

                elif len(self.getNodes(connections[0])) == 2 \
                and len(self.getNodes(connections[1])) == 2:
                    # If there are two outgoing edges from both,
                    # check that both connected nodes have two edges,
                    # one back at each side, and one continuing on each side to another node
                    middle_ctg = connections[0][:-1]
                    conn1_nodes = self.getNodes(connections[0])
                    conn2_nodes = self.getNodes(connections[1])

                    if len(conn1_nodes) == 2 and set(conn1_nodes) == set(conn2_nodes):
                        conn1_nodes.remove(node)
                        third_node_connections = self.getNodes(conn1_nodes[0])
                        if len(third_node_connections) == 2:
                            return Junction(node, conn1_nodes[0], [middle_ctg], GEMs[node]|GEMs[conn1_nodes[0]])

            # A third possibility exists: the same as above, but another edge also
            # exists between the starting and target nodes. There are thus 3
            # connections, out of which two are to the same contig (middle).
            elif len(connections) == 3 and len(set( [c[:-1] for c in connections] )) == 2:

                # Find which two are to opposite ends of the same contig
                if connections[0][:-1] == connections[1][:-1]:
                    middle_ctg = connections[0][:-1]
                    middle_ctg_node1 = connections[0]
                    middle_ctg_node2 = connections[1]
                    target_node = connections[2]
                elif connections[0][:-1] == connections[2][:-1]:
                    middle_ctg = connections[0][:-1]
                    middle_ctg_node1 = connections[0]
                    middle_ctg_node2 = connections[2]
                    target_node = connections[1]
                elif connections[1][:-1] == connections[2][:-1]:
                    middle_ctg = connections[1][:-1]
                    middle_ctg_node1 = connections[1]
                    middle_ctg_node2 = connections[2]
                    target_node = connections[0]

                # collect the connections from the connected nodes
                middle_ctg_node1_conn = self.getNodes(middle_ctg_node1)
                middle_ctg_node2_conn = self.getNodes(middle_ctg_node2)
                target_node_conn = self.getNodes(target_node)

                # We need to check these for other outgoing edges...
                # Only accept edges back to these same nodes
                if len(middle_ctg_node1_conn) == 2 \
                and len(middle_ctg_node2_conn) == 2 \
                and len(target_node_conn) == 3:

                    if node in middle_ctg_node1_conn \
                    and target_node in middle_ctg_node1_conn \
                    and node in middle_ctg_node2_conn \
                    and target_node in middle_ctg_node2_conn \
                    and node in target_node_conn \
                    and middle_ctg_node1 in target_node_conn \
                    and middle_ctg_node2 in target_node_conn:
                        return Junction(node, target_node, [middle_ctg], GEMs[node]|GEMs[target_node])

        return None

    def findPath(self,starting_node):
        '''
        From the given node, traverses the graph in both directions until
        it becomes impossible to continue
        Returns:
            list: crossed junctions in a list.
        '''
        visited = [starting_node] # Add starting node to list
        node = None
        path = []

        ext = self.extend(starting_node) #
        # Avoid extensions looping back to the same contig
        if ext != None and ext.target != None:
            if ext.target[:-1] != ext.start[:-1]:
                path.append(ext)
                node = ext.target

        # If there is a next node...
        while node != None:
            visited.append(node) # Add this node to list

            # Go to next node, which is always the opposite end of the newly
            # connected contig
            #ipdb.set_trace()
            node = self.opposite(node)
            visited.append(node) # Add this node to list
            ext = self.extend(node) # Go to next node

            if ext != None:
                node = ext.target
                path.append(ext)

            else:
                node = None

        # At this point we have reached the end of the right extension
        # Now do the same in the other direction
        node = None
        starting_node = self.opposite(starting_node)
        visited.insert(0, starting_node) # Insert at the beginning of the visited list
        ext = self.extend(starting_node) # Grab next nodes
        if ext != None and ext.target != None:
            if ext.target[:-1] != ext.start[:-1]:
                path.insert(0, Junction(ext.target,ext.start,ext.connections, GEMs[ext.target]|GEMs[ext.start]))
                node = ext.target

        while node != None:
            visited.insert(0, node) # Add this node to list

            # Go to next node, which is always the opposite end of the newly
            # connected contig
            node = self.opposite(node)
            visited.insert(0, node) # Add this node to list
            ext = self.extend(node) # Go to next node
            if ext != None:
                if ext.target != None:
                    node = ext.target
                    path.insert(0, Junction(ext.target,ext.start,ext.connections, GEMs[ext.target]|GEMs[ext.start]))
                else:
                    node = ext.target
                    path.insert(0, Junction(ext.target,ext.start,ext.connections, GEMs[ext.start]))
            else:
                node = None

        return path

    def unambiguousPaths(self):
        '''
        Append all unambiguous paths in the graph, i.e. paths with a single
        edge connecting every node in the path.
        '''
        visited_nodes = []
        for n in self.nodes:
            if n not in visited_nodes:
                path = self.findPath(n)
                if path != []:
                    self.paths.append(path)
                    visited_nodes = visited_nodes + \
                                    [ junc.start for junc in path ] + \
                                    [ junc.target for junc in path ]
                    # Also append starting and ending points
                    if path[0].start != None:
                        visited_nodes.append(self.opposite(path[0].start))
                    if path[-1].target != None:
                        visited_nodes.append(self.opposite(path[-1].target))

def makeNodes(contig_list):
    return [a+"s" for a in contig_list] + [a+"e" for a in contig_list]

def compareGEMlibs(lib1,lib2):
    '''
    Compares two sets of barcodes, collects all shared ones and
    counts them. Returns fraction of shared barcodes.
    '''
    # Find shared ones. Numpy arrays were tried here - way slower
    shared = lib1.intersection(lib2)
    totallength = len(lib1) + len(lib2) - len(shared) # Total number of unshared barcodes

    # Find the fraction of shared barcodes, avoid division by 0
    if totallength != 0:
        return len(shared) / totallength
    else:
        return 0

def pairwise_comparisons(GEMlist):
    '''
    Performs all pairwise comparisons between windows in GEMlist.
    '''
    # Compare the barcodes in every region to all other regions
    GEMcomparison = pd.DataFrame(np.zeros(( len(GEMlist), len(GEMlist) )), \
                                index=GEMlist.keys())
    GEMcomparison.columns = GEMcomparison.index

    # Iterate over rows in GEMcomparison
    # Index to keep track of position so we can skip calculating some fractions
    # twice
    idx = 0
    for idx, region1 in enumerate(GEMcomparison.index):
        lib1 = GEMlist[region1]

        # Report progress every 20 windows
        if idx in range(0,100000000,20):
            misc.printstatusFlush("[ BARCODE COMPARISON ]\t" + misc.reportProgress(idx+1, len(GEMlist)))

        fractions = [ compareGEMlibs(lib1,GEMlist[col]) for col in GEMcomparison.columns[idx:] ]
        GEMcomparison.loc[region1][idx:] = fractions # Update row values from idx
        GEMcomparison[region1][idx:] = fractions # Also update column values from idx

    misc.printstatus("[ BARCODE COMPARISON ]\t" + misc.reportProgress(idx+1, len(GEMlist)))

    return GEMcomparison

"""
def pairwise_comparisons(GEMlist):
    '''
    Performs all pairwise comparisons between windows in GEMlist.
    '''
    # Compare the barcodes in every region to all other regions
    GEMcomparison = {} # Dict of the structure { str(region): pd.Series(fractions, index = other_regions) }
    comp = {}

    # Iterate over GEMlist
    idx = 0

    import itertools
    for reg1, reg2 in itertools.combinations(GEMlist.items(), 2):
        if idx in range(0,100000000,2000):
            misc.printstatusFlush("[ BARCODE COMPARISON ]\t" + misc.reportProgress(idx+1, (len(GEMlist)*len(GEMlist))/2))
        idx += 1

        frac = compareGEMlibs(reg1[1],reg2[1])

        # Add both directions
        if reg1[0] in comp.keys():
            comp[reg1[0]].append((reg2[0], frac))
        else:
            comp[reg1[0]] = [(reg2[0], frac)]
        if reg2[0] in comp.keys():
            comp[reg2[0]].append((reg1[0], frac))
        else:
            comp[reg2[0]] = [(reg1[0], frac)]

    for k,v in comp.items():
        regions, fractions = zip(*[(i[0], i[1]) for i in v])
        GEMcomparison[k] = pd.Series(np.array(fractions), index = regions)

    misc.printstatus("[ BARCODE COMPARISON ]\t" + misc.reportProgress(idx+1, (len(GEMlist)*len(GEMlist))/2))

    return GEMcomparison


def compare(arg):
    reg1, reg2 = arg[0][0], arg[1][0]
    bc1, bc2 = arg[0][1], arg[1][1]
    frac = compareGEMlibs(bc1, bc2)
    GEMcomparison[reg1][reg2] = frac
    GEMcomparison[reg2][reg1] = frac

def pairwise_comparisons(GEMlist):
    '''
    Performs all pairwise comparisons between windows in GEMlist.
    '''
    # Compare the barcodes in every region to all other regions
    # Initiate matrix
    idx = list(GEMlist.keys())
    l = len(GEMlist)
    global GEMcomparison
    GEMcomparison = pd.DataFrame(np.zeros((len(idx), len(idx))), index = idx)
    GEMcomparison.columns = GEMcomparison.index

    # Iterate over GEMlist
    idx = 0

    # Create pool from columns of GEMcomparison
    import multiprocessing
    import itertools
    import tqdm

    combs = itertools.combinations(GEMlist.items(), 2)
    n_processes = 48
    #n_combs = len(list(combs)) do NOT use this: it wrecks the combs generator
    # for some reason
    # Calculate expected length of combs:
    import math
    n_combs = math.factorial(len(GEMlist.keys())) // math.factorial(2) // \
            math.factorial(len(GEMlist.keys())-2)


    #GEMcomparison = {k:{} for k in GEMlist.keys() } # Testing nested dict instead
    #GEMcomparison = {k:pd.Series(np.zeros(len(GEMlist.keys())), index = GEMlist.keys()) for k in GEMlist.keys() } # Dict of pd.Series
    print("starting mp")
    print("# comparisons to make: {}, split across {} processes".format(n_combs, n_processes))
    with multiprocessing.Pool(n_processes) as p:
        for _ in tqdm.tqdm(p.imap_unordered(compare, combs, chunksize = 50), total = n_combs):
            pass
    print("mp done")
    #179087275 16836

    return GEMcomparison
"""

def makeEdges(GEMcomparison, min_barcode_fraction):
    '''Create edges from the GEMcomparison dict.

    Args:
        GEMcomparison (dict): All-against-all comparison of the
            windows' barcodes.
        min_barcode_fraction (float): Minimum fraction of shared barcodes to
            draw an edge in the Linkgraph.
    Returns:
        list: Edges inferred from the fractions of shared barcodes.

    '''
    misc.printstatus("Number of windows: "+str(len(GEMcomparison.keys())))
    edges = []

    # Iterate over rows in GEMcomparison
    for idx, (region, fractions) in enumerate(GEMcomparison.items()):
        contig = region[:-1]
        window = region[-1]

        # Report progress every 100 windows
        if idx in range(0,10000000,100):
            misc.printstatusFlush("[ BARCODE LINKING ]\t" + misc.reportProgress(idx, len(GEMcomparison)))

        # Calculate outliers from the comparisons of window k to all other windows
        # outliers is a dict where each key is a connected window to region,
        # and value is the fraction is shared barcodes between region and window
        outliers = esd.getOutliers_QC(np.array(fractions),fractions.index,10)

        # If there are any outliers, i.e. edges to create, add them to the edges
        # list. Don't add edges for lower outliers (fractions < mean(fractions))
        # or where the fraction is less than
        # min_barcode_fraction (-f) and edges back to the same contig
        if len(outliers.keys()) > 0:
            edges = edges + [   (region, connected_window, fraction) \
                                for connected_window, fraction in outliers.items() \
                                if (connected_window[:-1] != region[:-1] ) \
                                and fraction > np.mean(fractions)]

    misc.printstatus("[ BARCODE LINKING ]\t" + misc.reportProgress(len(GEMcomparison), len(GEMcomparison)))

    return edges

"""
########### Deprecated
def makeEdges(GEMcomparison, min_barcode_fraction):
    '''Create edges from the GEMcomparison array.

    Args:
        GEMcomparison (pandas DataFrame): All-against-all comparison of the
            windows' barcodes. Columns = indices = window names. The dataframe
            contains the fraction of shared barcodes between every window.
        min_barcode_fraction (float): Minimum fraction of shared barcodes to
            draw an edge in the Linkgraph.
    Returns:
        list: Edges inferred from the fractions of shared barcodes.

    '''
    donewindows = 0
    misc.printstatus("Number of windows: "+str(len(GEMcomparison)))
    edges = []

    # Iterate over rows in GEMcomparison
    for region, fractions in GEMcomparison.iterrows():
        contig = region[:-1]
        window = region[-1]

        # Report progress every 100 windows
        if donewindows in range(0,10000000,100):
            misc.printstatusFlush("[ BARCODE LINKING ]\t" + misc.reportProgress(donewindows, len(GEMcomparison)))

        # Calculate outliers from the comparisons of window k to all other windows
        # outliers is a dict where each key is a connected window to region,
        # and value is the fraction is shared barcodes between region and window

        # Try to get rid of 0's,
        # could possibly speed up
        # FOr now in two steps
        columns = [ GEMcomparison.columns[i] for i, val in enumerate(fractions) if val != 0 ]
        fractions = list(filter(lambda a: a != 0, fractions))


        #outliers = esd.getOutliers_QC(np.array(fractions),GEMcomparison.columns,10)
        outliers = esd.getOutliers_QC(np.array(fractions),columns,10)

        # If there are any outliers, i.e. edges to create, add them to the edges
        # list. Don't add edges where the fraction is less than
        # min_barcode_fraction (-f) and edges back to the same contig
        if len(outliers.keys()) > 1:
            edges = edges + [   (region, connected_window, fraction) \
                                for connected_window, fraction in outliers.items() \
                                if (connected_window[:-1] != region[:-1] \
                                and fraction > min_barcode_fraction)]

        donewindows += 1
        #ipdb.set_trace()


    misc.printstatus("[ BARCODE LINKING ]\t" + misc.reportProgress(donewindows, len(GEMcomparison)))

    return edges
"""


def main(contig_lengths, GEMlist, min_barcode_number, min_barcode_fraction):
    '''Controller for graph_building.

    Args:
        contig_lengths (dict): Contig names (keys) and their lengths (values)
            from the input bam file.
        GEMlist (dict): Windows corresponding to start and end regions
            (size determined by -s) (keys) and the barcodes collected from
            these regions (values).
        min_barcode_number (int): Minimum number of shared barcodes to draw
            an edge in the Linkgraph.
        min_barcode_fraction (float): Minimum fraction of shared barcodes to
            draw an edge in the Linkgraph.
    Returns:
        Linkgraph: graph inferred from the input region barcodes.
        dict: dictionary of new scaffold names and which input contigs
            that should go into the new scaffold.
    '''

    global GEMs
    GEMs = GEMlist

    # Collect the fraction of shared barcodes in the all-against-all
    # comparison of windows
    GEMcomparison = pairwise_comparisons(GEMlist)

    # Infer linkage based on statistically significant outliers determined
    # by the ESD test to build the graph
    nodes = makeNodes(list(contig_lengths.keys()))
    edges = makeEdges(GEMcomparison, min_barcode_fraction)
    graph = Linkgraph(nodes, edges)

    return graph
