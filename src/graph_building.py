#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""graph_building

graph_building describes the Linkgraph class and implements
functions related to graph building during the ARBitR pipeline.

Copyright (c) 2020, Markus Hiltunen
Licensed under the GPL3 license. See LICENSE file.
"""

import mappy as mp
import numpy as np
import pandas as pd

import nuclseqTools as nt
import misc

class Junction:
    '''This class describes an ARBitR "junction". This is a part of a Linkgraph path
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
        '''Constructor for the Linkgraph class. Nodes are ends of the input contigs,
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
        '''Adds node to graph
        '''
        self.nodes.append(node)

    def addEdge(edge):
        '''Adds edge to graph. Edges have the format (node1, node2, weight),
        where nodes are directed, i.e. node1 is contig1s and node2 contig2s.
        This means that (node1, node2, weight) is the same as (node2, node1, weight).
        '''
        reverse_edge = (edge[1], edge[0], edge[2])
        self.edges.append(edge)
        self.edges.append(reverse_edge)

    def nodes():
        '''Returns all nodes in graph
        '''
        return self.nodes

    def edges():
        '''Returns all edges in graph
        '''
        return self.edges

    def opposite(self,node):
        '''Returns the opposite end of a given node.
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
        '''Returns all outgoing edges from the given node.
        '''
        outgoing_edges = []
        for e in self.edges:
            if e[0] == node:
                outgoing_edges.append(e)
        return outgoing_edges

    def getNodes(self,node):
        '''Returns:
            List: connected_nodes, all connected nodes to the given node.
        '''
        connected_nodes = []
        for e in self.edges:
            if e[0] == node:
                connected_nodes.append(e[1])
        return connected_nodes

    def extend(self,node):
        '''Returns a connected node to the given node, if there is
        an unambiguous connection. The connected node must have
        no other incoming edges.
        There are some possibilities:
        1. One outgoing edge.
        2. Two outgoing edges to s and e of the same contig.
        3. Three outgoing edges, two to s and e of the same contig and
        one to an additional node.

        Returns:
            Junction: If a single good connection is found, returns a junction
                describing this connection.
            None: if no good connection is found.
        '''

        # Edges are ordered, so we move from left to right along the edge
        connections = self.getNodes(node)

        # Check alternatives with one outgoing edge
        if len(connections) == 1:
            connected_node = connections[0]
            # Check the connected node for a single edge back to the original node
            if len(self.getNodes(connected_node)) == 1:
                return Junction(node, connected_node, [], GEMs[node]|GEMs[connected_node])
            # If the connected node has more than one outgoing edge
            elif len(self.getNodes(connected_node)) == 2:
                # We can get a junction in two cases from here:
                # 1. both new connections lead to connected_node, one on each
                # side.
                if self.getNodes(connected_node)[0][:-1] == self.getNodes(connected_node)[1][:-1]:
                    # Here we need also check if there is any other node involved.
                    # If not, we cannot confidently orient node.
                    if self.getNodes(self.opposite(node)) == 1:
                        return Junction(node, None, [connected_node[:-1]], GEMs[node]|GEMs[connected_node])
                    else:
                        return Junction(node, connected_node, [], GEMs[node]|GEMs[connected_node])
                # 2. one new connection leads elsewhere (third_node). Opposite end of connected_node
                # also leads here, and has only one outgoing edge.
                third_node = self.getNodes(connected_node)
                third_node.remove(node)
                third_node = third_node[0]
                # Next look at outgoing edges from opposite end of connections[0].
                opp_outg = self.getNodes(self.opposite(connected_node))
                # Test for 2.
                if len(opp_outg) == 1 and opp_outg[0] == third_node:
                    return Junction(node, connected_node, [], GEMs[node]|GEMs[connected_node])

        # Check alternatives with two outgoing edges
        elif len(connections) == 2:
            if connections[0][:-1] == connections[1][:-1]:
                # If to the opposite ends of the same contig without any connections,
                # return a Junction with None as the target node
                if len(self.getNodes(connections[0])) == 1 \
                and len(self.getNodes(connections[1])) == 1:
                    return Junction(node, None, [connections[0][:-1]], GEMs[node])

                # If there are connections everywhere
                elif len(self.getNodes(connections[0])) == 2 \
                and len(self.getNodes(connections[1])) == 2:
                    # If there are two outgoing edges from both,
                    # check how many outgoing edges they have
                    middle_ctg = connections[0][:-1]
                    conn1_nodes = self.getNodes(connections[0])
                    conn2_nodes = self.getNodes(connections[1])

                    if len(conn1_nodes) == 2 and set(conn1_nodes) == set(conn2_nodes):
                        # If both sides of the middle contig have the same outgoing
                        # edges, check that the connected node only leads back
                        conn1_nodes.remove(node)
                        third_node_connections = self.getNodes(conn1_nodes[0])
                        if len(third_node_connections) == 2:
                            return Junction(node, conn1_nodes[0], [middle_ctg], GEMs[node]|GEMs[conn1_nodes[0]])

                # Connections are uneven, but a path may still be found
                elif len(self.getNodes(connections[0])) == 1 \
                and (len(self.getNodes(connections[1])) > 1 \
                and node in set(self.getNodes(connections[1]))):
                    return Junction(node, connections[0], [], GEMs[node]|GEMs[connections[0]])
                elif len(self.getNodes(connections[1])) == 1 \
                and (len(self.getNodes(connections[0])) > 1 \
                and node in self.getNodes(connections[0])):
                    return Junction(node, connections[1], [], GEMs[node]|GEMs[connections[1]])

            # One connection may be back to a formerly connected node in the path
            elif len(self.getNodes(self.opposite(node))) == 1 and \
            self.getNodes(self.opposite(node))[0] in connections:
                # Remove this node from the connections
                connections.remove(self.getNodes(self.opposite(node))[0])
                if len(self.getNodes(connections[0])) == 1:
                    return Junction(node, connections[0], [], GEMs[node]|GEMs[connections[0]])

        # A third possibility exists: there are 3
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

    def findPath(self,node):
        '''From the given node, traverses the graph in both directions until
        it becomes impossible to continue.

        Description:
            Starting from an arbitrarily selected node, look for extensions
            first to the right and add these junctions to a path. Stop when
            a fork is encountered, or when a previously visited node
            is found. Then do the same to the left.

        Args:
            node (str): node to start from.

        Returns:
            list: path, crossed junctions in a list.
        '''
        starting_node = node
        path = [] # path to return, elements are Junctions
        visited = set() # keep track of visited contigs, elements are strings
        ext = self.extend(node) # Look for a possible extension

        while ext != None and ext.target and ext.target[:-1] not in visited:
            visited.add(node[:-1])
            path.append(ext)
            node = self.opposite(ext.target)
            ext = self.extend(node)

        # If loop was broken because there was no target, we still want to put
        # the junction in the path as there might be an interesting connection
        if ext and not ext.target:
            path.append(ext)

         # Need to add the last visited node to avoid adding too many junctions
         # to circular paths. It shouldn't affect linear paths
        visited.add(node[:-1])

        # At this point we have reached the end of the right extension
        # Now do the same in the other direction
        node = self.opposite(starting_node)
        ext = self.extend(node) # Look for a possible extension
        while ext != None and ext.target and ext.target[:-1] not in visited:
            visited.add(node[:-1])
            path.insert(0, Junction(ext.target,ext.start,ext.connections, GEMs[ext.target]|GEMs[ext.start]))
            node = self.opposite(ext.target)
            ext = self.extend(node)

        if ext and not ext.target:
            path.insert(0, Junction(ext.target,ext.start,ext.connections, GEMs[ext.start]))

        return path

    def unambiguousPaths(self):
        '''
        Append all unambiguous paths in the graph.
        '''
        visited_nodes = []
        for n in self.nodes:
            if n not in visited_nodes:
                path = self.findPath(n)
                if path != []:
                    self.paths.append(path)
                    # Append visisted nodes
                    for junc in path:
                        visited_nodes.append(junc.start)
                        visited_nodes.append(junc.target)
                        for conn in junc.connections:
                            visited_nodes.append(conn+"s")
                            visited_nodes.append(conn+"e")

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
    totallength = len(lib1|lib2) # Total number of barcodes

    # Find the fraction of shared barcodes, avoid division by 0
    return len(shared) / totallength if totallength != 0 else 0

def pairwise_comparisons(GEMlist):
    '''
    Performs all pairwise comparisons between windows in GEMlist.

    Returns:
        GEMcomparison (pd.DataFrame)
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

def calcOutliers(frac_series, factor = 3):
    '''Calculate the outliers of fracs. Default: major outlier.
    Factors > 10 usually give the best results.
    '''
    # Calculate quartiles and interquartile range
    quartiles = np.quantile(frac_series, [0.25,0.75])
    IQR = quartiles[1] - quartiles[0]
    upper_bound = quartiles[1] + (factor * IQR)
    return frac_series[frac_series > upper_bound]

def makeEdges(GEMcomparison, barcode_factor, min_barcode_fraction):
    '''Create edges from the GEMcomparison dataframe.

    Args:
        GEMcomparison (pd.DataFrame): All-against-all comparison of the
            windows' barcodes.
        barcode_factor (int): Factor for calculating outliers.
        min_barcode_fraction (float): Minimum fraction of shared barcodes to create
            an edge in the linkgraph.
    Returns:
        list: Edges inferred from the fractions of shared barcodes.
    '''

    misc.printstatus("Number of windows: "+str(len(GEMcomparison.keys())))
    edges = []

    with open("fractions.txt", "w") as out:
        for f in GEMcomparison.index:
            out.write("{}\t".format(f))
        out.write("\n")

        # Iterate over rows in GEMcomparison
        for idx, (region, fractions) in enumerate(GEMcomparison.iterrows()):
            contig = region[:-1]
            window = region[-1]

            out.write(region+"\t")
            for f in fractions:
                out.write("{}\t".format(f))
            out.write("\n")

            # Report progress every 100 windows
            if idx in range(0,10000000,100):
                misc.printstatusFlush("[ BARCODE LINKING ]\t" + misc.reportProgress(idx, len(GEMcomparison)))

            '''
            # Calculate outliers from the comparisons of window k to all other windows
            # outliers is a dict where each key is a connected window to region,
            # and value is the fraction of shared barcodes between region and window
            outliers = esd.getOutliers_QC(np.array(fractions),fractions.index,10)
            # Get rid of edges to the same contig.
            outliers = { k:v for k,v in outliers.items() if k[:-1] != region[:-1] \
                        and v > np.mean(fractions)}
            outliers = pd.Series(outliers)
            # If there are any outliers, i.e. edges to create, add them to the edges
            # list. Don't add edges for lower outliers (fractions < mean(fractions))
            # or where the fraction is less than
            # min_barcode_fraction (-f) and edges back to the same contig
            if len(outliers.keys()) > 1:
                sorted_outliers = outliers.sort_values(ascending = False)
                if sorted_outliers[0] > sorted_outliers[1] * barcode_factor:
                    outliers = outliers[outliers == sorted_outliers[0]]

            new_edges = [(region, connected_window, fraction) \
                        for connected_window, fraction in outliers.items()]

            # Let's try only writing single edges
            #if len(new_edges) == 1:
            for idx, mo in outliers.iteritems():
                edges.append( (region, idx, mo ) )

            '''

            # Ignore comparisons to the same contig and calculate outliers
            # In low coverage datasets the amount of 0's might cloud any
            # actual signal
            fractions = fractions.drop(labels = [contig+"s", contig+"e"], errors = "ignore")
            fractions = fractions[fractions > 0]
            if len(fractions) > 0:
                minor_outliers = calcOutliers(fractions, barcode_factor)
                minor_outliers = minor_outliers[minor_outliers > min_barcode_fraction]

                for ix, mo in minor_outliers.iteritems():
                    edges.append( (region, ix, mo ) )

        misc.printstatus("[ BARCODE LINKING ]\t" + misc.reportProgress(len(GEMcomparison), len(GEMcomparison)))

        return edges

def main(contig_lengths, GEMlist, barcode_factor, barcode_fraction):
    '''Controller for graph_building.

    Args:
        contig_lengths (dict): Contig names (keys) and their lengths (values)
            from the input bam file.
        GEMlist (dict): Windows corresponding to start and end regions
            (size determined by -s) (keys) and the barcodes collected from
            these regions (values).
        barcode_factor (int): Minimum fold difference between barcode fractions
            to create link.
        barcode_fraction (float): Minimum fraction of shared barcodes to create
            an edge in the linkgraph.
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
    edges = makeEdges(GEMcomparison, barcode_factor, barcode_fraction)
    graph = Linkgraph(nodes, edges)

    return graph
