import nuclseqTools as nt
import calcESD as esd
import misc
import mappy as mp

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
        # [(node1, node2, weight), (node1, node3, weight), ...] i.e. node1->node2 etc
        # Draw edges in both directions
        # self.edges = edges
        self.edges = list(set(edges + [(a[1], a[0], a[2]) for a in edges]))

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

    def findConnections(self,node):
        '''
        Returns all nodes connected to the given node.
        '''
        connected_nodes = []
        for e in self.edges:
            if e[0] == node:
                connected_nodes.append(e[1])
        return connected_nodes

    def getEdges(self,node):
        '''
        Returns all outgoing edges from the to the given node.
        '''
        outgoing_edges = []
        for e in self.edges:
            if e[0] == node:
                outgoing_edges.append(e)
        return outgoing_edges

    def bestOrientation(self, middle_tig, starting_node_outgoing_edges, target_node_outgoing_edges):
        '''
        Short contigs are often connected to both sides, i.e. both nodes
        of the contig have edges to two other nodes. In these cases,
        determine the orientation of the middle contig by the fraction weight
        to each connected node.
        '''
        path = []
        cutoff = 1.5
        #ipdb.set_trace()

        # We need to make four comparisons, between the right edges and the
        # left edges. Hopefully the weights support the same path.
        # Arbitrarily choose 50% difference as a cutoff, i.e. if one weight
        # is more than 50% higher, choose that one.
        starting_node = starting_node_outgoing_edges[0][0]
        target_node = target_node_outgoing_edges[0][0]

        starting_node_to_middle_tig_edges = \
        [a for a in starting_node_outgoing_edges if a[1][:-1] == middle_tig]
        target_node_to_middle_tig_edges = \
        [a for a in target_node_outgoing_edges if a[1][:-1] == middle_tig]

        left_edge1_node2 = starting_node_to_middle_tig_edges[0][1]
        left_edge1_weight = starting_node_to_middle_tig_edges[0][2]
        left_edge2_node2 = starting_node_to_middle_tig_edges[1][1]
        left_edge2_weight = starting_node_to_middle_tig_edges[1][2]

        right_edge1_node2 = target_node_to_middle_tig_edges[0][1]
        right_edge1_weight = target_node_to_middle_tig_edges[0][2]
        right_edge2_node2 = target_node_to_middle_tig_edges[1][1]
        right_edge2_weight = target_node_to_middle_tig_edges[1][2]

        # We have two possible choices for the path.
        # It always starts at starting_node (which we wont return as we already know this)
        # and ends at target_node. Question is the orientation of the middle
        # contig.
        if left_edge1_weight > cutoff * left_edge2_weight:
            if right_edge1_weight > cutoff * right_edge2_weight \
            and left_edge1_node2 != right_edge1_node2:
                path.append(left_edge1_node2)
                path.append(right_edge1_node2)
                path.append(target_node)

            elif right_edge2_weight > cutoff * right_edge1_weight \
            and left_edge1_node2 != right_edge2_node2:
                path.append(left_edge1_node2)
                path.append(right_edge2_node2)
                path.append(target_node)

        elif left_edge2_weight > cutoff * left_edge1_weight:
            if right_edge1_weight > cutoff * right_edge2_weight \
            and left_edge2_node2 != right_edge1_node2:
                path.append(left_edge2_node2)
                path.append(right_edge1_node2)
                path.append(target_node)

            elif right_edge2_weight > cutoff * right_edge1_weight \
            and left_edge2_node2 != right_edge2_node2:
                path.append(left_edge2_node2)
                path.append(right_edge2_node2)
                path.append(target_node)

        return path

    def extend(self,node):
        '''
        Returns a connected node to the given node, if there is
        an unambiguous connection.
        Ambiguous connections include:
            From starting node:
            1. More that one outgoing edge
            2. No outgoing edges
            If there is exactly one outgoing edge, the connected node has:
            3. More than one incoming edges
        '''

        if node not in self.nodes:
            raise Exception("Linkgraph error: node {} not in graph".format(node))

        else:
            # Edges are ordered, so we move from left to right along the edge
            connections = self.findConnections(node)

            # If only one outgoing edge, we continue to check the connected node
            # for a single edge back to the original node
            if len(connections) == 1:
                if len(self.findConnections(connections[0])) == 1:
                    # If so, we have found an unambigous connection and can return it
                    return connections

            # Now it starts getting complicated. Sometimes short contigs are connected
            # at both nodes to two different nodes. In these cases, use the weight
            # to determine the path.
            # Check that there are two connections, and that they are to the same contig
            elif len(connections) == 2 and connections[0][:-1] == connections[1][:-1]:
                # If true, also check that both connected nodes have four edges,
                # one back at each side, and one continuing on each side to another node
                middle_ctg = connections[0][:-1]
                conn1 = self.findConnections(connections[0])
                conn2 = self.findConnections(connections[1])
                if len(conn1) == 2 and len(conn2) == 2:
                    if node in conn1 and node in conn2:
                        conn1.remove(node)
                        conn2.remove(node)

                        # If the other connected node is the same,
                        # and has no other incoming edges
                        if conn1 == conn2:
                            third_node_connections = self.findConnections(conn1[0])
                            if len(third_node_connections) == 2:
                                # Success. Check orientation
                                ori = self.bestOrientation( middle_ctg, \
                                                            self.getEdges(node), \
                                                            self.getEdges(conn1[0]))
                                if ori != None:
                                    return ori
                                else:
                                    # If orientation cannot be found from the
                                    # weight, keep the contig in unknown
                                    # orientation for later to try to
                                    # orient it using overlaps
                                    return [connections[0][:-1]+"u", \
                                            connections[0][:-1]+"u", \
                                            conn1]


            # A third possibility exists: the same as above, but another edge also
            # exists between the starting and target nodes. There are thus 3
            # connections, out of which two are to the same contig (middle).
            elif len(connections) == 3 and len(set( [c[:-1] for c in connections] )) == 2:

                # If true, find which two are to opposite ends of the same contig
                # Three possibilities exist...
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
                middle_ctg_node1_conn = self.findConnections(middle_ctg_node1)
                middle_ctg_node2_conn = self.findConnections(middle_ctg_node2)
                target_node_conn = self.findConnections(target_node)

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

                        # Success. Try to find orientation
                        # Get rid of the edges between starting_node
                        # and target_node, then use the same function as before
                        ori = self.bestOrientation( middle_ctg, \
                                                    self.getEdges(node), \
                                                    self.getEdges(target_node))

                        if ori != []:
                            return ori
                        else:
                            # If orientation cannot be found from the
                            # weight, keep the contig in unknown
                            # orientation for later to try to
                            # orient it using overlaps
                            return [middle_ctg_node1[:-1]+"u", \
                                    middle_ctg_node1[:-1]+"u", \
                                    target_node]
        return None


    def findPath(self,starting_node):
        '''
        From the given node, traverses the graph in both directions until
        it becomes impossible.
        Returns every visited node in a list.
        '''
        visited = [starting_node] # Add starting node to list
        node = None

        def assessExtension(extension, side):
            '''
            Nested function to check list of path extension
            '''
            nonlocal visited

            if len(extension) == 1:
                node = extension[0]

            # The function extend should have already checked that
            # the 3 connections are traverseable
            elif len(extension) == 3:
                 # If extension to the right
                if side == "r":
                    visited = visited + extension[:2]
                    node = extension[2]

                # If extension to the left
                elif side == "l":
                    extension.reverse()
                    visited = extension[1:] + visited
                    node = extension[0]
            else:
                node = None
            return node

        ext = self.extend(starting_node) # Go to next node
        if ext != None:
            node = assessExtension(ext, "r")

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
                node = assessExtension(ext, "r")
            else:
                node = None

        # At this point we have reached the end of the right extension
        # Now do the same in the other direction
        node = self.opposite(starting_node)
        visited.insert(0, node) # Insert at the beginning of the visited list
        ext = self.extend(node) # Grab next nodes
        if ext != None:
            node = assessExtension(ext, "l")
        else:
            node = None

        while node != None:
            visited.insert(0, node) # Add this node to list

            # Go to next node, which is always the opposite end of the newly
            # connected contig
            node = self.opposite(node)
            visited.insert(0, node) # Add this node to list
            ext = self.extend(node) # Go to next node
            if ext != None:
                node = assessExtension(ext, "l")
            else:
                node = None

        return visited


def unambiguousPaths(linkgraph):
    '''
    Returns all unambiguous paths in the graph, i.e. paths with a single
    edge connecting every node in the path.
    '''
    visited_nodes = []
    paths = []
    for n in linkgraph.nodes:
        if n not in visited_nodes:
            path = linkgraph.findPath(n)
            paths.append(path)
            visited_nodes = visited_nodes + path

    return paths

def makeNodes(contig_list):
    return [a+"s" for a in contig_list] + [a+"e" for a in contig_list]


def compareGEMlibs(lib1,lib2):
    '''
    Compares two lists of barcodes, collects all shared ones and
    counts them. Returns fraction of shared barcodes.
    '''
    shared = set(lib1).intersection(lib2) # Find shared ones
    totallength = len(lib1) + len(lib2) - len(shared) # Total number of unshared barcodes

    # Find the fraction of shared barcodes, avoid division by 0
    if totallength != 0:
        fraction = len(shared) / totallength
    else:
        fraction = 0

    return fraction

def pairwise_comparisons(GEMlist):
    '''
    Performs all pairwise comparisons between windows in GEMlist.
    '''
    # Ccompare the barcodes in every region to all other regions
    donewindows = 0
    GEMcomparison = {}

    # Iterate over the GEMlist dict
    for region1, lib1 in GEMlist.items():
        contig = region1[:-1]

        # Report progress every 20 windows
        if donewindows in range(0,100000000,20):
            misc.printstatus(misc.reportProgress(donewindows, len(GEMlist)))

        nested_dict = {} # Each entry in the original dictionary is another dictionary

        # Iterate over the GEMlist dict inside the loop
        # This must be really inefficient...
        for region2, lib2 in GEMlist.items():
            fraction = compareGEMlibs(lib1,lib2) # Collect fraction of shared barcodes

            # Fill the nested dictionary
            nested_dict[region2] = fraction

        donewindows += 1
        GEMcomparison[region1] = nested_dict

    return GEMcomparison

def makeEdges(GEMcomparison, min_barcode_fraction):
    # Next step is to create a graph in a dict format where the
    # edges are based on outlying fraction values.
    # If a certain window has an outlier to another window, an edge is
    # created

    donewindows = 0 # Reuse this variable to keep track of progress
    misc.printstatus("Number of windows: "+str(len(GEMcomparison)))
    edges = []

    # Iterate over keys in GEMcomparison
    for region, fractions in GEMcomparison.items():
        contig = region[:-1]
        window = region[-1]

        # Report progress every 100 windows
        if donewindows in range(0,10000000,100):
            misc.printstatus(misc.reportProgress(donewindows, len(GEMcomparison)))

        # Calculate outliers from the comparisons of window k to all other windows
        outliers = esd.getOutliers(fractions)
        outliers_short = {}

        # Remove fractions less than min_barcode_fraction (-f), including lower outliers
        # and outliers with too few shared barcodes (-n)
        for k,v in outliers.items():
            if v > min_barcode_fraction:
                outliers_short[k] = v

        # Sort outliers_short dict by fraction value into the list sorted_outliers
        # skip forming edges within the same contig
        sorted_outliers = [(region, k, outliers_short[k]) \
                            for k in sorted(outliers_short, \
                                            key=outliers_short.get, \
                                            reverse=True) \
                            if region[:-1] != k[:-1]]

        if sorted_outliers != []:
            edges = edges + sorted_outliers
        donewindows += 1

    return edges

def makeScaffolds(paths):
    '''
    Given a list of paths, determines orientation of contigs.
    '''
    scaffolds = {}
    nr = 1
    for p in paths:
        linked_scaffold_name = "scaffold_{}".format(str(nr))
        # Check every other node to get contig and orientation
        for c in p[1::2]:
            tig = c[:-1]
            orientation = c[-1]

            if orientation == "e":
                if linked_scaffold_name not in scaffolds:
                    scaffolds[linked_scaffold_name] = [tig+"f"]
                else:
                    scaffolds[linked_scaffold_name].append(tig+"f")

            elif orientation == "s":
                if linked_scaffold_name not in scaffolds:
                    scaffolds[linked_scaffold_name] = [tig+"r"]
                else:
                    scaffolds[linked_scaffold_name].append(tig+"r")

            elif orientation == "u":
                if linked_scaffold_name not in scaffolds:
                    scaffolds[linked_scaffold_name] = [tig+"u"]
                else:
                    scaffolds[linked_scaffold_name].append(tig+"u")

        nr += 1
    return scaffolds

def main(contig_list, GEMlist, min_barcode_number, min_barcode_fraction):
    '''
    Controller for graph_building.
    '''
    # Collect the fraction of shared barcodes in the all-against-all
    # comparison of windows
    GEMcomparison = pairwise_comparisons(GEMlist)

    # Infer linkage based on statistically significant outliers determined
    # by the ESD test to build the graph
    nodes = makeNodes(contig_list)
    edges = makeEdges(GEMcomparison, min_barcode_fraction)
    graph = Linkgraph(nodes, edges)
    paths = unambiguousPaths(graph) # Determine unambiguous paths
    scaffolds = makeScaffolds(paths)

    #for k,v in scaffolds.items():
    #    print(k,v)

    return graph, scaffolds
