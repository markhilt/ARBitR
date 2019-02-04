import mummerTools as mt
import nuclseqTools as nt

def opposite(node):
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
    opposite_end = "".join([tig,opposite])
    return opposite_end

def findStartingNode(graph):
    '''
    ##### Find starting node #####
    Finding a good place to start a new contig is surprisingly challenging.
    For every given node, one must consider the number of edges from this node,
    the number of edges from every connected node,
    and the number of edges from the opposite side of the given node.
    There are a total of 5 different scenarios considered here.
    1.  Starting node has a single outgoing edge. Connected node has two edges, one back and one to
        the opposite end of starting node. Opposite end has only this edge.
    The rest of the scenarios have either no edges, one edge, two different, three different, or
    more than three edges at the opposite end. If one edge, this needs to be connected to a node
    with several connections.
    2.  Starting node has a single outgoing edge. Connected node has also only one (the same).
    3.  Starting node has two edges to the same node. Both connections only have edges back
        to the starting node.
    4.  Starting node has two edges to the same node. Both connections have two edges;
        one back and one to a new, third node. This third node has two connections, back
        to the middle node.
    5.  Starting node has three edges, two to the same node and one to a third node.
        Both connections of the middle node have two edges each, back to the starting
        and forward to the third node. Third node has three connections, back
        to both sides of the middle node and one to the starting node.

    "graph" is the graph in a dict format, and every node which has an edge should be
    represented in the keys
    '''

    linked_contig_ends = []
    # Go through every node
    for node in graph:

        # Collect incoming edges
        conn = findLink(node)

        # Collect opposite end of this node
        op = opposite(node)

        # And edge(s) to opposite end
        op_edge = findLink(op)

        # Check for 1
        if len(conn) == 1 \
        and conn == op_edge \
        and len(op_edge) == 1 \
        and len(findLink(conn[0])) == 2 \
        and findLink(conn[0])[0][:-1] == findLink(conn[0])[0][:-1]:
            node = node[:-1] + "u"
            linked_contig_ends.append(node)
            continue

        # Other scenarios all have shared qualities of opposite end, see above
        if len(op_edge) == 0 \
        or (len(op_edge) == 1 \
        and len(findLink(op_edge[0])) > 1) \
        or ( len(op_edge) == 2 \
        and op_edge[0][:-1] != op_edge[1][:-1] )\
        or ( len(op_edge) == 3 \
        and op_edge[0][:-1] != op_edge[1][:-1] \
        and op_edge[1][:-1] != op_edge[2][:-1] \
        and op_edge[0][:-1] != op_edge[2][:-1])\
        or len(op_edge) > 3:

            # Check for 2
            if len(conn) == 1 \
            and len(findLink(conn[0])) == 1:
                linked_contig_ends.append(node)
                continue

            # Check for 3,4
            if len(conn) == 2 \
            and conn[0][:-1] == conn[1][:-1]:

                # Check for 3
                if len(findLink(conn[0])) == 1 \
                and len(findLink(conn[1])) == 1:
                    linked_contig_ends.append(node)
                    continue

                # Check for 4
                if len(graph[conn[0]]) == 2 \
                and len(graph[conn[1]]) == 2 \
                and sorted(graph[conn[0]]) == sorted(graph[conn[1]]):
                    linked_contig_ends.append(node)
                    continue

            # Check for 5
            threew = find3way(node, conn)
            if len(conn) == 3 \
            and threew != None:
                linked_contig_ends.append(node)

    return linked_contig_ends

def findLink(region):
    '''
    Given a node in the format tigXXXs/e,
    returns all edges from input node
    '''
    # Collect all edges matching input region
    current_edges = []
    for (a,b) in edges_list:
        if a == region:
            current_edges.append(b)
        elif b == region:
            current_edges.append(a)
    return current_edges

def find3way(starting_node,list_of_edges):
    '''
    starting_node is a string of the format tigs/e, list_of_edges a list of len 3
    If successful, return two strings: same and diff
    same is the middle node which has two edges and
    diff the second node
    If unsuccessful return None
    '''
    if len(list_of_edges) != 3:
        return None
    else:
        # First find which 2 edges are to the same node on opposite ends
        node1 = list_of_edges[0][:-1]
        node2 = list_of_edges[1][:-1]
        node3 = list_of_edges[2][:-1]
        if node1 == node2:
            same = node1
            diff = list_of_edges[2]
        elif node1 == node3:
            same = node1
            diff = list_of_edges[1]
        elif node2 == node3:
            same = node2
            diff = list_of_edges[0]
        else:
            return None

        # If two of the three edges are in fact to the same node,
        # they would be expected to have two edges each:
        # 1 to starting_node and 1 to the next node (diff)
        if same != None:
            # Collect edges from both ends of the middle node (= same)
            links_same_start = findLink("".join( [same,"s"] ))
            links_same_end = findLink("".join( [same,"e"] ))

            # First check that both have len == 2,
            # i.e. there are no other nodes involved
            if not len(links_same_start) == 2 or not len(links_same_end) == 2:
                return None

            # Then check that i and diff are in both
            if starting_node in links_same_start and diff in links_same_start \
            and starting_node in links_same_end and diff in links_same_end:
                # If so, the middle node seems ok. Last thing
                # to do is to check diff for any other edges
                links_diff = findLink(diff)

                # Expected len == 3: both ends of middle node,
                # and starting_node
                if len(links_diff) != 3:
                    return None

                if starting_node in links_diff \
                and "".join( [same,"s"] ) in links_diff \
                and "".join( [same,"e"] ) in links_diff:
                    return (same,diff)

def create_linked_contig(starting_region):
    '''
    Links together contigs from the starting_region.
    '''
    tig = starting_region[:-1]
    side = starting_region[-1]
    i = findLink(starting_region)
    t = None

    if side == "e":
        orientation = "f"
    elif side == "s":
        orientation = "r"
    else:
        orientation = "u"
        i = findLink(tig+"e") # Collect link. Both ends are connected to the same node so "e" is arbitrarily chosen

    if orientation != None:
        # IF direction is still unknown, i.e. no good alignment was found,
        # don't add it into the [links]
        links = [ "".join([tig,orientation]) ]

    # If there is only one link, all is good
    if len(i) == 1:
        t = i[0]

    # Two links to the same node means mummerTools module must be used again
    # to determine direction
    elif len(i) == 2 \
    and i[0][:-1] == i[1][:-1]:
        #iseq = i[:-1]
        #idir = i[-1]
        #ref_fasta = fastafile.fetch(reference=tig, start=trimmed_fasta_coords[tig][0], end=trimmed_fasta_coords[tig][1])
        #query_fasta = fastafile.fetch(reference=iseq, start=trimmed_fasta_coords[iseq][0], end=trimmed_fasta_coords[iseq][1])
        #delta, reflen, querylen = mt.align(ref_fasta,query_fasta)
        #alignment, orientation = mt.findDirectionQuery(delta, reflen, querylen)

        if orientation != None:
            links.append( "".join([i[0][:-1],orientation]) )
        it = findLink(i[0])
        if len(it) == 1:
            t = None
        elif len(it) == 2:
            it.remove(starting_region)
            t = it[0]

    elif len(i) == 3:
        s = find3way(starting_region, i)
        if s != None:
            same, diff = s[0], s[1]
            links.append( "".join([same,"u"]) ) # Direction of middle node is unknown
            t = diff

    while t != None:
        tig = t[:-1]
        side = t[-1]
        if side == "s":
            orientation = "f"
            opposite = "e"
        elif side == "e":
            orientation = "r"
            opposite = "s"
        else:
            orientation = "u"
            opposite = "u"

        links.append( "".join([tig,orientation]) )

        # Look at opposite side of current node for unqiue edges
        next_edge = "".join([tig,opposite])
        l = findLink(next_edge)
        # If only one matching edge, check that connected node
        # has no other edges
        if len(l) == 1 and len(findLink(l[0])) == 1:
            t = l[0]

        # If two edges, check if they are to opposite sides
        # of the same node
        elif len(l) == 2 and l[0][:-1] == l[1][:-1]:
            links.append( "".join([l[0][:-1],"u"]) ) # Direction of middle node is unknown
            middle_edges1 = findLink(l[0])
            middle_edges2 = findLink(l[1])
            if middle_edges1 == middle_edges2 \
            and len(middle_edges1) == 2 and next_edge in middle_edges1:
                middle_edges1.remove(next_edge)
                t = middle_edges1[0]
            else:
                t = None

        # If 3 edges, there is a chance two of them are to the same
        # node, on opposite sides, and the last to another node
        # This means that the one where both ends are connected
        # belongs in between.
        elif len(l) == 3:
            s = find3way(next_edge,l)
            if s != None:
                same, diff = s[0], s[1]
                # Congrats, you have found a hit!
                links.append( "".join([same,"u"]) ) # Direction of middle node is unknown
                t = diff
            else:
                t = None

        else:
            t = None

    return links

def edges_as_list(graph):
    '''
    Converts a graph in dict format to a list of tuples of every edge
    '''
    global edges_list
    edges_list = []
    for i in graph:
        for a in graph[i]:
            # If edge was not already added in the opposite direction, add it
            if ( (a,i) ) not in edges_list:
                edges_list.append( (i, a) )
