"""formatting_tools

Functions to reformat objects for writing them to files.

Copyright (c) 2020, Markus Hiltunen
Licensed under the GPL3 license. See LICENSE file.
"""

def formatTable(dict1):
    '''
    To format a 2-dimensional dictionary to a human readable table
    '''
    firstrow = "Contig"
    remaining_table = ""
    for k in sorted(dict1):
        firstrow = firstrow + "\t" + k
        values = []

        for nk in sorted(dict1[k]):
            frac = dict1[k][nk][0]
            values.append(str(frac))

        nextrow = k + "\t" + "\t".join(values)
        remaining_table = remaining_table + nextrow + "\n"

    table = firstrow + "\n" + remaining_table
    return table

def formatGFA(linkgraph):
    '''
    Given linkgraph object, returns a string in gfa format.
    '''
    gfalist = []

    # Write link for every edge in the input graph
    for edge in linkgraph.edges:
        ref_tig, ref_side = edge[0][:-1], edge[0][-1]
        query_tig, query_side = edge[1][:-1], edge[1][-1]

        # Find out which contig ends are linked
        # Possibilities:
        # 1. start of ref to start of query
        # 2. start of ref to end of query
        # 3. end of ref to start of query
        # 4. end of ref to end of query

        # 1.
        if ref_side == "s" and query_side == "s":
            direction = "-"
            ldirection = "+"

        # 2.
        elif ref_side == "s" and query_side == "e":
            direction = "-"
            ldirection = "-"

        # 3.
        elif ref_side == "e" and query_side == "s":
            direction = "+"
            ldirection = "+"

        # 4.
        elif ref_side == "e" and query_side == "e":
            direction = "+"
            ldirection = "-"

        if direction and ldirection:
            gfalist.append("L\t{}\t{}\t{}\t{}\t*\n".format(ref_tig,direction,query_tig,ldirection) )

    return "".join(gfalist)


# Deprecated
def formatGFA_from_dict(dict1):
    '''
    Given a graph in a dict, returns a list of strings in gfa format
    '''
    gfalist = []

    # Write link for every edge in the input dict
    for k, val in dict1.items():
        contig = k[0:-1]
        side = k[-1]

        for v in val:
            ltig = v[0][0:-1]
            lside = v[0][-1]

            # Find out which contig ends are linked
            # Possibilities:
            # 1. start of ref to start of query
            # 2. start of ref to end of query
            # 3. end of ref to start of query
            # 4. end of ref to end of query

            # 1.
            if side == "s" and lside == "s":
                direction = "-"
                ldirection = "+"

            # 2.
            elif side == "s" and lside == "e":
                direction = "-"
                ldirection = "-"

            # 3.
            elif side == "e" and lside == "s":
                direction = "+"
                ldirection = "+"

            # 4.
            elif side == "e" and lside == "e":
                direction = "+"
                ldirection = "-"

            if direction and ldirection:
                gfalist.append("L\t{0}\t{1}\t{2}\t{3}\t*\n".format(contig,direction,ltig,ldirection) )

    return gfalist

def formatTSV(dict1, l):
    '''
    Given a dict, returns a string in tsv format
    '''
    tsvlist = []

    # Write link for every edge in the input dict
    for k in dict1:
        contig = k[:-1]
        side = k[-1]

        for v in dict1[k]:
            ltig = v.split("-")[0]
            lside = v.split["-"][1]

            # Orientation of contigs, to be filled later depending on which of 1-4 is true
            direction = ""
            ldirection = ""

            # 1.
            if side == "s" and lside == "s":
                direction = "r"
                ldirection = "f"

            # 2.
            elif side == "s" and lside == "e":
                direction = "r"
                ldirection = "r"

            # 3.
            elif side == "e" and lside == "s":
                direction = "f"
                ldirection = "f"

            # 4.
            elif side == "e" and lside == "e":
                direction = "f"
                ldirection = "r"

            bcnum = l[1][1] # Number of shared barcodes
            tsvlist.append("500\t{0}{1}\t{2}{3}\t{4}\t{5}\n".format(direction,contig,ldirection,ltig,bcnum,bcnum*10) )

    return tsvlist

def formatBed(dict):
    '''
    Given a dict where keys are ID's of linked contigs and values are lists
    containing tuples of the format ( str(feature), int(length) ), returns a string
    in the bed format.
    '''
    bed = []
    for key,value in dict.items():
        current_coordinate = 0 # Keep track of the coordinate to write new features to

        for feat in value:
            line = []
            line.append(key)
            line.append(str(current_coordinate))
            line.append(str(current_coordinate+feat[2]))
            line.append(feat[0])
            line.append("0")
            if feat[1] == "f":
                line.append("+")
            elif feat[1] == "r":
                line.append("-")
            else:
                line.append(".")
            current_coordinate = current_coordinate+feat[2]
            bed.append("\t".join(line))

    return "\n".join(bed)
