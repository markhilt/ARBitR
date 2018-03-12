# To format a 2-dimensional dictionary to a human readable table
def formatTable(dict1):
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

def formatGFA(dict1):
    gfalist = []

    # Write link for every edge in the input dict
    for k in dict1:
        contig = k[0:-1]
        side = k[-1]

        for v in dict1[k]:
            ltig = v[0:-1]
            lside = v[-1]

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
