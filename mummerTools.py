from subprocess import call

def align(seq1,seq2):
    '''
    Calls mummer externally to align two nucleotide sequences given in seq1 and seq2
    Writes a delta file in working directory, reads it and returns in as a list
    Requires mummer in $PATH
    '''
    # Write trimmed contigs to temp files for mummer alignment in working directory
    with open("ref.fa","w",encoding = "utf-8") as fastaout:
        fastaout.write(">seq1\n")
        fastaout.write(seq1)
    with open("query.fa","w",encoding = "utf-8") as fastaout:
        fastaout.write(">seq2\n")
        fastaout.write(seq2)

    # Call nucmer
    call(["nucmer", "ref.fa", "query.fa"])

    # Parse out.delta and reformat into dict
    with open("out.delta","r",encoding = "utf-8") as delta:
        passed = False # To check if inside a new alignment
        ali = [] # Reconstruct alignment in this list
        d = {} # Reconstruct delta in this dict
        reflen, querylen = None, None # To return in case there are no alignments

        for line in delta:
            line = line.strip()
            fields = line.split(" ")

            # Collect length of reference sequence
            if line.startswith(">") and len(fields) == 4:
                reflen = fields[2]
                querylen = fields[3]

            # Start collecting new alignment
            if len(fields) == 7:
                k = line
                ali = [] # Reset alignment
                passed = True

            # Current alignment stops at 0's in the delta file
            # Add current alignment to the dict d and go to next alignment
            if line == "0":
                d[k] = ali
                continue

            # Append as long as there is one int per line
            if passed == True and len(fields) == 1:
                ali.append(int(line))

    return d, reflen, querylen

def getCoords(alignment):
    '''
    Returns reference and query coordinates from given alignment
    '''
    fields = alignment.split(" ")
    alignment_ref = fields[0:2]
    alignment_query = fields[2:4]
    return alignment_ref, alignment_query

def findAlignment(delta, reflen):
    '''
    Expects to find an overlap between the end of seq1 and beginning of seq2
    in the given delta. If so, returns it, else return None.
    '''
    # Filter out correct alignment in the delta
    # New alignments have 7 fields
    if len(delta) > 0:
        for alignment in delta:
            alignment_ref, alignment_query = getCoords(alignment)

            # findAlignment assumes that sequences that went into alignment
            # were revcomped so that the end of ref should align to the beginning of
            # query
            if alignment_ref[1] == reflen and int(alignment_query[0]) == 1:
                # Correct alignment found. Return it as a list
                correct_alignment = [alignment] + delta[alignment]
                return correct_alignment
    return None

def findDirectionQuery(delta, reflen, querylen):
    '''
    Given a delta dict, attempts to find which direction the query contig should be merged
    into the linked contig based on which side of it that has the best alignment
    to the reference sequence. If both sides align, the direction is determined from
    the best alignment score (todo).
    '''
    if len(delta) > 0:
        for alignment in delta:
            alignment_ref, alignment_query = getCoords(alignment)

            # Alignments of interest end at the last reference coordinate
            if alignment_ref[1] == reflen:

                # Then look at the query start coordinate
                qstart = alignment_query[0]
                if qstart == 1:
                    # Query aligns in forward direction, return alignment and direction
                    correct_alignment = [alignment] + delta[alignment]
                    return correct_alignment,"f"

                if qstart == querylen:
                    # Query aligns in reverse direction, return alignment and direction
                    correct_alignment = [alignment] + delta[alignment]
                    return correct_alignment,"r"
    return None, None

def findDirectionRef(delta, reflen, querylen):
    '''
    Given a delta dict, attempts to find which direction the reference contig should be merged
    into the linked contig based on which side of it that has the best alignment
    to the reference sequence. If both sides align, the direction is determined from
    the best alignment score (todo). Assumes that query sequences that are linked
    in reverse orientation have been reverse complemented.
    '''
    if len(delta) > 0:
        for alignment in delta:
            alignment_ref, alignment_query = getCoords(alignment)
            # A bit more complex than findDirectionQuery, as
            # the delta format is based on the reference orientation
            # We expect an alignment that starts at the beginning of the query
            # sequence.
            if int(alignment_query[0]) == 1 and int(alignment_ref[1]) == reflen:
                # Referece aligns in reverse direction, return alignment and direction
                correct_alignment = [alignment] + delta[alignment]
                return correct_alignment,"f"

            # If reference aligns in reverse, this is at alignment_query[1]
            if int(alignment_query[1]) == 1 and int(alignment_ref[0]) == 1:

                # Referece aligns in reverse direction, return alignment and direction
                correct_alignment = [alignment] + delta[alignment]
                return correct_alignment,"r"

    return None, None
