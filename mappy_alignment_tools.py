import mappy as mp

def chopCigar(cigar):
    '''
    Given a cigar string, chops it up into a list.
    '''
    cig_list = []
    pos = 0
    for idx, c in enumerate(cigar):
        if c.isalpha():
            cig_list.append(cigar[pos:idx+1])
            pos  = idx + 1
    return cig_list

def findOverlap(seq1,seq2):
    '''
    Tries to find an alignment between the two given sequences
    '''

    # Write seq1 to fasta file
    with open("tmp.fasta", "w") as tmpfasta:
        tmpfasta.write(">tmp\n")
        tmpfasta.write(seq1)

    # Build index from seq1
    idx = mp.Aligner(fn_idx_in="tmp.fasta", preset="asm5")
    ref_len = idx.seq

    # Align
    alignments = idx.map(seq2)

    # Iterate over alignment and search for overlapping ends
    # Assumes that contigs have been oriented,
    # meaning that we search for an overlap between the suffix of
    # seq1 and the prefix of seq2
    if alignments:
        for aln in alignments:

            # Filter for alignments ending within the last 1 kb of seq1
            if aln.r_en > aln.ctg_len - 1000:

                # Then for alignments starting within the first 1 kb of seq2
                if aln.q_st < 1000:
                    return aln

    return None

def createConsensus(aln,string1,string2):
    '''
    Given two overlapping nucleotide strings and mappy alignment information,
    merges and returns the nucleotide strings.

    Assumes the overlap is string1 suffix to string2 prefix.
    '''
    # Track position in both strings
    ref_pos, query_pos = 0, 0

    # First add sequence of string1 up until alignment start
    output_string = [string1[:aln.r_st]]
    ref_pos += aln.r_st
    cig_list = chopCigar(aln.cigar_str) # Convert the cigar string to a list

    # Then walk through the aligned region in the cigar string,
    # gradually building the sequence. Because of the tendency of PacBio
    # to miss some bases, we will use the extra base at every indel position
    for cig in cig_list:

        # If strings match, there is no problem
        # Use whatever sequence, they are the same
        # Then move both position trackers
        if cig[-1] == "M":
            output_string.append( string1[ref_pos:ref_pos+int(cig[:-1])] )
            ref_pos += int(cig[:-1])
            query_pos += int(cig[:-1])

        # If insertion in query, add extra bases from query and increase position
        elif cig[-1] == "I":
            output_string.append( string2[query_pos:query_pos+int(cig[:-1])] )
            query_pos += int(cig[:-1])

        # If deletion in query, add extra bases from reference and increase position
        elif cig[-1] == "D":
            output_string.append( string1[ref_pos:ref_pos+int(cig[:-1])] )
            ref_pos += int(cig[:-1])

    # After iterating over the whole alignment, add remaining bases from
    # query sequence
    output_string.append( string2[query_pos:] )

    return ''.join(output_string)
