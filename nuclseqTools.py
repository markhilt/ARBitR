def reverse_complement(nuclstring):
    rev_comped = ""
    for l in reversed(nuclstring):
        if l == "A" or l == "a":
            rev_comped += "T"
        elif l == "T" or l == "t":
            rev_comped += "A"
        elif l == "C" or l == "c":
            rev_comped += "G"
        elif l == "G" or l == "g":
            rev_comped += "C"
        elif l == "N" or l == "n":
            rev_comped += "N"
    return rev_comped

def createConsensus(delta,string1,string2):
    '''
    Given two overlapping nucleotide strings (excluding overhangs) and alignment information,
    returns the merged sequence
    delta = [int,int,int, ...]
    strings = "ATCG..."
    '''
    new_string1 = string1
    new_string2 = string2
    start = 0
    for i in delta:
        if i > 0:
            start += i
            new_string2 = new_string2[:start-1] + "." + new_string2[start-1:]

        else:
            i = -i
            start += i
            new_string1 = new_string1[:start-1] + "." + new_string1[start-1:]

    # Write the new string. Insertions are favoured over deletions.
    # N's are put at mismatches
    cons = ""
    for i in range(0,len(new_string1),1):
        if new_string1[i] == new_string2[i]:
            cons += new_string1[i]
        elif new_string1[i] == ".":
            cons += new_string2[i]
        elif new_string2[i] == ".":
            cons += new_string1[i]
        else:
            cons += "N"

    return cons
