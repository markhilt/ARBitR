import numpy as np
import pysam
from scipy.stats import t
import misc

def collectGEMs(window, mapq):
    """
    Collects the barcodes from the given window
    """
    tig, start, stop = window[0], window[1], window[2]
    reads = samfile.fetch(tig, start, stop)
    BC_list = []
    for read in reads:
        if read.has_tag('BX') == True and read.mapping_quality >= mapq:
            BC = read.get_tag("BX")
            BC_list.append(BC)
    BC_list = set(BC_list)

    return BC_list

def getWindows(contig, region_size, contig_dict):
    '''
    To split up the contig into windows
    '''
    length = contig_dict[contig]
    start = (0,region_size)

    # If contig is shorter than region_size, just take all of it to "end" as well
    if (length - region_size) < 0:
        endstart = 0
    else:
        endstart = length - region_size
    end = (endstart,length)
    return [start,end]

def reportProgress(current,total):
    return "Completed: {0}% ({1} out of {2})".format( str(round( (current / total) * 100, 2)), current, total)


### CHECK WHERE IN USE
def countReads(contig,coords_to_check):
    '''
    To count the number of reads aligned to a region
    '''
    cov_arr = samfile.count_coverage(contig, coords_to_check[0],coords_to_check[1])
    cov = sum([sum(cov_arr[x]) for x in range(0,4,1)]) / 50
    return cov

def formatContigs(samfile):
    '''
    Creates a dict where keys are contig names and values their lengths
    '''
    i = 0
    gfa_header = []
    contig_dict = {}
    contig_list = samfile.references
    while i < len(contig_list):
        contig_dict[contig_list[i]] = samfile.lengths[i]
        gfa_header.append("S\t{0}\t*\tLN:i:{1}".format(contig_list[i],samfile.lengths[i]))
        i += 1

    return contig_dict, gfa_header

def main(input_bam, region_size, mapq):
    global samfile
    samfile = pysam.AlignmentFile(input_bam, "rb")
    contig_dict, gfa_header = formatContigs(samfile)
    GEMlist = {} # Inappropriately named "list"

    # First step is to collect all barcodes (passing -q cutoff) that are aligned
    # to each contigs first and last regions (-l)
    misc.printstatus("Starting barcode collection. Found {0} contigs.".format(len(contig_dict)))
    donecontigs = 0
    for contig in contig_dict:

        # Report progress every 20 windows
        if donecontigs in range(0,100000000,20):
            misc.printstatus("[ BARCODE COLLECTION ]\t"+reportProgress(donecontigs, len(contig_dict)))

        # Create windows from first and last X kb from each contig
        windowlist = getWindows(contig, region_size, contig_dict)

        # Collect barcodes from every window
        for window in windowlist:
            GEMs = collectGEMs( (contig, window[0], window[1]), mapq)

            if window == windowlist[0]:
                region = contig + "s"
            elif window == windowlist[1]:
                region = contig + "e"

            if GEMs:
                GEMlist[region] = GEMs

        donecontigs += 1

    samfile.close()

    return GEMlist, gfa_header
