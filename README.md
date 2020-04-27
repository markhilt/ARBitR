# AnVIL
## Assembly Validation and Improvement using Linked reads

Script collection for merging/scaffolding of genome assemblies using 10X Chromium linked reads.

## Dependencies
- Python3
- [numpy](https://numpy.org/)
- [scipy](https://www.scipy.org/)
- [pandas](https://pandas.pydata.org/)
- [Pysam](https://pysam.readthedocs.io/en/latest/api.html)
- [Mappy](https://pypi.org/project/mappy/)

## Summary
AnVIL is used for merging and scaffolding an existing genome assembly. It takes a position sorted read alignment file in the bam/sam format with 10X Chromium barcodes in the BX tag, and if additionally provided with the genome fasta file used for mapping, it will sort and merge the provided contigs into scaffolds. A key functionality is the possibility to resolve links between contigs by overlap-layout-consensus (OLC) instead of just inserting a gap between them.

AnVIL was developed for application on small genomes that have been assembled using long-read sequencing and an OLC method. Ideally, the genome size is <1Gb and the number of input contigs is <10000, though it may also be useful in other cases. E.g. it can additionally scaffold [Supernova](https://github.com/10XGenomics/supernova) assemblies of pure 10X Chromium data.

## Installation
Install dependencies using e.g. pip or conda. Then clone this repo and run AnVIL as below.

## Example usage:
```
./anvil.py <in.bam> # Builds the linkgraph and outputs it as a gfa file
./anvil.py -i <in.fasta> <in.bam> # Builds the linkgraph and uses it to merge the contigs in <in.fasta>
```

## Output
AnVIL will always produce a gfa file for visualization of the link graph with e.g. Bandage. It can also optionally be given a fasta file, in which case it will also produce a merged/scaffolded new fasta file. To merge and build scaffolds, the given fasta file must be the same as the one used for the read mapping.

## Algorithm overview
AnVIL works in three steps: (1) backbone graph creation, (2) filling of junctions with short contigs and (3) merging of contigs. During (1), GEM barcodes are collected from read alignments at the start and end regions of long contigs. The length is determined by -m, and should correspond approximately to the mean molecule size that went into 10X Chromium library creation. The region size to collect barcodes from is determined by -s and should be chosen so as high-quality read alignments are found (as contigs often end in repeat regions, reads mapping to these regions will have lower MAPQ - choose -s so that some anchoring region of the contigs is covered). Barcode profiles are then compared for every contig end and pairs with a significant fraction of shared barcodes are saved. A linkgraph object is then created, where nodes correspond to contig ends and edges to pairs with a significant fraction of shared barcodes. Paths through the linkgraph are then determined and described as lists of junction objects, where each junction has a start node, an end node, a barcode profile and a list of connected nodes. Each path will form a linked scaffold. The linkgraph is written to disk in gfa format. During the next step (2), the full lengths of the contigs shorter than -m are used for barcode collection. The barcode profiles in the path junctions are then used to pull in short contigs that have a significant fraction of shared barcodes to the junction. Significant hits are stored in the junction connections. During step (3), the paths are traversed, where the sequences in every junction are resolved by overlap-layout-consensus from the start to the end node. If no unambiguous overlap path through junction is found, a gap is instead introduced and connections that could not be merged into the junction are dropped to be written separately. A fasta file of the new scaffolds and remaining contigs is finally written to disk.

AnVIL is written and implemented in Python 3.

## Known limitations
AnVIL is highly dependent on high-quality read mappings. It is recommended to run [Pilon](https://github.com/broadinstitute/pilon) before AnVIL, and use the "align" program of [Longranger](https://github.com/10XGenomics/longranger) or [EMA](https://github.com/arshajii/ema) for mapping linked reads. As always, collapsed repeats are difficult to handle, and may appear to link together unrelated contigs. Short contigs made up entirely of repeats will be impossible to link because of the lack of seeding regions.
