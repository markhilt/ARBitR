# AnVIL
## Assembly Validation and Improvement using Linked reads

Script collection for quality control and merging/scaffolding of genome assemblies using 10X Chromium linked reads.

## Dependencies
- Python3
- [Pysam](https://github.com/pysam-developers/pysam)
- [Mummer](https://pypi.org/project/mappy/)

## Description
AnVIL is used for merging and scaffolding a genome assembly. It uses pysam to scan the ends of contigs in a given bam file, collects all GEM barcodes from reads that align to each contig end, and compares the barcodes in all contig ends. An object called a linkgraph is then built, where nodes (contig ends) with a higher-than-expected fraction of shared barcodes are connected with an edge. Unique paths through the graph can be used to build scaffolds. First, contig ends without read coverage are trimmed off. Mappy is then used to try to find an overlap between the contig ends. If an overlap is found, the contigs are fused together. If no overlap is found, the contigs are merged with an N*10 gap.

AnVIL requires as input **_only_** a bam file of mapped linked reads.

### Example usage:
```
./anvil.py <in.bam> # Builds the linkgraph and outputs it as a gfa file
./anvil.py -i <in.fasta> <in.bam> # Builds the linkgraph and uses it to merge the contigs in <in.fasta>
```

## Output
AnVIL will always produce a gfa file for visualization of the graph with e.g. Bandage. It can also optionally be given a fasta file, in which case it will also produce a merged/scaffolded new fasta file. To merge and scaffold, the given fasta file must be the same as the one used as a reference for the read mapping.

## Known limitations
AnVIL is highly dependent on high-quality read mappings. It is highly recommended to run [Pilon](https://github.com/broadinstitute/pilon) before AnVIL, and use [Longranger](https://github.com/10XGenomics/longranger) for mapping linked reads. As always, collapsed repeats are difficult to handle, and may appear to link together unrelated contigs. Short contigs made up entirely of repeats will be impossible to link because of the lack of seeding regions.
