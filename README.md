# AnVIL
## Assembly Validation and Improvement using Linked reads

Script collection for quality control and merging/scaffolding of genome assemblies using 10X Chromium linked reads. AnVIL contains the workhorse scripts chromQC and merge. chromQC is used for validation of a genome assembly and merge for merging and scaffolding the assembly. Tested on haploid long-read assemblies.

## Dependencies
- Python3
- Pysam
- Mummer (optional)

## Usage


## Output
merge will always produce a gfa file for visualization of the graph. It can also optionally be given a fasta file, in which case it will also produce a merged/scaffolded new fasta file and a bed file where the features are the merged contigs. To merge and scaffold the given fasta file must be the same as the one used as a reference for the read mapping, and Mummer is needed in $PATH.

chromQC will produce a plain text file containing a report of potential misassemblies, which can be used to manually split these contigs.

## Known limitations
AnVIL is highly dependent on high-quality read mappings. It is highly recommended to run Pilon before AnVIL. As always, collapsed repeats are difficult to handle, and may appear to link together unrelated contigs. Short contigs made up entirely of repeats will be impossible to link because of the lack of seeding regions.