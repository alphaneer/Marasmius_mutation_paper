This is the script collection for quality control and merging/scaffolding of genome assemblies using 10X Chromium linked reads used in the Marasmius mutation paper. For the development version of this pipeline, please see https://github.com/markhilt/AnVIL

![alt text](https://github.com/johannessonlab/Marasmius_mutation_paper/blob/master/genome_assembly/Marasmius_assembly_pipeline.png)
**Genome assembly pipeline developed for Marasmius oreades.** A. PacBio and Nanopore reads were assembled using [Canu v1.7](https://github.com/marbl/canu). To this initial assembly, 10X Chromium linked reads were mapped with [Longranger align](https://github.com/10XGenomics/longranger) and the assembly was polished with [Pilon](https://github.com/broadinstitute/pilon), before remapping of linked reads. B. [TigMint](https://github.com/bcgsc/tigmint) was used to create a bed file of inferred molecule spans of linked reads. The script convert_mol_cov.py was used to interpret this bed file. Probable misassemblies were broken manually within [Geneious](https://www.geneious.com/), and linked reads were remapped. C. We developed the script merge.py to infer linkage of contigs based on linked reads. This script was used in two iterations, with additional manual curation within Geneious and a final round of [Pilon](https://github.com/broadinstitute/pilon) polishing to obtain the final assembly. Raw sequencing reads are available at the Sequencing Read Archive under bioproject PRJNA525964 (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA525964/).

## Dependencies
- Python3
- [Pysam](https://github.com/pysam-developers/pysam)
- [Mummer](https://github.com/mummer4/mummer) (optional)

## Description
Contains the workhorse scripts **chromQC** and **merge**. 

chromQC is used for validation of a genome assembly. It works by taking a sliding window approach, where it collects the GEM barcodes in every window across the genome and compares them to every other window. Misassemblies can thus be spotted. 

merge is used for merging and scaffolding a genome assembly. It scans the ends of contigs in a given bam file, collects all GEM barcodes that align to each contig end, and compares the barcodes in all contig ends. A graph is then built, where contig ends with a higher-than-expected fraction of shared barcodes are connected with an edge. Unique paths through the graph form linked contigs, and are then merged or scaffolded. First, contig ends without read coverage are trimmed off. Mummer is then used to try to find an overlap between the contig ends. If an overlap is found, the contigs are fused together. If no overlap is found, the contigs are scaffolded by inserting N*10.

AnVIL requires as input **_only_** a bam file as produced by Longranger align.

### Example usage:
```
./merge.py <in.bam> # Creates the graph and outputs it as a gfa
./merge.py -i <in.fasta> <in.bam> # Creates the graph and uses it to merge the contigs in <in.fasta>
```

## Output
merge will always produce a gfa file for visualization of the graph with e.g. Bandage. It can also optionally be given a fasta file, in which case it will also produce a merged/scaffolded new fasta file and a bed file where the features are the merged contigs for visualization. To merge and scaffold, the given fasta file must be the same as the one used as a reference for the read mapping, and nucmer is needed in $PATH.

chromQC will produce a plain text file containing a report of potential misassemblies, which can be used to manually split these contigs.

## Known limitations
AnVIL is highly dependent on high-quality read mappings. It is highly recommended to run Pilon before AnVIL. As always, collapsed repeats are difficult to handle, and may appear to link together unrelated contigs. Short contigs made up entirely of repeats will be impossible to link because of the lack of seeding regions.