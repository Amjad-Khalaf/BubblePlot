# BubblePlot
Quick coverage vs GC content plots for genome assembles, with contigs coloured by likely taxonomic assignment based on BUSCO genes. Inspired by https://blobtoolkit.genomehubs.org/. **This is meant to be for rough and quick analyses**.

## About - what do I need?

<div align="justify">
  
Bubbleplot does not do any of the heavy lifting, you will unfortunately (for now) have to generate all the necessary files you need. These include: the `samtools index .fai file`, a `windowed coverage counts bed file` (example below in first code block), a tab-separated `gc_content` file (example below in second code block) and as many BUSCO `full_table.tsv` files as you want to plot. Briefly, BubblePlot will "assign" a contig to a taxon if it contains the most BUSCO genes from that taxon's `full_table.tsv` file.

```
ptg000001l      0       10000   38
ptg000001l      10000   20000   61
ptg000001l      20000   30000   72
```

```
ptg000001l      35.46
ptg000002l      35.36
```

</div>

## Installation

<div align="justify">

Installation is simple and straightforward, you need an environment with Python 3 and the following dependencies only:

```
argparse
sys
matplotlib
statistics
```

As for running, all you need is the python script in this repository.

</div>

## Quick - How do I use it?!

<div align="justify">
  
```
python bubbleplot.py [-h] --samtools_index_file $SAMTOOLS_FAI_FILE --coverage_bed_file $COVERAGE_BED_FILE --gc_content_file $GC_CONTENT_FILE --busco_full_tables "$BUSCO_TABLE_1,$BUSCO_TABLE_2" --busco_labels "$TAXON_1,$TAXON2"

```

For the full range of options, see:

```
usage: bubbleplot.py [-h] --samtools_index_file FILE
                                  --coverage_bed_file FILE --gc_content_file
                                  FILE --busco_full_tables LIST
                                  [--busco_labels LIST] [--plot_title STR]
                                  [--output_prefix STR] [--min_gc FLOAT]
                                  [--max_gc FLOAT] [--min_cov FLOAT]
                                  [--max_cov FLOAT]

options:
  -h, --help            show this help message and exit
  --samtools_index_file FILE
                        Samtools .fai file of fasta
  --coverage_bed_file FILE
                        BED file with average read depth per contig
  --gc_content_file FILE
                        GC content file generated using seqkit fx2tab -g -n
  --busco_full_tables LIST
                        Comma-separated list of BUSCO full table files
  --busco_labels LIST   Comma-separated labels for BUSCO tables
  --plot_title STR
  --output_prefix STR
  --min_gc FLOAT
  --max_gc FLOAT
  --min_cov FLOAT
  --max_cov FLOAT
```

</div>

## Example output

<div align="justify">

The figure output is a BlobToolKit style plot, with the rough taxonomic assignments relying on BUSCO genes.

<p align="center">
<img src="https://github.com/user-attachments/assets/f923b59b-5391-422c-bd24-def0343733b6" width="900">
</p>

The text output is a tab-separated file called `contig_list.tsv`, which contains columns in the following order: contigs, average coverage, average gc content, length, proposed assignment based on supplied BUSCO full tables.

```
ptg000001l      51.464233576642336      35.46   6846818.0       Insecta
ptg000002l      58.0204241948154        35.36   12721637.0      Insecta
ptg000043l      268.6969696969697       27.27   1317326.0       Microsporidia
```

</div>
