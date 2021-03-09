![Oxford Nanopore Technologies logo](https://github.com/epi2me-labs/mapula/raw/master/images/ONT_logo_590x106.png)

# Mapula
This package provides a command line tool that is able to parse alignments in `SAM` format and produce a range of useful stats.

Mapula provides several subcommands, use `--help` with each
one to find detailed usage instructions.

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/mapula/README.html)

## Installation
Count mapula can be installed following the usual Python tradition:
```
pip install mapula
```
Or via conda:
```
conda install mapula
```

## Usage

```
$ mapula -h
usage: mapping-stats <command> [<args>]"

Available subcommands are:
    count        Count mapping stats from a SAM/BAM file
    merge        Combine mapula count's json outputs

Collects stats from SAM/BAM files

positional arguments:
  command     Subcommand to run

optional arguments:
  -h, --help  show this help message and exit
```
The main subcommand command is `mapula count`. This command can accept one or several input `SAM` or `BAM` files and outputs mapping statistics.

Alignments are grouped by user-specifiable criteria, `-s`. These aggregations can then be output in two formats using `-f`. The default `.csv` format is the most easily iterpretable for a quick glance, or for onward programmatic analysis the `.json` output contains a more in-depth view of the data.

### Examples

*Output some stats in `.csv` format containing mapping stats:*
```
mapula count <paths_to_sam_or_bam> -r <path_to_a_reference_fasta>
```

*Split stats only by `read_group` and `barcode`:*
```
mapula count <paths_to_sam_or_bam> -r <path_to_a_reference_fasta> -s barcode read_group
```

*Output some stats in both `.csv` and `.json` format:*
```
mapula count
  <paths_to_sam_or_bam> -f all -r <path_to_a_reference_fasta> <path_to_a_reference_fasta>
```

*Accept multiple alignment and reference inputs*
```
mapula count
  <paths_to_sam_or_bam> <paths_to_sam_or_bam> -r <path_to_a_reference_fasta> <path_to_a_reference_fasta>
```

*Receive some `SAM` or `BAM` from stdin, output stats in `.csv`, and pipe the `SAM` records onwards:*
```
minimap2 -y -ax map-ont <path_to_a_reference_fasta> *_reads.fastq \
  | mapula -r <path_to_a_reference_fasta> -f csv -p \
  | samtools sort -o sorted.aligned.bam
```

### Important: tags

At present, for access to `barcode`, `run_id`, `read_group`, mapula depends on tags being available within the input `SAM` records, as follows:
- `barcode` = `bc`
- `run_id` = `rd`
- `read_group` = `rg`

If these are not available, Mapula will just provide a placeholder of `Unknown`. The minimap2 flag `-y` can carry information from the `.fastq`
header into the records it creates to faciliate this transfer of information.

---

Help
----

**Licence and Copyright**

© 2021- Oxford Nanopore Technologies Ltd.

`mapula` is distributed under the terms of the Mozilla Public License 2.0.

**Research Release**

Research releases are provided as technology demonstrators to provide early
access to features or stimulate Community development of tools. Support for
this software will be minimal and is only provided directly by the developers.
Feature requests, improvements, and discussions are welcome and can be
implemented by forking and pull requests. However much as we would
like to rectify every issue and piece of feedback users may have, the
developers may have limited resource for support of this software. Research
releases may be unstable and subject to rapid iteration by Oxford Nanopore
Technologies.
