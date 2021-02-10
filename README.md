![Oxford Nanopore Technologies logo](https://github.com/epi2me-labs/mapping-stats/raw/master/images/ONT_logo_590x106.png)


# Mapula

This package provides a command line tool that is able to parse alignments in
`SAM` format and produce a range of useful stats.

Mapula provides several subcommands, use `--help` with each one to find
detailed usage instructions.

## Installation

(Count) Mapula can be installed in the usual Python tradition:

```
pip install mapula
```


## Usage: count

```
$ mapula count -h
usage: mapula [-h] [-s SAM] [-o OUT] [-j JSON] -r [REFS [REFS ...]] [-e EXP]

Counts mapping stats from a SAM/BAM file

optional arguments:
  -h, --help            show this help message and exit
  -s SAM, --sam SAM     Input sam/bam file. (default: stdin)
  -o OUT, --out OUT     Output sam file. Use - for piping out. (default: None)
  -j JSON, --json JSON  Name of the output json (default: stats.mapula.json)
  -r [REFS [REFS ...]], --refs [REFS [REFS ...]]
                        List of references to which alignments have been made
  -e EXP, --exp EXP     A .CSV file containing expected counts by reference name
```
---


### Inputs

An example invocation is as follows:

```
mapula gather -s aligned.sam -r reference_1.fasta reference_2.fasta
```

An explanation of the available arguments follows:
- `-s` SAM/BAM
  - A path to a `SAM` or `BAM` file. Alignments will be parsed from this file and used to derive stats. Mapula also supports piping in, in which case this argument must be omitted e.g.:
    - ```minimap2 -y -ax map-ont ref.fasta *.fastq | mapula gather -r...etc```
- `-r` references
  - You must supply `mapula gather` reference `.fasta` files against which the alignments were made, this can take the form:
    - ```mapula gather -r reference_1.fasta reference_2.fasta```
- `-e` expected counts csv
  - In order to calculate correlations between observed and expected counts for a given set of reference sequences, you must use `-e` to provide a `.csv` containing two columns: 
    - `reference` (i.e. name of the reference sequence)
    - `expected_count`
- `-j` json_path
  - If provided, the `.json` output will be written at this location instead of `stats.mapula.json`
- `-o` sam_out
  - If provided, a `SAM` file constructed from the input records will be written to stdout. This argument is primarily used for piping, as in:
    - ```... | mapula gather -o | samtools sort - > sorted.bam```


### Stats & Terminology
For each alignment processed, `mapula count` updates various measurements.

#### Simple metrics
- alignment_count
- read_count
- primary_count
- secondary_count
- supplementary_count
- total_base_pairs

#### Distributions
- avg. alignment accuracy
- avg. read quality
- avg. read length
- reference coverage

#### Derived
- read n50
- median accuracy
- median quality
- cov80_count
- cov80_percent

Each of these stats are tracked at three levels:

1) **Global**: overall stats created from all alignments processed
2) **Group**: stats binned by group, i.e. run_id, barcode and reference file name
3) **Reference**: stats for every reference observed within a group

In addition, at the **global** and **group** levels, we also track correlations and their p_values:

- spearmans
- spearmans_p
- pearsons
- pearsons_p

By default these correlations will be 0, unless a `.csv` containing expected counts is provided using `-e`.



### Outputs
Mapula gather writes out several outputs by default.

#### JSON
By default, a `.json` file is produced, which has a nested structure, as per the levels described above:
```
# Top level
{
    ...global_stats
    children: {
      [group_name]: {
        ...group_stats,
        children {
          [reference_name]: {
            ...reference_stats
          },
          ...other_references
        }
      },
      ...other_groups
    }
}

```
The default filename of the `.json` file is `stats.mapula.json`.

The `.json` file is designed to support detailed downstream analysis. It is possible to disable creating it, however, if it is uneeded.

### CSV
Also by default, a set of `.csv` files are created which provide a more minimal representation of the stats collected at the 3 different levels.

By default, they are named:

1) `global.mapula.csv`
2) `groups.mapula.csv`
3) `reference.mapula.csv`

spearmans_rho_pvalue
> The p-value for the spearmans_rho statistic, it is the probability that the same correlation could have occurred by chance. Lower is better.

They contain the same overall stats as the `.json` file, but without the inclusion of the frequency distributions for accuracy, coverage, read length and read quality. However, the stats derived from these distributions, i.e. read n50, median accuracy, median quality and cov80 are retained.


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
