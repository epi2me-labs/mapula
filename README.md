![Oxford Nanopore Technologies logo](https://github.com/epi2me-labs/mapula/raw/master/images/ONT_logo_590x106.png)


# Mapula

This package provides a command line tool that is able to parse alignments in `SAM` format and produce a range of useful stats.

Mapula provides several subcommands, use `--help` with each
one to find detailed usage instructions.

## Installation
Count mapula can be installed following the usual Python tradition:
```
pip install mapula
```
## Usage: count
```
$ mapula count -h
usage: mapula [-h] [-s SAM] [-r [REFS [REFS ...]]] [-c [COUNTS [COUNTS ...]]] [-o SAM_OUT] [-j JSON_PATH]

Count mapping stats from a SAM/BAM file

optional arguments:
  -h, --help            show this help message and exit
  -s SAM, --sam SAM     Alignments in SAM format. By default, this script reads alignments from stdin. However, using this flag it is possible to pass in a file path.
  -r [REFS [REFS ...]], --refs [REFS [REFS ...]]
                        Provide reference .fasta files using the syntax: name=path_to_ref.
  -c [COUNTS [COUNTS ...]], --counts [COUNTS [COUNTS ...]]
                        Provide expected counts in csv format using the syntax: name=path_to_counts, where name should be equal to a name given to --refs.
  -o SAM_OUT, --sam_out SAM_OUT
                        Outputs a sam file from the parsed alignments. Use - for piping out. (default: None)
  -j JSON_PATH, --json_path JSON_PATH
                        Name of the output json (default: stats.mapula.json)
```

An example invocation is as follows:

```
mapula gather -s aligned.sam -r host=reference_1.fasta spikein=reference_2.fasta -c spikein=counts.csv
```

### Expected counts
The expected counts CSVs should have the following column headings: 

`reference`, `expected_count` 

The reference column should contain the ID of a sequence in the corresponding reference file. The expected_count column should equal the expected number of observations for that sequence.

---
### **Stats & Terminology**
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

Each of these stats are tracked at two levels:

2) **Group**: stats binned by group, i.e. run_id, barcode and reference file name
3) **Reference**: stats for every reference observed within a group

In addition, at the **group** levels, we also track correlations and their p_values:

- spearmans
- spearmans_p
- pearsons
- pearsons_p

By default these correlations will be 0, unless a `.csv` containing expected counts is provided using `-e`.

---
### **Outputs**
Mapula gather writes out several outputs by default.

#### JSON
By default, a `.json` file is produced, which has a nested structure, as per the levels described above:
```
# Top level
{
    [group_name]: {
      ...group_stats,
      references {
        [reference_name]: {
          ...reference_stats
        },
        ...other_references
      }
    },
    ...other_groups
}

```
The default filename of the `.json` file is `stats.mapula.json`.

The `.json` file is designed to support detailed downstream analysis. It is possible to disable creating it, however, if it is uneeded.

### CSV
Also by default, a set of `.csv` files are created which provide a more minimal representation of the stats collected at the 2 different levels.

By default, they are named:

1) `groups.mapula.csv`
2) `refs.mapula.csv`

They contain the same overall stats as the `.json` file, but without the inclusion of the frequency distributions for accuracy, coverage, read length and read quality. However, the stats derived from these distributions, i.e. read n50, median accuracy, median quality and cov80 are retained.

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
