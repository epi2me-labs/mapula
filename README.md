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
## Usage

At present there is only one subcommand: count (but with more to follow soon!).

```
$ mapula count -h
usage: mapula [-h] -r [[...]] [-c [[...]]] [-o] [-f] [-s  [...]] [-n]

Count mapping stats from a SAM/BAM file

positional arguments:
                Input alignments in SAM format. (Default: stdin).

optional arguments:
  -h, --help    show this help message and exit
  -r [ [ ...]]  Reference .fasta file(s). Format name=path_to_ref.
  -c [ [ ...]]  Expected counts CSV(s). Format name=path_to_counts. Expected columns: reference,expected_count.
  -o            Outputs a sam file from the parsed alignments. Use - for piping out. (Default: None).
  -f            Sets the format(s) in which to output results. [Choices: csv, json, all] (Default: csv).
  -s  [ ...]    Split by these criteria, comma separated. [Choices: group,run_id,barcode,read_group,reference] (Default: group,run_id,barcode).
  -n            Prefix of the output files, if there are any.
```

An example invocation is as follows:

```
mapula count <path_to_sam_or_bam> -r <name>=<path_to_a_reference_fasta>
```

In the above example, name should be a short-hand nickname for the group of reference sequences in your fasta. E.g. `host`. You may provide multiple values for `-r`.

You may also split the alignments by different criteria using `-s`, e.g. `-s read_group reference` to
group the alignments into rows based on their read_group and aligned reference.

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
