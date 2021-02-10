# Mapping Stats

This package is designed to provide a command line tool that is able to
parse alignments in `SAM` format and from them produce a range of useful
stats.

As they are parsed, alignments are binned into groups based on their `run_id`,
the name of the file containing the reference sequence they aligned to (or not, as in
unmapped reads) and their assigned barcode (if that information is available).

The statistics are then calculated per-group, as well as per-reference within
each group, and are output in `CSV` as well as JSON format.

In addition, the `ercc` subcommand provides the same functionality, but with
the added requirement that a `CSV` of expected counts will be provided. With this,
correlations between the expected and observed amounts of each ERCC reference
can be determined, which is useful in control experiments.


## Usage
---
Mapping Stats provides several subcommands, use `--help` with each
one to find detailed usage instructions.

```
$ mapping_stats -h
usage: mapping-stats <command> [<args>]

Available subcommands are:
   gather       Gather mapping stats from a SAM/BAM file
   aggregate    Combine mapping stats .JSON outputs

Collects stats from SAM/BAM files

positional arguments:
  command     Subcommand to run

optional arguments:
  -h, --help  show this help message and exit
```