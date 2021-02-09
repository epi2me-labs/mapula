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
   ercc         Same as gather, but also calculate ERCC correlations
   aggregate    Combine mapping stats .JSON outputs

Collects stats from SAM/BAM files

positional arguments:
  command     Subcommand to run

optional arguments:
  -h, --help  show this help message and exit
```


## Stats
---
The `gather` and `ercc` commands both collect broadly the same set of stats from the alignments they parse. Here is a brief description of what each one means.
```
group
> ERCC, Host or unmapped, this specifies which references we are talking about

run_id
> Can be a valid ONT run_id or unknown, if that information is unavailable

barcode
> The barcode detected during demultiplexing

base_pairs
> The sum of the read lengths which have been aligned to this group

observed_count
> The number of alignments made to this group. This should be high for ERCC and Host and the barcode you used, and low for all others.

unique_observed_references
> The number of references within a group (e.g. ERCC contains 92 references) which have attracted alignments. This should be high for ERCC, we expect in most cases (except very low throughput runs) to recover most of the transcripts.

coverage80_count
> The number of alignments which have covered 80% or more of the length of their alignment target. This should be high for ERCC, since the transcripts are short, and practically non-existent for Human since the targets are chromosome length.

coverage80_percent
> The same as the above, but shown as a proportion of alignments in that group. E.g. 50% of alignments have spanned 80% or more of the target.

median_accuracy
> This is the average alignment accuracy for alignments made within the group. Given ONTs error profile, you should see values in the range of 91% - 95%. Lower than 91% should be a red flag and require further investigation.

median_quality
> This is the average per-read Phred quality score for reads which have been aligned within the group. This should be similar across groups.

spearmans_rho
> The correlation coefficient, i.e. measuring the strength of association between the observed counts and the expected counts for this group. Since only the ERCC group will have expected counts attached, this measure should be ignored for all non-ERCCs. Higher is better, and we expect to see a correlation coefficient of 99.5+ under normal circumstances.

spearmans_rho_pvalue
> The p-value for the spearmans_rho statistic, it is the probability that the same correlation could have occurred by chance. Lower is better.