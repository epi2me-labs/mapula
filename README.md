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
The main command is `mapula count`. This command can accept one or several input `SAM` or `BAM` files and output various useful statistics.

For available options, see:
```
mapula count -h
```

**Per-alignment (outputs to stdout)**

If you have no need of aggregated stats, you can opt to pipe
per-alignment stats to stdout (WIP: currently only in json format). Or for instance if mapula is part of a chain of piping operators, you may relay the `SAM` records.

Pipe per-alignment stats to stdout in `.json` format:
```
mapula count <paths_to_sam_or_bam> -r <path_to_a_reference_fasta> -p json
```
Pipe `SAM` records from multiple files to stdout:
```
mapula count <paths_to_sam_or_bam> <paths_to_sam_or_bam> -p sam
```

**Aggregation (outputs to file)**

Aggregated stats can be calculated by enabling the `-a` flag. Alignments are grouped by user-specified criteria, `-s`. These aggregations can then be output 
in two formats using `-f`. The `.csv` format is the most easily iterpretable for a quick glance, or for onward programmatic analysis the `.json` output contains a more in-depth view of the data.

Output some stats in `.csv` format containing mapping stats:
```
mapula count <paths_to_sam_or_bam> -r <path_to_a_reference_fasta> -a
```

Output some stats in `.csv` format split by `read_group` and `barcode`:
```
mapula count <paths_to_sam_or_bam> -r <path_to_a_reference_fasta> -a -s barcode read_group
```

Output some stats in `.csv` and `.json` format split by `SAM` and reference `fasta`:
```
mapula count \
  <paths_to_sam_or_bam> <paths_to_sam_or_bam> \
  -r <path_to_a_reference_fasta> <path_to_a_reference_fasta> \
  -a \
  -s source fasta
```

Receive some `SAM` or `BAM` from stdin, calculate some aggregated stats and output in `.csv`, and pipe the `SAM` records onwards:
```
minimap2 -y -ax map-ont <path_to_a_reference_fasta> *_reads.fastq \
  | mapula -r <path_to_a_reference_fasta> -a source fasta run_id barcode -p sam \
  | samtools sort -o sorted.aligned.bam
```

## Important: tags

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
