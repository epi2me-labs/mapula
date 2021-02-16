import os
import sys
import tqdm
import pysam
import argparse
import numpy as np
import pandas as pd
from pysam import AlignmentFile
from dataclasses import dataclass
from typing import Union, Dict, List
from scipy.stats import pearsonr, spearmanr
from mapula.lib.const import (
    UNKNOWN,
    UNMAPPED,
    UNCLASSIFIED,
    DISTS,
    Groupers,
    Format
)
from mapula.lib.utils import (
    get_total_alignments,
    parse_cli_key_value_pairs,
    errprint,
    write_data,
)
from mapula.lib.bio import (
    BASE_ACCS,
    BASE_QUALS,
    get_alignment_tag,
    get_alignment_mean_qscore,
    get_median_from_frequency_dist,
    get_alignment_accuracy,
    get_n50_from_frequency_dist,
    get_alignment_coverage,
)


@dataclass
class TrackedReference(object):
    group: str
    filename: str
    length: int
    expected_count: int


class ObservationGroup():

    def __init__(
        self,
        # Ident
        name,
        identity,
        # Tracked
        base_pairs=0,
        observations=0,
        primary_count=0,
        secondary_count=0,
        supplementary_count=0,
        read_n50=0,
        cov80_count=0,
        cov80_percent=0,
        median_accuracy=0,
        median_quality=0,
        alignment_accuracies=None,
        alignment_coverages=None,
        aligned_qualities=None,
        read_lengths=None,
        # Optional (ident dependent)
        has_counts: bool = False,
        tracked_references: Dict[str, TrackedReference] = {},
        observed_references: Dict[str, int] = {}
    ) -> None:
        """
        Represents an instance of a reference sequence to
        which reads have been aligned, and tracks many
        with respect to the alignments.
        """
        self.name = name
        self.identity = identity

        self.base_pairs = base_pairs
        self.observations = observations
        self.primary_count = primary_count
        self.secondary_count = secondary_count
        self.supplementary_count = supplementary_count

        self.read_n50 = read_n50
        self.cov80_count = cov80_count
        self.cov80_percent = cov80_percent
        self.median_accuracy = median_accuracy
        self.median_quality = median_quality

        self.alignment_accuracies = (
            alignment_accuracies or list([0] * 1001))
        self.alignment_coverages = (
            alignment_coverages or list([0] * 101))
        self.aligned_qualities = (
            aligned_qualities or list([0] * 600))
        self.read_lengths = (
            read_lengths or list([0] * 1000))

        self.has_counts = has_counts
        self.spearmans_rho: float = 0
        self.spearmans_rho_pval: float = 0
        self.pearson: float = 0
        self.pearson_pval: float = 0

        self.observed_references = observed_references
        self.tracked_references = tracked_references

    #
    # Speedups
    #
    __slots__ = (
        'name', 'base_pairs', 'observations',
        'primary_count', 'secondary_count', 'supplementary_count',
        'cov80_count', 'cov80_percent', 'median_accuracy',
        'median_quality', 'alignment_accuracies', 'read_n50',
        'alignment_coverages', 'aligned_qualities', 'read_lengths',
        '__dict__', 'identity'
    )

    def update(
        self,
        aln: pysam.AlignedSegment,
    ) -> None:
        if aln.is_supplementary:
            self.supplementary_count += 1
            return

        if aln.is_secondary:
            self.secondary_count += 1
            return

        length = aln.query_length

        self.base_pairs += length
        self.observations += 1

        try:
            # E.g. LEN: 1000 bins, 50 bin width, 50,000 max length
            self.read_lengths[int(length / 50)] += 1
            self._update_read_n50()
        except IndexError:
            pass

        if aln.is_unmapped:
            return

        self.primary_count += 1
        reference = aln.reference_name

        try:
            self.observed_references[reference] += 1
        except KeyError:
            self.observed_references[reference] = 1

        try:
            quality = get_alignment_mean_qscore(
                aln.query_alignment_qualities)
            self.aligned_qualities[int(quality / 0.1)] += 1
            self._update_median_quality()
        except (IndexError, TypeError):
            pass

        try:
            accuracy = get_alignment_accuracy(aln) or 0
            self.alignment_accuracies[int(accuracy / 0.1)] += 1
            self._update_median_accuracy()
        except (IndexError, TypeError):
            pass

        try:
            coverage = get_alignment_coverage(
                aln.reference_length, length
            ) or 0
            self.alignment_coverages[int(coverage)] += 1

            if coverage >= 80:
                self.cov80_count += 1

            self._update_cov80_percent()
        except IndexError:
            pass

        if self.has_counts:
            self._update_correlations()

    def _update_read_n50(self):
        self.read_n50 = get_n50_from_frequency_dist(
            self.read_lengths, 50, self.base_pairs
        )

    def _update_median_accuracy(self):
        self.median_accuracy = get_median_from_frequency_dist(
            BASE_ACCS, np.array(self.alignment_accuracies)
        )

    def _update_cov80_percent(self):
        self.cov80_percent = (
            100 * self.cov80_count / self.observations
        )

    def _update_median_quality(self):
        self.median_quality = get_median_from_frequency_dist(
            BASE_QUALS, np.array(self.aligned_qualities)
        )

    def _update_correlations(self):
        observed = []
        expected = []
        for key, val in self.tracked_references.items():
            expected.append(val.expected_count)
            observed.append(self.observed_references.get(key, 0))

        coef, p = spearmanr(observed, expected)
        self.spearmans_rho = coef
        self.spearmans_rho_pval = p

        coef2, p2 = pearsonr(observed, expected)
        self.pearson = coef2
        self.pearson_pval = p2

    def __add__(self, new):
        if not isinstance(new, self.__class__):
            raise TypeError

        for attr in (
            "alignment_count",
            "read_count",
            "total_base_pairs",
            "primary_count",
            "secondary_count",
            "supplementary_count",
            "cov80_count"
        ):
            old = getattr(self, attr)
            new = getattr(self, attr)
            setattr(self, attr, old + new)

        self.alignment_accuracies = list(map(sum, zip(
            self.alignment_accuracies, new.alignment_accuracies)))
        self.alignment_coverages = list(map(sum, zip(
            self.alignment_coverages, new.alignment_coverages)))
        self.aligned_qualities = list(map(sum, zip(
            self.aligned_qualities, new.aligned_qualities)))
        self.read_lengths = list(map(sum, zip(
            self.read_lengths, new.read_lengths)))

        self._update_read_n50()
        self._update_median_quality()
        self._update_median_accuracy()
        self._update_cov80_percent()

        if self.has_counts:
            self.tracked_references.update(
                new.tracked_references)

            for k, v in new.observed_references.items():
                if self.observed_references.get(k):
                    self.tracked_references[k] += v
                    continue
                self.observed_references[k] = v

            self._update_correlations

        return self

    def to_dict(self, output_corrs=False):
        """
        Serialises this object to a dictionary
        format.
        """
        corrs = {
            "spearmans_rho": self.spearmans_rho,
            "spearmans_rho_pval": self.spearmans_rho_pval,
            "pearson": self.pearson,
            "pearson_pval": self.pearson_pval,
        } if output_corrs else {}

        return {
            **self.identity,
            "base_pairs": self.base_pairs,
            "observations": self.observations,
            "primary_count": self.primary_count,
            "secondary_count": self.secondary_count,
            "supplementary_count": self.supplementary_count,
            "read_n50": self.read_n50,
            "cov80_count": self.cov80_count,
            "cov80_percent": self.cov80_percent,
            "median_accuracy": self.median_accuracy,
            "median_quality": self.median_quality,
            "alignment_accuracies": self.alignment_accuracies,
            "alignment_coverages": self.alignment_coverages,
            "aligned_qualities": self.aligned_qualities,
            "read_lengths": self.read_lengths,
            **corrs
        }


class CountMappingStats(object):

    def __init__(
        self,
        sam: str,
        refs: Dict[str, str],
        counts: Dict[str, str],
        splitby: List[str],
        output_sam: Union[str, None],
        output_name: str,
        output_format: str
    ) -> None:
        """
        A subcommand that runs a process designed to
        scan alignments made in SAM format and accumulate
        many useful statistics which are binned into groups
        and reported in CSV or JSON format.
        """
        errprint("Running: Mapula (count)")

        self.sam = sam
        self.splitby = splitby
        self.output_sam = output_sam
        self.output_name = output_name
        self.output_format = output_format
        self.output_corrs = bool(counts)
        self.total_records = None

        self.reference_files = parse_cli_key_value_pairs(refs)
        self.counts_files = parse_cli_key_value_pairs(counts)

        try:
            if self.sam not in [sys.stdin]:
                self.total_records = get_total_alignments(self.sam)
            self.records = AlignmentFile(sam, "r")
        except OSError:
            errprint(
                "[Error]: Can't find SAM: {}".format(sam)
            )
            sys.exit(1)

        outfile = None
        if output_sam:
            outfile = AlignmentFile(
                output_sam, "w", template=self.records
            )

        errprint("[1/4] Loading references")
        self.tracked = self.get_tracked_references(
            self.reference_files,
            self.counts_files
        )

        errprint("[2/4] Parsing alignments")
        self.observed = self.get_observed_references(
            self.splitby,
            self.counts_files,
            self.total_records,
            outfile,
            self.tracked,
        )

        errprint("[3/4] Writing data")
        if self.output_format in [Format.ALL, Format.JSON]:
            self.write_stats_to_json(
                '{}.json'.format(output_name),
                self.observed,
                output_corrs=self.output_corrs
            )

        if self.output_format in [Format.ALL, Format.CSV]:
            self.write_stats_to_csv(
                '{}.csv'.format(output_name),
                self.observed,
                mask=DISTS,
                sort=['observations', 'base_pairs'],
                round_fields={
                    "cov80_percent": 2, "spearmans_rho": 2,
                    "spearmans_rho_pval": 2, "pearson": 2,
                    "pearson_pval": 2
                },
                output_corrs=self.output_corrs
            )

    def get_tracked_references(
        self,
        reference_files: Dict[str, str],
        counts_files: Dict[str, str]
    ) -> Dict:
        """
        """

        tracked = {}
        for name, path in reference_files.items():
            basename = os.path.basename(path)

            # See if we have counts for this fasta
            df_counts = None
            if matching := counts_files.get(name):
                df_counts = pd.read_csv(matching, index_col=0)

            # Open the fasta and read refs into a dict
            open_ref = pysam.FastaFile(path)
            for ref in open_ref.references:

                try:
                    # Collect expected count data if possible
                    count = df_counts.at[ref, 'expected_count']
                except (KeyError, AttributeError):
                    count = 0

                tracked[ref] = TrackedReference(
                    group=name,
                    filename=basename,
                    length=open_ref.get_reference_length(ref),
                    expected_count=count
                )

        # Also create an unmapped
        # row for binning unaligned reads
        tracked[UNMAPPED] = TrackedReference(
            group=UNMAPPED,
            filename=UNMAPPED,
            length=0,
            expected_count=0
        )

        return tracked

    def get_observed_references(
        self,
        mask: List[str],
        counts_files: Dict[str, str],
        total_records: Union[int, None],
        outfile: Union[None, pysam.AlignmentFile],
        tracked_references: Dict[str, TrackedReference]
    ) -> Dict:
        """
        """

        ticks = tqdm.tqdm(total=total_records, leave=False) or None

        observations = {}
        for aln in self.records.fetch(until_eof=True):

            if total_records:
                ticks.update(1)

            if outfile:
                outfile.write(aln)

            reference = aln.reference_name or UNMAPPED
            group = tracked_references[reference].group
            run_id = get_alignment_tag(aln, 'RD', default=UNKNOWN)
            read_group = get_alignment_tag(aln, 'RG', default=UNKNOWN)
            barcode = get_alignment_tag(aln, 'BC', default=UNCLASSIFIED)

            identity = {
                'group': group,
                'run_id': run_id,
                'barcode': barcode,
                'read_group': read_group,
                'reference': reference
            }

            obs_name = "-".join(identity[i] for i in mask)
            obs_ident = {m: identity[m] for m in mask}

            matching_obs = observations.get(obs_name)
            if not matching_obs:
                matching_obs = ObservationGroup(obs_name, obs_ident)
                observations[obs_name] = matching_obs

                # At the moment, the only use for tracking
                # references within Observations is to be
                # able to calculate correlations when counts
                # are provided via CSV. This is only possible
                # in the situaton where the mask includes "group"
                # i.e. the reference .fasta in which the current
                # observed reference sequence resides.
                if "group" in mask:
                    matching_obs.has_counts = bool(
                        "reference" not in mask
                        and group in counts_files
                    )

                    if matching_obs.has_counts:
                        matching_obs.tracked_references = ({
                            trkname: trkref
                            for trkname, trkref
                            in tracked_references.items()
                            if trkref.group == group
                        })

            matching_obs.update(aln)

        return observations

    def write_stats_to_json(
        self,
        path: str,
        data: Dict,
        **kwargs
    ) -> None:
        """
        """

        output = {
            key: val.to_dict(**kwargs)
            for key, val in data.items()
        }
        write_data(path, output)

    def write_stats_to_csv(
        self,
        path: str,
        data: Dict,
        mask: List[str],
        sort: List[str],
        round_fields: Dict[str, int],
        **kwargs
    ) -> None:
        """
        """

        output = []
        for val in data.values():
            output.append({
                i: j for i, j in val.to_dict(**kwargs).items()
                if i not in mask
            })

        df = pd.DataFrame(output)
        df = df.sort_values(sort, ascending=False)
        df = df.reset_index(drop=True)
        df = df.round(round_fields)
        df.to_csv(path, index=False)

    @ classmethod
    def execute(cls, argv) -> None:
        """
        Parses command line arguments and
        initialises a CountMappingStats object.
        """
        parser = argparse.ArgumentParser(
            description="Count mapping stats from a SAM/BAM file",
        )

        parser.add_argument(
            help=(
                "Input alignments in SAM format. (Default: stdin)."
            ),
            dest="sam",
            default=sys.stdin,
            metavar='',
            nargs='?'
        )

        parser.add_argument(
            '-r',
            help=(
                "Reference .fasta file(s). Format name=path_to_ref."
            ),
            dest="refs",
            required=True,
            metavar='',
            nargs='*',
        )

        parser.add_argument(
            '-c',
            help=(
                "Expected counts CSV(s). Format name=path_to_counts. "
                "Expected columns: reference,expected_count."
            ),
            dest="counts",
            default=[],
            required=False,
            metavar='',
            nargs='*'
        )

        parser.add_argument(
            '-o',
            dest="output_sam",
            default=None,
            required=False,
            metavar='',
            help=(
                "Outputs a sam file from the parsed alignments. "
                "Use - for piping out. (Default: None)."
            ),
        )

        parser.add_argument(
            "-f",
            dest="format",
            required=False,
            default=Format.CSV,
            choices=Format.choices,
            metavar='',
            help=(
                "Sets the format(s) in which to output results. "
                "[Choices: csv, json, all] (Default: csv)."
            )
        )

        parser.add_argument(
            "-s",
            dest="splitby",
            required=False,
            default=[Groupers.GROUP, Groupers.RUN_ID, Groupers.BARCODE],
            choices=Groupers.choices,
            metavar='',
            nargs="+",
            help=(
                "Split by these criteria, space separated. "
                "[Choices: {}] (Default: group run_id "
                "barcode).".format(' '.join(Groupers.choices))
            )
        )

        parser.add_argument(
            "-n",
            dest="name",
            required=False,
            default="mapula",
            metavar='',
            help=(
                "Prefix of the output files, if there are any."
            )
        )

        args = parser.parse_args(argv)

        if args.counts and not (Groupers.GROUP in args.splitby):
            errprint(
                "If you want to provide expected counts, you "
                "are required to split by group. ( -s group )"
            )
            sys.exit(1)

        cls(
            sam=args.sam,
            refs=args.refs,
            counts=args.counts,
            splitby=args.splitby,
            output_sam=args.output_sam,
            output_name=args.name,
            output_format=args.format,
        )
