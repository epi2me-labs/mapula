import os
import sys
import tqdm
import pysam
import argparse
import pandas as pd
from pysam import AlignmentFile
from typing import Union, Dict, List
from mapula.lib.bio import get_alignment_tag
from mapula.lib.observation import ObservationGroup, TrackedReference
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


class CountMappingStats(object):

    def __init__(
        self,
        sams: List[str],
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

        self.sams = sams
        self.splitby = splitby
        self.output_sam = output_sam
        self.output_name = output_name
        self.output_format = output_format
        self.output_corrs = bool(counts)

        self.reference_files = parse_cli_key_value_pairs(refs)
        self.counts_files = parse_cli_key_value_pairs(counts)

        self.total_alignments = self.get_total_count(sams)
        self.alignment_files = self.get_alignment_files(sams)

        outfile = None
        if output_sam:
            outfile = AlignmentFile(output_sam, "w", 
                template=self.alignment_files[0]
            )

        errprint("[1/3] Loading references")
        self.references = self.get_references(
            self.reference_files,
            self.counts_files
        )

        errprint("[2/3] Parsing alignments")
        self.observations = self.get_observations(
            self.splitby,
            self.counts_files,
            self.total_alignments,
            self.alignment_files,
            outfile,
            self.references,
        )

        if not self.observations:
            errprint(
                "No observations made, is the "
                "BAM file empty?")
            sys.exit(0)

        errprint("[3/3] Writing data")
        self.write_observations(
            self.output_name,
            self.output_format,
            self.observations
        )

    @staticmethod
    def get_total_count(
        paths: List[str]
    ) -> Union[int, None]:
        """
        Iterates over the supplied SAM/BAM paths, and
        totals their counts, if and only if none of
        the sources are from stdin, in which case
        returns None.
        """
        if sys.stdin in paths:
            return None

        total = 0
        for p in paths:
            total += get_total_alignments(p)

        return total

    @staticmethod
    def get_alignment_files(
        paths: List[str]
    ) -> List[AlignmentFile]:
        """
        Iterates over the paths to SAM/BAM files given
        and created pysam AlignmentFile objects for each,
        returning the objects as a list.
        """
        alignment_files = []
        for path in paths:
            try:
                alignment_files.append(AlignmentFile(path, "r"))
            except OSError:
                errprint(
                    "[Error]: File not found: {}".format(path))
                sys.exit(1)

        return alignment_files

    @staticmethod
    def get_references(
        reference_files: Dict[str, str],
        counts_files: Dict[str, str]
    ) -> Dict[str, TrackedReference]:
        """
        Iterates over the input reference and counts 
        files, and collects all of the contigs (a.k.a. 
        sequences, chromosomes) into a dictionary 
        structure keyed by contig name, with values as
        TrackedReference objects, containing useful 
        information such as containing filename, length 
        and expected_count if available.
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

    def get_observations(
        self,
        splitby: List[str],
        counts_files: Dict[str, str],
        total_records: Union[int, None],
        alignment_files: List[AlignmentFile],
        outfile: Union[None, pysam.AlignmentFile],
        tracked_references: Dict[str, TrackedReference]
    ) -> Dict[str, ObservationGroup]:
        """
        Iterates over the input alignments and
        groups them by the splitby criteria, creating
        for each group an ObservationGroup object,
        which calculates and stores alignment stats. 
        """
        ticks = tqdm.tqdm(
            total=total_records, leave=False) or None

        observations = {}
        for alignment_file in alignment_files:
            for aln in alignment_file.fetch(until_eof=True):

                if total_records:
                    ticks.update(1)

                if outfile:
                    outfile.write(aln)

                self._update_observations(
                    aln, splitby, observations,
                    counts_files, tracked_references
                )

        return observations

    @staticmethod
    def _update_observations(
        aln: pysam.AlignedSegment,
        splitby: List[str],
        observations: Dict[str, ObservationGroup],
        counts_files: Dict[str, str],
        tracked_references: Dict[str, TrackedReference]
    ) -> Dict[str, ObservationGroup]:
        """
        Accepts a dictionary containing ObservationsGroups
        and updates it with new information from the given
        AlignedSegment. May add new observations to the dict,
        or update an existing one.
        """
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

        obs_name = "-".join(identity[i] for i in splitby)
        obs_ident = {m: identity[m] for m in splitby}

        matching_obs = observations.get(obs_name)
        if not matching_obs:
            matching_obs = ObservationGroup(obs_name, obs_ident)
            observations[obs_name] = matching_obs

            # TODO: Improve this approach
            # At the moment, the only use for tracking
            # references within Observations is to be
            # able to calculate correlations when counts
            # are provided via CSV. This is only possible
            # in the situaton where the mask includes "group"
            # i.e. the reference .fasta in which the current
            # observed reference sequence resides.
            if "group" in splitby:
                matching_obs.has_counts = bool(
                    "reference" not in splitby
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

    def write_observations(
        self,
        name: str,
        formats: str,
        observations: Dict[str, ObservationGroup]
    ) -> None:
        """
        Determines the requested output format, fields 
        and formatting options and then produces the 
        output from the input observations.
        """
        if formats in [Format.ALL, Format.JSON]:
            self.write_stats_to_json(
                '{}.json'.format(name),
                observations,
                output_corrs=self.output_corrs
            )

        if formats in [Format.ALL, Format.CSV]:
            self.write_stats_to_csv(
                '{}.csv'.format(name),
                observations,
                mask=DISTS,
                sort=['observations', 'base_pairs'],
                round_fields={
                    "cov80_percent": 2, "spearmans_rho": 2,
                    "spearmans_rho_pval": 2, "pearson": 2,
                    "pearson_pval": 2
                },
                output_corrs=self.output_corrs
            )

    @staticmethod
    def write_stats_to_json(
        path: str,
        observations: Dict[str, ObservationGroup],
        **kwargs
    ) -> None:
        """
        Iterates over the input ObservationGroup
        objects and transforms into them to a 
        serialisable format, combines the data, 
        and then writes it out to file. Supplied
        **kwargs are passed through to the to_dict
        method of ObservationGroup.
        """
        output = {
            key: val.to_dict(**kwargs)
            for key, val in observations.items()
        }
        write_data(path, output)

    @staticmethod
    def write_stats_to_csv(
        path: str,
        observations: Dict[str, ObservationGroup],
        mask: List[str],
        sort: List[str],
        round_fields: Dict[str, int],
        **kwargs
    ) -> None:
        """
        Iterates over the input ObservationGroup
        objects and forms a row in a dataframe for
        each one. Omitting whichever fields are in
        the mask list, e.g. frequency dists. The
        dataframe is then sorted, columns rounded, 
        and finally written out to file in .csv 
        format. Supplied **kwargs are passed through 
        to the to_dict method of ObservationGroup.
        """
        output = []
        for val in observations.values():
            output.append({
                i: j for i, j in val.to_dict(**kwargs).items()
                if i not in mask
            })

        df = pd.DataFrame(output)
        df = df.sort_values(sort, ascending=False)
        df = df.reset_index(drop=True)
        df = df.round(round_fields)
        df.to_csv(path, index=False)

    @classmethod
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
                "Input alignments in SAM format. (Default: [stdin])."
            ),
            dest="sams",
            default=[sys.stdin],
            metavar='',
            nargs='*',
        )

        parser.add_argument(
            '-r',
            help=(
                "Reference .fasta file(s). Format name=path_to_ref."
            ),
            dest="refs",
            required=True,
            metavar='',
            nargs='*'
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
            sams=args.sams,
            refs=args.refs,
            counts=args.counts,
            splitby=args.splitby,
            output_sam=args.output_sam,
            output_name=args.name,
            output_format=args.format,
        )
