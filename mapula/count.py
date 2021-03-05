import os
import sys
import tqdm
import json
import pysam
import argparse
import pandas as pd
from pysam import AlignmentFile
from typing import Union, Dict, List, TextIO, Tuple
from mapula.observation import (
    ObservationGroup, TrackedReference, Observation)
from mapula.const import (
    UNMAPPED,
    Pipers,
    Groupers,
    Format
)
from mapula.utils import (
    get_total_alignments,
    errprint,
    write_data,
)


class CountMappingStats(object):

    def __init__(
        self,
        sams: List[str],
        refs: List[str],
        counts: str,
        pipe: Union[str, None],
        aggregate: bool,
        splitby: List[str],
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
        self.pipe = pipe
        self.aggregate = aggregate
        self.splitby = splitby
        self.output_name = output_name
        self.output_format = output_format
        steps = 3 if self.aggregate else 2

        self.reference_files = refs
        self.counts_file = counts

        self.total_alignments = self.get_total_count(sams)
        self.alignment_files = self.get_alignment_files(sams)

        #
        # Todo: Refactor into function call, de-uglify, fix headers
        #
        self.outfile = None
        if pipe == Pipers.SAM:
            self.outfile = AlignmentFile(
                '-', "w", template=self.alignment_files[0][1])
        elif pipe == Pipers.JSON:
            self.outfile = sys.stdout

        errprint(f"[1/{steps}] Loading references")
        self.references = self.get_references(self.reference_files)
        self.counts = self.get_expected_counts(self.counts_file)

        #
        # Todo: Refactor to make use of coroutines so
        # complexity is distributed across functions
        #
        errprint(f"[2/{steps}] Parsing alignments")
        self.aggregations = self.parse_alignments(
            self.alignment_files,
            self.pipe,
            self.aggregate,
            self.splitby,
            self.total_alignments,
            self.outfile,
            self.references,
            self.counts
        )

        if not self.aggregate:
            return

        if self.aggregate and not self.aggregations:
            errprint(
                "No aggregations made, is the "
                "BAM file empty?")
            sys.exit(0)

        errprint(f"[3/{steps}] Writing aggregations")
        self.write_aggregations(
            self.output_name,
            self.output_format,
            self.aggregations
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

        for p in paths:
            if os.stat(p).st_size > 1000000000:
                return None

        total = 0
        for p in paths:
            total += get_total_alignments(p)

        return total

    @staticmethod
    def get_alignment_files(
        paths: List[str]
    ) -> List[Tuple[str, AlignmentFile]]:
        """
        Iterates over the paths to SAM/BAM files given
        and created pysam AlignmentFile objects for each,
        returning the objects as a list.
        """
        alignment_files = []
        for path in paths:
            try:
                name = (
                    os.path.basename(path) 
                    if path != sys.stdin else "stdin"
                )
                alignment_files.append((
                    name,
                    AlignmentFile(path, "r")
                ))
            except OSError:
                errprint(
                    f"[Error]: File not found: {path}")
                sys.exit(1)

        return alignment_files

    @staticmethod
    def get_references(
        reference_files: List[str],
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
        tracked = {
            UNMAPPED: TrackedReference(
                name=UNMAPPED,
                fasta=UNMAPPED,
                length=0
            )
        }

        for path in reference_files:
            basename = os.path.basename(path)

            try:
                open_ref = pysam.FastaFile(path)
            except FileNotFoundError:
                errprint(f'[Error]: {path} does not exist.')
                sys.exit(1)

            for ref in open_ref.references:
                tracked[ref] = TrackedReference(
                    name=ref,
                    fasta=basename,
                    length=open_ref.get_reference_length(ref),
                )

        return tracked

    @staticmethod
    def get_expected_counts(
        counts_file: Union[None, str]
    ) -> Union[None, Dict[str, int]]:
        """
        """
        if not counts_file:
            return None

        if not os.path.exists(counts_file):
            errprint(
                "[Error]: Supplied counts file does not "
                "exist. Exiting."
            )
            sys.exit(1)

        counts = pd.read_csv(counts_file)
        counts.columns = map(str.lower, counts.columns)

        if not all(
            col in counts.columns
            for col in ['reference', 'expected_count']
        ):
            errprint(
                "[Error]: Supplied counts file does not "
                "contain the required columns, "
                "'reference,expected_count'"
            )
            sys.exit(1)

        return counts.to_dict()['reference']

    def parse_alignments(
        self,
        alignments_files: List[Tuple[str, AlignmentFile]],
        pipe: Union[str, None],
        aggregate: bool,
        splitby: List[str],
        total_records: Union[int, None],
        outfile: Union[None, TextIO, pysam.AlignmentFile],
        references: Dict[str, TrackedReference],
        counts: Union[None, Dict[str, int]]
    ) -> Dict[str, ObservationGroup]:
        """
        Iterates over the input alignments and
        groups them by the splitby criteria, creating
        for each group an ObservationGroup object,
        which calculates and stores alignment stats. 
        """
        ticks = tqdm.tqdm(total=total_records or None, leave=False)
        stream_sam = bool(pipe == Pipers.SAM)
        stream_json = bool(pipe == Pipers.JSON)

        aggregations = {}
        for alname, alhandle in alignments_files:
            for aln in alhandle.fetch(until_eof=True):

                obs = Observation.from_alignment(
                    aln, alname, references, stream_json)

                if stream_sam:
                    outfile.write(aln)
                elif stream_json:
                    outfile.write(f'{json.dumps(obs.__dict__)}\n')

                if aggregate:
                    self._update_aggregations(
                        obs, splitby, aggregations, counts)

                ticks.update(1)

        if aggregate:
            for grp in aggregations.values():
                grp._update_summary_stats()

        return aggregations

    @staticmethod
    def _update_aggregations(
        obs: Observation,
        splitby: List[str],
        aggregations: Dict[str, ObservationGroup],
        counts: Union[None, Dict[str, int]]
    ) -> Dict[str, ObservationGroup]:
        """
        Accepts a dictionary containing ObservationsGroups
        and updates it with new information from the given
        AlignedSegment. May add new observations to the dict,
        or update an existing one.
        """
        identity = {
            'source': obs.source,
            'fasta': obs.fasta,
            'run_id': obs.run_id,
            'barcode': obs.barcode,
            'read_group': obs.read_group,
            'reference': obs.reference
        }

        obs_name = "-".join(identity[i] for i in splitby)
        obs_ident = {m: identity[m] for m in splitby}

        matching_obs = aggregations.get(obs_name)
        if not matching_obs:
            matching_obs = ObservationGroup(
                obs_name, obs_ident, counts=counts)
            aggregations[obs_name] = matching_obs

        matching_obs.update(obs, update_summary_stats=False)
        return aggregations

    def write_aggregations(
        self,
        name: str,
        formats: str,
        aggregations: Dict[str, ObservationGroup]
    ) -> None:
        """
        Determines the requested output format, fields 
        and formatting options and then produces the 
        output from the input observations.
        """
        if formats in [Format.ALL, Format.JSON]:
            self.write_stats_to_json(
                '{}.json'.format(name),
                aggregations,
            )

        if formats in [Format.ALL, Format.CSV]:
            self.write_stats_to_csv(
                '{}.csv'.format(name),
                aggregations,
                sort=['observations', 'base_pairs'],
                round_fields={
                    "cov80_percent": 2, "spearmans_rho": 2,
                    "spearmans_rho_pval": 2, "pearson": 2,
                    "pearson_pval": 2
                },
            )

    @staticmethod
    def write_stats_to_json(
        path: str,
        observations: Dict[str, ObservationGroup],
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
            key: val.to_dict(json=True)
            for key, val in observations.items()
        }
        write_data(path, output)

    @staticmethod
    def write_stats_to_csv(
        path: str,
        observations: Dict[str, ObservationGroup],
        sort: List[str],
        round_fields: Dict[str, int],
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
            output.append(val.to_dict())

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
            help="Reference .fasta file(s).",
            dest="refs",
            required=True,
            metavar='',
            nargs='*'
        )

        parser.add_argument(
            '-c',
            help=(
                "Expected counts CSV. "
                "Required columns: reference,expected_count."
            ),
            dest="counts",
            required=False,
            metavar=''
        )

        parser.add_argument(
            '-p',
            dest="pipe",
            default=None,
            choices=Pipers.choices,
            required=False,
            metavar='',
            help=(
                "Choose what to pipe to stdout, if anything. "
                "[Choices: sam, json] (Default: None)."
            ),
        )

        parser.add_argument(
            "-a",
            dest="aggregate",
            required=False,
            default=False,
            action="store_true",
            help=(
                "Switch on aggregation, which will group alignments "
                "and produce summary statistics for each group."
            )
        )

        parser.add_argument(
            "-f",
            dest="format",
            required=False,
            default=Format.CSV,
            choices=Format.choices,
            metavar='',
            help=(
                "If aggregating [-a], output results in this format. "
                "[Choices: csv, json, all] (Default: csv)."
            )
        )

        parser.add_argument(
            "-s",
            dest="splitby",
            required=False,
            default=[Groupers.FASTA],
            choices=Groupers.choices,
            metavar='',
            nargs="+",
            help=(
                "If aggregating [-a], split by these criteria, space separated. "
                "[Choices: {}] (Default: fasta).".format(
                    ' '.join(Groupers.choices))
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

        if not (args.pipe or args.aggregate):
            errprint(
                "[Error]: Neither piping nor aggregating is enabled, "
                "and 'just spin the wheels' is not a valid option. Choose "
                "-a or -p <arg>."
            )
            sys.exit(1)

        if args.counts and (Groupers.REFERENCE in args.splitby):
            errprint(
                "[Error]: If you want to provide expected counts, you "
                "cannot also split by reference. ( try -s fasta )"
            )
            sys.exit(1)

        cls(
            sams=args.sams,
            refs=args.refs,
            counts=args.counts,
            pipe=args.pipe,
            aggregate=args.aggregate,
            splitby=args.splitby,
            output_name=args.name,
            output_format=args.format,
        )
