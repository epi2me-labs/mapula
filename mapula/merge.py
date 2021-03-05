import os
import sys
import argparse
from mapula.const import Format
from typing import List, Dict, Union
from mapula.count import CountMappingStats
from mapula.utils import errprint, load_data
from mapula.observation import ObservationGroup


class MergeMappingStats(CountMappingStats):

    def __init__(
        self,
        json: List[str],
        counts: str,
        output_name: str,
        output_format: str
    ) -> None:
        """
        A subcommand that runs a process designed
        to combine the JSON outputs from several independent
        runs of the count subcommand into a single file.
        """
        errprint("Running: Mapula (merge)")

        self.json_paths = json
        self.counts_file = counts
        self.output_name = output_name
        self.output_format = output_format
        self.output_corrs = bool(counts)

        self.counts = self.get_expected_counts(self.counts_file)

        errprint("[1/2] Merging aggregations")
        self.aggregations = self.merge_aggregations(
            self.json_paths,
            self.counts
        )

        errprint("[2/2] Writing data")
        self.write_aggregations(
            self.output_name,
            self.output_format,
            self.aggregations
        )

    def merge_aggregations(
        self,
        json_paths: List[str],
        counts: Union[None, Dict[str, int]]
    ) -> Dict[str, ObservationGroup]:
        """

        """
        aggregations = {}
        splitby = []

        for path in json_paths:

            try:
                existing_data = load_data(path)
            except FileNotFoundError:
                errprint(f'[Error]: {path} does not exist.')
                sys.exit(1)

            for obsname, obsdata in existing_data.items():
                obsgrp = ObservationGroup.from_dict(obsdata, counts)

                # Consistency check: All input observation groups should
                # have the same identifying/grouping attributes in order
                # to be able to be aggregated.
                if not splitby:
                    splitby = list(obsgrp.identity.keys())
                elif splitby != list(obsgrp.identity.keys()):
                    errprint(
                        '[Error]: Cannot aggregate. The json input appears '
                        'to have been created using different splitby criteria. '
                        'I.e. ({}) vs ({})'.format(
                            ','.join(splitby),
                            ','.join(obsgrp.identity.keys())
                        )
                    )
                    sys.exit(1)

                if aggregations.get(obsname):
                    aggregations[obsname] += obsgrp
                    continue
                aggregations[obsname] = obsgrp

        return aggregations

    @classmethod
    def execute(cls, argv) -> None:
        """
        Parses command line arguments and 
        initialises a MergeMappingStats object.
        """
        parser = argparse.ArgumentParser(
            description="Combine .json outputs from mapula count"
        )

        parser.add_argument(
            help=(
                "Input .json files from mapula count. "
                "(Default: [stdin])."
            ),
            dest="json",
            default=[sys.stdin],
            metavar='',
            nargs='*',
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
        cls(
            json=args.json,
            counts=args.counts,
            output_name=args.name,
            output_format=args.format
        )
