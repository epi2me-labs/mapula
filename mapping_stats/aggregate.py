import os
import argparse
import dataclasses
from typing import List
from mapping_stats.count import Alignments
from mapping_stats.lib.refmap import RefMap
from mapping_stats.lib.core import BaseSubcommand
from mapping_stats.lib.utils import load_data, write_data


class AggregateMappingStats(BaseSubcommand):
    """
    A subcommand that runs a process designed
    to combine the JSON outputs from runs of the
    Gather or ERCC subcommands into a single file.
    """

    def __init__(
        self,
        json_paths: List[str],
        fasta_paths: List[str],
        exp_counts_path: str = None,
    ) -> None:
        self.observed = Alignments()
        self.json_paths = json_paths
        self.refmap = RefMap(fasta_paths, exp_counts_path)

        self.load()
        self.write()

    def load(self) -> None:
        for path in self.json_paths:
            abspath = os.path.abspath(path)
            data = load_data(abspath)

            new = Alignments.fromdict(**data)
            self.observed.add(new, self.refmap)

    def write(self) -> None:
        data = dataclasses.asdict(self.observed)
        write_data("merged.stats.mapula.json", data)

    @classmethod
    def execute(cls, argv) -> None:
        parser = argparse.ArgumentParser(
            description="Combine mapping stats .JSON outputs"
        )

        parser.add_argument(
            "-j", "--JSON", nargs="*", required=False, default=["mapping-stats.json"]
        )

        parser.add_argument(
            "-r",
            "--refs",
            nargs="*",
            required=True,
            help="List of references to which alignments have been made",
        )

        parser.add_argument(
            "-e",
            "--exp",
            required=False,
            default=None,
            help="A .CSV file containing expected counts by reference name",
        )

        args = parser.parse_args(argv)
        cls(
            json_paths=args.JSON,
            fasta_paths=args.refs,
            exp_counts_path=args.exp,
        )
