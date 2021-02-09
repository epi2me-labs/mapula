import os
import argparse
import dataclasses
from typing import List
from mapping_stats.gather import Alignments
from mapping_stats.lib.core import BaseSubcommand
from mapping_stats.lib.utils import load_data, write_data


class AggregateMappingStats(BaseSubcommand):
    """
    A subcommand that runs a process designed
    to combine the JSON outputs from runs of the
    Gather or ERCC subcommands into a single file.
    """

    def __init__(self, out_path: str, json_paths: List[str], use_ercc=False) -> None:
        self.observed = Alignments()
        self.out_path = out_path
        self.json_paths = json_paths
        self.use_ercc = use_ercc
        self.load()
        self.write(out_path)

    def load(self) -> None:
        for path in self.json_paths:
            abspath = os.path.abspath(path)
            data = load_data(abspath)

            new = Alignments.fromdict(data)
            self.observed += new

    def write(self, json_path) -> None:
        data = dataclasses.asdict(self.observed)
        write_data(json_path, data)

    @classmethod
    def execute(cls, argv) -> None:
        parser = argparse.ArgumentParser(
            description="Combine mapping stats .JSON outputs"
        )

        parser.add_argument(
            "-o",
            "--OUT",
            default="merged.mapping-stats.json",
            help="Output merged json file. (default: merged.mapping-stats.json)",
        )

        parser.add_argument(
            "-j", "--JSON", nargs="*", required=False, default=["mapping-stats.json"]
        )

        args = parser.parse_args(argv)
        cls(out_path=args.OUT, json_paths=args.JSON)
