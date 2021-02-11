import os
import sys
import argparse
from typing import List, Dict
from mapula.lib.refmap import RefMap
from mapula.count import AlignmentGroup
from mapula.lib.const import UNMAPPED
from mapula.lib.utils import (
    errprint,
    load_data,
    write_data,
    parse_cli_key_value_pairs,
)


class AggregateMappingStats:

    def __init__(
        self,
        json_paths: List[str],
        refs: Dict[str, str],
        counts: Dict[str, str],
    ) -> None:
        """
        A subcommand that runs a process designed
        to combine the JSON outputs from several independent
        runs of the count subcommand into a single file.
        """
        errprint("--- Running: Mapula (aggregate)")

        self.groups = {}
        self.json_paths = json_paths
        self.use_corrs = bool(counts)
        self.refs = parse_cli_key_value_pairs(refs)
        self.counts = parse_cli_key_value_pairs(counts)

        self.refmaps = self.load_refmaps()
        self.aggregate()
        self.write_stats_to_json()

        errprint("--- Operation successful")

    def aggregate(self) -> None:
        """

        """
        for path in self.json_paths:
            existing_data = load_data(path)
            for key, val in existing_data.items():
                name = val['name']
                refmap = self.refmaps[name]
                group = AlignmentGroup.fromdict(
                    refmap, self.use_corrs, **val)

                if not self.groups.get(key):
                    self.groups[key] = group
                else:
                    self.groups[key] += group

    def load_refmaps(self):
        """
        Given the reference and counts paths provided, 
        iterates over each pair and builds RefMap objects. 
        In addition, makes a RefMap for unmapped.
        """
        refmaps = {}
        for name, ref in self.refs.items():
            counts = self.counts.get(name)
            refmaps[name] = RefMap(
                name, ref, counts,
            )
        refmaps[UNMAPPED] = RefMap(UNMAPPED)
        return refmaps

    def _load_group(
        self,
        path: str
    ) -> dict:
        """
        Instantiate a working state from a .json file, 
        if path exists, otherwise returns an empty dict.
        """
        if not os.path.exists(path):
            errprint("[Error]: {} does not exist.".format(path))
            sys.exit(1)

        existing_data = load_data(path)
        for key, val in existing_data.items():
            try:
                name = val['name']
                refmap = self.refmaps[name]
                existing_data[key] = AlignmentGroup.fromdict(
                    refmap, **val)
            except KeyError:
                errprint("[Error]: existing .json data is malformed")
                sys.exit(1)
        return existing_data

    def write_stats_to_json(self) -> None:
        """
        Gathers the stats data together for each
        AlignedGroup in self.groups, and dumps it to file.
        """
        data = {k: v.todict() for k, v in self.groups.items()}
        write_data("merged.mapula.json", data)

    @classmethod
    def execute(cls, argv) -> None:
        """
        Parses command line arguments and 
        initialises a AggregateMappingStats object.
        """
        parser = argparse.ArgumentParser(
            description="Combine mapping stats .JSON outputs"
        )

        parser.add_argument(
            "-j", "--JSON", nargs="*", required=False, default=["mapping-stats.json"]
        )

        parser.add_argument(
            '-r',
            '--refs',
            help='Provide reference .fasta files using the syntax: '
            'name=path_to_ref.',
            nargs='*',
            default=[]
        )

        parser.add_argument(
            '-c',
            '--counts',
            help='Provide expected counts in csv format using the syntax: '
            'name=path_to_counts, where name should be equal to a name given '
            'to --references. Expected column headings: reference,expected_count. '
            'The reference column should contain the ID of a sequence in the '
            'corresponding reference file. The expected counts column should '
            'equal the expected number of observations for that sequence.',
            nargs='*',
            default=[]
        )

        args = parser.parse_args(argv)
        cls(
            json_paths=args.JSON,
            refs=args.refs,
            counts=args.counts,
        )
