import sys
import os
import csv
import pysam
from mapula.lib.utils import errprint
from typing import TypedDict, Union


class RefMap(dict):
    class RefMapping(TypedDict):
        filename: str
        length: int
        expected_count: float

    def __init__(
        self,
        name: str,
        path: str = None,
        exp_counts_path: str = None,
    ) -> None:
        """
        RefMap is designed to store basic information about the
        references found in .fasta files. Each reference is stored
        by name, along with their length, and if present, 
        their expected count.
        """
        super(RefMap, self).__init__()
        self.name = name
        self.has_counts = bool(exp_counts_path)

        if path:
            self.path = path
            self.basename = os.path.basename(path)
            self.fasta = pysam.FastaFile(path)
            self.exp_counts_path = exp_counts_path
            self._load_mappings(self.exp_counts_path)

    def _load_mappings(
        self,
        exp_counts_path: str = None
    ) -> dict:
        """
        Iterates over the provided fasta_paths and creates
        a record of each reference found.

        fasta_paths:
        List of paths to files in .fasta format

        exp_counts_path:
        A path to a .csv file containing two columns (with a header):
        "reference" and "expected_count"
        """
        if exp_counts_path:
            exp_counts = self._load_exp_counts(exp_counts_path)
        else:
            exp_counts = {}

        for reference in self.fasta.references:
            self[reference] = self.RefMapping(
                length=self.fasta.get_reference_length(reference),
                expected_count=exp_counts.get(reference),
            )

    def _load_exp_counts(self, path: str) -> dict:
        """
        Assuming path is a csv containing two columns (with a header):
        "reference" and "expected_count", loads it into a dict.
        """
        if not os.path.exists(path):
            errprint("Error opening expected counts path, {}.".format(path))
            return

        reader=csv.DictReader(open(path))
        try:
            counts={
                line["reference"]: float(line["expected_count"])
                for line in reader
            }
            return counts
        except KeyError:
            errprint(
                "Expected counts CSV is incorrectly formatted: "
                "Two columns and a header are expected: "
                "reference, expected_count."
            )

    def _get_ref_key(self, ref: str, key: str, default=None) -> dict:
        """
        Utility for acquiring the value of a key for a given ref without
        having to use a nested .get().
        """
        try:
            return self[ref][key]
        except KeyError:
            return default

    def get_ref_expected_count(self, ref: str) -> Union[str, None]:
        """
        Gets the expected count for a ref, or None
        """
        return self._get_ref_key(ref, "expected_count")

    def get_ref_length(self, ref: str) -> Union[int, None]:
        """
        Gets the length of a ref, or None
        """
        return self._get_ref_key(ref, "length")

    def check_counts(self, filename: str) -> bool:
        # Do we have counts for all the references
        # in the refmap that belong to this filename
        has_counts=False
        for k, v in self.filter_by_filename(
            filename
        ).items():
            if v.get('expected_count') != None:
                has_counts=True
            else:
                if not has_counts:
                    continue
                errprint(
                    "Error: some references from {} do "
                    "not have counts in the expected counts "
                    "csv, but others do".format(filename)
                )
                sys.exit(1)
        return has_counts
