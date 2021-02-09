import os
import csv
import pysam
from typing import TypedDict, Union, List


class RefMap(dict):
    """
    Accepts a list of fasta files and creates a mapping of
    each reference sequence name to the file they originate
    from and their length.
    """

    class RefMapping(TypedDict):
        filename: str
        length: int
        expected_count: float

    def __init__(self, fasta_paths: List[str], exp_counts_path: str) -> None:
        super(RefMap, self).__init__()

        self.fasta_paths = fasta_paths
        self.exp_counts_path = exp_counts_path

        self.load_mappings(self.fasta_paths, self.exp_counts_path)

    def load_mappings(self, fasta_paths: List[str], exp_counts_path: str) -> dict:
        exp_counts = self._load_exp_counts(exp_counts_path)

        for path in fasta_paths:
            fasta = pysam.FastaFile(path)
            basename = os.path.basename(path)
            for reference in fasta.references:
                self[reference] = self.RefMapping(
                    filename=basename,
                    length=fasta.get_reference_length(reference),
                    expected_count=exp_counts.get(reference, None),
                )

    def _load_exp_counts(self, path: str) -> dict:
        counts = {}
        if os.path.exists(path):
            reader = csv.DictReader(open(path))
            for line in reader:
                counts[line["reference"]] = float(line["expected_count"])
        return counts

    def _get_ref_key(self, ref: str, key: str, default=None) -> dict:
        try:
            return self[ref][key]
        except KeyError:
            return default

    def get_ref_filename(self, ref: str) -> Union[str, None]:
        return self._get_ref_key(ref, "filename")

    def get_ref_expected_count(self, ref: str) -> Union[str, None]:
        return self._get_ref_key(ref, "expected_count")

    def get_ref_length(self, ref: str) -> Union[int, None]:
        return self._get_ref_key(ref, "length")

    def filter_by_filename(self, filename: str) -> dict:
        mappings = {}
        for k, v in self.items():
            if v["filename"] == filename:
                mappings[k] = v
        return mappings
