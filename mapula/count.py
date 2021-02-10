import os
import sys
import pysam
import argparse
from typing import List, Union
from pysam import AlignmentFile
from dataclasses import dataclass
from argparse import RawTextHelpFormatter
from mapula.lib.refmap import RefMap
from mapula.lib.bio import get_alignment_tag
from mapula.lib.stats import CoreStats, CorrelationStats
from mapula.lib.const import UNKNOWN, UNMAPPED, UNCLASSIFIED
from mapula.lib.core import BaseSubcommand, MappingStatsContainer
from mapula.lib.utils import (
    write_data,
    load_data,
    get_data_slots,
    add_attrs,
    add_dists,
    get_group_name,
)


@dataclass
class AlignedReferenceAttributes(object):
    """
    A base dataclass providing descriptive and required fields for
    ObservedReference.

    Note:
    These fields exist in a separate class to enable them to be placed
    at the beginning of the inheritance chain, which is necessary
    because they do not have default values.
    """

    name: str
    length: int


@dataclass
class AlignedReference(CoreStats, MappingStatsContainer, AlignedReferenceAttributes):
    """
    Represents an instance of a reference sequence to which reads from
    a given group (where group means binned by barcode, run_id and reference
    file name) have aligned. Tracks many statistics of the alignments to this
    reference.
    """

    # Optional optimisation
    __slots__ = get_data_slots(AlignedReferenceAttributes, CoreStats)

    def update(self, aln: pysam.AlignedSegment, refmap: RefMap) -> None:
        self._update_core_stats(self, aln, refmap)

    @staticmethod
    def _add(old, new):
        add_attrs(
            old,
            new,
            "alignment_count",
            "read_count",
            "total_base_pairs",
            "primary_count",
            "secondary_count",
            "supplementary_count",
        )

        add_dists(old, new, "alignment_accuracies")
        add_dists(old, new, "read_qualities")
        add_dists(old, new, "read_lengths")
        add_dists(old, new, "alignment_coverages")

        old.update_read_n50(old.total_base_pairs)
        old.update_median_quality()
        old.update_median_accuracy()
        old.update_alignment_cov80(old.read_count)

        return old

    def add(self, new, refmap: RefMap):
        return self._add(self, new)


@dataclass
class AlignmentGroupAttributes(object):
    """
    A base dataclass providing descriptive and required fields for
    ObservedGroup.

    Note:
    These fields exist in a separate class to enable them to be placed at
    the beginning of the inheritance chain, which is necessary because they
    do not have default values.
    """

    name: str
    run_id: str
    barcode: str


@dataclass
class AlignmentGroup(
    CoreStats, CorrelationStats, MappingStatsContainer, AlignmentGroupAttributes
):
    """
    Represents a set of binned alignments by reference filename, run_id
    and barcode if available. Tracks both summary statistics (i.e. representing
    the whole group) as well as per-refefence statistics (i.e. broken down
    by individual aligned reference sequence).
    """

    _child_type = AlignedReference

    # Optional optimisation
    __slots__ = get_data_slots(AlignmentGroupAttributes, CorrelationStats, CoreStats)

    def update(self, aln: pysam.AlignedSegment, refmap: RefMap):
        ref_name = aln.reference_name or UNMAPPED
        reference = self.children.get(ref_name)
        if not reference:
            length = refmap.get_ref_length(ref_name) or 0
            reference = AlignedReference(name=ref_name, length=length)
            self.children[ref_name] = reference
        reference.update(aln, refmap)

        AlignedReference._update_core_stats(self, aln, refmap)
        self.update_correlations(self.children, refmap)

    def add(self, new, refmap: RefMap):
        for rn, rv in new.children.items():
            if not self.children.get(rn):
                self.children[rn] = rv
            else:
                self.children[rn].add(rv, refmap)
        self._child_type._add(self, new)
        self.update_correlations(self.children, refmap)
        return self


@dataclass
class Alignments(CoreStats, CorrelationStats, MappingStatsContainer):
    """
    Represents a set of alignments made to some group of reference
    sequences contained by any number of reference files. Bins each
    alignment into a group based on its target reference sequence,
    run_id and barcode (if found), and calculates useful statistics for
    both the entire group and for every reference sequenced within
    each group.
    """

    _child_type = AlignmentGroup

    # Optional optimisation
    __slots__ = get_data_slots(CoreStats, CorrelationStats, MappingStatsContainer)

    def update(self, aln: pysam.AlignedSegment, refmap: RefMap) -> None:
        AlignedReference._update_core_stats(self, aln, refmap)
        self.update_correlations(self.children, refmap)

        group = self._assign_group(aln, refmap)
        group.update(aln, refmap)

    def _assign_group(self, aln: pysam.AlignedSegment, refmap: RefMap):
        reference = aln.reference_name
        filename = refmap.get_ref_filename(reference) or UNMAPPED
        run_id = get_alignment_tag(aln, "RD", UNKNOWN)
        barcode = get_alignment_tag(aln, "BC", UNCLASSIFIED)
        group_name = get_group_name(filename, run_id, barcode)

        group = self.children.get(group_name)
        if not group:
            # If we have not already constructed a group
            # object do so now
            group = AlignmentGroup(name=filename, run_id=run_id, barcode=barcode)
            self.children[group_name] = group

        return group

    def add(self, new, refmap: RefMap):
        for gn, gv in new.children.items():
            if not self.children.get(gn):
                self.children[gn] = gv
            else:
                self.children[gn].add(gv, refmap)
        return self


class CountMappingStats(BaseSubcommand):
    """
    A subcommand that runs a process designed to scan alignments
    made in SAM format and accumulate many useful statistics which are
    binned into groups and reported in JSON format.
    """

    def __init__(
        self,
        sam_path: str,
        out_path: Union[str, None],
        json_path: str,
        fasta_paths: List[str],
        exp_counts_path: str = None,
    ) -> None:
        self.json_path = json_path
        self.out_path = out_path
        self.sam_path = sam_path

        self.records = AlignmentFile(sam_path, "r")
        if out_path:
            self.outfile = AlignmentFile(out_path, "w", template=self.records)
        self.refmap = RefMap(fasta_paths, exp_counts_path)

        self.observed = self.load()
        self.update()
        self.write()

    def load(
        self,
    ) -> Alignments:
        if not os.path.exists(self.json_path):
            return Alignments()
        data = load_data(self.json_path)
        return Alignments.fromdict(**data)

    def update(
        self,
    ) -> None:
        for aln in self.records.fetch(until_eof=True):
            if self.out_path:
                self.outfile.write(aln)
            self.observed.update(aln, self.refmap)

    def write(
        self,
    ) -> None:
        data = self.observed.asdict()
        write_data(self.json_path, data)

    @classmethod
    def execute(cls, argv) -> None:
        parser = argparse.ArgumentParser(
            description="Count mapping stats from a SAM/BAM file",
            formatter_class=RawTextHelpFormatter,
        )

        parser.add_argument(
            "-s",
            "--sam",
            default=sys.stdin,
            help="Input sam/bam file. (default: stdin)",
        )

        parser.add_argument(
            "-o",
            "--out",
            default=None,
            help="Output sam file. Use - for piping out. (default: None)",
        )

        parser.add_argument(
            "-j",
            "--json",
            required=False,
            default="stats.mapula.json",
            help="Name of the output json (default: stats.mapula.json)",
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
            sam_path=args.sam,
            out_path=args.out,
            json_path=args.json,
            fasta_paths=args.refs,
            exp_counts_path=args.exp,
        )
