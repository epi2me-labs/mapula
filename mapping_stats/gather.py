import os
import sys
import csv
import pysam
import argparse
import dataclasses
from typing import List
from pysam import AlignmentFile
from dataclasses import dataclass, field
from mapping_stats.lib.refmap import RefMap
from mapping_stats.lib.bio import get_alignment_tag
from mapping_stats.lib.stats import CoreStats, CorrelationStats
from mapping_stats.lib.const import UNKNOWN, UNMAPPED, UNCLASSIFIED
from mapping_stats.lib.core import BaseSubcommand, UpdatingStatsItem
from mapping_stats.lib.utils import (
    write_data,
    load_data,
    get_data_slots,
    add_attrs,
    add_dists,
    get_group_name,
)


@dataclass
class BaseAlignedReference(object):
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
class AlignedReference(CoreStats, BaseAlignedReference, UpdatingStatsItem):
    """
    Represents an instance of a reference sequence to which reads from
    a given group (where group means binned by barcode, run_id and reference
    file name) have aligned. Tracks many statistics of the alignments to this
    reference.
    """

    # Optional optimisation
    __slots__ = get_data_slots(BaseAlignedReference, CoreStats)

    def update(
        self, aln: pysam.AlignedSegment, refmap: RefMap, *args, **kwargs
    ) -> None:
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

    def __add__(self, new):
        return self._add(self, new)

    @classmethod
    def fromdict(cls, data: dict):
        return cls(**data)


@dataclass
class BaseAlignmentGroup(object):
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
    references: dict = field(default_factory=lambda: {})
    _reference_class = AlignedReference


@dataclass
class AlignmentGroup(
    CoreStats, CorrelationStats, BaseAlignmentGroup, UpdatingStatsItem
):
    """
    Represents a set of binned alignments by reference filename, run_id
    and barcode if available. Tracks both summary statistics (i.e. representing
    the whole group) as well as per-refefence statistics (i.e. broken down
    by individual aligned reference sequence).
    """

    # Optional optimisation
    __slots__ = get_data_slots(BaseAlignmentGroup, CorrelationStats, CoreStats)

    def update(self, aln: pysam.AlignedSegment, refmap: RefMap, *args, **kwargs):
        ref_name = aln.reference_name or UNMAPPED

        # Get or create the reference
        reference = self.references.get(ref_name)
        if not reference:
            length = refmap.get_ref_length(ref_name) or 0
            reference = self._reference_class(name=ref_name, length=length)
            self.references[ref_name] = reference
        reference.update(aln, refmap)

        self._reference_class._update_core_stats(self, aln, refmap)
        self.update_correlations(self.references, refmap)

    def __add__(self, new):
        for rn, rv in new.references.items():
            if not self.references.get(rn):
                self.references[rn] = rv
            else:
                self.references[rn] += rv
        self._reference_class._add(self, new)
        return self

    @classmethod
    def fromdict(cls, data: dict):
        for rn, rv in data.get("references", {}).items():
            data["references"][rn] = cls._reference_class.fromdict(rv)
        return cls(**data)


@dataclass
class BaseAlignments(object):
    """
    A base dataclass providing descriptive and required fields for
    Observed objects.

    Note:
    These fields exist in a separate class to enable them to be
    placed at the beginning of the inheritance chain, which is
    necessary because they do not have default values.
    """

    groups: dict = field(default_factory=lambda: {})
    _group_class = AlignmentGroup


@dataclass
class Alignments(BaseAlignments, UpdatingStatsItem):
    """
    Represents a set of alignments made to some group of reference
    sequences contained by any number of reference files. Bins each
    alignment into a group based on its target reference sequence,
    run_id and barcode (if found), and calculates useful statistics for
    both the entire group and for every reference sequenced within
    each group.
    """

    # Optional optimisation
    __slots__ = get_data_slots(BaseAlignments)

    def update(
        self, aln: pysam.AlignedSegment, refmap: RefMap, *args, **kwargs
    ) -> None:
        group = self._assign_group(aln, refmap=refmap)
        group.update(aln, refmap, *args, **kwargs)

    def _assign_group(self, aln: pysam.AlignedSegment, refmap: RefMap):
        reference = aln.reference_name
        filename = refmap.get_ref_filename(reference) or UNMAPPED
        run_id = get_alignment_tag(aln, "RD", UNKNOWN)
        barcode = get_alignment_tag(aln, "BC", UNCLASSIFIED)
        group_name = get_group_name(filename, run_id, barcode)

        group = self.groups.get(group_name)
        if not group:
            # If we have not already constructed a group
            # object do so now
            group = self._group_class(
                name=filename,
                run_id=run_id,
                barcode=barcode,
            )
            self.groups[group_name] = group

        return group

    def __add__(self, new):
        for gn, gv in new.groups.items():
            if not self.groups.get(gn):
                self.groups[gn] = gv
            else:
                self.groups[gn] += gv
        return self

    @classmethod
    def fromdict(cls, data: dict):
        for gn, gv in data.get("groups", {}).items():
            data["groups"][gn] = cls._group_class.fromdict(gv)
        return cls(**data)


class GatherMappingStats(BaseSubcommand):
    """
    A subcommand that runs a process designed to scan alignments
    made in SAM format and accumulate many useful statistics which are
    binned into groups and reported in JSON format.
    """

    def __init__(
        self,
        sam_path: str,
        out_path: str,
        json_path: str,
        fasta_paths: List[str],
        exp_counts_path: str = None,
    ) -> None:
        self.json_path = json_path
        self.out_path = out_path
        self.sam_path = sam_path
        self.fasta_paths = fasta_paths
        self.refmap = RefMap(fasta_paths, exp_counts_path)

        self.records = AlignmentFile(sam_path, "r")
        self.outfile = AlignmentFile(out_path, "w", template=self.records)

        self.observed = self.load()
        self.update()
        self.write()

    def load(
        self,
    ) -> Alignments:
        if not os.path.exists(self.json_path):
            return Alignments()
        data = load_data(self.json_path)
        return Alignments.fromdict(data)

    def update(
        self,
    ) -> None:
        for aln in self.records.fetch(until_eof=True):
            if self.outfile:
                self.outfile.write(aln)

            self.observed.update(aln, refmap=self.refmap)

    def write(
        self,
    ) -> None:
        data = dataclasses.asdict(self.observed)
        write_data(self.json_path, data)

    @classmethod
    def execute(cls, argv) -> None:
        parser = argparse.ArgumentParser(
            description="Gather mapping stats from a SAM/BAM file"
        )

        parser.add_argument(
            "-s",
            "--SAM",
            default=sys.stdin,
            help="Input sam/bam file. (default: stdin)",
        )

        parser.add_argument(
            "-o", "--OUT", default="-", help="Output sam file. (default: stdout)"
        )

        parser.add_argument(
            "-j", "--JSON", required=False, default="mapping-stats.json"
        )

        parser.add_argument(
            "-f",
            "--FASTA",
            nargs="*",
            required=True,
            help="List of references to which alignments have been made",
        )

        parser.add_argument(
            "-e",
            "--EXP",
            required=False,
            default=None,
            help="A .CSV file containing expected counts by reference name",
        )

        args = parser.parse_args(argv)
        cls(
            sam_path=args.SAM,
            out_path=args.OUT,
            json_path=args.JSON,
            fasta_paths=args.FASTA,
            exp_counts_path=args.EXP,
        )
