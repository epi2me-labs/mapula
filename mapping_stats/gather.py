import os
import sys
import pysam
import argparse
import dataclasses
from typing import List
from pysam import AlignmentFile
from dataclasses import dataclass, field
from mapping_stats.lib.stats import CoreStats
from mapping_stats.lib.const import UNKNOWN, UNMAPPED, UNCLASSIFIED
from mapping_stats.lib.core import BaseSubcommand, UpdatingStatsItem
from mapping_stats.lib.bio import ReferenceFastas, FastaFile, get_alignment_tag
from mapping_stats.lib.utils import (write_data, load_data, get_data_slots, 
    add_attrs, add_dists, get_group_name)


@dataclass
class BaseAlignedReference(object):
    """
    A base dataclass providing descriptive and
    required fields for ObservedReference. 
    
    Note:
    These fields exist in a separate class to 
    enable them to be placed at the beginning of 
    the inheritance chain, which is necessary 
    because they do not have default values.
    """
    name: str
    length: int


@dataclass
class AlignedReference(
        CoreStats,
        BaseAlignedReference,
        UpdatingStatsItem
    ):
    """
    Represents an instance of a reference 
    sequence to which reads from a given group
    (where group means binned by barcode, run_id
    and reference file name) have aligned. Tracks 
    many statistics of the alignments to this 
    reference. 
    """

    # Optional optimisation
    __slots__ = get_data_slots(
        BaseAlignedReference,
        CoreStats
    )

    def update(
        self,
        aln: pysam.AlignedSegment,
    ) -> None:
        self._update_item(self, aln, self.length)

    @staticmethod
    def _update_item(
        item,
        aln: pysam.AlignedSegment,
        reference_length: int
    ) -> dict:
        item.alignment_count += 1

        if aln.is_supplementary:
            item.supplementary_count += 1
            return

        if aln.is_secondary:
            item.secondary_count += 1
            return

        item.read_count += 1
        item.update_total_base_pairs(aln)

        item.update_read_length_dist(aln)
        item.update_read_n50(
            item.total_base_pairs
        )

        item.update_read_quality_dist(aln)
        item.update_median_quality()

        if aln.is_unmapped:
            return

        item.primary_count += 1

        item.update_alignment_accuracy_dist(aln)
        item.update_median_accuracy()

        item.update_alignment_coverage_dist(
            aln, reference_length)
        item.update_alignment_cov80(
            item.read_count)

    @staticmethod
    def _add(
        old, 
        new
    ):
        add_attrs(
            old,
            new,
            'alignment_count',
            'read_count',
            'total_base_pairs',
            'primary_count',
            'secondary_count',
            'supplementary_count',
        )

        add_dists(
            old,
            new,
            'alignment_accuracies'
        )

        add_dists(
            old,
            new,
            'read_qualities'
        )
       
        add_dists(
            old,
            new,
            'read_lengths'
        )

        add_dists(
            old,
            new,
            'alignment_coverages'
        )

        old.update_read_n50(
            old.total_base_pairs
        )
        old.update_median_quality()
        old.update_median_accuracy()
        old.update_alignment_cov80(
            old.read_count
        )

        return old

    def __add__(
        self, 
        new
    ):
        return self._add(self, new)

    @classmethod
    def fromdict(
        cls,
        data: dict
    ):
        return cls(**data)


@dataclass
class BaseAlignmentGroup(object):
    """
    A base dataclass providing descriptive and
    required fields for ObservedGroup. 
    
    Note:
    These fields exist in a separate class to 
    enable them to be placed at the beginning of 
    the inheritance chain, which is necessary 
    because they do not have default values.
    """
    name: str
    run_id: str
    barcode: str
    references: dict = field(default_factory=lambda: {})
    _reference_class = AlignedReference


@dataclass
class AlignmentGroup(
        CoreStats,
        BaseAlignmentGroup,
        UpdatingStatsItem
    ):
    """
    Represents a set of binned alignments by 
    reference filename, run_id and barcode if
    available. Tracks both summary statistics 
    (i.e. representing the whole group) as well
    as per-refefence statistics (i.e. broken down
    by individual aligned reference sequence).
    """

    # Optional optimisation
    __slots__ = get_data_slots(
        BaseAlignmentGroup,
        CoreStats
    )

    def update(
        self,
        aln: pysam.AlignedSegment,
        fasta: FastaFile = None
    ):
        name = aln.reference_name or UNMAPPED

        # Get or create the reference
        reference = self.references.get(name)
        if not reference:
            length = (
                fasta.get_reference_length(name) 
                if fasta else 0
            )
            reference = self._reference_class(
                name=name,
                length=length
            )
            self.references[name] = reference
        reference.update(aln)

        # Ensure we're also updating the aggregate
        # stats.
        self._reference_class._update_item(
            self, aln, reference.length)

    def __add__(
        self,
        new
    ):
        for rn, rv in new.references.items():
            if not self.references.get(rn):
                self.references[rn] = rv
            else:
                self.references[rn] += rv
        self._reference_class._add(self, new)
        return self

    @classmethod
    def fromdict(
        cls,
        data: dict
    ):
        for rn, rv in data.get('references', {}).items():
            data['references'][rn] = cls._reference_class.fromdict(rv)
        return cls(**data)


@dataclass
class BaseAlignments(object):
    """
    A base dataclass providing descriptive and
    required fields for Observed objects. 
    
    Note:
    These fields exist in a separate class to 
    enable them to be placed at the beginning of 
    the inheritance chain, which is necessary 
    because they do not have default values.
    """
    groups: dict = field(default_factory=lambda: {})
    _group_class = AlignmentGroup


@dataclass
class Alignments(
        BaseAlignments,
        UpdatingStatsItem
    ):
    """
    Represents a set of alignments made to some
    group of reference sequences contained by
    any number of reference files. Bins each
    alignment into a group based on its target
    reference sequence, run_id and barcode (if 
    found), and calculates useful statistics for
    bohh the entire group and for every reference 
    sequenced within each group.
    """

    # Optional optimisation
    __slots__ = get_data_slots(
        BaseAlignments
    )

    def update(
        self,
        aln: pysam.AlignedSegment,
        fastas: ReferenceFastas
    ) -> None:
        group = self._assign_group(aln, fastas)
        fasta = fastas.get_fasta_file(group.name)
        group.update(aln, fasta=fasta)

    def _assign_group(
        self,
        aln: pysam.AlignedSegment,
        fastas: ReferenceFastas
    ):
        reference = aln.reference_name
        name = fastas.get_fasta_name(reference) or UNMAPPED
        run_id = get_alignment_tag(aln, 'RD', UNKNOWN)
        barcode = get_alignment_tag(aln, 'BC', UNCLASSIFIED)
        group_name = get_group_name(name, run_id, barcode)

        group = self.groups.get(group_name)
        if not group:
            # If we have not already constructed a group 
            # object do so now
            group = self._group_class(
                name=name,
                run_id=run_id, 
                barcode=barcode,
            )
            self.groups[group_name] = group

        return group

    def __add__(
        self,
        new
    ):
        for gn, gv in new.groups.items():
            if not self.groups.get(gn):
                self.groups[gn] = gv
            else:
                self.groups[gn] += gv
        return self

    @classmethod
    def fromdict(
        cls,
        data: dict
    ):
        for gn, gv in data.get('groups', {}).items():
            data['groups'][gn] = cls._group_class.fromdict(gv)
        return cls(**data)


class GatherMappingStats(BaseSubcommand):
    """
    A subcommand that runs a process designed
    to scan alignments made in SAM format and
    accumulate many useful statistics which are
    binned into groups and reported in JSON format.
    """
    def __init__(
        self,
        sam_path: str,
        out_path: str,
        json_path: str,
        fasta_paths: List[str]
    ) -> None:
        self.json_path = json_path
        self.out_path = out_path
        self.sam_path = sam_path
        self.fasta_paths = fasta_paths
        self.observed = self.load(json_path)

        self.fastas = ReferenceFastas(fasta_paths)
        self.records = AlignmentFile(sam_path, "r")
        self.outfile = AlignmentFile(out_path, "w", 
            template=self.records)

        self.update(self.outfile)
        self.write(json_path)

    def write(
        self,
        json_path
    ) -> None:
        data = dataclasses.asdict(self.observed)
        write_data(json_path, data)

    @staticmethod
    def load(
        path: str = None
    ) -> Alignments:
        if not os.path.exists(path):
            return Alignments()
        data = load_data(path)
        return Alignments.fromdict(data)

    def update(
        self,
        outfile: AlignmentFile = None
    ) -> None:
        for aln in self.records.fetch(
            until_eof=True
        ):
            if outfile:
                outfile.write(aln)

            self.observed.update(
                aln, 
                fastas=self.fastas
            )

    @classmethod
    def execute(
        cls,
        argv
    ) -> None:
        parser = argparse.ArgumentParser(
            description='Gather mapping stats from a SAM/BAM file')
        
        parser.add_argument(
            '-s',
            '--SAM',
            default=sys.stdin,
            help="Input sam/bam file. (default: stdin)"
        )

        parser.add_argument(
            "-o", 
            "--OUT",
            default='-',
            help="Output sam file. (default: stdout)"
        )

        parser.add_argument(
            '-j',
            '--JSON',
            required=False,
            default='mapping-stats.json'
        )

        parser.add_argument(
            "-f",
            "--FASTA",
            nargs="*",
            required=True,
            help="List of references to which alignments can be counted"
        )

        args = parser.parse_args(argv)
        cls(
            sam_path=args.SAM, 
            out_path=args.OUT, 
            json_path=args.JSON, 
            fasta_paths=args.FASTA
        )