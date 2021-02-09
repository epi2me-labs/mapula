import os
import sys
import csv
import pysam
import argparse
from pysam import AlignmentFile
from dataclasses import dataclass, field
from typing import List, Dict, TypedDict
from mapping_stats.lib.core import UpdatingStatsItem
from mapping_stats.lib.bio import ReferenceFastas, FastaFile
from mapping_stats.lib.utils import load_data, get_data_slots
from mapping_stats.lib.stats import CoreStats, ERCCCorrelationStats
from mapping_stats.gather import (GatherMappingStats, BaseAlignmentGroup,
    Alignments, AlignmentGroup)


@dataclass
class AlignmentGroupERCC(
        AlignmentGroup,
        ERCCCorrelationStats
    ):
    """
    Represents a set of binned alignments by 
    reference filename, run_id and barcode if
    available. Tracks both summary statistics 
    (i.e. representing the whole group) as well
    as per-refefence statistics (i.e. broken down
    by individual aligned reference sequence).

    Note:
    Although this class is named `ObservedGroupERCC`
    it does not mean that each group will contain
    alignments of ERCC reads to ERCC references.
    The name of this class reflects that it has
    slightly modified behaviour to accomodate ERCC
    specific counting.
    """

    # Optional optimisation
    __slots__ = get_data_slots(
        BaseAlignmentGroup,
        CoreStats,
        ERCCCorrelationStats
    )

    def update(
        self,
        aln: pysam.AlignedSegment,
        fasta: FastaFile = None,
        ercc_counts: Dict[str, float] = None,
        ercc_reference: str = None
    ) -> None:
        super(AlignmentGroupERCC, self).update(
            aln, fasta=fasta
        )

        if ercc_reference == self.name:
            self.update_correlations(
                self.references,
                ercc_counts
            )


@dataclass
class BaseAlignmentsERCC(object):
    """
    A base dataclass providing descriptive and
    required fields for ObservedERCC objects. 
    
    Note:
    These fields exist in a separate class to 
    enable them to be placed at the beginning of 
    the inheritance chain, which is necessary 
    because they do not have default values.
    """
    ercc_reference: str
    ercc_counts: Dict[str, int]
    groups: dict = field(default_factory=lambda: {})


@dataclass
class AlignmentsERCC(
        Alignments,
        BaseAlignmentsERCC,
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
    _group_class = AlignmentGroupERCC

    # Optional optimisation
    __slots__ = get_data_slots(
        Alignments,
        BaseAlignmentsERCC
    )

    def update(
        self,
        aln: pysam.AlignedSegment,
        fastas: ReferenceFastas,
    ) -> None:
        group = self._assign_group(aln, fastas)
        fasta = fastas.get_fasta_file(group.name)
        group.update(
            aln, 
            fasta=fasta,
            ercc_counts=self.ercc_counts,
            ercc_reference=self.ercc_reference
        )

    def __add__(
        self,
        new
    ):
        super(AlignmentsERCC, self).__add__(new)
        for gv in self.groups.values():
            gv.update_correlations(
                gv.references, 
                self.ercc_counts
            )
        return self

    @classmethod
    def fromdict(
        cls,
        data: dict,
    ):
        for gn, gv in data.get('groups', {}).items():
            data['groups'][gn] = cls._group_class.fromdict(gv)
        return cls(**data)


class ERCCMappingStats(GatherMappingStats):
    """
    A subcommand that runs a process designed
    to scan alignments made in SAM format and
    accumulate many useful statistics which are
    binned into groups and reported in JSON format.

    In addition to the Gather subcommand, this
    process collects additional information specific
    to spike-in control experiments, such as
    correlations of counts with expected counts. 
    """
    def __init__(
        self,
        sam_path: str,
        out_path: str,
        json_path: str,
        ercc_counts: str,
        ercc_reference: str,
        fasta_paths: List[str]
    ) -> None:
        self.out_path = out_path
        self.sam_path = sam_path
        self.json_path = json_path

        self.ercc_reference = os.path.basename(
            ercc_reference
        )
        self.ercc_counts = self.load_counts(
            ercc_counts
        )
        self.fasta_paths = set(
            fasta_paths + [ercc_reference]
        )

        self.observed = self.load(json_path)
        self.fastas = ReferenceFastas(self.fasta_paths)
        self.records = AlignmentFile(sam_path, "r")
        self.outfile = AlignmentFile(out_path, "w", 
            template=self.records)

        self.update(self.outfile)
        self.write(json_path)

    @staticmethod
    def load_counts(
        ercc_counts: str
    ) -> dict:
        counts = {}
        reader = csv.DictReader(open(ercc_counts))
        for line in reader:
            counts[line['Reference']] = float(
                line['expected_count'])
        return counts

    def load(
        self,
        path: str = None
    ) -> AlignmentsERCC:
        if not os.path.exists(path):
            return AlignmentsERCC(
                ercc_counts=self.ercc_counts,
                ercc_reference=self.ercc_reference
            )
        data = load_data(path)
        return AlignmentsERCC.fromdict(data)

    @classmethod
    def execute(
        cls,
        argv
    ) -> None:
        parser = argparse.ArgumentParser(
            description='Gather mapping stats from a SAM/BAM file'
        )
        
        parser.add_argument(
            '-s',
            '--SAM',
            default=sys.stdin,
            help="Input sam/bam file. (default: stdin)"
        )

        parser.add_argument(
            '-r',
            '--REF',
            required=True,
            help="The reference .FASTA for ERCC transcripts"
        )

        parser.add_argument(
            '-e',
            '--EXP',
            required=True,
            help="A .CSV file containing expected ERCC counts"
        )

        parser.add_argument(
            "-o", 
            "--OUT",
            default='-',
            help="Output SAM file. (default: stdout)"
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
            required=False,
            default=[],
            help="Other reference .FASTA against which to count alignments"
        )

        args = parser.parse_args(argv)
        cls(
            sam_path=args.SAM,
            out_path=args.OUT,
            json_path=args.JSON, 
            ercc_reference=args.REF,
            ercc_counts=args.EXP, 
            fasta_paths=args.FASTA
        )