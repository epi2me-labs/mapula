import abc
import pysam
from functools import partial
from typing import Union, Dict
from mapping_stats.lib.core import U
from dataclasses import dataclass, field
from scipy.stats import pearsonr, spearmanr
from mapping_stats.lib.bio import (get_alignment_mean_qscore, 
    get_median_from_frequency_dist, get_alignment_accuracy, 
    get_n50_from_frequency_dist, get_alignment_coverage)


@dataclass
class BaseAlignmentStats(object):
    """
    Dataclass providing basic statistics 
    alignment counting statistics which can
    be incremented on each update.
    """
    alignment_count: int = 0
    total_base_pairs: int = 0
    read_count: int = 0
    primary_count: int = 0
    secondary_count: int = 0
    supplementary_count: int = 0

    def update_total_base_pairs(
        self, 
        aln: pysam.AlignedSegment
    ) -> None:
        self.total_base_pairs += aln.query_length


@dataclass
class AlignedReadQualityStats(object):
    """
    Dataclass providing statistics related
    to measurements of underlying read quality.
    """
    read_qualities: list = field(default_factory=lambda: [0] * 600)
    median_quality: float = 0

    def _update_read_quality_dist(
        self, 
        quality: Union[float, int, None]
    ) -> None:
        self.read_qualities[int(quality / 0.1)] += 1
    
    def update_read_quality_dist(
        self, 
        aln: pysam.AlignedSegment
    ) -> None:
        quality = get_alignment_mean_qscore(
            aln.query_qualities
        )
        self._update_read_quality_dist(quality)
    
    def update_median_quality(
        self
    ) -> None:
        self.median_quality = get_median_from_frequency_dist(
            self.read_qualities, 0.1
        )


@dataclass
class AlignmentAccuracyStats(object):
    """
    Dataclass providing statistics related
    to measurements of alignment accuracy.
    """
    alignment_accuracies: list = field(default_factory=lambda: [0] * 1001)
    median_accuracy: float = 0

    def _update_alignment_accuracy_dist(
        self, 
        accuracy: Union[float, int]
    ) -> None:
        self.alignment_accuracies[int(accuracy / 0.1)] += 1

    def update_alignment_accuracy_dist(
        self, 
        aln: pysam.AlignedSegment
    ) -> None:
        accuracy = get_alignment_accuracy(aln) or 0
        self._update_alignment_accuracy_dist(accuracy)

    def update_median_accuracy(
        self
    ) -> None:
        self.median_accuracy = get_median_from_frequency_dist(
            self.alignment_accuracies, 0.1
        )


@dataclass
class AlignedReadLengthStats(object):
    """
    Dataclass providing stats pertaining to
    the length of the reads (query sequences)
    with which the alignments have been made.

    read n50 is the read length at which 50% of
    all bases exist within reads of that length
    or greater.
    """
    read_lengths: list = field(default_factory=lambda: [0] * 1000)
    read_n50: int = 0
    
    def _update_read_length_dist(
        self, 
        length: int
    ) -> None:
        self.read_lengths[int(length / 50)] += 1

    def update_read_length_dist(
        self, 
        aln: pysam.AlignedSegment
    ) -> None:
        self._update_read_length_dist(aln.query_length)

    def update_read_n50(
        self,
        total_base_pairs: int
    ) -> None:
        self.read_n50 = get_n50_from_frequency_dist(
            self.read_lengths, 50, total_base_pairs
        )


@dataclass
class AlignedCoverageStats(object):
    """
    Dataclass providing stats which measure
    how much of the target (or reference) 
    sequence the alignments of the queries
    are covering.

    cov80 count and percent report the
    proportion of alignments that cover 80% or
    more of their targets.
    """
    alignment_coverages: list = field(default_factory=lambda: [0] * 101)
    cov80_percent: int = 0
    cov80_count: int = 0

    def _update_alignment_coverage_dist(
        self,
        coverage: float
    ) -> None:
        self.alignment_coverages[int(coverage)] += 1

    def update_alignment_coverage_dist(
        self, 
        aln: pysam.AlignedSegment,
        reference_length: int
    ) -> None:
        coverage = get_alignment_coverage(
            aln.query_alignment_length, reference_length
        ) or 0
    
        self._update_alignment_coverage_dist(coverage)

    def update_alignment_cov80(
        self,
        total_count: int
    ) -> None:
        if total_count:
            total_count = sum(self.alignment_coverages)
            self.cov80_count = sum(self.alignment_coverages[80:])
            if self.cov80_count:
                self.cov80_percent = (100 * self.cov80_count / total_count)


@dataclass
class CoreStats(
    AlignedCoverageStats,
    AlignedReadQualityStats,
    AlignmentAccuracyStats,
    AlignedReadLengthStats,
    BaseAlignmentStats,
):
    """
    Collects all the core stats into one easy,
    importable dataclass.
    """


@dataclass
class ERCCCorrelationStats(object):
    """
    Dataclass providing stats which measure
    the correlation between the expected and
    observed counts of alignments made to given
    reference sequences.

    Calculates spearmans_rho and pearsons, with
    associated p_values.
    """
    spearman = 0
    spearman_p = 0
    pearson = 0
    pearson_p = 0

    def _read_count_or_0(
        self,
        references: Dict[str, U],
        key: str
    ) -> int:
        result = references.get(key, 0)
        if result:
            result = result.read_count
        return result

    def update_correlations(
        self,
        references: Dict[str, U],
        ercc_counts: Dict[str, int]
    ):  
        par = partial(self._read_count_or_0, references)

        obs = [par(k) for k, v in ercc_counts.items()]
        exp = [v for k, v in ercc_counts.items()]

        if sum(obs) > 1:
            self.spearman, self.spearman_p = spearmanr(obs, exp)
            self.pearson, self.pearson_p = pearsonr(obs, exp)