import pysam
from typing import Union
from dataclasses import dataclass, field
from scipy.stats import pearsonr, spearmanr
from mapula.lib.refmap import RefMap
from mapula.lib.bio import (
    get_alignment_mean_qscore,
    get_median_from_frequency_dist,
    get_alignment_accuracy,
    get_n50_from_frequency_dist,
    get_alignment_coverage,
)


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

    def update_total_base_pairs(self, aln: pysam.AlignedSegment) -> None:
        self.total_base_pairs += aln.query_length


@dataclass
class AlignedReadQualityStats(object):
    """
    Dataclass providing statistics related
    to measurements of underlying read quality.
    """

    read_qualities: list = field(default_factory=lambda: [0] * 600)
    median_quality: float = 0

    def _update_read_quality_dist(self, quality: Union[float, int, None]) -> None:
        self.read_qualities[int(quality / 0.1)] += 1

    def update_read_quality_dist(self, aln: pysam.AlignedSegment) -> None:
        quality = get_alignment_mean_qscore(aln.query_qualities)
        self._update_read_quality_dist(quality)

    def update_median_quality(self) -> None:
        self.median_quality = get_median_from_frequency_dist(self.read_qualities, 0.1)


@dataclass
class AlignmentAccuracyStats(object):
    """
    Dataclass providing statistics related
    to measurements of alignment accuracy.
    """

    alignment_accuracies: list = field(default_factory=lambda: [0] * 1001)
    median_accuracy: float = 0

    def _update_alignment_accuracy_dist(self, accuracy: Union[float, int]) -> None:
        self.alignment_accuracies[int(accuracy / 0.1)] += 1

    def update_alignment_accuracy_dist(self, aln: pysam.AlignedSegment) -> None:
        accuracy = get_alignment_accuracy(aln) or 0
        self._update_alignment_accuracy_dist(accuracy)

    def update_median_accuracy(self) -> None:
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

    def _update_read_length_dist(self, length: int) -> None:
        self.read_lengths[int(length / 50)] += 1

    def update_read_length_dist(self, aln: pysam.AlignedSegment) -> None:
        self._update_read_length_dist(aln.query_length)

    def update_read_n50(self, total_base_pairs: int) -> None:
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

    def _update_alignment_coverage_dist(self, coverage: float) -> None:
        self.alignment_coverages[int(coverage)] += 1

    def update_alignment_coverage_dist(
        self, aln: pysam.AlignedSegment, reference_length: int
    ) -> None:
        coverage = (
            get_alignment_coverage(aln.query_alignment_length, reference_length) or 0
        )
        self._update_alignment_coverage_dist(coverage)

    def update_alignment_cov80(self, total_count: int) -> None:
        if total_count:
            total_count = sum(self.alignment_coverages)
            self.cov80_count = sum(self.alignment_coverages[80:])
            if self.cov80_count:
                self.cov80_percent = 100 * self.cov80_count / total_count


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
    importable dataclass, with a handy method
    to update all of the core stats on an object
    which inherits from this class.
    """

    @staticmethod
    def _update_core_stats(item, aln: pysam.AlignedSegment, refmap: RefMap) -> dict:
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
        item.update_read_n50(item.total_base_pairs)

        item.update_read_quality_dist(aln)
        item.update_median_quality()

        if aln.is_unmapped:
            return

        item.primary_count += 1

        item.update_alignment_accuracy_dist(aln)
        item.update_median_accuracy()

        length = refmap.get_ref_length(aln.reference_name)
        item.update_alignment_coverage_dist(aln, length)
        item.update_alignment_cov80(item.read_count)

        return item


@dataclass
class CorrelationStats(object):
    """
    Dataclass providing stats which measure
    the correlation between the expected and
    observed counts of alignments made to given
    reference sequences.

    Calculates spearmans_rho and pearsons, with
    associated p_values.
    """

    spearman: float = 0
    spearman_p: float = 0
    pearson: float = 0
    pearson_p: float = 0

    def update_correlations(self, references: dict, refmap: RefMap):
        expected_references = []
        for k, v in references.items():
            if exp_count := refmap.get_ref_expected_count(k):
                expected_references.append((v.read_count, exp_count))

        obs = [v[0] for v in expected_references]
        exp = [v[1] for v in expected_references]

        if len(obs) > 2:
            self.spearman, self.spearman_p = spearmanr(obs, exp)
            self.pearson, self.pearson_p = pearsonr(obs, exp)