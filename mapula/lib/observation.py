import pysam
import numpy as np
from typing import Dict
from dataclasses import dataclass
from scipy.stats import pearsonr, spearmanr
from mapula.lib.bio import (
    BASE_ACCS,
    BASE_QUALS,
    get_alignment_tag,
    get_alignment_mean_qscore,
    get_median_from_frequency_dist,
    get_alignment_accuracy,
    get_n50_from_frequency_dist,
    get_alignment_coverage,
)


@dataclass
class TrackedReference(object):
    group: str
    filename: str
    length: int
    expected_count: int


class ObservationGroup():

    def __init__(
        self,
        # Ident
        name,
        identity,
        # Tracked
        base_pairs=0,
        observations=0,
        primary_count=0,
        secondary_count=0,
        supplementary_count=0,
        read_n50=0,
        cov80_count=0,
        cov80_percent=0,
        median_accuracy=0,
        median_quality=0,
        alignment_accuracies=None,
        alignment_coverages=None,
        aligned_qualities=None,
        read_lengths=None,
        # Optional (ident dependent)
        has_counts: bool = False,
        tracked_references: Dict[str, TrackedReference] = {},
        observed_references: Dict[str, int] = {}
    ) -> None:
        """
        Represents an instance of a reference sequence to
        which reads have been aligned, and tracks many
        with respect to the alignments.
        """
        self.name = name
        self.identity = identity

        self.base_pairs = base_pairs
        self.observations = observations
        self.primary_count = primary_count
        self.secondary_count = secondary_count
        self.supplementary_count = supplementary_count

        self.read_n50 = read_n50
        self.cov80_count = cov80_count
        self.cov80_percent = cov80_percent
        self.median_accuracy = median_accuracy
        self.median_quality = median_quality

        self.alignment_accuracies = (
            alignment_accuracies or list([0] * 1001))
        self.alignment_coverages = (
            alignment_coverages or list([0] * 101))
        self.aligned_qualities = (
            aligned_qualities or list([0] * 600))
        self.read_lengths = (
            read_lengths or list([0] * 1000))

        self.has_counts = has_counts
        self.spearmans_rho: float = 0
        self.spearmans_rho_pval: float = 0
        self.pearson: float = 0
        self.pearson_pval: float = 0

        self.observed_references = observed_references
        self.tracked_references = tracked_references

    #
    # Speedups
    #
    __slots__ = (
        'name', 'base_pairs', 'observations',
        'primary_count', 'secondary_count', 'supplementary_count',
        'cov80_count', 'cov80_percent', 'median_accuracy',
        'median_quality', 'alignment_accuracies', 'read_n50',
        'alignment_coverages', 'aligned_qualities', 'read_lengths',
        '__dict__', 'identity'
    )

    def update(
        self,
        aln: pysam.AlignedSegment,
    ) -> None:
        if aln.is_supplementary:
            self.supplementary_count += 1
            return

        if aln.is_secondary:
            self.secondary_count += 1
            return

        length = aln.query_length

        self.base_pairs += length
        self.observations += 1

        try:
            # E.g. LEN: 1000 bins, 50 bin width, 50,000 max length
            self.read_lengths[int(length / 50)] += 1
            self._update_read_n50()
        except IndexError:
            pass

        if aln.is_unmapped:
            return

        self.primary_count += 1
        reference = aln.reference_name

        try:
            self.observed_references[reference] += 1
        except KeyError:
            self.observed_references[reference] = 1

        try:
            quality = get_alignment_mean_qscore(
                aln.query_alignment_qualities)
            self.aligned_qualities[int(quality / 0.1)] += 1
            self._update_median_quality()
        except (IndexError, TypeError):
            pass

        try:
            accuracy = get_alignment_accuracy(aln) or 0
            self.alignment_accuracies[int(accuracy / 0.1)] += 1
            self._update_median_accuracy()
        except (IndexError, TypeError):
            pass

        try:
            coverage = get_alignment_coverage(
                aln.reference_length, length
            ) or 0
            self.alignment_coverages[int(coverage)] += 1

            if coverage >= 80:
                self.cov80_count += 1

            self._update_cov80_percent()
        except IndexError:
            pass

        if self.has_counts:
            self._update_correlations()

    def _update_read_n50(self):
        self.read_n50 = get_n50_from_frequency_dist(
            self.read_lengths, 50, self.base_pairs
        )

    def _update_median_accuracy(self):
        self.median_accuracy = get_median_from_frequency_dist(
            BASE_ACCS, np.array(self.alignment_accuracies)
        )

    def _update_cov80_percent(self):
        self.cov80_percent = (
            100 * self.cov80_count / self.observations
        )

    def _update_median_quality(self):
        self.median_quality = get_median_from_frequency_dist(
            BASE_QUALS, np.array(self.aligned_qualities)
        )

    def _update_correlations(self):
        observed = []
        expected = []
        for key, val in self.tracked_references.items():
            expected.append(val.expected_count)
            observed.append(self.observed_references.get(key, 0))

        coef, p = spearmanr(observed, expected)
        self.spearmans_rho = coef
        self.spearmans_rho_pval = p

        coef2, p2 = pearsonr(observed, expected)
        self.pearson = coef2
        self.pearson_pval = p2

    def __add__(self, new):
        if not isinstance(new, self.__class__):
            raise TypeError

        for attr in (
            "alignment_count",
            "read_count",
            "total_base_pairs",
            "primary_count",
            "secondary_count",
            "supplementary_count",
            "cov80_count"
        ):
            old = getattr(self, attr)
            new = getattr(self, attr)
            setattr(self, attr, old + new)

        self.alignment_accuracies = list(map(sum, zip(
            self.alignment_accuracies, new.alignment_accuracies)))
        self.alignment_coverages = list(map(sum, zip(
            self.alignment_coverages, new.alignment_coverages)))
        self.aligned_qualities = list(map(sum, zip(
            self.aligned_qualities, new.aligned_qualities)))
        self.read_lengths = list(map(sum, zip(
            self.read_lengths, new.read_lengths)))

        self._update_read_n50()
        self._update_median_quality()
        self._update_median_accuracy()
        self._update_cov80_percent()

        if self.has_counts:
            self.tracked_references.update(
                new.tracked_references)

            for k, v in new.observed_references.items():
                if self.observed_references.get(k):
                    self.tracked_references[k] += v
                    continue
                self.observed_references[k] = v

            self._update_correlations

        return self

    def to_dict(self, output_corrs=False):
        """
        Serialises this object to a dictionary
        format.
        """
        corrs = {
            "spearmans_rho": self.spearmans_rho,
            "spearmans_rho_pval": self.spearmans_rho_pval,
            "pearson": self.pearson,
            "pearson_pval": self.pearson_pval,
        } if output_corrs else {}

        return {
            **self.identity,
            "base_pairs": self.base_pairs,
            "observations": self.observations,
            "primary_count": self.primary_count,
            "secondary_count": self.secondary_count,
            "supplementary_count": self.supplementary_count,
            "read_n50": self.read_n50,
            "cov80_count": self.cov80_count,
            "cov80_percent": self.cov80_percent,
            "median_accuracy": self.median_accuracy,
            "median_quality": self.median_quality,
            "alignment_accuracies": self.alignment_accuracies,
            "alignment_coverages": self.alignment_coverages,
            "aligned_qualities": self.aligned_qualities,
            "read_lengths": self.read_lengths,
            **corrs
        }