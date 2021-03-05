import pysam
import numpy as np
from typing import Dict, List
from dataclasses import dataclass
from typing import TypedDict, Union
from scipy.stats import pearsonr, spearmanr
from mapula.const import UNKNOWN, UNMAPPED, UNCLASSIFIED
from mapula.bio import (
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
    name: str
    fasta: str
    length: int


@dataclass
class Observation(object):
    reference: str
    run_id: str
    read_group: str
    barcode: str
    fasta: str
    source: str
    length: int
    alignment_type: str
    quality: Union[None, float] = None
    coverage: Union[None, float] = None
    accuracy: Union[None, float] = None

    @classmethod
    def from_alignment(
        cls,
        aln: pysam.AlignedSegment,
        alignment_filename: str,
        references: Dict[str, TrackedReference],
        non_primary_stats=False
    ):
        reference = aln.reference_name or UNMAPPED
        run_id = get_alignment_tag(aln, 'RD', default=UNKNOWN)
        read_group = get_alignment_tag(aln, 'RG', default=UNKNOWN)
        barcode = get_alignment_tag(aln, 'BC', default=UNCLASSIFIED)
        fasta = getattr(references.get(reference, None), 'fasta', UNKNOWN)

        length = aln.query_length

        if aln.is_supplementary:
            alignment_type = "supplementary"
        elif aln.is_secondary:
            alignment_type = "secondary"
        elif aln.is_unmapped:
            alignment_type = "unmapped"
        else:
            alignment_type = "primary"

        accuracy = None
        quality = None
        coverage = None
        if alignment_type == "primary" or non_primary_stats:
            accuracy = get_alignment_accuracy(aln) or 0
            coverage = get_alignment_coverage(
                aln.reference_length, length) or 0
            quality = get_alignment_mean_qscore(
                aln.query_alignment_qualities) or 0

        return cls(
            reference, run_id, read_group, barcode, fasta,
            alignment_filename, length, alignment_type,
            quality, coverage, accuracy
        )


class ObservationGroup():

    def __init__(
        self,
        name: str,
        identity: dict,
        counts: Union[None, Dict[str, int]] = None,
        **defaults
    ) -> None:
        """
        Represents an instance of a reference sequence to
        which reads have been aligned, and tracks many
        with respect to the alignments.
        """
        self.name = name
        self.identity = identity
        self.counts = counts

        if self.identity.get('reference') and self.counts:
            raise ValueError(
                "Counts have been supplied but cannot be "
                "determined when splitting by reference."
            )

        self.base_pairs: int = defaults.get(
            'base_pairs', 0)
        self.observations: int = defaults.get(
            'observations', 0)
        self.primary_count: int = defaults.get(
            'primary_count', 0)
        self.secondary_count: int = defaults.get(
            'secondary_count', 0)
        self.supplementary_count: int = defaults.get(
            'supplementary_count', 0)

        self.read_n50: float = defaults.get(
            'read_n50', 0)
        self.cov80_count: int = defaults.get(
            'cov80_count', 0)
        self.cov80_percent: float = defaults.get(
            'cov80_percent', 0)
        self.median_accuracy: int = defaults.get(
            'median_accuracy', 0)
        self.median_quality: int = defaults.get(
            'median_quality', 0)

        self.alignment_accuracies: List[int] = defaults.get(
            'alignment_accuracies', [0] * 1001)
        self.alignment_coverages: List[int] = defaults.get(
            'alignment_coverages', [0] * 101)
        self.aligned_qualities: List[int] = defaults.get(
            'aligned_qualities', [0] * 600)
        self.read_lengths: List[int] = defaults.get(
            'aligned_qualities', [0] * 1000)

        if self.counts:
            self.spearmans_rho: float = defaults.get(
                'spearmans_rho', 0)
            self.spearmans_rho_pval: float = defaults.get(
                'spearmans_rho_pval', 0)
            self.pearson: float = defaults.get(
                'pearson', 0)
            self.pearson_pval: float = defaults.get(
                'pearson_pval', 0)

        if not self.identity.get('reference'):
            self.observed_references: dict[str, int] = defaults.get(
                'observed_references', {})
            self.observed_reference_count: int = len(
                self.observed_references)

    #
    # Speedups
    #
    __slots__ = (
        'name', 'base_pairs', 'observations',
        'primary_count', 'secondary_count', 'supplementary_count',
        'cov80_count', 'cov80_percent', 'median_accuracy',
        'median_quality', 'read_n50', 'alignment_accuracies',
        'alignment_coverages', 'aligned_qualities', 'read_lengths',
        'identity', 'counts', '__dict__'
    )

    def update(
        self,
        obs: Observation,
        update_summary_stats: bool = True
    ) -> None:

        if obs.alignment_type == "supplementary":
            self.supplementary_count += 1
            return

        if obs.alignment_type == "secondary":
            self.secondary_count += 1
            return

        self.base_pairs += obs.length
        self.observations += 1

        try:
            # E.g. LEN: 1000 bins, 50 bin width, 50,000 max length
            self.read_lengths[int(obs.length / 50)] += 1
        except IndexError:
            pass

        if obs.alignment_type == "unmapped":
            return

        self.primary_count += 1
        reference = obs.reference
        if not self.identity.get('reference'):
            try:
                self.observed_references[reference] += 1
            except KeyError:
                self.observed_references[reference] = 1

        try:
            self.aligned_qualities[int(obs.quality / 0.1)] += 1
        except (IndexError, TypeError):
            pass

        try:
            self.alignment_accuracies[int(obs.accuracy / 0.1)] += 1
        except (IndexError, TypeError):
            pass

        try:
            self.alignment_coverages[int(obs.coverage or 0)] += 1
            if obs.coverage >= 80:
                self.cov80_count += 1
        except (IndexError, TypeError):
            pass

        if update_summary_stats:
            self._update_summary_stats()

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
        for key, val in self.counts.items():
            expected.append(val)
            observed.append(self.observed_references.get(key, 0))

        if not sum(observed):
            return

        coef, p = spearmanr(observed, expected)
        self.spearmans_rho = coef
        self.spearmans_rho_pval = p

        coef2, p2 = pearsonr(observed, expected)
        self.pearson = coef2
        self.pearson_pval = p2

    def _update_summary_stats(self):
        self._update_cov80_percent()
        self._update_read_n50()
        self._update_median_quality()
        self._update_median_accuracy()

        if self.counts:
            self._update_correlations()

        if not self.identity.get('reference'):
            self.observed_reference_count = len(
                self.observed_references)

    def __add__(self, new):
        if not isinstance(new, self.__class__):
            raise TypeError

        self.base_pairs += new.base_pairs
        self.observations += new.observations
        self.primary_count += new.primary_count
        self.secondary_count += new.secondary_count
        self.supplementary_count += new.supplementary_count
        self.cov80_count += new.cov80_count

        self.alignment_accuracies = self._add_elementwise(
            self.alignment_accuracies, new.alignment_accuracies)
        self.alignment_coverages = self._add_elementwise(
            self.alignment_coverages, new.alignment_coverages)
        self.aligned_qualities = self._add_elementwise(
            self.aligned_qualities, new.aligned_qualities)
        self.read_lengths = self._add_elementwise(
            self.read_lengths, new.read_lengths)

        for k, v in new.observed_references.items():
            if self.observed_references.get(k):
                self.observed_references[k] += v
                continue
            self.observed_references[k] = v

        self._update_summary_stats()

        return self

    def _add_elementwise(self, arr1: List, arr2: List):
        return [a + b for a, b in zip(arr1, arr2)]

    @classmethod
    def from_dict(
        cls,
        data: dict,
        counts
    ):
        """
        Creates an AlignedReference object from
        a dictionary.
        """
        identity = {
            'fasta': data.pop('fasta', None),
            'run_id': data.pop('run_id', None),
            'barcode': data.pop('barcode', None),
            'read_group': data.pop('read_group', None),
            'reference': data.pop('reference', None)
        }

        obs_ident = {k: v for k, v in identity.items() if v}
        obs_name = "-".join(i for i in obs_ident.values())
        return cls(obs_name, obs_ident, counts=counts, **data)

    def to_dict(
        self,
        json=False
    ):
        """
        Serialises this object to a dictionary
        format.
        """
        optional_fields = {}

        optional_fields.update({
            "spearmans_rho": self.spearmans_rho,
            "spearmans_rho_pval": self.spearmans_rho_pval,
            "pearson": self.pearson,
            "pearson_pval": self.pearson_pval,
        } if self.counts else {})

        optional_fields.update({
            "alignment_accuracies": self.alignment_accuracies,
            "alignment_coverages": self.alignment_coverages,
            "aligned_qualities": self.aligned_qualities,
            "read_lengths": self.read_lengths,
        } if json else {})

        optional_fields.update({
            "observed_references": self.observed_references,
        } if json and not self.identity.get('reference') else {})

        optional_fields.update({
            "observed_reference_count": self.observed_reference_count,
        } if not self.identity.get('reference') else {})

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
            **optional_fields
        }
