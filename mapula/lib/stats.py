import pysam
from typing import Union, Dict
from mapula.lib.refmap import RefMap
from scipy.stats import pearsonr, spearmanr
from mapula.lib.utils import errprint, add_attrs, add_dists
from mapula.lib.bio import (
    get_alignment_mean_qscore,
    get_median_from_frequency_dist,
    get_alignment_accuracy,
    get_n50_from_frequency_dist,
    get_alignment_coverage,
)


class DataLayer(object):
    """
    A minimal object for holding key-value pairs
    in a property called data.
    """
    DEFAULTS = {}

    def __init__(self, **defaults) -> None:
        super().__init__()
        self._data = defaults

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, value):
        errprint("Do not set data directly")


class BaseAlignmentStats(DataLayer):
    """
    Dataclass providing basic statistics
    alignment counting statistics which can
    be incremented on each update.
    """
    ALN_COUNT = 'alignment_count'
    READ_COUNT = 'read_count'
    TOTAL_BASE_PAIRS = 'total_base_pairs'

    PRIMARY_COUNT = 'primary_count'
    SECONDARY_COUNT = 'secondary_count'
    SUPPLEMENTARY_COUNT = 'supplementary_count'

    ALN_ACCS = 'alignment_accuracies'
    MED_ACC = 'median_accuracy'
    READ_QUALS = 'read_qualities'
    MED_QUAL = 'median_quality'
    READ_LENS = 'read_lengths'
    READ_N50 = 'read_n50'

    DEFAULTS = {
        # Basic stats
        ALN_COUNT: 0,
        READ_COUNT: 0,
        TOTAL_BASE_PAIRS: 0,
        PRIMARY_COUNT: 0,
        SECONDARY_COUNT: 0,
        SUPPLEMENTARY_COUNT: 0,
        # Distributions
        ALN_ACCS: [0] * 1001,
        READ_QUALS: [0] * 600,
        READ_LENS: [0] * 1000,
        # Derived
        MED_ACC: 0,
        MED_QUAL: 0,
        READ_N50: 0
    }

    BASIC = (
        ALN_COUNT,
        READ_COUNT,
        TOTAL_BASE_PAIRS,
        PRIMARY_COUNT,
        SECONDARY_COUNT,
        SUPPLEMENTARY_COUNT
    )

    DISTS = (
        ALN_ACCS,
        READ_QUALS,
        MED_QUAL,
        READ_LENS
    )

    DERIVED = (
        MED_ACC,
        MED_QUAL,
        READ_N50
    )

    #
    # Alignment count
    #
    @property
    def alignment_count(self):
        return self.data[self.ALN_COUNT]

    @alignment_count.setter
    def alignment_count(self, value):
        self.data[self.ALN_COUNT] = value

    #
    # Read count
    #
    @property
    def read_count(self):
        return self.data[self.READ_COUNT]

    @read_count.setter
    def read_count(self, value):
        self.data[self.READ_COUNT] = value

    #
    # Base pairs
    #
    @property
    def total_base_pairs(self):
        return self.data[self.TOTAL_BASE_PAIRS]

    @total_base_pairs.setter
    def total_base_pairs(self, value):
        self.data[self.TOTAL_BASE_PAIRS] = value

    def update_total_base_pairs(
        self,
        aln: pysam.AlignedSegment
    ) -> None:
        self.total_base_pairs += aln.query_length

    #
    # Primary alignments
    #
    @property
    def primary_count(self):
        return self.data[self.PRIMARY_COUNT]

    @primary_count.setter
    def primary_count(self, value):
        self.data[self.PRIMARY_COUNT] = value

    #
    # Secondary alignments
    #
    @property
    def secondary_count(self):
        return self.data[self.SECONDARY_COUNT]

    @secondary_count.setter
    def secondary_count(self, value):
        self.data[self.SECONDARY_COUNT] = value

    #
    # Supplementary alignments
    #
    @property
    def supplementary_count(self):
        return self.data[self.SUPPLEMENTARY_COUNT]

    @supplementary_count.setter
    def supplementary_count(self, value):
        self.data[self.SUPPLEMENTARY_COUNT] = value

    #
    # Read qualities
    #
    @property
    def read_qualities(self):
        return self.data[self.READ_QUALS]

    @read_qualities.setter
    def read_qualities_setter(self, value):
        self.data[self.READ_QUALS] = value

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

    #
    # Median quality
    #
    @property
    def median_quality(self):
        return self.data[self.MED_QUAL]

    @median_quality.setter
    def median_quality(self, value):
        self.data[self.MED_QUAL] = value

    def update_median_quality(self):
        self.median_quality = get_median_from_frequency_dist(
            self.read_qualities, 0.1
        )

    #
    # Alignment accuracies
    #
    @property
    def alignment_accuracies(self):
        return self.data[self.ALN_ACCS]

    @alignment_accuracies.setter
    def alignment_accuracies_setter(self, value):
        self.data[self.ALN_ACCS] = value

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

    #
    # Median accuracy
    #
    @property
    def median_accuracy(self):
        return self.data[self.MED_ACC]

    @median_accuracy.setter
    def median_accuracy(self, value):
        self.data[self.MED_ACC] = value

    def update_median_accuracy(self):
        self.median_accuracy = get_median_from_frequency_dist(
            self.alignment_accuracies, 0.1
        )

    #
    # Read lengths
    #
    @property
    def read_lengths(self):
        return self.data[self.READ_LENS]

    @read_lengths.setter
    def read_lengths(self, value):
        self.data[self.READ_LENS] = value

    def _update_read_length_dist(
        self,
        length: int
    ):
        self.read_lengths[int(length / 50)] += 1

    def update_read_length_dist(
        self,
        aln: pysam.AlignedSegment
    ) -> None:
        self._update_read_length_dist(aln.query_length)

    #
    # Read N50
    #
    @property
    def read_n50(self):
        return self.data[self.READ_N50]

    @read_n50.setter
    def read_n50(self, value):
        self.data[self.READ_N50] = value

    def update_read_n50(self):
        self.read_n50 = get_n50_from_frequency_dist(
            self.read_lengths, 50, self.total_base_pairs
        )

    #
    # Helper methods for updating the
    # above properties.
    #
    def update_basic_stats(
        self,
        aln: pysam.AlignedSegment
    ) -> dict:
        self.alignment_count += 1

        if aln.is_supplementary:
            self.supplementary_count += 1
            return

        if aln.is_secondary:
            self.secondary_count += 1
            return

        self.read_count += 1
        self.update_total_base_pairs(aln)

        self.update_read_length_dist(aln)
        self.update_read_n50()

        self.update_read_quality_dist(aln)
        self.update_median_quality()

        if aln.is_unmapped:
            return

        self.primary_count += 1

        self.update_alignment_accuracy_dist(aln)
        self.update_median_accuracy()

        return self

    def add_basic_stats(self, new):
        add_attrs(
            self,
            new,
            "alignment_count",
            "read_count",
            "total_base_pairs",
            "primary_count",
            "secondary_count",
            "supplementary_count",
        )

        add_dists(self, new, "alignment_accuracies")
        add_dists(self, new, "read_qualities")
        add_dists(self, new, "read_lengths")

        self.update_read_n50()
        self.update_median_quality()
        self.update_median_accuracy()

        return self


class AlignedCoverageStats(DataLayer):
    """
    Dataclass providing stats which measure
    how much of the target (or reference)
    sequence the alignments of the queries
    are covering.

    cov80 count and percent report the proportion 
    of alignments that cover 80% or more of their 
    targets.
    """
    ALN_COVS = 'alignment_coverages'
    COV80_COUNT = 'cov80_count'
    COV80_PERCENT = 'cov80_percent'

    DEFAULTS = {
        ALN_COVS: [0] * 101,
        COV80_COUNT: 0,
        COV80_PERCENT: 0
    }

    @property
    def alignment_coverages(self):
        return self.data[self.ALN_COVS]

    @alignment_coverages.setter
    def alignment_coverages(self, value):
        self.data[self.ALN_COVS] = value

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
        coverage = (
            get_alignment_coverage(
                aln.reference_length,
                reference_length
            ) or 0
        )
        self._update_alignment_coverage_dist(coverage)

    @property
    def cov80_count(self):
        return self.data[self.COV80_COUNT]

    @cov80_count.setter
    def cov80_count(self, value):
        self.data[self.COV80_COUNT] = value

    @property
    def cov80_percent(self):
        return self.data[self.COV80_PERCENT]

    @cov80_percent.setter
    def cov80_percent(self, value):
        self.data[self.COV80_PERCENT] = value

    def update_alignment_cov80(self, total_count: int) -> None:
        self.cov80_count = sum(self.alignment_coverages[80:])
        self.cov80_percent = 100 * self.cov80_count / total_count

    #
    # Helper methods for updating the
    # above properties.
    #
    def update_coverage_stats(
        self,
        aln: pysam.AlignedSegment,
        refmap: RefMap,
        read_count: int
    ) -> dict:
        if (
            aln.is_supplementary 
            or aln.is_secondary 
            or aln.is_unmapped
        ):
            return

        length = refmap.get_ref_length(aln.reference_name)
        self.update_alignment_coverage_dist(aln, length)
        self.update_alignment_cov80(read_count)

        return self

    def add_coverage_stats(self, new, read_count):
        add_dists(self, new, "alignment_coverages")
        self.update_alignment_cov80(read_count)

        return self


class CorrelationStats(DataLayer):
    """
    Dataclass providing stats which measure
    the correlation between the expected and
    observed counts of alignments made to given
    reference sequences.

    Calculates spearmans_rho and pearsons, with
    associated p_values.
    """
    SPEARMAN = 'spearman'
    SPEARMAN_P = 'spearman_p'
    PEARSON = 'pearson'
    PEARSON_P = 'pearson_p'

    DEFAULTS = {
        SPEARMAN: 0,
        SPEARMAN_P: 0,
        PEARSON: 0,
        PEARSON_P: 0
    }

    CORRS = (SPEARMAN, SPEARMAN_P, PEARSON, PEARSON_P)

    @property
    def spearman(self):
        return self.data[self.SPEARMAN]

    @spearman.setter
    def spearman(self, value):
        self.data[self.SPEARMAN] = value

    @property
    def spearman_p(self):
        return self.data[self.SPEARMAN_P]

    @spearman_p.setter
    def spearman_p(self, value):
        self.data[self.SPEARMAN_P] = value

    @property
    def pearson(self):
        return self.data[self.PEARSON]

    @pearson.setter
    def pearson(self, value):
        self.data[self.PEARSON] = value

    @property
    def pearson_p(self):
        return self.data[self.PEARSON_P]

    @pearson_p.setter
    def pearson_p(self, value):
        self.data[self.PEARSON_P] = value

    #
    # Helper methods for updating the
    # above properties.
    #
    def update_correlation_stats(
        self,
        observed: Dict[str, int],
        expected: Dict[str, int]
    ):
        exp = [v for k, v in expected.items()]
        obs = [
            observed[k] if observed.get(k) else 0
            for k, v in expected.items()
        ]

        if not (exp and obs):
            return

        self.spearman, self.spearman_p = spearmanr(obs, exp)
        self.pearson, self.pearson_p = pearsonr(obs, exp)
