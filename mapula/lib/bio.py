import math
import pysam
from typing import Union, Iterable, List
from mapula.lib.const import UNKNOWN

LOOKUP = [pow(10, -0.1 * q) for q in range(100)]


def get_alignment_tag(
    alignment: pysam.AlignedSegment, tag: str, default: str = UNKNOWN
):
    """
    Inspects the tags of the input AlignedSegment
    and checks for the presence of a barcode tag.
    If it exists, returns the barcode labeled within,
    else returns a string of 'Unknown'.
    """
    return alignment.get_tag(tag) if alignment.has_tag(tag) else default


def get_median_from_frequency_dist(arr: Iterable, width: Union[int, float]):
    """
    Returns the median value from an array whose
    positions represent frequency counts of values
    at that index in a given range.
    """
    arr_sum: int = sum(arr)
    half_way_pos = arr_sum / 2
    is_odd = bool(arr_sum % 2)

    lower = None

    accumulator = 0
    for idx, count in enumerate(arr):
        accumulator += count

        if accumulator > half_way_pos:
            if lower is not None:
                avg = ((lower * width) + (idx * width)) / 2
                return float(format(avg, ".2f"))

            return float(format(idx * width, ".2f"))

        if accumulator == half_way_pos:
            if is_odd:
                return float(format(idx * width, ".2f"))

            if lower is None:
                lower = idx
                continue


def get_alignment_accuracy(alignment: pysam.AlignedSegment):
    """
    Returns the percentage accuracy of a given aligned
    segment as a float.
    """
    if alignment.is_unmapped:
        return None

    el_count = [0] * 10
    for el in alignment.cigartuples:
        el_count[el[0]] += el[1]

    # Number of ambiguous bases
    nn = 0
    if alignment.has_tag("nn"):
        nn = alignment.get_tag("nn")

    # NM = #mismatches + #I + #D + #ambiguous_bases
    nm = alignment.get_tag("NM")
    dels = el_count[2]
    ins = el_count[1]

    mismatches = nm - dels - ins - nn
    matches = el_count[0] - mismatches

    accuracy = float(matches) / (matches + nm) * 100

    return accuracy


def get_n50_from_frequency_dist(
    arr: Iterable, width: int, total: int
) -> Union[float, int]:
    """
    Calculates the N50 from an array whose positions represent
    frequency counts of values at that index in a given range.

    N50 is defined as:

    Length N for which 50% of all bases in the sequences
    are in a sequence of length L < N

    Returns the approximate 'size' of the value at which
    the N50 is achieved.
    """
    n50 = 0
    cumulative_value = 0
    half_total = total / 2

    for bin_number, bin_count in enumerate(arr):
        bin_value = (bin_number * width) + (width / 2)

        cumulative_value += bin_value * bin_count
        if cumulative_value >= half_total:
            n50 = bin_value
            break

    return n50


def get_alignment_mean_qscore(scores: List[int]) -> Union[float, None]:
    """
    Returns the phred score corresponding to the mean of
    the probabilities associated with the phred scores
    provided.
    """
    if scores is None:
        return None

    if not scores:
        return 0.0

    sum_prob = 0.0
    for val in scores:
        sum_prob += LOOKUP[val]

    mean_prob = sum_prob / len(scores)

    return -10.0 * math.log10(mean_prob)


def get_alignment_coverage(query_alignment_length: int, reference_length: int):
    """
    Computes the percentage coverage of a given alignment
    length of the reference length.

    query_alignment_length: length of the aligned segment
    reference_length: length of the reference
    """
    if (query_alignment_length is None) or (reference_length is None):
        return None

    return 100 * float(query_alignment_length) / float(reference_length)
