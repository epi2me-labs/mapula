UNMAPPED = "unmapped"
UNKNOWN = "unknown"
UNCLASSIFIED = "unclassified"

DISTS = [
    'alignment_accuracies', 'alignment_coverages',
    'aligned_qualities', 'read_lengths'
]

CORRS = [
    'spearmans_rho', 'spearmans_rho_pval',
    'pearson', 'pearson_pval'
]


class Groupers(object):
    GROUP = 'group'
    RUN_ID = 'run_id'
    BARCODE = 'barcode'
    READ_GROUP = 'read_group'
    REFERENCE = 'reference'

    choices = [
        GROUP, RUN_ID, BARCODE,
        READ_GROUP, REFERENCE
    ]


class Format(object):
    CSV = 'csv'
    JSON = 'json'
    ALL = 'all'

    choices = [CSV, JSON, ALL]
