UNMAPPED = "unmapped"
UNKNOWN = "unknown"
UNCLASSIFIED = "unclassified"


class Pipers(object):
    SAM = 'sam'
    JSON = 'json'

    choices = [SAM, JSON]


class Groupers(object):
    SOURCE = 'source'
    FASTA = 'fasta'
    RUN_ID = 'run_id'
    BARCODE = 'barcode'
    READ_GROUP = 'read_group'
    REFERENCE = 'reference'

    choices = [
        SOURCE, FASTA, RUN_ID, BARCODE,
        READ_GROUP, REFERENCE
    ]


class Format(object):
    CSV = 'csv'
    JSON = 'json'
    ALL = 'all'

    choices = [CSV, JSON, ALL]
