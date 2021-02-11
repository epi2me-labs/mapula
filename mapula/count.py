import os
import sys
import pysam
import argparse
import pandas as pd
from typing import Union, Dict
from pysam import AlignmentFile
from mapula.lib.refmap import RefMap
from mapula.lib.bio import get_alignment_tag
from mapula.lib.const import UNKNOWN, UNMAPPED, UNCLASSIFIED
from mapula.lib.stats import (
    BaseAlignmentStats,
    CorrelationStats,
    AlignedCoverageStats
)
from mapula.lib.utils import (
    parse_cli_key_value_pairs,
    errprint,
    write_data,
    load_data,
    get_group_name
)


class AlignedReference(BaseAlignmentStats):
    NAME = 'name'
    LENGTH = 'length'

    IDENT = (
        NAME,
        LENGTH
    )

    DEFAULTS = {
        **BaseAlignmentStats.DEFAULTS,
    }

    def __init__(
        self,
        name: str,
        length: int,
        refmap: RefMap,
        **defaults
    ) -> None:
        """
        Represents an instance of a reference sequence to 
        which reads have been aligned, and tracks many 
        with respect to the alignments.
        """
        self.refmap = refmap

        super(AlignedReference, self).__init__(
            name=name,
            length=length,
            **(defaults or self.DEFAULTS)
        )

    #
    # Basic attributes
    #
    @property
    def name(self):
        return self._data[self.NAME]

    @name.setter
    def name(self, value):
        self._data[self.NAME] = value

    @property
    def length(self):
        return self._data[self.LENGTH]

    @length.setter
    def length(self, value):
        self._data[self.LENGTH] = value

    #
    # Methods
    #
    def update(self, aln: pysam.AlignedSegment) -> None:
        """
        Given the input alignment, re-calculates
        basic statistics.
        """
        self.update_basic_stats(aln)

    @classmethod
    def fromdict(cls, refmap: RefMap, **data: dict):
        """
        Creates an AlignedReference object from
        a dictionary.
        """
        name = data.pop(cls.NAME)
        length = data.pop(cls.LENGTH)
        return cls(name, length, refmap, **data)

    def todict(self, output_dists=True):
        """
        Serialises this object to a dictionary
        format.
        """
        data = {**self.data}
        if not output_dists:
            data = {
                k: v for k, v in data.items()
                if k not in self.DISTS
            }

        return data

    def __add__(self, new):
        """
        Given another AlignedReference, adds the
        values of its stats properties to self.
        """
        return self.add_basic_stats(new)


class AlignmentGroup(
    BaseAlignmentStats,
    AlignedCoverageStats,
    CorrelationStats
):
    NAME = 'name'
    RUN_ID = 'run_id'
    BARCODE = 'barcode'
    REFERENCES = 'references'

    DEFAULTS = {
        **AlignedCoverageStats.DEFAULTS,
        **BaseAlignmentStats.DEFAULTS,
        **CorrelationStats.DEFAULTS
    }

    IDENT = (
        NAME, RUN_ID, BARCODE
    )

    def __init__(
        self,
        name: str,
        run_id: str,
        barcode: str,
        refmap: RefMap,
        use_corrs: bool = False,
        references={},
        **data
    ) -> None:
        """
        Represents a set of binned alignments by reference 
        name, run_id and barcode if available. Tracks 
        both summary statistics (i.e. representing the whole 
        group) as well as per-refefence statistics (i.e. 
        broken down by individual aligned reference sequences).
        """
        self.refmap = refmap
        self.references = references
        self.use_corrs = use_corrs

        if not data:
            data = {
                **AlignedCoverageStats.DEFAULTS,
                **BaseAlignmentStats.DEFAULTS,
            }

            if use_corrs:
                data.update(
                    CorrelationStats.DEFAULTS
                )

        super(AlignmentGroup, self).__init__(
            name=name, run_id=run_id, barcode=barcode, **data
        )

    #
    # Basic attributes
    #
    @property
    def name(self):
        return self._data[self.NAME]

    @name.setter
    def name(self, value):
        self._data[self.NAME] = value

    @property
    def run_id(self):
        return self._data[self.RUN_ID]

    @run_id.setter
    def run_id(self, value):
        self._data[self.RUN_ID] = value

    @property
    def barcode(self):
        return self._data[self.BARCODE]

    @barcode.setter
    def barcode(self, value):
        self._data[self.BARCODE] = value

    #
    # Methods
    #
    def update(self, aln: pysam.AlignedSegment):
        """
        Given the input alignment, re-calculates
        group and per-reference statistics.
        """
        self.update_references(aln)
        self.update_basic_stats(aln)
        self.update_coverage_stats(
            aln, self.refmap, self.read_count
        )
        if self.use_corrs and self.refmap.has_counts:
            self.update_correlations()

    def update_references(self, aln: pysam.AlignedSegment):
        """
        Given an alignment, finds or creates the 
        correct reference instance within this group
        and updates it.
        """
        ref_name = aln.reference_name or UNMAPPED
        reference = self.references.get(ref_name)
        if not reference:
            length = self.refmap.get_ref_length(
                ref_name) or 0
            reference = AlignedReference(
                name=ref_name,
                length=length,
                refmap=self.refmap
            )
            self.references[ref_name] = reference
        reference.update(aln)

    def update_correlations(self):
        """
        Handles updating correlation specific
        statistics.
        """
        obs = {
            k: v.read_count
            for k, v in self.references.items()
        }
        exp = {
            k: v['expected_count']
            for k, v in self.refmap.items()
        }

        self.update_correlation_stats(obs, exp)

    @classmethod
    def fromdict(
        cls,
        refmap: RefMap,
        use_corrs,
        **data: dict
    ):
        """
        Creates an AlignedGroup object and
        nested AlignedReference objects from
        a dictionary.
        """
        name = data.pop('name')
        run_id = data.pop('run_id')
        barcode = data.pop('barcode')
        references = data.pop('references')

        for key, val in references.items():
            rname = val.pop('name')
            rlength = val.pop('length')
            references[key] = AlignedReference(
                rname, rlength, refmap, **val
            )

        if not use_corrs:
            errprint(
                "[Warning]: expected counts have not been "
                "provided but the existing data contains "
                "correlation metrics. Quitting."
            )
            sys.exit(1)

        return cls(
            name, run_id, barcode, refmap,
            use_corrs, references, **data
        )

    def todict(self):
        """
        Serialises this object to a dictionary
        format.
        """
        data = {**self.data, self.REFERENCES: {}}
        for key, val in self.references.items():
            data[self.REFERENCES][key] = val.todict()

        return data

    def __add__(self, new):
        """
        Given another AlignedGroup, adds the
        values of its stats properties to self.
        """
        for rn, rv in new.references.items():
            if not self.references.get(rn):
                self.references[rn] = rv
            else:
                self.references[rn] += rv
        self.add_basic_stats(new)
        self.add_coverage_stats(new, self.read_count)
        if self.use_corrs:
            self.update_correlations()
        return self


class CountMappingStats(object):

    def __init__(
        self,
        sam: str,
        refs: Dict[str, str],
        counts: Dict[str, str],
        sam_out: Union[str, None],
        json_path: str,
    ) -> None:
        """
        A subcommand that runs a process designed to 
        scan alignments made in SAM format and accumulate 
        many useful statistics which are binned into groups 
        and reported in JSON format.
        """
        errprint("--- Running: Mapula (count)")

        self.sam = sam
        self.refs = parse_cli_key_value_pairs(refs)
        self.counts = parse_cli_key_value_pairs(counts)
        self.use_corrs = bool(counts)
        self.sam_out = sam_out
        self.json_path = json_path

        try:
            self.records = AlignmentFile(sam, "r")
        except OSError:
            errprint(
                '[Error]: Could not find SAM: {}'.format(
                    sam
                )
            )
            sys.exit(1)

        self.refmaps = self.load_refmaps()
        self.groups = self.load_groups()

        if sam_out:
            self.outfile = AlignmentFile(
                sam_out, "w", template=self.records
            )

        errprint("--- Generating stats")

        self.update_groups()

        errprint("--- Writing data")

        self.write_stats_to_json()
        self.write_groups_to_csv()
        self.write_refs_to_csv()

        errprint("--- Operation successful")

    def load_refmaps(self):
        """
        Given the reference and counts paths provided, 
        iterates over each pair and builds RefMap objects. 
        In addition, makes a RefMap for unmapped.
        """
        refmaps = {}
        for name, ref in self.refs.items():
            counts = self.counts.get(name)
            refmaps[name] = RefMap(
                name, ref, counts,
            )
        refmaps[UNMAPPED] = RefMap(UNMAPPED)
        return refmaps

    def load_groups(self) -> dict:
        """
        Instantiate a working state from a .json file, 
        if path exists, otherwise returns an empty dict.
        """
        if not os.path.exists(self.json_path):
            return {}

        errprint('[Notice]: Adding to existing data {}.'.format(
            self.json_path
        ))

        existing_data = load_data(self.json_path)
        for key, val in existing_data.items():
            try:
                name = val['name']
                refmap = self.refmaps[name]
                existing_data[key] = AlignmentGroup.fromdict(
                    refmap, self.use_corrs, **val)
            except KeyError:
                errprint(
                    "[Error]: existing .json data is malformed"
                )
                sys.exit(1)
        return existing_data

    def update_groups(
        self,
    ) -> None:
        """
        Iterates over the alignments contained within
        the input sam file or stream, and for each one
        finds the appropriate AlignmentGroup, and calls
        .update on that group to increment its stats.
        """
        for aln in self.records.fetch(until_eof=True):
            if self.sam_out:
                self.outfile.write(aln)

            self._get_or_create_group(aln).update(aln)

    def _get_or_create_group(self, aln: pysam.AlignedSegment):
        """
        Uses the alignment to discover the name of the
        aligned reference, and from that the group of
        references to which this alignment belongs. In
        addition, finds the run_id and barcode information
        if available for this alignment. Altogether, this
        data is used to match the correct AlignmentGroup
        from self.groups if we have already created it, or
        otherwise creates it on the spot.
        """
        reference = aln.reference_name
        run_id = get_alignment_tag(aln, "RD", UNKNOWN)
        barcode = get_alignment_tag(aln, "BC", UNCLASSIFIED)

        name = UNMAPPED
        for name, refmap in self.refmaps.items():
            if reference in refmap:
                break

        group_name = get_group_name(name, run_id, barcode)
        group = self.groups.get(group_name)

        if not group:
            group = AlignmentGroup(
                name=name, run_id=run_id, barcode=barcode,
                refmap=self.refmaps.get(name),
                use_corrs=self.use_corrs
            )
            self.groups[group_name] = group

        return group

    def write_stats_to_json(self) -> None:
        """
        Gathers the stats data together for each
        AlignedGroup in self.groups, and dumps it to file.
        """
        data = {
            key: val.todict()
            for key, val in self.groups.items()
        }
        write_data(self.json_path, data)

    def write_groups_to_csv(self) -> None:
        """
        Gathers the stats data together for each
        AlignedGroup in self.groups, and dumps 
        it to file.
        """
        mask = [
            *AlignmentGroup.IDENT,
            *AlignmentGroup.BASIC,
            *AlignmentGroup.DERIVED,
            *AlignmentGroup.CORRS
        ]

        data = []
        for grp in self.groups.values():
            data.append({
                i: j for i, j in grp.todict().items()
                if i in mask
            })

        df = pd.DataFrame(data)
        df = df.sort_values(
            ['read_count', 'total_base_pairs'], ascending=False)
        df = df.reset_index(drop=True)
        df = df.round({"cov80_percent": 2, "spearman": 2,
                       "spearman_p": 2, "pearson": 2, "pearson_p": 2})
        df.to_csv('groups.mapula.csv', index=False)

    def write_refs_to_csv(self) -> None:
        """
        Gathers the stats data together for each
        AlignedReference within each group in 
        self.groups, and dumps it to file.
        """
        mask = [
            *AlignedReference.IDENT,
            *AlignedReference.BASIC,
            *AlignedReference.DERIVED,
        ]

        data = []
        for grp in self.groups.values():
            for ref in grp.references.values():
                data.append({
                    i: j for i, j in ref.todict().items()
                    if i in mask
                })

        df = pd.DataFrame(data)
        df.sort_values(['read_count', 'total_base_pairs'], ascending=False)
        df.reset_index(drop=True)
        df.to_csv('refs.mapula.csv', index=False)

    @classmethod
    def execute(cls, argv) -> None:
        """
        Parses command line arguments and 
        initialises a CountMappingStats object.
        """
        parser = argparse.ArgumentParser(
            description="Count mapping stats from a SAM/BAM file",
            formatter_class=argparse.RawDescriptionHelpFormatter
        )

        parser.add_argument(
            "-s",
            "--sam",
            help='Alignments in SAM format. By default, this script '
            'reads alignments from stdin. However, using this flag '
            'it is possible to pass in a file path.',
            default=sys.stdin,
        )

        parser.add_argument(
            '-r',
            '--refs',
            help='Provide reference .fasta files using the syntax: '
            'name=path_to_ref.',
            default=[],
            nargs='*'
        )

        parser.add_argument(
            '-c',
            '--counts',
            help='Provide expected counts in csv format using the syntax: '
            'name=path_to_counts, where name should be equal to a name given '
            'to --references. Expected column headings: reference,expected_count. '
            'The reference column should contain the ID of a sequence in the '
            'corresponding reference file. The expected_count column should '
            'equal the expected number of observations for that sequence.',
            default=[],
            nargs='*'
        )

        parser.add_argument(
            "-o",
            "--sam_out",
            default=None,
            help='Outputs a sam file from the parsed alignments. Use - for '
            'piping out. (default: None)',
        )

        parser.add_argument(
            "-j",
            "--json_path",
            required=False,
            default="stats.mapula.json",
            help='Name of the output json (default: stats.mapula.json)',
        )

        args = parser.parse_args(argv)
        cls(
            sam=args.sam,
            refs=args.refs,
            counts=args.counts,
            sam_out=args.sam_out,
            json_path=args.json_path,
        )
