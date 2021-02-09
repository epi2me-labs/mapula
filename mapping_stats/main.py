import sys
import argparse
from mapping_stats.ercc import ERCCMappingStats
from mapping_stats.gather import GatherMappingStats
from mapping_stats.aggregate import AggregateMappingStats


class MappingStats(object):

    def __init__(self):
        parser = argparse.ArgumentParser(
            description='Collects stats from SAM/BAM files',
            usage='''mapping-stats <command> [<args>]

Available subcommands are:
   gather       Gather mapping stats from a SAM/BAM file
   ercc         Same as gather, but also calculate ERCC correlations
   aggregate    Combine mapping stats .JSON outputs
'''
        )

        parser.add_argument('command', help='Subcommand to run')

        args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.command):
            print('Unrecognized command')
            parser.print_help()
            exit(1)

        getattr(self, args.command)(sys.argv[2:])

    def gather(self, argv):
        GatherMappingStats.execute(argv)

    def aggregate(self, argv):
        AggregateMappingStats.execute(argv)

    def ercc(self, argv):
        ERCCMappingStats.execute(argv)


def run_main():
    MappingStats()


if __name__ == '__main__':
    run_main()