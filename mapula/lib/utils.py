import sys
import json
import pysam
from operator import add


def errprint(*args, **kwargs):
    """
    Prints, but to stderr
    """
    sys.stderr.write(*args, **kwargs)
    sys.stderr.write("\n")


def load_data(
    path: str,
) -> dict:
    """
    Attempts to load json data from the path given.
    """
    with open(path) as data:
        try:
            return json.load(data)
        except json.decoder.JSONDecodeError:
            errprint("Error loading data file {}.".format(path))
            raise


def write_data(path: str, data: dict) -> None:
    """
    Writes self._data out to a file located
    at self.json_path.
    """
    with open(path, "w") as out:
        json.dump(data, out)


def parse_cli_key_value_pairs(arg):
    """
    Does what it says on the tin!
    """
    try:
        return {
            key: value for key, value in (
                item.split('=') for item in arg
            )
        }
    except ValueError:
        errprint('Error, invalid format: {}'.format(arg))
        errprint('Format of extra references must be key=value')
        sys.exit(1)


def get_total_alignments(
    samfile_path: str,
):
    alignments = pysam.AlignmentFile(samfile_path, "r")
    total = alignments.count(until_eof=True)
    alignments.close()
    return total