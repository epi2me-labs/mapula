import sys
import json
from operator import add


def errprint(*args, **kwargs):
    """
    Prints, but to stderr
    """
    sys.stderr.write(*args, **kwargs)
    sys.stderr.write("\n")


def get_data_slots(*classes):
    """
    Method useful for extracting fields from dataclasses
    in order to apply them for speedups within __slots__.
    """
    fields = []
    fields_attr = "__dataclass_fields__"
    for clss in classes:
        if not hasattr(clss, fields_attr):
            continue
        fields = fields + list(getattr(clss, fields_attr).keys())
    return tuple(set(fields))


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


def add_dists(old, new, attr):
    """
    Given two lists, map over them and sum their elements
    together in a pairwise fashion and update old.attr with
    the result
    """
    result = list(map(add, getattr(old, attr), getattr(new, attr)))
    old.data[attr] = result


def add_attrs(old, new, *attrs):
    """
    Given a list of attribute names, old.attr and new.attr and
    update old.attr with the result.
    """
    for attr in attrs:
        result = getattr(old, attr) + getattr(new, attr)
        old.data[attr] = result


def get_group_name(
    name: str,
    run_id: str,
    barcode: str,
) -> str:
    """
    Returns a formatted string
    """
    return "{}-{}-{}".format(name, run_id, barcode)


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
