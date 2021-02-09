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
    with open(path) as data:
        try:
            return json.load(data)
        except json.decoder.JSONDecodeError:
            print("Error loading data file {}.".format(path))
            raise


def write_data(path: str, data: dict) -> None:
    """
    Writes self._data out to a file located
    at self.json_path.
    """
    with open(path, "w") as out:
        json.dump(data, out)


def add_dists(old, new, attr):
    result = list(map(add, getattr(old, attr), getattr(new, attr)))
    setattr(old, attr, result)


def add_attrs(old, new, *attrs):
    for attr in attrs:
        result = getattr(old, attr) + getattr(new, attr)
        setattr(old, attr, result)


def get_group_name(
    name: str,
    run_id: str,
    barcode: str,
) -> str:
    return "{}-{}-{}".format(name, run_id, barcode)
