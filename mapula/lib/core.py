import abc
import pysam
import dataclasses
from dataclasses import dataclass, field
from mapula.lib.refmap import RefMap


class BaseSubcommand(abc.ABC):
    @abc.abstractclassmethod
    def execute(cls, argv):
        raise NotImplementedError


@dataclass
class BaseMappingStatsContainer(abc.ABC):
    @abc.abstractmethod
    def update(self, aln: pysam.AlignedSegment, refmap: RefMap):
        raise NotImplementedError

    @abc.abstractclassmethod
    def fromdict(cls, **data: dict):
        raise NotImplementedError

    @abc.abstractmethod
    def asdict(self):
        raise NotImplementedError

    @abc.abstractmethod
    def add(self, new, refmap):
        raise NotImplementedError


@dataclass
class MappingStatsContainer(BaseMappingStatsContainer):
    _child_type = None
    children: dict = field(default_factory=lambda: {})

    @classmethod
    def fromdict(cls, **data: dict):
        for key, val in data.get("children", {}).items():
            data["children"][key] = cls._child_type.fromdict(**val)
        return cls(**data)

    def asdict(self):
        data = dataclasses.asdict(self)
        return {k: v for k, v in data.items() if not (k.startswith("_") or v is None)}
