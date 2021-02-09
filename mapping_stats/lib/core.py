import abc
import pysam
from typing import TypeVar


class BaseSubcommand(abc.ABC):

    @abc.abstractclassmethod
    def execute():
        raise NotImplementedError

    @abc.abstractmethod
    def load():
        raise NotImplementedError
    
    @abc.abstractmethod
    def write():
        raise NotImplementedError


class UpdatingStatsItem(abc.ABC):

    @abc.abstractmethod
    def update(
        self,
        aln: pysam.AlignedSegment,
        **kwargs
    ):
        raise NotImplementedError

    @abc.abstractclassmethod
    def fromdict(
        cls,
        data: dict
    ):
        raise NotImplementedError

    @abc.abstractmethod
    def __add__(
        self, 
        new
    ):
        raise NotImplementedError


U = TypeVar('U', bound=UpdatingStatsItem)