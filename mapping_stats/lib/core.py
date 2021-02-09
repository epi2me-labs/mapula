import abc
from typing import TypeVar


class BaseSubcommand(abc.ABC):
    @abc.abstractclassmethod
    def execute():
        raise NotImplementedError


class UpdatingStatsItem(abc.ABC):

    # The signature of this method is
    # expected to change when implemented by
    # a derived class, hence follows the
    # Liskov substitution principle.
    @abc.abstractmethod
    def update(self, *args, **kwargs):
        raise NotImplementedError

    @abc.abstractclassmethod
    def fromdict(cls, data: dict):
        raise NotImplementedError

    @abc.abstractmethod
    def __add__(self, new):
        raise NotImplementedError


U = TypeVar("U", bound=UpdatingStatsItem)
