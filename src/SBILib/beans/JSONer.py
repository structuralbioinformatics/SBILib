import json
from abc import ABCMeta, abstractmethod


class JSONer(object, metaclass=ABCMeta):
    @abstractmethod
    def as_dict(self):
        pass

    def json(self, pretty=False):
        if pretty:
            return json.dumps(self.as_dict(), indent=2, separators=(',', ':'))
        return json.dumps(self.as_dict(), separators=(',', ':'))
