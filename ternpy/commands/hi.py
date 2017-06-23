"""The hi command."""


from json import dumps

from .base import Base


class Hi(Base):
    """Say hello, world!"""

    def run(self):
        print('Hello, world!')
        print('You supplied the following options:', dumps(self.options, indent=2, sort_keys=True))
