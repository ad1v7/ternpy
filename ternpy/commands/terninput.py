"""The terninput command."""

from .base import Base
from ternpy.inputgenerators import input_generator


class TernInput(Base):
    """An tt command."""

    def run(self):
        configfile = self.options.get('FILE')
        ternary = self.options.get('TERNFILE')
        IG = input_generator.InputGenerator(ternary, configfile)
        IG.generate_files()
