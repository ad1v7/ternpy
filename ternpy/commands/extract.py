"""The extract command."""

from .base import Base
from ternpy.configcreator import phaseconfig


class Extract(Base):
    """An extract command."""

    def run(self):
        configfile = self.options.get('FILE')
        jobdir = phaseconfig.read_dftdir(configfile)
        phaseconfig.create_datafiles(configfile, jobdir)
