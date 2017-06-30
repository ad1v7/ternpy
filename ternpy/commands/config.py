"""The config command."""


from .base import Base
from ternpy.configcreator import phaseconfig


class Config(Base):
    """Prepare config file"""

    def run(self):
        # phaselist is a list created by user
        if self.options['-p'] or self.options['--phaselist']:
            phaselist = self.options.get('FILE')
        else:
            phaselist = 'phaselist.conf'
        if self.options['-j'] or self.options['--dftdir']:
            dftdir = self.options.get('PATH')

        phaseconfig.extract(phaselist, dftdir)
