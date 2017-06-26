"""The config command."""


from .base import Base
from ternpy.configcreator import phaseconfig


class Config(Base):
    """Prepare config file"""

    def run(self):
        if self.options['-l']:
            phaselist = self.options.get('PHASELIST')
        else:
            phaselist = 'phaselist.conf'

        print('phaselist')
        #phaseconfig.create_config(phaselist, "joboutput", "configs")
        pass
