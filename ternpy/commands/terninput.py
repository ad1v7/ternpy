"""The terninput command."""

from .base import Base
from ternpy.inputgenerators import input_generator
from ternpy.configcreator import phaseconfig

try:
    input = raw_input
except NameError:
    pass


class TernInput(Base):
    """An tt command."""

    def run(self):
        configfile = self.options.get('FILE')
        allphases = phaseconfig.read_config(configfile)
        phasenames = [ph for ph in allphases]

        for idx, name in enumerate(phasenames):
            print(idx, name)

        print('       c       ')
        print('      / \      ')
        print('     a---b     ')
        idx = []
        while len(idx) != 3:
            idx = input('Select ternary corners (a b c): ').split()
            if len(idx) != 3:
                print('Select 3 numbers from the list separated by spaces')
        ternary = [phasenames[int(i)] for i in idx]
        print(ternary)
        ternary_name = input('Name ternary: ')

        IG = input_generator.InputGenerator(ternary, configfile, ternary_name)
        IG.generate_files()


