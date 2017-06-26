"""
ternpy

Usage:
  ternpy hi
  ternpy config [((-p | --phaselist) FILE)] [((-j | --dftdir) PATH)]
            [((-d | --projectdir) PROJECTDIR)]
  ternpy extract
  ternpy -h | --help
  ternpy --version

Options:
  -h --help                         Show this screen.
  --version                         Show version.

Examples:
  ternpy config -p PhaseListFile -j ~/MyData -d ~/MyProject

Help:
  For help using this tool, please open an issue on the Github repository:
  https://github.com/ad1v7/ternpy
"""


from inspect import getmembers, isclass

from docopt import docopt

from . import __version__ as VERSION


def main():
    """Main CLI entrypoint."""
    import ternpy.commands
    options = docopt(__doc__, version=VERSION)

    # Here we'll try to dynamically match the command the user is trying to run
    # with a pre-defined command class we've already created.
    for (k, v) in options.items(): 
        if hasattr(ternpy.commands, k) and v:
            module = getattr(ternpy.commands, k)
            ternpy.commands = getmembers(module, isclass)
            command = [cmd[1] for cmd in ternpy.commands if cmd[0] != 'Base'][0]
            command = command(options)
            command.run()
