import ternpy
import sys
import numpy as np


def main():
    print(ternpy.joke())


def gen_input_files():
    data = np.genfromtxt(sys.argv[1], names=True)
    config = np.genfromtxt(sys.argv[2], usecols=(1, 2, 3))
    IG = ternpy.InputGenerator(config, data)
    print(hasattr(IG, 'dataa'))
    print('finito')
    IG.generate_files()
