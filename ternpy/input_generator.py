import numpy as np
from Plots.phase import balance
from Plots.phase import balance_interface
from Plots.phase import get_balance_matrix
from collections import OrderedDict
#import sys

allphases = {'Brucite': {'atoms': ['Mg', 'O', 'H'], 'natoms': [1, 2, 2],
                      'Structures': {'plotname': ['P4', 'P3']}},
          'Enstatite': {'atoms': ['Mg', 'Si', 'O'], 'natoms': [1, 1, 3],
                        'Structures': {'plotname': ['Bridgmanite',
                                                    'Clinoenstatite',
                                                    'HP-clino',
                                                    'ilmenite',
                                                    'ortho']}},
          'PhaseB': {'atoms': ['Mg', 'Si', 'O', 'H'], 'natoms': [12, 4, 21, 2],
                     'Structures': {'plotname': ['Phase B_struct']}},
          'SiO2': {'atoms': ['Si', 'O'], 'natoms': [1, 2], 'Structures':
                   {'plotname': ['alpha', 'coe', 'stishovite']}},
          'MgO': {'atoms': ['Mg', 'O'], 'natoms': [1, 1], 'Structures':
                  {'plotname': ['MgO']}},
          'H2O': {'atoms': ['O', 'H'], 'natoms': [1, 2], 'Structures':
                  {'plotname': ['H2O']}}}

ternary = ['MgO', 'SiO2', 'H2O']


class InputGenerator:

    def __init__(self, allphases, ternary):
        self.allphases = allphases
        self.tern_phases, self.compositions = self.get_phases(ternary)
        print self.compositions
        print '++++++++++++++++'
        print self.tern_phases.keys()
        print '++++++++++++++++'
        print self.get_coords()
        #mydata = self.get_phase_data(['SiO2', 'MgO', 'H2O', 'Brucite'])
        #print balance_interface(mydata, 3)

    # returns all phases within ternary diagram
    def get_phases(self, ternary):
        phases = []
        composition = []
        for phase in self.allphases.keys():
            test, arr = self.is_phase_in_ternary(phase, ternary)
            if test:
                phases.append(phase)
                composition.append(arr)
        return self.get_phase_data(phases), composition

    # returns True if phase is found in ternary list
    # and list representing its decompostion reaction
    # phase is a string
    # ternary is a list of strings
    def is_phase_in_ternary(self, phase, ternary):
        if phase in ternary:
            idx = ternary.index(phase)
            if idx == 0:
                arr = [1, 0, 0]
            elif idx == 1:
                arr = [0, 1, 0]
            elif idx == 2:
                arr = [0, 0, 1]
            return True, arr
        else:
            temp = ternary[:]
            temp.append(phase)
            data = self.get_phase_data(temp)
            test, arr = balance_interface(data, 3)
            return test, arr[:3]

    # returns 'full' phase record given list of phase names
    def get_phase_data(self, phaselist):
        phase_dict = OrderedDict()
        for corner in phaselist:
            try:
                phase_dict[corner] = self.allphases[corner]
            except KeyError:
                print('No key in main dict: '+corner)
        return phase_dict

    # returns list of (x,y) pairs
    # order of pairs match order in config file
    def get_coords(self):
        coords = []
        for c in self.compositions:
            x = 0.5 * (2. * c[1] + c[2]) / (c[0] + c[1] + c[2])
            y = np.sqrt(3) / 2. * c[2] / (c[0] + c[1] + c[2])
            coords.append((x, y))
        return coords

    #####################################################
    #                                                   #
    #  This is it so far, above should work             #
    #                                                   #
    #     Below still needs changes                     #
    #      to proceed further you will need to create   #
    #                                                   #
    #       simple database file and load it with       #
    #       genfromtext                                 #
    #                                                   #
    #####################################################

    # returns enthalpy of formation per molecule
    # p - index of pressure point, i row of phase in config file
    # a,b,c are h of constituents times number of molecules/atoms
    def enthalpy_of_formation(self, p, i):
        num_of_molecules = (
            self.config[i][0] + self.config[i][1] + self.config[i][2])
        a = self.config[i][0]*self.data[p][1]
        b = self.config[i][1]*self.data[p][2]
        c = self.config[i][2]*self.data[p][3]
        phase = self.data[p][i+1]
        enthalpy = (phase-a-b-c)/num_of_molecules
        # Below is needed to avoid later issues when generating
        # ternary graphs. The problem is that tricontour plot
        # goes wild when enthalpy is close to zero so better remove it
        #if enthalpy < -0.001:
        #    return enthalpy
        #elif enthalpy == 0:
        #    return 0.0
        #else:
        #    return 0.1
        return enthalpy

    # returns array containing enthalpies
    # of formation for all phases
    def get_all_enthalpies(self, p):
        rows, cols = self.config.shape
        arr = np.zeros((rows, 3))
        for i in range(rows):
            x, y = self.coords[i]
            arr[i][0] = x
            arr[i][1] = y
            arr[i][2] = self.enthalpy_of_formation(p, i)
        return arr

    # generate input files for each pressure point
    def generate_files(self):
        for idx, line in enumerate(self.data):
            p = str(int(line[0]))
            arr = self.get_all_enthalpies(idx)
            # remove rows with positive enthalpies
            rows, cols = self.config.shape
            pos_h = []
            for i in range(rows):
                if arr[i][2] > 0:
                    pos_h.append(i)
            arr = np.delete(arr, pos_h, axis=0)
            # save only files which contains at least one meta(?)stable phase
            rows, cols = arr.shape
            if rows > 3:
                header = str(rows)+' '+p+'\n'
                np.savetxt(p+'.in', arr, header=header, fmt=['%.15f', '%.15f', '%.15f'])

if __name__ == '__main__':
    print 'Test'
    IG = InputGenerator(allphases, ternary)
#data = np.genfromtxt(sys.argv[1], names=True)
#config = np.genfromtxt(sys.argv[2], usecols=(1, 2, 3))
#IG = InputGenerator(config, data)
#print(hasattr(IG, 'dataa'))
#IG.generate_files()
