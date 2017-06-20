import numpy as np
from Plots.phase import balance
from Plots.phase import balance_interface
from Plots.phase import get_balance_matrix
from collections import OrderedDict
#import sys

allphases = {'Brucite': {'atoms': ['Mg', 'O', 'H'], 'natoms': [1, 2, 2],
                      'Structures': {'name': ['P4', 'P3']}},
          'Enstatite': {'atoms': ['Mg', 'Si', 'O'], 'natoms': [1, 1, 3],
                        'Structures': {'name': ['Bridgmanite',
                                                    'clino',
                                                    'ilmenite']}},
          'PhaseB': {'atoms': ['Mg', 'Si', 'O', 'H'], 'natoms': [12, 4, 21, 2],
                     'Structures': {'name': ['PhaseB']}},
          'SiO2': {'atoms': ['Si', 'O'], 'natoms': [1, 2], 'Structures':
                   {'name': ['alpha',  'stishovite']}},
          'MgO': {'atoms': ['Mg', 'O'], 'natoms': [1, 1], 'Structures':
                  {'name': ['MgO']}},
          'H2O': {'atoms': ['O', 'H'], 'natoms': [1, 2], 'Structures':
                  {'name': ['Ice8']}}}

ternary = ['MgO', 'SiO2', 'H2O']
p4 = np.genfromtxt('P4.dat', names=True)
data = {'Brucite': {'P3': p4, 'P4': []}}

print data['Brucite']['P3']['H']

class InputGenerator:

    # allphases = dict created during data extraction from the VASP run
    # ternary = list of names for ternary corners (e.g. 'Mgo')
    # names must match names in allphases dict
    # list order is as follows
    # ['a', 'b', 'c'] corresponds to:
    #                c
    #               a b
    def __init__(self, allphases, ternary):
        self.allphases = allphases
        self.ternary = ternary
        self.tern_phases, self.compositions = self.get_phases(ternary)
        coords = self.get_coords()
        for phase, comp, xy in zip(self.tern_phases, self.compositions, coords):
            self.tern_phases[phase]['comp'] = comp
            self.tern_phases[phase]['coord'] =xy
        print '++++++++++++++++'
        print self.tern_phases.keys()
        #print self.tern_phases
        print '++++++++++++++++'
        #mydata = self.get_phase_data(['SiO2', 'MgO', 'H2O', 'Brucite'])
        #print balance_interface(mydata, 3)
        self.load_data()

        
    def load_data(self):
        newdata = OrderedDict() # just in case
        for phase in self.tern_phases.keys():
            for poly in self.tern_phases[phase]['Structures']['name']:
                try:
                    filename = poly+'.dat'
                    newdata[phase] = {poly: np.genfromtxt(filename, names=True)}
                except IOError:
                    print 'no file for', poly
        return newdata

    # returns all phases within ternary diagram
    # AND phase decomposition wrt to ternary corners
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
    # AND integer list representing its decompostion reaction
    # within given ternary diagram
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
        for phase in phaselist:
            try:
                phase_dict[phase] = self.allphases[phase]
            except KeyError:
                print('No key in main dict: ' + phase)
        return phase_dict

    # returns list of (x,y) pairs
    # same order as phases, compositions, etc
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
    # p - index of pressure point
    # a,b,c are h of constituents times number of molecules(or atoms)
    def enthalpy_of_formation(self, p, phase, poly):
        num_of_molecules = sum(self.tern_phases[phase]['comp'])
        abc = self.tern_phases[phase]['comp']
        a = self.get_most_stable(self.ternary[0], p)
        b = self.get_most_stable(self.ternary[1], p)
        c = self.get_most_stable(self.ternary[2], p)
        a = abc[0]*a
        b = abc[1]*b
        c = abc[2]*c
        phase = self.data[phase][poly]['H'][p]
        enthalpy = (phase-a-b-c)/num_of_molecules
        return enthalpy

    def get_most_stable(self, phase, p):
        return min([h for h in self.data[phase]])

    '''
    def enthalpy_of_formation(self, p, i):
        num_of_molecules = (
            self.config[i][0] + self.config[i][1] + self.config[i][2])
        a = self.config[i][0]*self.data[p][1]
        b = self.config[i][1]*self.data[p][2]
        c = self.config[i][2]*self.data[p][3]
        phase = self.data[p][i+1]
        enthalpy = (phase-a-b-c)/num_of_molecules
        return enthalpy
    '''

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
