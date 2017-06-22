import numpy as np
from collections import OrderedDict
import json
from ternpy.utils.stoichbalancer import balance_interface
from ternpy.configcreator import phaseconfig

#from ...utils.stoichbalancer import balance_interface
#from ...configcreator import phaseconfig

#from .. import utils
#from utils import stoichbalancer
#from stoichbalancer import balance_interface

#from .. import configcreator
#from configcreator import phaseconfig

ternary = ['MgO', 'SiO2', 'Ice']


class InputGenerator:

    # allphases = dict created during data extraction from the VASP run
    # ternary = list of names for ternary corners (e.g. 'Mgo')
    # ternary corner names must match names in allphases dict
    # list order is as follows
    # ['a', 'b', 'c'] corresponds to:
    #                c
    #               a b
    def __init__(self, ternary):
        # self.allphases = allphases
        self.allphases = phaseconfig.read_config("configs")
        self.ternary = ternary
        self.tern_phases, self.compositions = self.get_phases(ternary)
        coords = self.get_coords()

        # and 'comp' (decomposition) and 'coord' (x,y coordinates)
        # keys to the dictionary and asign corresponding values
        for ph, comp, xy in zip(self.tern_phases, self.compositions, coords):
            self.tern_phases[ph]['comp'] = comp
            self.tern_phases[ph]['coords'] = xy
        self.data = self.load_data()
        self.pressures = self.data.itervalues().next().itervalues().next()['P']
        # save config file which contains only phases within ternary diagram
        self.save_config()
        # print(self.tern_phases)

    # load energy data from files and returns dictionary
    def load_data(self):
        d = {}
        direc = 'energies/'
        for ph in self.tern_phases.keys():
            try:
                d[ph] = {poly: np.genfromtxt(direc+poly+'.dat', names=True) for
                         poly in self.tern_phases[ph]['structures']}
            except IOError:
                print('no file for', ph)
        return d

    def save_config(self):
        fname = (self.ternary[0] + '-' +
                 self.ternary[1] + '-' +
                 self.ternary[2] + '.json')
        with open(fname, "w") as f:
            json.dump(self.tern_phases, f, indent=3)

    # returns all phases within ternary diagram
    # AND phases decomposition wrt to ternary corners
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

    # returns enthalpy of formation per molecule
    # p - index of pressure point
    # a,b,c are h of constituents times number of molecules(or atoms)
    #       c
    #      a b
    #
    def enthalpy_of_formation(self, phase, poly, p):
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

    # returns lowest value of enthalpy
    # for given phase and its polymorphs
    def get_most_stable(self, phase, p):
        return min([self.data[phase][poly]['H'][p] for poly in
                    self.data[phase]])

    # returns array containing enthalpies
    # of formation for all phases
    def get_all_enthalpies(self, p):
        #rows = sum(len(v['structures']['name']) for v in
        #           self.tern_phases.itervalues())
        #arr = np.zeros((rows, 3))
        arr = []
        for phase in self.tern_phases:
            for poly in self.tern_phases[phase]['structures']:
                x, y = self.tern_phases[phase]['coords']
                z = self.enthalpy_of_formation(phase, poly, p)
                arr.append([x, y, z])
        return arr

    # This function should be called
    # to generate files ready for convex hull plotter
    def generate_files(self):
        for p, press in enumerate(self.pressures):
            arr = self.get_all_enthalpies(p)
            pos_h = []
            for i, row in enumerate(arr):
                if row[2] > 0.:
                    pos_h.append(i)
            arr = np.delete(arr, pos_h, axis=0)
            fname = str(press)+'.in'
            rows, cols = arr.shape
            if rows > 3:
                with open(fname, 'w') as f:
                    f.write('3      Pressure: '+str(press)+' kbar\n')
                    f.write(str(rows)+'\n')
                    np.savetxt(f, arr, fmt=['%.15f', '%.15f', '%.15f'])

if __name__ == '__main__':
    IG = InputGenerator(ternary)
    # data = np.genfromtxt(sys.argv[1], names=True)
    # config = np.genfromtxt(sys.argv[2], usecols=(1, 2, 3))
    # IG = InputGenerator(config, data)
    IG.generate_files()
