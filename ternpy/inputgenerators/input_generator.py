import os
import numpy as np
from collections import OrderedDict
import json
from ternpy.utils.stoichbalancer import balance
from ternpy.configcreator import phaseconfig


class InputGenerator:

    # configfile: path to the config file
    # ternary: path to file with ternary corners
    # ternary corner names must match names in allphases dict
    # list order is as follows
    # ['a', 'b', 'c'] corresponds to:
    #                c
    #               a b
    def __init__(self, ternary, configfile, ternary_name):
        # self.allphases = allphases
        self.allphases = phaseconfig.read_config(configfile)
        #self.ternary = self.read_ternary(ternary)
        self.ternary = ternary
        self.tern_phases, self.compositions = self.get_phases(self.ternary)
        coords = self.get_coords()
        self.projectdir = os.path.dirname(os.path.abspath(configfile))
        self.projectdir += '/'+ternary_name

        # and 'comp' (decomposition) and 'coord' (x,y coordinates)
        # keys to the dictionary and asign corresponding values
        for ph, comp, xy in zip(self.tern_phases, self.compositions, coords):
            self.tern_phases[ph]['comp'] = comp
            self.tern_phases[ph]['coords'] = xy

        self.data = self.load_data()
        # self.pressures is a  list of floats
        self.pressures = self.get_press_range()
        # save config file which contains only phases within ternary diagram
        self.save_config()
        # keep only data for pressures in self.pressures
        for phase in self.data:
            for poly in self.data[phase]:
                matrix = self.data[phase][poly]
                newmatrix = []
                for line in matrix:
                    if line[0] in self.pressures:
                        newmatrix.append(line)
                self.data[phase][poly] = np.array(newmatrix)

    # Returns list of pressure points which exists for every structure
    def get_press_range(self):
        pressures = self.data.itervalues().next().itervalues().next()['P']
        for data in self.data.itervalues():
            for entry in data.itervalues():
                pressures = list(set(pressures).intersection(entry['P']))
        return sorted(pressures)

    # read ternary file
    # returns list of ternary corners
    # NOT USED
    def read_ternary(self, ternary):
        with open(ternary, 'r') as f:
            return f.readline().split()

    # load energy data from files and returns dictionary
    # TODO load only data in the pressure range
    def load_data(self):
        d = {}
        direc = 'energies/'
        for ph in self.tern_phases.keys():
            try:
                d[ph] = {poly: np.genfromtxt(direc+ph+'-'+poly+'.dat', names=True) for
                         poly in self.tern_phases[ph]['structures']}
            except IOError:
                print('no file for', ph)
        return d

    def save_config(self):
        fname = (self.ternary[0] + '-' +
                 self.ternary[1] + '-' +
                 self.ternary[2] + '.json')
        fname = self.projectdir + '/' + fname
        print(fname)
        if not os.path.exists(self.projectdir):
            os.makedirs(self.projectdir)
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
                arr = [1, 0, 0, 1]
            elif idx == 1:
                arr = [0, 1, 0, 1]
            elif idx == 2:
                arr = [0, 0, 1, 1]
            return True, arr
        else:
            temp = ternary[:]
            temp.append(phase)
            data = self.get_phase_data(temp)
            test, arr = balance(data, 3)
            return test, arr

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
        num_of_molecules = sum(self.tern_phases[phase]['comp'][:3])
        pabc = self.tern_phases[phase]['comp']
        a = self.get_most_stable(self.ternary[0], p)
        b = self.get_most_stable(self.ternary[1], p)
        c = self.get_most_stable(self.ternary[2], p)
        a = pabc[0]*a
        b = pabc[1]*b
        c = pabc[2]*c
        phase = self.data[phase][poly]['H'][p]
        enthalpy = (phase-a-b-c) / (num_of_molecules*pabc[3])
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
            # p index is not always right
            # TODO
            arr = self.get_all_enthalpies(p)
            pos_h = []  # list to store indices of phases with rel_h > 0
            for i, row in enumerate(arr):
                if row[2] > 0.:
                    pos_h.append(i)
            arr = np.delete(arr, pos_h, axis=0)
            fname = self.projectdir + '/inputfiles/' + str(press)+'.in'
            if not os.path.exists(os.path.dirname(fname)):
                os.makedirs(os.path.dirname(fname))
            rows, cols = arr.shape
            if rows > 3:
                with open(fname, 'w') as f:
                    f.write('3      Pressure: '+str(press)+' kbar\n')
                    f.write(str(rows)+'\n')
                    np.savetxt(f, arr, fmt=['%.15f', '%.15f', '%.15f'])

if __name__ == '__main__':
    #IG = InputGenerator(ternary, configfile)
    # data = np.genfromtxt(sys.argv[1], names=True)
    # config = np.genfromtxt(sys.argv[2], usecols=(1, 2, 3))
    # IG = InputGenerator(config, data)
    #IG.generate_files()
    pass
