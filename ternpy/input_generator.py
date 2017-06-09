import numpy as np
import sys


class InputGenerator:
    #print data['O']
    #print data.dtype.names

    def __init__(self, config, data):
        self.config = config
        self.data = data
        self.coords = self.get_coords()

    # returns list of (x,y) pairs
    # order of pairs match order in config file
    def get_coords(self):
        coords = []
        for line in self.config:
            x = 0.5 * (2. * line[1] + line[2]) / (line[0] + line[1] + line[2])
            y = np.sqrt(3) / 2. * line[2] / (line[0] + line[1] + line[2])
            coords.append((x, y))
        return coords

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

#data = np.genfromtxt(sys.argv[1], names=True)
#config = np.genfromtxt(sys.argv[2], usecols=(1, 2, 3))
#IG = InputGenerator(config, data)
#print(hasattr(IG, 'dataa'))
#IG.generate_files()
