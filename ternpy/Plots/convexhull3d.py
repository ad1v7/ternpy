import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import time
import math
from scipy.spatial import ConvexHull
from matplotlib.tri import Triangulation
import re

# consider adding title
# where corner labels are comming from
# sort colorbar issues
# docstrings and unit test...


class ConvexHullData:
    def __init__(self, file_list):
        self.data, self.pressures, self.filenames = self.load_data(file_list)
        self.convex_hulls = self.get_all_convex_hulls(self.data)
        self.minz = self.min_z(self.data, self.convex_hulls)

    # returns array containg xyz data and pressures
    def load_data(self, file_list):
        new_file_list = []
        data = []
        pressures = []
        for entry in file_list:
            try:
                data.append(np.genfromtxt(entry, delimiter="", skip_header=2))
                with open(entry, 'r') as f:
                    press = f.readline().split()[3]
                    f.close()
                pressures.append(press)
                new_file_list.append(entry)
            except ValueError:
                print("Oops! Looks like a wrong file format... " +
                      os.path.basename(str(entry)))
            except IOError:
                print("Ommiting directory... " +
                      os.path.basename(entry))
        return data, self.kbar_to_gpa(pressures), new_file_list

    def kbar_to_gpa(self, pressures):
        return [str(float(item)/10) for item in pressures]

    # returns array containing convex hull objects
    def get_all_convex_hulls(self, data):
        hulls = []
        for entry in data:
            hulls.append(ConvexHull(entry))
        return hulls

    # returns the value of the lowest point (z coordinate)
    # from the convex hull list
    def min_z(self, data, hulls):
        minz = 0.
        for idx, hull in enumerate(hulls):
            if minz > min(data[idx][hull.vertices, 2]):
                minz = min(data[idx][hull.vertices, 2])
        return minz


class PlotConvexHull:
    def __init__(self, ConvexHullData, **kwargs):
        self.CHD = ConvexHullData
        self.steps = kwargs.get('steps', 15)
        self.showaxis = kwargs.get('showaxis', 'off')
        self.showpressure = kwargs.get('showpressure', True)
        self.savefig = kwargs.get('savefig', True)
        self.showfig = kwargs.get('showfig', False)
        self.fontsize = kwargs.get('fontsize', 20)
        self.showbar = kwargs.get('showbar', False)
        self.contour_step = abs(self.CHD.minz)/self.steps
        if self.savefig:
            self.newdir = time.strftime("%d-%m-%y_%H:%M:%S")
            os.mkdir(self.newdir, 0o755)
        # appending 2 small values to avoid issues (white background)
        # when convex hull z value is close to 0
        self.levels = np.append(np.arange(self.CHD.minz, -0, 5e-4),
                                [-1e-8, -1e-9, -0])
        self.labels = self.get_nice_chem_formulas(['MgO', 'SiO2', 'H2O'])

    def savefigure(self, filename):
        savefile = os.path.basename(filename)
        savefile = os.path.splitext(savefile)[0] + '.png'
        plt.savefig(self.newdir+"/"+savefile, bbox_inches='tight', dpi=300)
        print (os.path.basename(filename) +
               " output saved to " + self.newdir+"/"+savefile)

    def colorbar(self, cax):
        # colorbar label
        label_cb = '$\Delta$h of formation per molecule / eV'
        if self.showbar:
            cb = plt.colorbar(cax, format='%.2f', orientation='vertical')
            cb.set_label(label_cb, fontsize=int(self.fontsize/1.3))
            cb.ax.tick_params(labelsize=14)
            cb.ax.invert_xaxis()

    def get_nice_chem_formulas(self, labels):
        labels = [re.sub("([0-9])", "_\\1", k) for k in labels]
        return ['$\mathregular{'+k+'}$' for k in labels]

    def plot_all(self, ConvexHullData):
        for i in range(len(ConvexHullData.data)):
            data = ConvexHullData.data[i]
            hull_data = ConvexHullData.convex_hulls[i]
            pressure = ConvexHullData.pressures[i]
            filename = ConvexHullData.filenames[i]
            self.plot(data, hull_data, pressure, filename)

    def plot(self, data, hull_data, pressure, filename):
        x, y, z = data.T
        simplices = hull_data.simplices
        # spacing between contours based on min value of z
        contour_spacing = 1 + int(abs(min(z))/self.contour_step)
        # triangulation based on the convex hull
        tri = Triangulation(x, y, triangles=simplices)
        fig, ax = plt.subplots()
        # plot lines (triangles)
        ax.triplot(tri, color='k', linewidth=1.)
        # set diagram window size
        ax.axis([-0.1, 1.1, -0.05, 0.95])
        # show x-y axis
        ax.axis(self.showaxis)
        # plot contour lines based on convex hull facets
        ax.tricontour(x, y, z, contour_spacing, linewidths=0.5,
                      colors='k', triangles=simplices)
        # plot colour map
        cax = ax.tricontourf(x, y, z, cmap=plt.get_cmap('plasma'),
                             triangles=simplices, levels=self.levels)
        # plot convex hull vertices.
        ax.scatter(data[hull_data.vertices, 0],
                   data[hull_data.vertices, 1],
                   marker='o', c='k', s=30, zorder=10)
        # add ternary diagram labels
        ax.text(0, -0.01, self.labels[0], fontsize=self.fontsize,
                ha='center', va='top')
        ax.text(1, -0.01, self.labels[1], fontsize=self.fontsize,
                ha='center', va='top')
        ax.text(0.5, math.sqrt(3)/2+0.02, self.labels[2],
                fontsize=self.fontsize, ha='center', va='bottom')
        # add colorbar
        self.colorbar(cax)
        # add presure label
        if self.showpressure:
            ax.text(0.85, 0.8, 'P='+pressure+' GPa',
                    fontsize=self.fontsize, ha='center', va='center')
        if self.savefig:
            self.savefigure(filename)
        if self.showfig:
            plt.show()


if __name__ == '__main__':
    hull = ConvexHullData(sys.argv[1:])
    hull_plotter = PlotConvexHull(hull, steps=12)
    hull_plotter.plot_all(hull)
