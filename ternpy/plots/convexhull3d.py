import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import time
import math
import re

from scipy.spatial import ConvexHull
from matplotlib.tri import Triangulation

# consider adding title
# where corner labels are comming from
# sort colorbar issues
# docstrings and unit test...

# consider:
# self.ternary = ['Mgo', 'SiO2', 'Phase H']
# same as for input generator, names must correspond
# to phase names in original config file + match names of *.dat files
#
# How to identify phases in convex hull plot?
# one way is too look at x,y and possibly z coordinates
# and then get name from InputGenerator.tern_phases


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
                with open(entry, 'r') as f:
                    press = f.readline().split()[2]
                    # This is an idiot test to verify input files
                    teststring = f.readline().split()[1]
                    if teststring != 'Pressure:':
                        raise IndexError
                    f.close()
                data.append(np.genfromtxt(entry, delimiter="", skip_header=2))
                pressures.append(press)
                new_file_list.append(entry)
            except ValueError:
                print("Oops! Looks like a wrong file format... " +
                      os.path.basename(str(entry)))
            except IOError:
                print("Ommiting directory... " +
                      os.path.basename(entry))
            except IndexError:
                print('Something wrong with ' +
                      os.path.basename(str(entry)))
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
    def __init__(self, CHD, **kwargs):
        self.CHD = CHD
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

# Class to find and hold metastable points
# together with relevant functions
# such as finding decomposition reaction and its strength
class FindMetastable:
    def __init__(self, CHD):
        self.hulls = CHD.convex_hulls
        self.data = CHD.data
        self.pressures = CHD.pressures
        # below all have the same order and
        # should have the same number of elements
        self.vertices = self.get_all_vertices()
        self.points = self.get_all_points()
        self.metastable = self.get_all_metastable()
        self.triangles = self.get_all_triangles()

    # returns array of convex hull vertices
    def get_all_vertices(self):
        vertices = []
        for i, entry in enumerate(self.data):
            vertices.append(entry[self.hulls[i].vertices])
        return vertices

    # returns array of all points used to compute convex hull
    def get_all_points(self):
        points = []
        for hull in self.hulls:
            points.append(hull.points)
        return points

    # returns array of all metastable points 
    def get_all_metastable(self):
        metastable = []
        for i, vertices in enumerate(self.vertices):
            metastable.append(self.get_metastable(vertices, self.points[i]))
        return np.array(metastable)

    # returns array of metastable points for a given configuration
    def get_metastable(self, vertices, points):
        metastable = []
        for p in points:
            if not self.is_point_in_arr(p, vertices):
                self.is_smallest(p, metastable)
        return metastable

    # 2D test is pont in array; only x,y coords are used
    def is_point_in_arr(self, p, arr):
        for entry in arr:
            if (entry[0] == p[0]) and (entry[1] == p[1]):
                return True
        return False

    # test is point already in metastable array
    # by comparing x and y coords
    # if found it tests for value of z components
    # the point with lowest z component is saved in metastable
    def is_smallest(self, p, metastable):
        found_one = False
        for meta in metastable:
            if (p[0] == meta[0]) and (p[1] == meta[1]):
                if p[2] < meta[2]:
                    meta[2] = p[2]
                found_one = True
        if not found_one:
            metastable.append(p)

    # returns distance between point and plane ABC
    def point_distance_to_plane(self, point, tri):
        nhat = self.norm_to_plane(tri)
        return abs(np.dot(point, nhat))

    # returns unit vector normal to plane
    def norm_to_plane(self, tri):
        A, B, C = tri
        v1 = A - B
        v2 = A - C
        nhat = np.cross(v1, v2)
        return self.vec_to_unitvec(nhat)

    # point distance to line in 3D
    # can also return position on the line
    # http://www.fundza.com/vectors/point2line/index.html
    def point_distance_to_line(self, point, line):
        A, B = line
        line_vec = B - A
        point_vec = point - A
        linelength = np.linalg.norm(line_vec)
        line_unitvec = self.vec_to_unitvec(line_vec)
        point_vec_scaled = 1/linelength * point_vec
        test = np.dot(line_unitvec, point_vec_scaled)
        if test < 0.:
            test = 0.
        elif test > 1.:
            test = 1.
        nearest = test * line_vec
        dist = np.linalg.norm(nearest-point_vec)
        nearest = nearest + A
        return (dist, nearest)
        #return dist

    # converts given vector to unit vector
    def vec_to_unitvec(self, vec):
        return vec / np.sqrt(np.dot(vec, vec))

    # return true if points in 2D are collinear
    # ignore z coordinate
    def is_collinear2D(self, A, B, C):
        x1, y1 = B[0] - A[0], B[1] - A[1]
        x2, y2 = C[0] - A[0], C[1] - A[1]
        return abs(x1 * y2 - x2 * y1) == 0

    # returns true if point is on AB line segment
    # in 2D (ignore z coordinate)
    def is_point_between(self, point, A, B):
        dotproduct = np.dot(B-A, point-A)
        sqlengthAB = (B[0]-A[0])*(B[0]-A[0]) + (B[1]-A[1])*(B[1]-A[1])
        return (abs(np.cross(A-point, B-point)) == 0 and
                dotproduct > 0 and sqlengthAB > dotproduct)

    # returns true if point is found on one of triangle
    # edges and relevant vertices in 2D (ignore z coord)
    def is_point_on_tri_edge(self, point, tri):
        A, B, C = tri
        if self.is_point_between(point[:2], A[:2], B[:2]):
            return True, [point, A, B]
        elif self.is_point_between(point[:2], B[:2], C[:2]):
            return True, [point, B, C]
        elif self.is_point_between(point[:2], A[:2], C[:2]):
            return True, [point, A, C]
        else:
            return False, []

    # returns true either if point is inside triangle OR on triangle side
    def is_point_in_triangle(self, point, tri):
        A, B, C = tri
        if (self.same_side(point, A, B, C) and self.same_side(point, B, A, C)
                and self.same_side(point, C, A, B)):
            return True, [point, A, B, C]
        else:
            return False, []

    # return true only if point is inside triangle but NOT on its edges
    def is_point_inside_triangle(self, point, tri):
        A, B, C = tri
        intribool, temp = self.is_point_in_triangle(point, tri)
        ontriedgebool, temp2 = self.is_point_on_tri_edge(point, tri)
        if intribool and not ontriedgebool:
            return True, [point, A, B, C]
        else:
            return False, []

    # tests are points p1 and p2 on the same side
    # of the line segment ab
    def same_side(self, p1, p2, a, b):
        cp1 = np.cross(b[:2]-a[:2], p1[:2]-a[:2])
        cp2 = np.cross(b[:2]-a[:2], p2[:2]-a[:2])
        if np.dot(cp1, cp2) >= 0:
            return True
        else:
            return False

    # returns array of convex hull triangles (simplices)
    # triangles which are perpendicular to the z=0 plane are removed
    # in other words only triangles from the projection onto
    # z plane are saved.
    # Also in case where # of triangles > 1 the triangle which corresponds
    # to ternary diagram vertices is removed to simplify later calculations
    def get_all_triangles(self):
        all_triangles = []
        for data, hull in zip(self.data, self.hulls):
            triangles = []
            counter = 0
            for simplex in hull.simplices:
                A, B, C = data[simplex]
                if not self.is_collinear2D(A, B, C):
                    triangles.append([A, B, C])
                    counter += 1
            if counter > 1:
                triangles = np.delete(triangles,obj=triangles[0], axis=0)
            all_triangles.append(triangles)
        return np.array(all_triangles)

    # find relevant decomposition reaction
    # aka find triangle or triangle edge to which point belongs to
    # from the set of all triangles
    #
    # THIS FUNCTION NEEDS OUTPUT TO BE FORMATTED
    #
    # ISSUE: if decomposition line is shared by two triangles
    # function will print the same line twice
    def find_decomposition(self, point, tri_set):
        for tri in tri_set:
            tribool, triangle = self.is_point_inside_triangle(point, tri)
            linebool, line = self.is_point_on_tri_edge(point, tri)
            if tribool:
                print('Triangle:')
                print(triangle)
                print('Distance:')
                print(self.point_distance_to_plane(point, triangle[1:]))
                print('--------------')
            elif linebool:
                print('Line:')
                print(line)
                print('Distance:')
                print(self.point_distance_to_line(point, line[1:]))
                print('--------------')
                # this is because line can be shared by more than one triangle
                break

    # currently prints all metastable points at a given pressure
    # and relevant decomopostion reaction
    # (and loop is over pressures/files)
    #
    # THIS FUNCTION NEEDS OUTPUT TO BE FORMATTED
    #
    def find_all_decomposition(self):
        for points, press, tri_set in zip(self.metastable, self.pressures,
                                          self.triangles):
            print('-------')
            print(press)
            print('-------')
            for point in points:
                self.find_decomposition(point, tri_set)

    # function finds phase name given
    # numpy array (x,y) coordinate
    #
    # ? where are phase names comming from ?
    #
    # in fact this helper function is likely to be
    # usefull for many other situations so it should go to utils
    def xy_to_name(self, point):
        return 'Hmmmpff'



if __name__ == '__main__':
    hull = ConvexHullData(sys.argv[1:])
    hull_plotter = PlotConvexHull(hull, steps=12)
    hull_plotter.plot_all(hull)
    MS = FindMetastable(hull)
    point = np.array([1, .5, 1])
    A = np.array([0.33, 0.33, 0.33])
    B = np.array([1., 1., 1.])
    C = np.array([0.5, 1, 0])

    #point = MS.metastable[0][0]
    #point = MS.vertices[0][2]
    #A = MS.triangles[0][0][0]
    #B = MS.triangles[0][0][1]
    #C = MS.triangles[0][0][2]

    #tri = MS.triangles[0][0]
    #print(MS.is_point_on_tri_edge(point, tri))
    #print(MS.is_point_in_triangle(point, A, B, C))
    #print(MS.is_point_inside_triangle(point, tri))
    #print("----------------")
    #print(MS.is_point_between(point[:2], B[:2], C[:2]))
    #for point in MS.metastable[1]:
    #    print 'Point is: ', point, '\n'
    #    MS.find_decomposition(point)
    print('xxxxxxxx')
    print('MetaStable')
    print(MS.metastable)
    print(MS.find_all_decomposition())
