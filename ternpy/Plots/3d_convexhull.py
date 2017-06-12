import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import time
import math
from scipy.spatial import ConvexHull
from matplotlib.tri import Triangulation


# returns array containg data for ConvexHull()
def load_data(file_list):
    data = []
    pressures = []
    for entry in file_list:
        try:
            data.append(np.genfromtxt(entry, delimiter="", skip_header=2))
            with open(entry, 'r') as f:
                press = f.readline().split()[3]
                f.close()
            pressures.append(press)
        except ValueError:
            print(os.path.basename(entry) +
                  " Oops! Looks like wrong file format...")
    return data, pressures



# returns array containing convex hull objects
def get_all_convex_hulls(data):
    hulls = []
    for entry in data:
        hulls.append(ConvexHull(entry))
    return hulls


# returns the value of the lowest point (z coordinate)
# from the convex hull list
def min_z(data, hulls):
    minz = 0.
    for idx, hull in enumerate(hulls):
        if minz > min(data[idx][hull.vertices, 2]):
            minz = min(data[idx][hull.vertices, 2])
    return minz

if __name__ == '__main__':
    print len(sys.argv)
    #convex_hull_3d(sys.argv[1])
    #print(len(get_all_convex_hulls(load_data(sys.argv[1:]))))
    data, pressures = load_data(sys.argv[1:])
    chulls = get_all_convex_hulls(data)
    #print(data[chulls.vertices, 2])
    print(min_z(data, chulls))
    print len(sys.argv)
    print pressures
