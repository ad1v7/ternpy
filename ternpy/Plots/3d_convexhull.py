import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import time
import math
from scipy.spatial import ConvexHull
from matplotlib.tri import Triangulation

# generate 3D convex hull for a given file
def convex_hull_3d(input_file):
    try:
        X = np.genfromtxt(input_file, delimiter="", skip_header=2)
        return ConvexHull(X)

    except ValueError:
        print(os.path.basename(input_file) +
              " Oops! Looks like wrong file format...")


# returns array containing 3D convex hull objects
def get_all_convex_hulls(files):
    convex_hulls = []
    for file in files:
        convex_hulls.append(convex_hull_3d(file))

convex_hull_3d(sys.argv[1])
