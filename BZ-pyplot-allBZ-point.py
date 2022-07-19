#!/usr/bin/env python3

import sys
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

# For debuging with ipython (optionnal)
# from IPython import embed
# embed()


a = 5.43e-10
pi = np.pi

DESCRIPTION_SCRIPT="""
This script aims to plot with matplotlib the convex hull of the BZ. It is thus
just a 3d scatter plot of all W dots of the BZ. Then all lines of the convex
hull are drawed.

The approach here is pragmatic and dummy. We list all W dots and all lines, and
plot everything.
"""


def main_BZ_pyplot_allBZ_point():
    """
    See DESCRIPTION_SCRIPT variable. 
    """

    y_max_xp = np.array([1/4,  1/2,  0])
    y_max_xm = np.array([-1/4, 1/2,  0])
    y_max_zp = np.array([0,    1/2,  1/4])
    y_max_zm = np.array([0,    1/2,  -1/4])

    y_min_xp = np.array([1/4,  -1/2,  0])
    y_min_xm = np.array([-1/4, -1/2,  0])
    y_min_zp = np.array([0,    -1/2,  1/4])
    y_min_zm = np.array([0,    -1/2,  -1/4])


    x_max_yp = np.array([1/2, 1/4,   0])
    x_max_ym = np.array([1/2, -1/4,  0])
    x_max_zp = np.array([1/2, 0,     1/4])
    x_max_zm = np.array([1/2, 0,     -1/4])

    x_min_yp = np.array([-1/2, 1/4,   0])
    x_min_ym = np.array([-1/2, -1/4,  0])
    x_min_zp = np.array([-1/2, 0,     1/4])
    x_min_zm = np.array([-1/2, 0,     -1/4])


    z_max_yp = np.array([0,    1/4,   1/2])
    z_max_ym = np.array([0,    -1/4,  1/2])
    z_max_xp = np.array([1/4,  0,     1/2])
    z_max_xm = np.array([-1/4, 0,     1/2])

    z_min_yp = np.array([0,    1/4,   -1/2])
    z_min_ym = np.array([0,    -1/4,  -1/2])
    z_min_xp = np.array([1/4,  0,     -1/2])
    z_min_xm = np.array([-1/4, 0,     -1/2])

    all_lines = [
            [ y_max_zm, y_max_xm], # ymax square
            [ y_max_xm, y_max_zp],
            [ y_max_zp, y_max_xp],
            [ y_max_xp, y_max_zm],
            [ z_max_ym, z_max_xm], # zmax square
            [ z_max_xm, z_max_yp],
            [ z_max_yp, z_max_xp],
            [ z_max_xp, z_max_ym],
            [ x_max_ym, x_max_zm], # xmax square
            [ x_max_zm, x_max_yp],
            [ x_max_yp, x_max_zp],
            [ x_max_zp, x_max_ym],
            [ y_min_zm, y_min_xm], # ymin square
            [ y_min_xm, y_min_zp],
            [ y_min_zp, y_min_xp],
            [ y_min_xp, y_min_zm],
            [ z_min_ym, z_min_xm], # zmin square
            [ z_min_xm, z_min_yp],
            [ z_min_yp, z_min_xp],
            [ z_min_xp, z_min_ym],
            [ x_min_ym, x_min_zm], # xmin square
            [ x_min_zm, x_min_yp],
            [ x_min_yp, x_min_zp],
            [ x_min_zp, x_min_ym],
            [ x_max_zp, z_max_xp], # Here start the link between square, from xmax
            [ x_max_zm, z_min_xp],
            [ x_max_yp, y_max_xp],
            [ x_max_ym, y_min_xp],
            [ x_min_zp, z_max_xm], # From xmin
            [ x_min_zm, z_min_xm],
            [ x_min_yp, y_max_xm],
            [ x_min_ym, y_min_xm],
            [ z_max_yp, y_max_zp], # The others
            [ z_max_ym, y_min_zp],
            [ z_min_yp, y_max_zm],
            [ z_min_ym, y_min_zm]
    ]

    fig = plt.figure()
    ax = plt.axes(projection='3d')

    for lines in all_lines:
        plot_segment(lines[0], lines[1], ax)

    plt.show()


def plot_segment(dot_1, dot_2, ax):
    ax.plot3D([dot_1[0], dot_2[0]], [dot_1[1], dot_2[1]], [dot_1[2], dot_2[2]], c='b')

if __name__ == "__main__":
     main_BZ_pyplot_allBZ_point()

