# Pyntegrals - Copyright (C) 2022, Lucia Giorgi (quacksort)
# https://github.com/nekobanana/pyntegrals
# Contact me at quacksort@gmail.com

import numbers
from typing import Callable, Union

import numpy as np
import tripy
from scipy import integrate
from shapely import geometry


def integrate_over_polygon(function: Union[Callable, float], polygon: list[tuple]):
    f = function
    if isinstance(function, numbers.Number):
        f = lambda x, y: function
    triangles = tripy.earclip(polygon)
    r = 0
    for triangle in triangles:
        int_x, int_y, A, t = __lintrasf(geometry.Polygon(triangle))
        tmp = integrate.dblquad(
            lambda y, x: __get_dblquad_function(f, A, t, x, y),
            0, int_x,
            lambda x: 0, lambda x: __line_from_axis_intersection(int_x, int_y, x))
        r += tmp[0]
    return r


def __lintrasf(triangle: geometry.Polygon):
    coords = np.array(geometry.polygon.orient(triangle, sign=1).exterior.coords)
    transl_vect = coords[0]
    w1 = np.array(coords[1] - transl_vect)
    vers_w1 = w1 / (w1 ** 2).sum() ** 0.5
    w2 = np.array(coords[2] - transl_vect)
    vers_w2 = w2 / (w2 ** 2).sum() ** 0.5
    A = np.array([vers_w1, vers_w2]).swapaxes(0, 1)
    new_coord1 = np.linalg.solve(A, coords[1] - transl_vect)
    new_coord2 = np.linalg.solve(A, coords[2] - transl_vect)
    int_x = new_coord1[0]
    int_y = new_coord2[1]
    return int_x, int_y, A, transl_vect


def __invtrasf(v, A, t):
    ret = np.matmul(A, v) + t
    return ret


def __line_from_axis_intersection(int_x, int_y, x):
    return -int_y / int_x * x + int_y


def __get_dblquad_function(function, A, t, x, y):
    f_point = __invtrasf(np.array([x, y]), A, t)
    return function(f_point[0], f_point[1]) * np.linalg.det(A)
