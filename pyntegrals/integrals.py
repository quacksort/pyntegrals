# Pyntegrals - Copyright (C) 2022, Lucia Giorgi (quacksort)
# https://github.com/nekobanana/pyntegrals
# Contact me at quacksort@gmail.com

import numbers
from typing import Callable, Union

import numpy as np
import tripy
from scipy import integrate
from shapely import geometry


def integrate_over_polygon(function: Union[Callable, float], polygon: list[tuple]) -> float:
    """
        Computes the integral of a two-variables function over a polygonal domain in :math:`R^2`.

        Returns the definite integral of ``function(x, y)`` = :math:`f(x,y)` over
        ``polygon`` = :math:`\\Sigma`, that is:

        :math:`\\int_{\\Sigma} f(x, y) \\, d\\Sigma`.

        Parameters
        ----------
        function : Union[Callable, float]
            A scalar or a Python function with two arguments: the first one is x and the
            second one is y in R^2 space.

        polygon: list[tuple]
            A polygon defined by a list of its vertices. Each vertex is
            described by a tuple containing its coordinated in the form (x, y).
            The vertices in the list must be in clockwise or
            counter-clockwise order (it makes no difference)

        Returns
        -------
        r : float
            The definite integral of ``function(x, y)`` over ``polygon``.

        Example
        --------
        Compute the integral of ``x**2 + y`` over a 50 x 40 rectangle :math:`\\Sigma` with a
        vertex in (0, 0) and another in (50, 40)

        :math:`\\int_{\\Sigma} (x^2 + y) \\, d\\Sigma`.

        >>> from pyntegrals.integrals import integrate_over_polygon
        >>> polygon = [(0, 0), (50, 0), (50, 40), (0, 40)]
        >>> rho = lambda x, y: x**2 + y
        >>> mass = integrate_over_polygon(function=rho, polygon=polygon)
            1706666.666666667

        For more example see examples/example.py

    """
    f = function
    if isinstance(function, numbers.Number):
        f = lambda x, y: function
    triangles = tripy.earclip(polygon)
    r: float = 0
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
