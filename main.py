import numpy as np
import tripy
from scipy import integrate
from shapely.geometry import Polygon
from shapely.geometry.polygon import orient


def lintrasf(triangle: Polygon):
    coords = np.array(orient(triangle, sign=-1).exterior.coords)
    # w1 = coords[0]
    transl_vect = coords[0]
    w1 = np.array(coords[1] - transl_vect)
    vers_w1 = w1 / (w1 ** 2).sum() ** 0.5
    # w2 = coords[1]
    w2 = np.array(coords[2] - transl_vect)
    vers_w2 = w2 / (w2 ** 2).sum() ** 0.5
    # w3 = np.array([coords[1, 0] - coords[2, 0], coords[1, 1] - coords[2, 1]])
    A = np.array([vers_w1, vers_w2]).swapaxes(0, 1)
    new_coord1 = np.linalg.solve(A, coords[0] - transl_vect)
    new_coord2 = np.linalg.solve(A, coords[1] - transl_vect)
    new_coord3 = np.linalg.solve(A, coords[2] - transl_vect)
    C = np.array([new_coord1, new_coord2, new_coord3])
    return C, A, transl_vect


def invtrasf(v, A, t):
    ret = np.matmul(A, v) + t
    return ret


class Parabole_concrete():
    def __init__(self, a, b):
        self.a = a
        self.b = b

    def __call__(self, y):
        return self.a * y ** 2 + self.b * y


def line_from_axis_intersection(int_x, int_y, x):
    return -int_y / int_x * x + int_y


def f(x):
    a = np.ones(x[0].shape)
    return a


def f1(x):
    a = -5000000 * (x[1] * 0.0001) ** 2 + 20000 * (x[1] * 0.0001)
    return a


def main():
    concreteSSL = ConcreteSSL(20, 0.002, 0.0048)
    highest_y = 500
    h = 48
    h_y = highest_y - h
    parabole = Parabole_concrete(-5000000, 20000)
    polygon = [(0, 452), (300, 452), (300, 472), (0, 472)]
    triangles = tripy.earclip(polygon)
    sub1_parabole = lambda y: parabole((y - h_y) * concreteSSL.eps_cu / h)

    r = 0
    for triangle in triangles:
        C, A, transl_vect = lintrasf(Polygon(triangle))
        tmp = integrate.dblquad(lambda y, x: sub1_parabole(invtrasf(np.array([x, y]), A, transl_vect)[1]) * np.linalg.det(A),
        # tmp = integrate.dblquad(lambda y, x: np.linalg.det(A),
                                0, C[1, 0],
                                lambda x: 0, lambda x: line_from_axis_intersection(C[1, 0], C[2, 1], x))
        r += tmp[0]
    # tri_list = np.array(triangles).swapaxes(0, 1)
    # get a "good" scheme of degree 10
    # scheme = quadpy.t2.get_good_scheme(10)
    # val = scheme.integrate(lambda x: parabole((x - h_y) * concreteSSL.eps_cu / h, 1), tri_list)
    r = -r
    print(r)
    print("OK")


if __name__ == '__main__':
    main()
