from pyntegrals.integrals import integrate_over_polygon


def main():

    polygon = [(0, 0), (300, 0), (300, 500), (0, 500)]
    area = integrate_over_polygon(function=1, polygon=polygon)
    print(f"Area = {area}")
    ry = integrate_over_polygon(function=lambda x, y: y, polygon=polygon) / area
    print(f"Momentum axis (gx) = {ry}")
    I = integrate_over_polygon(function=lambda x, y: (y - ry)**2, polygon=polygon)
    print(f"Inertia wrt gx = {I}")


if __name__ == '__main__':
    main()
