#!/usr/bin/env python3

import geometry
import sys


def main():
    ct = geometry.CoordinateTree()
    ct.load_from_file(sys.argv[1])

    print(ct.get_projection_matrix("system3", "system1"))

    f = open("3.pov", "w");
    f.write(ct.to_pov_ray(sys.argv[2]))

main()
