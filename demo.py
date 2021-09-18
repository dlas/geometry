#!/usr/bin/env python3

import geometry
import sys
import os


def main():
    ct = geometry.CoordinateTree()
    ct.load_from_file(sys.argv[1])

    print(ct.get_projection_matrix("root", "system1"))
    print(ct.get_projection_matrix("system1", "root"))
    print(ct.get_projection_matrix("system2", "root"))
    print(ct.get_projection_matrix("root", "system2"))
    print(ct.get_projection_matrix("system3", "system2"))
    print(ct.get_projection_matrix("system3", "system1"))

    f = open("demo.pov", "w");
    f.write(ct.to_pov_ray(sys.argv[2]))
    f.close()

    os.system("povray demo.pov -H2000 -W2000")
    os.system("eog demo.png")

main()
