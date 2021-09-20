#!/usr/bin/env python3

import geometry
import sys
import os
import numpy

def tsolve2(ct):
    u =  {"system1":["x"], "system2":["y"], "system3":["z"]};

    v = numpy.matrix([[5],[5],[5],[1]])
    z = numpy.matrix([[0],[0],[0],[1]])
    r = ct.solve("system4", u, v);
    print (r);

    ct.moveall(r);

    m = ct.get_projection_matrix("root", "system4")
    print(m)
    print(m*v)


def tsolve(ct):
    u =  {"systemx2": ["x","y", "z"]};

    v = numpy.matrix([[5],[5],[5],[1]])
    z = numpy.matrix([[0],[0],[0],[1]])
    r = ct.solve("systemx2", u, v);
    print (r);

    ct.moveall(r);

    m = ct.get_projection_matrix("root", "systemx2")
    print(m)
    print(m*v)

def main():
    ct = geometry.CoordinateTree()
    ct.load_from_file(sys.argv[1])

    print(ct.get_projection_matrix("root", "system1"))
    print(ct.get_projection_matrix("system1", "root"))
    print(ct.get_projection_matrix("simple", "root"))
    print(ct.get_projection_matrix("root", "simple"))
    print(ct.get_projection_matrix("system3", "system2"))
    print(ct.get_projection_matrix("system3", "system1"))


    tsolve2(ct)
    return
    f = open("demo.pov", "w");
    f.write(ct.to_pov_ray(sys.argv[2]))
    f.close()

    os.system("povray demo.pov -H2000 -W2000")
    os.system("eog demo.png")

main()
