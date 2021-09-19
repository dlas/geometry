#!/usr/bin/env python3

import geometry
import sys
import os
import numpy

def tsolve(ct):
    u =  {"systemx": ["x","y", "z"]};

    v = numpy.matrix([[5],[5],[5],[1]])
    z = numpy.matrix([[0],[0],[0],[1]])
    r = ct.solve("systemx", u, v);
    print (r);

    ct.vecmove(r);

    print (ct.get_projection_matrix("systemx", "root") * z)

def main():
    ct = geometry.CoordinateTree()
    ct.load_from_file(sys.argv[1])

    print(ct.get_projection_matrix("root", "system1"))
    print(ct.get_projection_matrix("system1", "root"))
    print(ct.get_projection_matrix("system2", "root"))
    print(ct.get_projection_matrix("root", "system2"))
    print(ct.get_projection_matrix("system3", "system2"))
    print(ct.get_projection_matrix("system3", "system1"))


    tsolve(ct)

    return
    f = open("demo.pov", "w");
    f.write(ct.to_pov_ray(sys.argv[2]))
    f.close()

    os.system("povray demo.pov -H2000 -W2000")
    os.system("eog demo.png")

main()
