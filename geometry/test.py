#!/usr/bin/env python3

import geometry
import numpy
import sys
import pytest


def test():
    ct = geometry.CoordinateTree()
    ct.load_from_file("test_1.json")

    m = ct.get_projection_matrix("root", "system1")

    assert(m[0,0]==pytest.approx(0))
    assert(m[1,1]==pytest.approx(0))
    assert(m[1,0]==pytest.approx(1))
    assert(m[0,1]==pytest.approx(-1))

def test2():
    ct = geometry.CoordinateTree()
    ct.load_from_file("test_1.json")

    m = ct.get_projection_matrix("system1", "system3")
    n = ct.get_projection_matrix("system3", "system1")

    v = numpy.matrix([[1],[2],[3],[1]]);

    assert((m)*(n * v) == pytest.approx(v))

def test_solver1():
    ct = geometry.CoordinateTree();

    m = numpy.identity(4)

    v1 = numpy.matrix([[1],[2],[3],[1]]);
    v2 = numpy.matrix([[4],[5],[6],[1]]);

    r = ct.solve([m,m,m,m], {0:{"row":0, "var":0}, 1:{"row":2, "var":1}, 2:{"row":1, "var":2}}, v1, v2);

    print(r)

