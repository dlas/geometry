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

