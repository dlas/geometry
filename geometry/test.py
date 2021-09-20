#!/usr/bin/env python3

import geometry
import numpy
import sys
import pytest


def test():
    ct = geometry.CoordinateTree()
    ct.load_from_file("test_1.json")

    m = ct.get_projection_matrix("root", "system1")

    assert(m[0, 0] == pytest.approx(0))
    assert(m[1, 1] == pytest.approx(0))
    assert(m[1, 0] == pytest.approx(1))
    assert(m[0, 1] == pytest.approx(-1))


def test2():
    ct = geometry.CoordinateTree()
    ct.load_from_file("test_1.json")

    m = ct.get_projection_matrix("system1", "system3")
    n = ct.get_projection_matrix("system3", "system1")

    v = numpy.matrix([[1], [2], [3], [1]])

    assert((m)*(n * v) == pytest.approx(v))


def test_solver1():

    ct = geometry.CoordinateTree()
    ct.load_from_file("test_2.json")
    u = {"system1": ["x"], "system2": ["y"], "system3": ["z"]}

    v = numpy.matrix([[5], [5], [5], [1]])
    z = numpy.matrix([[0], [0], [0], [1]])
    r = ct.solve("system4", u, v)

    ct.moveall(r)

    m = ct.get_projection_matrix("root", "system4")

    assert(m*v == pytest.approx(z))
