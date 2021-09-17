#!/usr/bin/env python3

import geometry
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

