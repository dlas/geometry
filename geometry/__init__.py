#!/usr/bin/env python3


import numpy
import sys
import json
import math


def deg2rad(x: float):
    """ Convert from degreese to radians"""
    return math.pi*x/180


class CoordinateSystem:
    """Define a single coordinate system, possibly related to some root. 
    A coordinate system is defined by an angular and translational offset
    from its parent.
    """

    offset_x: float
    offset_y: float
    offset_z: float
    angle_x:  float
    angle_y:  float
    angle_z:  float
    parent_name: str
    my_name: str

    def get_active_matrix(self):
        """ Compute a numpy matrix representig the active translation of coordiantes
        in this coordinate system to coordinates in the parent coordinate system
        """
        s1 = math.sin(deg2rad(self.angle_x))
        s2 = math.sin(deg2rad(self.angle_y))
        s3 = math.sin(deg2rad(self.angle_z))
        c1 = math.cos(deg2rad(self.angle_x))
        c2 = math.cos(deg2rad(self.angle_y))
        c3 = math.cos(deg2rad(self.angle_z))

        v1 = [c2*c3, -s2, c2*s3, self.offset_x]
        v2 = [s1*s3 + c1*c3*s2, c1*c2, c1*s2*s3-c3*s1, self.offset_y]
        v3 = [c3*s1*s2 - c1*s3, c2*s1, c1*c3 + s1*s2*s3, self.offset_z]
        v4 = [0, 0, 0, 1]

        return numpy.matrix([v1, v2, v3, v4])

    def get_passive_matrix(self):
        """ Compute a numpy matrix representing the passive translation of
        coordinates in the parent coordinate system into this coordinate 
        system.
        """

        # Given the symmetries of the rotate/tanslate matrix, there is probably an easier
        # way to do this... But i'm lazy!
        return numpy.linalg.inv(self.get_active_matrix())

    def __init__(self, name, parent, vec):
        """
          name    - the name of this coordinate system
          parent  - the name of the parent coordinate system
          vec     - the parameters defining this coordinate system
                    angles are in degreese and offsets are in millimeters. 
                    Translational offsets are applied before rotation and the 
                    order of rotation is X, Y, Z.
        """
        self.offset_x, \
            self.offset_y, \
            self.offset_z, \
            self.angle_x,  \
            self.angle_y,  \
            self.angle_z = vec

        self.name = name
        self.parent_name = parent


class CoordinateTreeNode:
    cs: CoordinateSystem
    parent: 'CoordinateTreeNode'
    children: 'CoordinateTreeNode'
    level: int

    def __init__(self, cs: CoordinateSystem, parent: 'CoordinateTreeNode', level: int):
        self.cs = cs
        self.parent = parent
        self.level = level
        self.children = []


class CoordinateTree:
    """Describe a tree of coordinates systems wherein we can project a point
    from one coordinate system to another
    """

    root: CoordinateSystem      # Special Root coordinate system
    rootT: CoordinateTreeNode   # Tree Node fro the root coordinate system
    NameToCS = {}               # Mapping of coordinate system names

    def get_projection_matrix(self, frm: str, to: str):
        """Compute a matrix that translates points in the coordinate
        system of "frm" to the coordinate system "to".
          to   -  the name of the new coordinate system.
          from -  the name of the old coordinate system.
        this will yield a 4x4 matrix where the UL 3x3 matrix
        represents the rotation, the right row the translation. 
        If you multiply it by a vector [[X],[Y],[Z],[1]] you will
        transform it from the 'to' system to the 'frm' system
        """

        f = self.NameToCS[frm]
        t = self.NameToCS[to]

        m = numpy.identity(4)

        # Walk up the tree, multiplying by the appropriate matrix
        # for each coordinate system we pass through

        while (not (t is f)):
            if t.level >= f.level:
                m = m * t.cs.get_active_matrix()
                t = t.parent
            elif f.level >= t.level:
                #m = m * f.cs.get_passive_matrix();
                m = f.cs.get_passive_matrix() * m
                f = f.parent
        return m

    def load_from_file(self, fn: str):
        """Load coordinate systems from a file
          - fn a JSON file containing coordinate system information
        """

        # BUG -- failed open
        fp = open(fn)

        # Bug bad JSON
        j = json.load(fp)

        for entry in j:
            # BUG missing or incorrect fields
            newcs = CoordinateSystem(
                entry["name"], entry["parent"], entry["relative"])

            self.NameToCS[newcs.name] = CoordinateTreeNode(newcs, None, 0)

        # Build a root node.
        # BUG -- we assume everything is reachable from root
        self.root = CoordinateSystem("root", None, [0, 0, 0, 0, 0, 0])
        self.rootT = CoordinateTreeNode(self.root, None, 0)
        self.NameToCS["root"] = self.rootT

        # fill out the parent and child relationships for each node
        for csts in self.NameToCS:
            cst = self.NameToCS[csts]
            if (cst.cs.parent_name != None):
                cst.parent = self.NameToCS[cst.cs.parent_name]
                cst.parent.children.append(cst)

        # Fill out the level for each node
        self.resolve_levels(self.rootT, 0)

    def resolve_levels(self, cur: CoordinateTreeNode, level: int):
        cur.level = level
        for n in cur.children:
            self.resolve_levels(n, level+1)

    def to_pov_ray(self, r="root"):
        """ Write out a snippet of pov-ray code that shows every coordinate system"""

        output_str = ""
        xp_vec = numpy.matrix([[1], [0], [0], [1]])
        yp_vec = numpy.matrix([[0], [1], [0], [1]])
        zp_vec = numpy.matrix([[0], [0], [1], [1]])

        xn_vec = numpy.matrix([[-1], [0], [0], [1]])
        yn_vec = numpy.matrix([[0], [-1], [0], [1]])
        zn_vec = numpy.matrix([[0], [0], [-1], [1]])

        for csts in self.NameToCS:
            tx = self.get_projection_matrix(r, csts)
            xp_proj = tx*xp_vec
            yp_proj = tx*yp_vec
            zp_proj = tx*zp_vec
            xn_proj = tx*xn_vec
            yn_proj = tx*yn_vec
            zn_proj = tx*zn_vec

            output_str += self.povray_cylinder(xp_proj, xn_proj, "Red")
            output_str += self.povray_cylinder(yp_proj, yn_proj, "Green")
            output_str += self.povray_cylinder(zp_proj, zn_proj, "Blue")
            output_str +=\
                'text {{ internal 1 "{0}" 0.1, 0  scale .2 pigment {{Yellow}} translate <{1}, {2}, {3}>}}\n'\
                .format(csts, xp_proj[0, 0], xp_proj[1, 0]+.1, xp_proj[2, 0])
        return self.pov_ray_preamble() + output_str

    def pov_ray_preamble(self):
        p = """
        #include "colors.inc"
        background {color Grey}
        camera {
            location <5, 5, -20>
            look_at  <0, 1, 2>
        }
        light_source { <-1, -1, -30> color White }
        """
        return p

    def povray_cylinder(self, v1, v2, color):
        v_str =\
            'cylinder {{ <{0}, {1}, {2}>, <{3}, {4}, {5}>, 0.1 open texture {{ pigment {{ color {6}}}}}}}\n'\
            .format(v1[0, 0], v1[1, 0], v1[2, 0], v2[0, 0], v2[1, 0], v2[2, 0], color)

        return v_str
