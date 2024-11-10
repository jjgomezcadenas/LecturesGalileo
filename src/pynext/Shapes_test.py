from . TpcEL import TpcEL
from pynext.system_of_units import *
from . Shapes import Sphere
from . Shapes import SphereShell
from . Shapes import Cylinder
from . Shapes import CylinderShell
from . Shapes import Disk
from . Shapes import FlatPlate
from . Shapes import Brick

from pytest import approx
from pytest import fixture
from math import pi
from collections import namedtuple

ShapeParams = namedtuple('ShapeParams',
                         """inner_volume shell_volume
                         inner_surface outer_surface thickness_surface thickness""")
@fixture(scope='module')
def Next100PV():
    L    = 130 * cm
    R    =  52 * cm
    ts   =   1 * cm
    tc   =  10 * cm
    Rin  = R
    Rout = Rin + ts
    pv = CylinderShell(Rin, Rout, L)
    fp = FlatPlate(R, tc)
    return pv, fp

def assert_shape(sp, SP):
    assert sp.inner_volume()      == SP.inner_volume
    assert sp.shell_volume()      == SP.shell_volume
    assert sp.inner_surface()     == SP.inner_surface
    assert sp.outer_surface()     == SP.outer_surface
    assert sp.thickness_surface() == SP.thickness_surface
    assert sp.thickness()         == SP.thickness


def test_sphere():
    R=1
    sphere = ShapeParams(
    inner_volume      = (4/3) * pi * R**3,
    shell_volume      = 0,
    inner_surface     = 4 * pi * R**2,
    outer_surface     = 4 * pi * R**2,
    thickness_surface = 0,
    thickness         = 0
    )
    sp = Sphere(R = 1)
    assert_shape(sp, sphere)


def test_sphere_shell():
    Rin=1
    Rout = 2
    sphereShell = ShapeParams(
    inner_volume      = (4/3) * pi * Rin**3,
    shell_volume      = (4/3) * pi * (Rout**3 - Rin**3),
    inner_surface     = 4 * pi * Rin**2,
    outer_surface     = 4 * pi * Rout**2,
    thickness_surface = 0,
    thickness         = Rout - Rin
    )
    sp = SphereShell(Rin, Rout)
    assert_shape(sp, sphereShell)


def test_cylinder():
    R=1
    L = 1
    cyl = ShapeParams(
    inner_volume      = pi * R**2 * L,
    shell_volume      = 0,
    inner_surface     = 2 * pi * R * L,
    outer_surface     = 2 * pi * R * L,
    thickness_surface = 0,
    thickness         = 0
    )
    sp = Cylinder(R, L)
    assert_shape(sp, cyl)


def test_cylinder_shell():
    Rin  =1
    Rout = 2
    L    = 1
    scyl = ShapeParams(
    inner_volume      = pi * Rin**2 * L,
    shell_volume      = pi * (Rout**2 - Rin**2) * L,
    inner_surface     = 2 * pi * Rin * L,
    outer_surface     = 2 * pi * Rout * L,
    thickness_surface = pi * (Rout**2 - Rin**2),
    thickness         = Rout - Rin
    )
    sp = CylinderShell(Rin, Rout, L)
    assert_shape(sp, scyl)


def test_disk():
    R = 1
    t = 1
    disk = ShapeParams(
    inner_volume      = pi * R**2 * t,
    shell_volume      = 0,
    inner_surface     = pi * R**2,
    outer_surface     = pi * R**2,
    thickness_surface = 2 * pi *  R * t,
    thickness         = t
    )
    sp = Disk(R, t)
    assert_shape(sp, disk)


def test_brick():
    width  = 1
    heigth = 1
    length =1
    brick = ShapeParams(
    inner_volume      = width * heigth * length,
    shell_volume      = 0,
    inner_surface     = 2 * (width * heigth + width * length + length * heigth),
    outer_surface     = 2 * (width * heigth + width * length + length * heigth),
    thickness_surface = 0,
    thickness         = 0
    )
    sp = Brick(width, heigth, length)
    assert_shape(sp, brick)


def test_cylinder_next(Next100PV):
    pv, _ = Next100PV


    assert pv.Rin             / cm  == approx(52, rel=1e-3)
    assert pv.Rout            / cm  == approx(53, rel=1e-3)
    assert pv.L               / cm  == approx(130, rel=1e-3)
    assert pv.inner_volume()  / m3  == approx(1.10, rel=1e-2)
    assert pv.shell_volume()  / m3  == approx(4.29e-2, rel=1e-3)
    assert pv.inner_surface() / m2  == approx(4.25, rel=1e-3)
    assert pv.outer_surface() / m2  == approx(4.33, rel=1e-3)


def test_flat_plate(Next100PV):
    _, fp = Next100PV

    assert fp.R / m  == approx(0.52, rel=1e-3)
    assert fp.t / cm == approx(10, rel=1e-3)
    assert fp.S / mm2 == approx(8.49e+05, rel=1e-2)
    assert fp.V / mm3 == approx(8.49e+07, rel=1e-2)
