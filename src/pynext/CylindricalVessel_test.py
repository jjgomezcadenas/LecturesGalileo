from math import pi, exp, log
from . system_of_units import *
from . PhysicalVolume import PhysicalVolume
from . Material import RadioactiveMaterial
from . Material import SelfAtt
from . Material import vacuum, ti316, cu12, cu03, pb

from pytest import approx
from pytest import fixture
from operator import itemgetter, attrgetter
from math import exp
from collections import namedtuple

from . CylindricalVessel import CylindricalVessel
from . CylindricalVessel import  CVD
from . activity_functions import CVA
from . NextData import NextPVData
from scipy.integrate import quad
from . math_functions import attenuation_factor


@fixture(scope='module')
def next100pv():
    npvd = NextPVData()
    cvd_pv = CVD(name    = 'Next100PV',
                 R       = npvd.pv_inner_radius,
                 th_body = npvd.pv_body_thickness,
                 L       = npvd.pv_length,
                 th_head = npvd.pv_head_thickness)

    pv = CylindricalVessel(name=cvd_pv.name, material=ti316, cvd=cvd_pv)
    return  pv


def test_att(next100pv):
    pv = next100pv
    att = SelfAtt(mu=pv.cv.body.material.mu, L=pv.body_thickness)
    tf, _ = quad(att.f, -pi/2 + 0.0001, pi/2 -  0.0001)
    attf = tf / (2 * pi)
    assert pv.body_self_shield_activity_bi214  == approx(pv.body_activity_bi214 * attf, rel=1e-5)


def test_pv_body(next100pv):
    pv = next100pv
    assert pv.cv.body.shape.inner_volume()  == approx(
    pi * pv.cv.body.radius**2 * pv.length, rel=1e-5)
    Rout = pv.radius + pv.body_thickness
    assert pv.cv.body.shape.shell_volume()  == approx(
    pi * (Rout**2 - pv.radius**2) * pv.length, rel=1e-5)
    assert pv.cv.body.shape.inner_surface()  == approx(
    2 * pi * pv.cv.body.radius * pv.length, rel=1e-5)
    assert pv.cv.body.shape.outer_surface()  == approx(
    2 * pi * Rout * pv.length, rel=1e-5)


def test_pv_head(next100pv):
    pv = next100pv
    assert pv.cv.head.shape.inner_volume()  == approx(
    pi * pv.radius**2 * pv.head_thickness, rel=1e-5)
    assert pv.cv.head.shape.shell_volume()  == 0
    assert pv.cv.head.shape.inner_surface()  == approx(
    pi * pv.radius**2, rel=1e-5)
    assert pv.cv.head.shape.outer_surface()  == approx(
    pi * pv.radius**2, rel=1e-5)


def test_mass_pv_body(next100pv):
    pv = next100pv
    mass   =  pv.cv.body.shape.shell_volume() * ti316.rho
    assert pv.body_mass  == approx(mass, rel=1e-5)


def test_mass_pv_head(next100pv):
    pv = next100pv
    mass   =  pv.cv.head.shape.inner_volume() * ti316.rho
    assert pv.head_mass  == approx(mass, rel=1e-5)


def test_activity_pv_body(next100pv):
    pv = next100pv
    activity_bi214     =  pv.body_mass * pv.cv.body.material.mass_activity_bi214
    activity_tl208     =  pv.body_mass * pv.cv.body.material.mass_activity_tl208

    assert pv.body_activity_bi214 == approx(activity_bi214, rel=1e-5)
    assert pv.body_activity_tl208 == approx(activity_tl208, rel=1e-5)


def test_activity_pv_head(next100pv):
    pv = next100pv
    activity_bi214     =  pv.head_mass * pv.cv.head.material.mass_activity_bi214
    activity_tl208     =  pv.head_mass * pv.cv.head.material.mass_activity_tl208

    assert pv.head_activity_bi214 == approx(activity_bi214, rel=1e-5)
    assert pv.head_activity_tl208 == approx(activity_tl208, rel=1e-5)


def test_activity_self_shield_pv_body(next100pv):
    pv = next100pv
    activity_bi214_shielded = pv.cv.body.activity_bi214_self_shield(pv.cvd.th_body)
    activity_tl208_shielded = pv.cv.body.activity_tl208_self_shield(pv.cvd.th_body)

    assert pv.body_self_shield_activity_bi214  == approx(activity_bi214_shielded, rel=1e-5)
    assert pv.body_self_shield_activity_tl208  == approx(activity_tl208_shielded, rel=1e-5)

    att = attenuation_factor(pv.cv.body.material.mu, pv.cvd.th_body)
    assert pv.body_self_shield_activity_bi214  == approx(pv.body_activity_bi214 * att, rel=1e-5)
    assert pv.body_self_shield_activity_tl208  == approx(pv.body_activity_tl208 * att, rel=1e-5)


def test_activity_self_shield_pv_head(next100pv):
    pv = next100pv
    activity_bi214_shielded = pv.cv.head.activity_bi214_self_shield(pv.cvd.th_head)
    activity_tl208_shielded = pv.cv.head.activity_tl208_self_shield(pv.cvd.th_head)

    assert pv.head_self_shield_activity_bi214  == approx(activity_bi214_shielded, rel=1e-5)
    assert pv.head_self_shield_activity_tl208  == approx(activity_tl208_shielded, rel=1e-5)

    att = attenuation_factor(pv.cv.head.material.mu, pv.cvd.th_head)
    assert pv.head_self_shield_activity_bi214  == approx(pv.head_activity_bi214 * att, rel=1e-5)
    assert pv.head_self_shield_activity_tl208  == approx(pv.head_activity_tl208 * att, rel=1e-5)


def test_pv_properties(next100pv):
    pv = next100pv
    assert pv.name                           == 'Next100PV'
    assert pv.material_name                  == ti316.name
    assert pv.radius / mm                     == approx(1360 / 2, rel=1e-3)
    assert pv.body_thickness / mm            == approx(10, rel=1e-3)
    assert pv.head_thickness / mm            == approx(12, rel=1e-3)
