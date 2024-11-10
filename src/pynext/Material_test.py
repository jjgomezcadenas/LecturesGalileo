from . Material import PhysicalMaterial
from . Material import PVMaterial
from pynext.system_of_units import *

from pytest import approx
from pytest import fixture
from operator import itemgetter, attrgetter
from math import exp

@fixture(scope='module')
def latt():

    rho_2020 = 124.3 * kg/m3
    rho_1020 = 58 * kg/m3
    rho_0520 = 30 * kg/m3
    rho_list = [
    ('rho_2020' , rho_2020),
    ('rho_1020' , rho_1020),
    ('rho_0520' , rho_0520),
    ]
    ordered = sorted([ (name, r) for (name, r) in rho_list ],
                     key=itemgetter(1))
    LATT = {}
    for name, rho in ordered:
        xe = PhysicalMaterial(name='GXe', rho=rho, mu_over_rho=0.039 * cm2/g)
        LATT[name] = xe.Latt /m
    return LATT


@fixture(scope='module')
def materials():

    ti316 = PVMaterial(name='316Ti',
                                rho = 7.87 * g/cm3,
                                mu_over_rho = 0.039 * cm2/g,
                                a_bi214 = 1.0 * mBq/kg,
                                a_tl208 = 0.4 * mBq/kg,
                                Sm      = 139 * MPa )

    cu =   PVMaterial(name='Cu',
                                rho = 8.96 * g/cm3,
                                mu_over_rho = 0.039 * cm2/g,
                                a_bi214 = 3 * muBq/kg,
                                a_tl208 = 3 * muBq/kg,
                                Sm      = 10 * MPa)
    return ti316, cu

def test_latt(latt):

    assert latt['rho_0520']   == approx(8.55, rel=1e-2)
    assert latt['rho_1020']   == approx(4.42, rel=1e-2)
    assert latt['rho_2020']   == approx(2.06, rel=1e-2)

def test_material(materials):
    ti316, _ = materials
    g_cm3 = g / cm3
    cm2_g = cm2 / g
    icm   = 1 / cm
    mbq_kg = mBq / kg


    assert ti316.name                                  == '316Ti'
    assert ti316.density / g_cm3                       == approx(7.87, rel=1e-2)
    assert ti316.mass_attenuation_coefficient / cm2_g  == approx(0.039, rel=1e-2)
    assert ti316.attenuation_coefficient / icm         == approx(0.31, rel=1e-2)
    assert ti316.attenuation_length / cm               == approx(3.26, rel=1e-2)
    assert ti316.maximum_allowable_strength / MPa      == approx(139, rel=1e-2)
    assert ti316.mass_activity_bi214 / mbq_kg          == approx(1, rel=1e-3)
    assert ti316.mass_activity_tl208 / mbq_kg          == approx(0.4, rel=1e-3)
