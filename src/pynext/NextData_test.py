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
from . NextData import next100_lead_shield
from . NextData import next100_PV
from . NextData import next100_copper_shield
from . NextData import next100_envelop


def test_lead_shield():
    npvd = NextPVData()
    n100_pb = next100_lead_shield()
    assert npvd.pb_inner_radius  == approx(npvd.pv_outer_radius, rel=1e-5)
    Rout = npvd.pb_inner_radius + npvd.pb_body_thickness
    assert npvd.pb_outer_radius  == approx(Rout, rel=1e-5)
    V = pi * (npvd.pb_outer_radius**2 - npvd.pb_inner_radius**2) * npvd.pb_length
    rho = n100_pb.cv.body.material.rho

    n100_pb.body_mass
    assert n100_pb.body_mass  == approx(rho * V, rel=1e-5)
