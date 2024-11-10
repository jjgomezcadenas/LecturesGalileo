from math import pi, exp, log
from . system_of_units import *
from . PhysicalVolume import PhysicalVolume
from . Material import RadioactiveMaterial
from . Shapes import Brick

from pytest import approx
from pytest import fixture
from operator import itemgetter, attrgetter
from collections import namedtuple


@fixture(scope='module')
def leadBrick():

    pb =   RadioactiveMaterial(name='Pb',
                              rho = 11.33 * g/cm3,
                              mu_over_rho = 0.044 * cm2/g,
                              a_bi214 = 370 * muBq/kg,
                              a_tl208 = 73 * muBq/kg )

    bs =   Brick             (width= 10 * cm, heigth= 10 * cm, length=10 * cm)

    return pb, bs, PhysicalVolume('leadBrick', pb, bs)


def test_brick(leadBrick):
    pb, bs, lb = leadBrick

    assert lb.shape.V    == bs.width * bs.heigth * bs.length
