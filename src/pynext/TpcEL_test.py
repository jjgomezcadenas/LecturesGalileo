from . TpcEL import TpcEL
from pynext.system_of_units import *
from pytest import approx

def test_tpcel():
    tpc = TpcEL()
    kv_cm_b = kilovolt / (cm * bar)
    kv_cm   = kilovolt / cm
    cm_bar = 1 / (cm * bar)


    tpc.EoverP / kv_cm_b              == approx(3.5, rel=1e-5)
    tpc.DriftVoltage / kv_cm          == approx(0.5, rel=1e-5)
    tpc.Pressure / bar                == approx(15, rel=1e-5)
    tpc.GridDistance / mm             == approx(5, rel=1e-5)
    tpc.DriftLength / m               == approx(1.20, rel=1e-5)
    tpc.GridVoltage / kilovolt        == approx(26.25, rel=1e-5)
    tpc.CathodeVoltage / kilovolt     == approx(86.25, rel=1e-5)
    tpc.OpticalGain                   == approx(2800, rel=1e-5)
    tpc.ScintillationPhotons(1 * MeV) == approx(41667, rel=1e-5)
    tpc.IonizationElectrons(1 * MeV)  == approx(62500, rel=1e-5)
    tpc.ELPhotons(1 * MeV)            == approx(1.7531e+08, rel=1e-5)
