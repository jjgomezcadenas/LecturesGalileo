from . system_of_units import *
from . CylindricalVessel import  CVD
from . activity_functions import CVA
from . CylindricalVessel import CylindricalVessel
from pynext.Material import vacuum, ti316, cu12, cu03, pb

class RFlux:
    def __init__(self,
                 name  = 'LSC2010',
                 U238  = 0.55 * Bq / cm2,
                 Th232 = 0.36 * Bq / cm2
                 ):
        self.name  = name
        self.U238  = U238
        self.Th232 = Th232
        self.Bi214 = self.U238 * (0.7 / 100)
        self.Tl208 = self.Th232 * (13.65 / 100)
        self.energy_gamma_bi214 = 2447.9 * keV
        self.energy_gamma_tl208 = 2614.5 * keV

    def __str__(self):

        s = """

        {:s}
        ------------------

        U238  flux  = {:7.2f} mBq / cm2
        Th232 flux  = {:7.2f} mBq / cm2
        Bi214 flux  = {:7.2f} mBq / cm2
        Tl208 flux  = {:7.2f} mBq / cm2

    """.format(self.name,
               self.U238 / (mBq / cm2),
               self.Th232 / (mBq / cm2),
               self.Bi214 / (mBq / cm2),
               self.Tl208 / (mBq / cm2))
        return s

    __repr__ = __str__


class NextPVData:
    def __init__(self,
                 name='Next100PVData',
                 pv_inner_diameter = 1360 * mm,
                 pv_length         = 1600 * mm,
                 pv_body_thickness =   10 * mm,
                 pv_head_thickness =   12 * mm,
                 cs_body_thickness =  120 * mm,
                 cs_head_thickness =  120 * mm,
                 pb_body_thickness =  200 * mm,
                 pb_head_thickness =  200 * mm,
                 ):
        self.name              = name
        self.pv_inner_diameter = pv_inner_diameter
        self.pv_length         = pv_length
        self.pv_body_thickness = pv_body_thickness
        self.pv_head_thickness = pv_head_thickness
        self.pv_inner_radius   = self.pv_inner_diameter / 2.
        self.pv_outer_diameter = self.pv_inner_diameter + 2 * self.pv_body_thickness
        self.pv_outer_radius   = self.pv_outer_diameter / 2.

        self.cs_head_thickness = cs_head_thickness
        self.cs_body_thickness = cs_body_thickness
        self.cs_inner_diameter = self.pv_inner_diameter - 2 * self.cs_body_thickness
        self.cs_length         = self.pv_length
        self.cs_inner_radius   = self.cs_inner_diameter / 2.
        self.cs_outer_diameter = self.pv_inner_diameter
        self.cs_outer_radius   = self.cs_outer_diameter / 2.

        self.pb_head_thickness = pb_head_thickness
        self.pb_body_thickness = pb_body_thickness
        self.pb_inner_diameter = self.pv_outer_diameter
        self.pb_length         = self.pv_length
        self.pb_inner_radius   = self.pb_inner_diameter / 2.
        self.pb_outer_diameter = self.pb_inner_diameter + 2 * pb_body_thickness
        self.pb_outer_radius   = self.pb_outer_diameter / 2.

    def __str__(self):

        s = """

        {:s}
        ------------------

        PV :
        inner diameter  = {:7.2f} mm
        inner radius    = {:7.2f} mm
        outer diameter  = {:7.2f} mm
        outer radius    = {:7.2f} mm
        body thickness  = {:7.2f} mm
        head thickness  = {:7.2f} mm
        length          = {:7.2f} mm

        CS :
        inner diameter  = {:7.2f} mm
        inner radius    = {:7.2f} mm
        outer diameter  = {:7.2f} mm
        outer radius    = {:7.2f} mm
        body thickness  = {:7.2f} mm
        head thickness  = {:7.2f} mm
        length          = {:7.2f} mm

        PB :
        inner diameter  = {:7.2f} mm
        inner radius    = {:7.2f} mm
        outer diameter  = {:7.2f} mm
        outer radius    = {:7.2f} mm
        body thickness  = {:7.2f} mm
        head thickness  = {:7.2f} mm
        length          = {:7.2f} mm

    """.format(self.name,
               self.pv_inner_diameter / mm,
               self.pv_inner_radius / mm,
               self.pv_outer_diameter / mm,
               self.pv_outer_radius / mm,
               self.pv_body_thickness / mm,
               self.pv_head_thickness / mm,
               self.pv_length,
               self.cs_inner_diameter / mm,
               self.cs_inner_radius / mm,
               self.cs_outer_diameter / mm,
               self.cs_outer_radius / mm,
               self.cs_body_thickness / mm,
               self.cs_head_thickness / mm,
               self.cs_length,
               self.pb_inner_diameter / mm,
               self.pb_inner_radius / mm,
               self.pb_outer_diameter / mm,
               self.pb_outer_radius / mm,
               self.pb_body_thickness / mm,
               self.pb_head_thickness / mm,
               self.pb_length)
        return s

    __repr__ = __str__


def next100_lead_shield():
    n100d = NextPVData()
    cvd_pb = CVD(name    ='PBShield',
                 R       =  n100d.pb_inner_radius,
                 th_body = n100d.pb_body_thickness,
                 L       = n100d.pb_length,
                 th_head = n100d.pb_head_thickness)
    n100_pb = CylindricalVessel(name='Next100Pb', material=pb, cvd=cvd_pb)
    return n100_pb


def next100_PV():
    n100d = NextPVData()
    cvd_pv = CVD(name    = 'Next100PV',
                 R       = n100d.pv_inner_radius,
                 th_body = n100d.pv_body_thickness,
                 L       = n100d.pv_length,
                 th_head = n100d.pv_head_thickness)
    # Pressure Vessel
    n100_pv = CylindricalVessel(name=cvd_pv.name, material=ti316, cvd=cvd_pv)
    return n100_pv


def next100_copper_shield():
    n100d = NextPVData()
    cvd_cu = CVD(name    ='CUShield',
                 R       = n100d.cs_inner_radius,
                 th_body = n100d.cs_body_thickness,
                 L       = n100d.cs_length,
                 th_head = n100d.cs_head_thickness)
    n100_cu = CylindricalVessel(name='Next100CU', material=cu12, cvd=cvd_cu)
    return n100_cu


def next100_envelop():
    n100d = NextPVData()
    cvd_pv = CVD(name    = 'Next100PV',
                 R       = n100d.pv_inner_radius,
                 th_body = n100d.pv_body_thickness,
                 L       = n100d.pv_length,
                 th_head = n100d.pv_head_thickness)
    # Pressure Vessel
    n100_envelop = CylindricalVessel(name='Next100Envelop', material=vacuum, cvd=cvd_pv)
    return n100_envelop
