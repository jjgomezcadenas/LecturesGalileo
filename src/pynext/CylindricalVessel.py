"""
Physical Volume
Defines a material medium + Shape
"""

from math import pi, exp, log
from . system_of_units import *
from . PhysicalVolume import PhysicalVolume
from . Shapes import CylinderShell
from . Shapes import SphereShell
from . Shapes import Disk
from . Material import PVMaterial
from . Material import RadioactiveMaterial
from . activity_functions import Activity
from . activity_functions import CVA
from pynext.activity_functions import str_activity
from collections import namedtuple

# Cylindrical Vessel Dimensions (CVD)
CVD = namedtuple('CVD', 'name R th_body L th_head')



class CylindricalDetector:
    """Base class for detectors with a CylinderShell shape"""

    def __init__(self, name, inner_diameter, length, thickness, material):

        self.name              = name
        self.inner_diameter    = inner_diameter
        self.inner_radius      = inner_diameter / 2
        self.outer_diameter    = inner_diameter + 2 * thickness
        self.outer_radius      = self.outer_diameter / 2
        self.length            = length
        self.thickness         = thickness
        self.material          = material
        cs                     = CylinderShell(Rin=self.inner_radius,
                                               Rout=self.outer_radius,
                                               L=self.length)
        self.detector         = PhysicalVolume(name, material, cs)

    def __str__(self):

        s = """

        {:s}
        ------------------
        material         = {:s}
        inner diameter   = {:7.2f} mm
        inner radius     = {:7.2f} mm
        outer diameter   = {:7.2f} mm
        outer radius     = {:7.2f} mm
        thickness        = {:7.2f} mm
        length           = {:7.2f} mm
        mass             = {:7.2f} kg
        activity Bi214   = {:7.2e} mBq
        activity Tl208   = {:7.2e} mBq


    """.format(self.name,
               self.material.name,
               self.inner_diameter / mm,
               self.inner_radius / mm,
               self.outer_diameter / mm,
               self.outer_radius / mm,
               self.thickness / mm,
               self.length,
               self.detector.mass / kg,
               self.detector.activity_bi214 / mBq,
               self.detector.activity_tl208 / mBq)
        return s

    __repr__ = __str__


class NextFieldCage(CylindricalDetector):
    def __init__(self,
                 name='Next100FieldCage',
                 inner_diameter  = 1050 * mm,
                 length          = 1300 * mm,
                 thickness       =   25 * mm,
                 electrode_pitch =   12 * mm,
                 material        =      None,
                 electrode       =      None,
                 resitstorActivityFC =  None
                 ):

        super().__init__(name, inner_diameter, length, thickness, material)
        self.electrode         = electrode
        self.electrode_pitch   = electrode_pitch
        self.nof_electrodes    = self.length / self.electrode_pitch
        self.nof_resistors     = 2 * self.nof_electrodes
        self.resitstorActivity = resitstorActivityFC

        self.activityElectrodes = Activity(name = 'ActivityElectrodesFC',
                    bi214 = self.electrode.detector.activity_bi214 * self.nof_electrodes / 2,
                    tl208 = self.electrode.detector.activity_tl208 * self.nof_electrodes / 2)

        self.activityResistors = Activity(name = 'ActivityResistorsFC',
                    bi214 = self.resitstorActivity.bi214 * self.nof_resistors / 2,
                    tl208 = self.resitstorActivity.tl208 * self.nof_resistors / 2)

        self.activityPoly = Activity(name = 'ActivityPoly',
                    bi214 = self.detector.activity_bi214 / 2,
                    tl208 = self.detector.activity_tl208 / 2)
    @property
    def activity_electrodes(self):
        return self.activityElectrodes

    @property
    def activity_resistors(self):
        return self.activityResistors

    @property
    def activity_poly(self):
        return self.activityPoly

    def __str__(self):
        s2 = super().__str__()
        s = """
        ------------------
        electrode pitch  = {:7.2f} mm
        nof_electrodes   = {:7.2f}
        nof_resistors    = {:7.2f}

        ---
        Activity
        {:s}
        {:s}
        {:s}

    """.format(
               self.electrode_pitch,
               self.nof_electrodes,
               self.nof_resistors,
               str_activity('ActivityElectrodesFC', self.activityElectrodes, unit='mBq'),
               str_activity('ActivityResistorsFC',  self.activityResistors,  unit='mBq'),
               str_activity('ActivityPoly',         self.activityPoly,       unit='mBq')
               )
        return s2 + s

    __repr__ = __str__



class CylindricalVessel:
    def __init__(self, name, material, cvd):
        """material fill the Vessel
           body is a cylindrical shell
           head is a disk
           cvd is namedtuple (CVD = Cylindrical Vessel Dimensions)
           """

        CV = namedtuple('CV', 'name material body head')

        self.cvd = cvd
        Rout = cvd.R + cvd.th_body
        cs =   CylinderShell(Rin=cvd.R, Rout=Rout, L=cvd.L)
        ch =   Disk         (R=cvd.R, t=cvd.th_head)

        self.cv = CV(name = name,
                     material = material,
                     body     = PhysicalVolume(name, material, cs),
                     head     = PhysicalVolume(name, material, ch))

    @property
    def name(self):
        return self.cvd.name

    @property
    def material_name(self):
        return self.cv.material.name

    @property
    def radius(self):
        return self.cvd.R

    @property
    def length(self):
        return self.cvd.L

    @property
    def body_thickness(self):
        return self.cvd.th_body

    @property
    def head_thickness(self):
        return self.cvd.th_head

    @property
    def body_surface(self):
        return self.cv.body.S

    @property
    def body_volume(self):
        return self.cv.body.volume

    @property
    def body_mass(self):
        return self.cv.body.mass

    @property
    def head_surface(self):
        return self.cv.head.S

    @property
    def head_volume(self):
        return self.cv.head.volume

    @property
    def head_mass(self):
        return self.cv.head.mass

    @property
    def body_activity_bi214(self):
        return self.cv.body.activity_bi214

    @property
    def body_activity_tl208(self):
        return self.cv.body.activity_tl208

    @property
    def head_activity_bi214(self):
        return self.cv.head.activity_bi214

    @property
    def head_activity_tl208(self):
        return self.cv.head.activity_tl208

    @property
    def body_self_shield_activity_bi214(self):
        return self.cv.body.activity_bi214_self_shield(self.cvd.th_body)

    @property
    def body_self_shield_activity_tl208(self):
        return self.cv.body.activity_tl208_self_shield(self.cvd.th_body)

    @property
    def head_self_shield_activity_bi214(self):
        return self.cv.head.activity_bi214_self_shield(self.cvd.th_head)

    @property
    def head_self_shield_activity_tl208(self):
        return self.cv.head.activity_tl208_self_shield(self.cvd.th_head)

    @property
    def body_transmittance(self):
        return self.cv.body.transmittance(self.cvd.th_body)

    @property
    def head_transmittance(self):
        return self.cv.head.transmittance(self.cvd.th_head)

    @property
    def body_absorption(self):
        return self.cv.absorption_at_qbb(self.cvd.th_body)

    @property
    def head_absorption(self):
        return self.cv.absorption_at_qbb(self.cvd.th_head)

    def __str__(self):

        s = """
        Cylindrical Vessel:

        ----------------
        name      = {:s}
        material  = {:s}

        specific activity of material:
        Bi-214    = {:7.2e} mBq/kg
        Tl-208    = {:7.2e} mBq/kg

        body:
        R              = {:7.2f} mm
        body thickness = {:7.2f} mm
        head thickness = {:7.2f} mm
        length         = {:7.2f} mm
        surface        = {:7.2e} mm2
        volume         = {:7.2e} mm3
        mass           = {:7.2f} kg
        activity Bi-214 = {:7.2f} mBq, self-shielded ={:7.2f} mBq
        activity Tl-208 = {:7.2f} mBq, self-shielded ={:7.2f} mBq
        transmittance   = {:7.2e}

        heads:
        thickness = {:7.2f} mm
        surface   = {:7.2e} mm2
        volume    = {:7.2e} mm3
        mass      = {:7.2f} kg
        activity Bi-214 = {:7.2f} mBq, self-shielded ={:7.2f} mBq
        activity Tl-208 = {:7.2f} mBq, self-shielded ={:7.2f} mBq
        transmittance   = {:7.2e}

        """.format(self.name, self.material_name,
                   self.cv.material.mass_activity_bi214 / (mBq/kg),
                   self.cv.material.mass_activity_tl208 / (mBq/kg),
                   self.radius / mm,
                   self.body_thickness / mm,
                   self.head_thickness / mm,
                   self.length / mm,
                   self.body_surface / mm2,
                   self.body_volume / mm3,
                   self.body_mass / kg,
                   self.body_activity_bi214 / mBq,
                   self.body_self_shield_activity_bi214 / mBq,
                   self.body_activity_tl208 / mBq,
                   self.body_self_shield_activity_tl208 / mBq ,
                   self.body_transmittance,
                   self.head_thickness / mm,
                   2 * self.head_surface / mm2,
                   2 * self.head_volume / mm3,
                   2 * self.head_mass / kg,
                   self.head_activity_bi214 / mBq,
                   self.head_self_shield_activity_bi214 / mBq,
                   self.head_activity_tl208 / mBq,
                   self.head_self_shield_activity_tl208 / mBq,
                   self.head_transmittance
                  )
        return s

    __repr__ = __str__
