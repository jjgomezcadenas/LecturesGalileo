from . system_of_units import *
from math import pi

class PMT:
    def __init__(self,
                 name ='R1141',
                 D    = 3 * 2.5 * cm,
                 QE   = 0.3,
                 a_pmt_bi214 =  0.35 * mBq, a_pmt_tl208  = 0.19 * mBq,
                 a_base_bi214 = 1.22 * mBq, a_base_tl208 = 0.33 * mBq,
                 a_window_bi214 = 0.61 * mBq, a_window_tl208 = 0.38 * mBq
                 ):
        """
        Defines a PMT of diameter D, quantum efficiency QE and Activity A
        """
        self.name           = name
        self.D              = D
        self.R              = D / 2
        self.QE             = QE
        self.a_pmt_bi214    = a_pmt_bi214    / 2 #solid angle
        self.a_pmt_tl208    = a_pmt_tl208    / 2 # assume only forward radiation
        self.a_base_bi214   = a_base_bi214   / 2
        self.a_base_tl208   = a_base_tl208   / 2
        self.a_window_bi214 = a_window_bi214 / 2
        self.a_window_tl208 = a_window_tl208 / 2
        self.a_bi214        = self.a_pmt_bi214 + self.a_base_bi214 + self.a_window_bi214
        self.a_tl208        = self.a_pmt_tl208 + self.a_base_tl208 + self.a_window_tl208
        self.S            = pi * self.R**2

    def __str__(self):
        s="""
        PMT name                = %s
        PMT Diameter            = %7.2f cm
        PMT Surface             = %7.2f cm2
        PMT QE                  = %7.2f
        total activity:
        bi214 activity          = %7.2f mBq
        tl208 activity          = %7.2f mBq

        partial activities
        bi214 activity (PMT)    = %7.2f mBq
        tl208 activity (PMT)    = %7.2f mBq
        bi214 activity (Base)   = %7.2f mBq
        tl208 activity (Base)   = %7.2f mBq
        bi214 activity (Window) = %7.2f mBq
        tl208 activity (Window) = %7.2f mBq
        """%(self.name, self.D / cm, self.S / cm2, self.QE,
             self.a_bi214 / mBq, self.a_tl208 / mBq,
             self.a_pmt_bi214  / mBq,
             self.a_pmt_tl208  / mBq,
             self.a_base_bi214 / mBq,
             self.a_base_tl208 / mBq,
             self.a_window_bi214 / mBq,
             self.a_window_tl208 / mBq)
        return s

    __repr__ = __str__


class SiPM:

    def __init__(self,
                 name  = 'SENSL',
                 L     = 1.0 * mm,
                 QE    = 0.5,
                 TPB   = 0.5,
                 a_bi214 = 0.1 * muBq, a_tl208= 0.1 * muBq):
        """
        Defines a SiPM of size L, quantum efficiency QE and TPB efficiency TPB.
        """

        self.name    = name
        self.L       = L
        self.QE      = QE
        self.TPB     = TPB
        self.PDE     = QE*TPB
        self.S       = L**2
        self.a_bi214 = a_bi214 / 2
        self.a_tl208 = a_tl208 / 2

    def __str__(self):
        s="""
        SiPM name          = %s
        SiPM size          = %7.2f mm
        SiPM Surface       = %7.2f mm2
        SiPM QE            = %7.2f
        SiPM TPB eff       = %7.2f
        SiPM global PDE    = %7.2f
        bi214 activity     = %7.2f muBq
        bi214 activity     = %7.2f muBq
        """%(self.name, self.L / mm,
             self.S / mm2,
             self.QE, self.TPB,
             self.PDE, self.a_bi214/muBq,
             self.a_tl208 / muBq)
        return s

    __repr__ = __str__

class KDB:
    """A Kapton Dice Board"""

    def __init__(self,
                 L        = 100  * mm,
                 pitch    = 10   * mm,
                 nof_sipm = 64,
                 a_bi214  = 31 * muBq, a_tl208= 15 * muBq):
        """
        Defines a dice board (e.g, a substrate) where the SiPMs are mounted at a given pitch.
        """
        self.nof_sipm = nof_sipm
        self.L        = L
        self.pitch    = pitch
        self.S        = L**2
        self.a_bi214  = a_bi214 / 2
        self.a_tl208  = a_tl208 / 2

    def __str__(self):
        s="""
        nof_sipm          = %d
        KDB size          = %7.2f mm
        KDB pitch         = %7.2f mm
        KDB Surface       = %7.2f cm2
        bi214 activity    = %7.2f muBq
        bi214 activity    = %7.2f muBq
        """%(self.nof_sipm, self.L / mm, self.pitch / mm,
             self.S / cm2,
             self.a_bi214 / muBq,
             self.a_tl208 / muBq)
        return s

    __repr__ = __str__
