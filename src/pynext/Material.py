"""
Material
Defines a material medium
"""
from . system_of_units import *
from math import pi, exp, log
import numpy as np
from scipy.integrate import quad


class SelfAtt:
    def __init__(self, mu, L):
        self.mu = mu
        self.L  = L

    def f(self, theta):
        a1 = np.cos(theta) / (self.mu * self.L)
        a2 = np.exp(-self.mu * self.L / np.cos(theta))
        return a1 * (1 - a2)

class PhysicalMaterial:
    """mu_over_rho is the mass attenuation coefficient at 2.5 MeV"""
    def __init__(self, name, rho, mu_over_rho):

        self.name        = name
        self.rho         = rho
        self.mu_over_rho = mu_over_rho
        self.mu          = mu_over_rho * rho
        self.Latt        = 1 / self.mu

    @property
    def density(self):
        return self.rho

    @property
    def mass_attenuation_coefficient(self):
        return self.mu_over_rho

    @property
    def attenuation_coefficient(self):
        return self.mu

    @property
    def attenuation_length(self):
        return self.Latt

    def transmittance_at_qbb(self, z):
        return exp(-z*self.mu)

    def absorption_at_qbb(self, z):
        return 1 - exp(-z*self.mu)

    def __str__(self):
        g_cm3 = g / cm3
        cm2_g = cm2 / g
        icm   = 1 / cm

        s= """
        material                                   = %s
        density (rho)                              = %7.2f g/cm3
        mass attenuation coefficient (mu_over_rho) = %7.2f cm2/g
        attenuation coefficient (mu)               = %7.2f cm^-1
        attenuation length (Latt)                  = %7.2f cm
    """%(self.name,
         self.density / g_cm3,
         self.mass_attenuation_coefficient / cm2_g,
         self.attenuation_coefficient / icm,
         self.attenuation_length / cm
         )

        return s

    __repr__ = __str__


class RadioactiveMaterial(PhysicalMaterial):
    def __init__(self, name, rho, mu_over_rho, a_bi214, a_tl208):

        super().__init__(name, rho, mu_over_rho)
        self.a_bi214        = a_bi214
        self.a_tl208        = a_tl208
        self.C = 1 / 3

    @property
    def mass_activity_bi214(self):
        return self.a_bi214

    @property
    def mass_activity_tl208(self):
        return self.a_tl208

    # def surface_activity(self, z, isotope='Bi214'):
    #     """Assuming that the material is an infinte slab
    #        with an specific activity A0 (Bq/[M]) and a thickness z ([L])
    #        then the number of decays (N) is:
    #
    #        N = A0 * M
    #
    #        The flux Phi (N /[S]) that escapes the material (not self-shielded) trhough
    #        one of the surfaces is given by:
    #
    #        Phi  = rho * A0 * z * I / (2pi),
    #        where:
    #        I = Integral (-pi/2, pi/2) { (cos(theta)/mu L) (1 - exp(-mu L/cos(theta))) }
    #
    #     """
    #
    #     self.SA =  self.rho * self.a_bi214 * z # Bi 214 by default
    #     if isotope == 'Tl208':
    #         self.SA =  self.rho * self.a_tl208 * z
    #
    #     att = SelfAtt(mu=self.mu, L=z)
    #     tf, _ = quad(att.f, -np.pi/2 + 0.0001, np.pi/2 -  0.0001)
    #
    #     return self.SA  * (tf / (2 * np.pi))
    #
    def __str__(self):
        bq_kg = Bq / kg

        s = super().__str__() + """
        activity Bi-214                = %7.2e Bq /kg
        activity Tl-208                = %7.2e Bq /kg
    """%(self.mass_activity_bi214 / bq_kg,
         self.mass_activity_tl208 / bq_kg)
    
        return s

    __repr__ = __str__

class PVMaterial(RadioactiveMaterial):
    """Material used for construction of PV
    Sm is the maximum allowable strength of the material
    """

    def __init__(self, name, rho, mu_over_rho, a_bi214, a_tl208, Sm):

        super().__init__(name, rho, mu_over_rho, a_bi214, a_tl208)
        self.Sm           = Sm

    @property
    def maximum_allowable_strength(self):
            return self.Sm


    def __str__(self):

        s = super().__str__() + """
        maximum_allowable_strength = %7.2e MPa
    """%(self.maximum_allowable_strength / MPa)

        return s

    __repr__ = __str__

class GXe:
    def __init__(self, rho=89.9 * kg/m3):

        self.rho = {'rho_2020': 124.3 * kg/m3,
                    'rho_3020': 203.35 * kg/m3,
                    'rho_1520':   89.9 * kg/m3,
                    'rho_1020':    58 * kg/m3,
                    'rho_0520':    30 * kg/m3,
                    'rho_0720':    40 * kg/m3,
                    'rho_0920':    50 * kg/m3}

        self.xed = {}
        self.xe = PhysicalMaterial(name='GXe', rho=rho, mu_over_rho=0.039 * cm2/g)
        for name, r in self.rho.items():
            self.xed[name] = PhysicalMaterial(name='GXe + %s'%name,
                                              rho=r, mu_over_rho=0.039 * cm2/g)



A_BI214_316Ti   =    1    * mBq/kg
A_TL208_316Ti   =    0.15 * mBq/kg
A_BI214_CU_LIM  =   12    * muBq/kg
A_TL208_CU_LIM  =    1.4  * muBq/kg
A_BI214_CU_BEST =    3    * muBq/kg
A_TL208_CU_BEST =    1.4  * muBq/kg
A_BI214_PB      =  370    * muBq/kg
A_TL208_PB      =   73    * muBq/kg
A_BI214_Poly    =   62    * muBq/kg
A_TL208_Poly    =    8    * muBq/kg

vacuum = RadioactiveMaterial(name='vacuum', rho=1e-25 * g/cm3, mu_over_rho=1e-25 * cm2/g,
                             a_bi214=0, a_tl208=0)
ti316  = RadioactiveMaterial(name='316ti', rho=7.87 * g/cm3, mu_over_rho=0.039 * cm2/g,
                            a_bi214=A_BI214_316Ti, a_tl208=A_TL208_316Ti)

cu12 = RadioactiveMaterial(name='CuUpperLimits',
                           rho = 8.96 * g/cm3,
                           mu_over_rho = 0.039 * cm2/g,
                           a_bi214 = A_BI214_CU_LIM,
                           a_tl208 = A_TL208_CU_LIM )
cu03 = RadioactiveMaterial(name='CuBest',
                           rho = 8.96 * g/cm3,
                           mu_over_rho = 0.039 * cm2/g,
                           a_bi214 = A_BI214_CU_BEST,
                           a_tl208 = A_TL208_CU_BEST )

pb =   RadioactiveMaterial(name='Pb',
                           rho = 11.33 * g/cm3,
                           mu_over_rho = 0.044 * cm2/g,
                           a_bi214 = A_BI214_PB,
                           a_tl208 = A_TL208_PB )

poly =   RadioactiveMaterial(name='Poly',
                           rho = 0.97 * g/cm3,
                           mu_over_rho = 1e-6 * cm2/g,
                           a_bi214 = A_BI214_Poly,
                           a_tl208 = A_TL208_Poly )
