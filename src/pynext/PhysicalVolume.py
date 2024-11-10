"""
Physical Volume
Defines a material medium + Shape
"""

from . system_of_units import *
from . math_functions import attenuation_factor

class PhysicalVolume:
    def __init__(self,name, material,shape):
        """
       Defines a physical volume

       """
        self.name     = name
        self.material = material
        self.shape    = shape
        self.M        = self.shape.V * self.material.rho

    @property
    def V(self):
        return self.shape.V

    @property
    def S(self):
        return self.shape.S

    @property
    def radius(self):
        return self.shape.radius()

    @property
    def thickness(self):
        return self.shape.thickness()

    @property
    def volume(self):
        return self.shape.V

    @property
    def surface(self):
        return self.shape.S

    @property
    def mass(self):
        return self.M

    @property
    def activity_bi214(self):
        return self.M * self.material.mass_activity_bi214

    @property
    def activity_tl208(self):
        return self.M * self.material.mass_activity_tl208

    @property
    def activity_tl208(self):
        return self.M * self.material.mass_activity_tl208

    def transmittance(self, z):
        return self.material.transmittance_at_qbb(z)

    def absorption(self, z):
        return self.material.absorption_at_qbb(z)

    def activity_bi214_self_shield(self, z):
        return self.activity_bi214 * attenuation_factor(self.material.mu, z)

    def activity_tl208_self_shield(self, z):
        return self.activity_tl208 * attenuation_factor(self.material.mu, z)


    def __str__(self):

        s =  """ Physical Volume:
        Shape           = %s
        Material        = %s
        radius          = %7.2e  mm
        thickness       = %7.2e  mm
        Volume          = %7.2e  m3
        Surface         = %7.2e  m2
        Mass            = %7.2e  kg
        activity Bi-214 = %7.2e  mBq
        activity Tl-208 = %7.2e  mBq
    """%(self.shape,
         self.material,
         self.radius / mm,
         self.thickness / mm,
         self.V / m3,
         self.S / m2,
         self.M / kg,
         self.activity_bi214 / mBq,
         self.activity_tl208 / mBq)

        return s

    __repr__ = __str__
