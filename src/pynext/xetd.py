from . system_of_units import *
from . math_functions import attenuation_factor
from . Material import PhysicalMaterial
from . PhysicalVolume import PhysicalVolume
from . Shapes import Cylinder, CylinderShell, FlatPlate

class Xenon:
    def __init__(self, gxe):

        self.rho = {'xe_2020': 124.3 * kg/m3,
                    'xe_3020': 203.35 * kg/m3,
                    'xe_1520':   89.9 * kg/m3,
                    'xe_1020':    58 * kg/m3,
                    'xe_0520':    30 * kg/m3,
                    'xe_0720':    40 * kg/m3,
                    'xe_0920':    50 * kg/m3,
                    'lxe'    :    3  * g/cm3}

        self.xe = PhysicalMaterial(name=gxe, rho=self.rho[gxe], mu_over_rho=0.039 * cm2/g)
        
    def __str__(self):
        g_cm3 = g / cm3
        cm2_g = cm2 / g
        icm   = 1 / cm

        s= """
        material                                   = %s
        density (rho)                              = %7.2g g/cm3
        mass attenuation coefficient (mu_over_rho) = %7.2g cm2/g
        attenuation coefficient (mu)               = %7.2g cm^-1
        attenuation length (Latt)                  = %7.2g cm
    """%(self.xe.name,
         self.xe.density / g_cm3,
         self.xe.mass_attenuation_coefficient / cm2_g,
         self.xe.attenuation_coefficient / icm,
         self.xe.attenuation_length / cm
         )

        return s

    __repr__ = __str__


class CXe:
    """
    Represents a cylinder filled up with xenon (can be gas or liquid)
    
    Xe -->xe (gas at P, T or Lxe)
    R --> fiducial volume radius
    L --> fiducial length
    Rb ---> radius of buffer gas
    Cb --> z of buffer gas in cathode
    Ab --> z of buffer gas in anode
    """
    def __init__(self, Xe, R, L, Rb, Cb, Ab):
        self.R  =R
        self.L  = L
        self.Xe = Xe 
        self.Rb = Rb
        self.Cb = Cb
        self.Ab = Ab
        self.xeShape  = Cylinder(R + Rb, L + Cb + Ab)
        self.fidShape = Cylinder(R,L)
        self.xe       = PhysicalVolume("XeVolume", Xe.xe, self.xeShape)
        self.fid      = PhysicalVolume("FiducialVolume", Xe.xe, self.fidShape)
       
        self.LbShape  = CylinderShell(R, R + Rb, L + Cb + Ab)
        self.CbShape  = FlatPlate(R + Rb, Cb)
        self.AbShape  = FlatPlate(R + Rb, Ab)

        self.xeLb     = PhysicalVolume("LongBuffer", Xe.xe, self.LbShape)
        self.xeCb     = PhysicalVolume("CathodeBuffer", Xe.xe, self.CbShape)
        self.xeAb     = PhysicalVolume("AnodeBuffer", Xe.xe,self.AbShape)

    @property
    def fiducialRadius(self):
        return self.R

    @property
    def xeRadius(self):
        return self.R+self.Rb

    @property
    def fiducialLength(self):
        return self.L
        
    @property
    def xeLength(self):
        return self.L + self.Cb + self.Ab

    @property
    def fiducialVolume(self):
        return self.fid.volume

    @property
    def xeVolume(self):
        return self.xe.volume

    @property
    def fiducialMass(self):
        return self.fid.mass
        
    @property
    def xeMass(self):
        return self.xe.mass
        
    @property
    def bufferMass(self):
        return self.xeLb.mass + self.xeCb.mass + self.xeAb.mass

    @property
    def readoutSurface(self):
        return self.xeCb.surface

    def __str__(self):
        s= """
        %s
        """%(self.Xe)

        s+= """
        Fiducial Radius = %7.2f cm  xe Radius = %7.2f cm
        Fiducial Length = %7.2f cm  Xe Length = %7.2f cm
        """%(self.fiducialRadius/cm,self.xeRadius/cm,
             self.fiducialLength/cm,self.xeLength/cm)

        s+="""
        Fiducial Volume = %7.2f m3 Xe Volume = %7.2f m3
        Fiducial mass = %7.2f kg Xe mass = %7.2f kg
        """%(self.fiducialVolume/m3, self.xeVolume/m3,
             self.fiducialMass/kg,self.xeMass/kg)

        s+= """
        Radial buffer radius = %7.2f cm (E/P ~ 2)
        Cathode buffer thickness = %7.2f cm   (E/P~1)
        Anode buffer thickness = %7.2f cm
        """%(self.Rb/cm, self.Cb/cm,self.Ab/cm)

        s+= """
        Long buffer Volume = %7.2f m3
        Long buffer mass = %7.2f kg
        """%(self.xeLb.volume/m3,self.xeLb.mass/kg)

        s+= """
        Cathode buffer Volume = %7.2f m3
        Cathode buffer mass = %7.2f kg
        """%(self.xeCb.volume/m3,self.xeCb.mass/kg)

        s+= """
        Anode buffer Volume = %7.2f m3
        Anode buffer mass = %7.2f kg
        """%(self.xeAb.volume/m3,self.xeAb.mass/kg)

        s+= """
        Total buffer Mass = %7.2f kg
        Cathode Readout Surface = %7.2f m2
        """%(self.bufferMass/kg,self.readoutSurface/m2)

        return s

    __repr__ = __str__