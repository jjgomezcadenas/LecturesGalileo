from dataclasses import dataclass
import numpy as np
import pandas as pd

from  . system_of_units import *





@dataclass
class WLS:
    name : str
    qeff : float

    def emission_tpb_(self, l):
        A     = 0.782
        alpha = 3.7e-2
        s1    = 15.43
        mu1   = 418.1
        s2    = 9.72
        mu2   = 411.2

        t1 = A * (alpha/2) * np.exp((alpha/2) * (2 * mu1 + alpha * s1**2 - 2 * l))
        t2  = erfc((mu1 + alpha * s1**2 - l)/(s1 * np.sqrt(2)))
        t3 = (1 - A) * (1/np.sqrt(2 * s2**2 * np.pi)) * np.exp(-(l - mu2)**2/(2*s2**2))
        return t1 * t2 + t3

    def __post_init__(self):
        # TIR in core to clad 1
        L = np.arange(350, 550, 1)
        e = self.emission_tpb_(L)
        self.emx = np.max(e)


    @property
    def quantum_efficiency(self):
        return self.qeff

    def emission_tpb(self, lamda):
        return self.emission_tpb_(lamda/nm) / self.emx


@dataclass
class FiberWLS:
    d      : float  #diameter for cyl, lateral size for square
    wls    : WLS
    qfib   : float  # quantum efficiency
    qptfe  : float  # reflectivity of PTFE around fiber
    ncore  : float
    nclad1 : float
    nclad2 : float
    latt   : float  # attenuation length  d
    path   : str = '/Users/jjgomezcadenas/Projects/Development/pynextsw/pynext/'

    def __post_init__(self):
        # TIR in core to clad 1
        self.thetac1 = np.arcsin(self.nclad1/self.ncore)
        self.thetat1  = 0.5 * np.pi - self.thetac1
        self.ptir1 = (1 - np.cos(self.thetat1))  # 2 x 2 pi (1 - cos(theta)) /4 pi (forward and backward)
        #self.ptir1 = np.cos(self.thetac1)  # 2 x 2 pi (cos(theta) /4 pi (forward and backward)

        # refracted to clad2 : critical
        self.theta2 = np.arcsin((self.nclad1 / self.nclad2) * np.sin(self.thetac1))

        # TIR in clad1 to clad2
        self.thetac2 = np.arcsin(self.nclad2/self.nclad1)
        self.thetat2  = 0.5 * np.pi - self.thetac2
        self.ptir2 = (1 - np.cos(self.thetat2))  # fraction between two claddings
        #self.ptir2 = np.cos(self.thetac2)

        self.teff1 = self.qfib * self.ptir1
        self.teff2 = self.qfib * self.ptir2
        self.fabs = pd.read_csv('stGobainAbs.csv', delimiter=',')
        self.fem = pd.read_csv('stGobainEm.csv', delimiter=',')
        self.fmu = pd.read_csv('elijanMu.csv', delimiter=',')
        L        = np.arange(350, 600, 1) * nm
        num      = np.sum(self.wls.emission_tpb(L) *  self.blue_absorption(L) )
        den      = np.sum(self.wls.emission_tpb(L))
        self.Pa  = num / den

    @property
    def quantum_efficiency(self)->float:
        return self.qfib

    @property
    def diameter(self)->float:
        return self.d

    @property
    def trapping_efficiency_c1(self)->float:
        return self.teff1

    @property
    def trapping_efficiency_c2(self)->float:
        return self.teff2

    @property
    def trapping_efficiency(self)->float:
        return self.teff1 + self.teff2

    def transmittance(self, d)->float:
        return np.exp(-d / self.latt)

    def absorption(self, d)->float:
        return 1 - self.transmittance(d)

    def blue_transmittance(self, lamda)->float:
        mu = self.wls_mu(lamda)
        dcm = self.d/cm
        #print(lamda/nm, dcm, mu, np.exp(-dcm * mu))
        return np.exp(-dcm * mu)

    def blue_absorption(self, lamda)->float:
        return 1 - self.blue_transmittance(lamda)

    @property
    def blue_absorption_probability(self)->float:
        return  self.Pa

    def wls_absorption(self, lamda : float)->float:
        return np.interp(lamda/nm, self.fabs.WL.values, self.fabs.A.values)

    def wls_emission(self, lamda : float)->float:
        return np.interp(lamda/nm, self.fem.WL.values, self.fem.A.values)

    def wls_mu(self, lamda : float)->float:
        return np.interp(lamda/nm, self.fmu.WL.values, self.fmu.mu.values)

    def __str__(self):
        s= f"""
        diameter ={self.d/mm} mm, Q = {self.qfib}, PTFE refl = {self.qptfe}
        ncore = {self.ncore}, nclad1 ={self.nclad1}, nclad2 ={self.nclad2}
        Absoprtion prob at 450 nm     = {self.blue_absorption(450*nm)}
        Trapping efficieny            = {self.trapping_efficiency}
        Fiber coated with WLS         = {self.wls.name}
        WLS QE                        = {self.wls.quantum_efficiency}

    """
        return s

    __repr__ = __str__



@dataclass
class SiPM:
    name   : str
    xsize  : float
    PDE    : float
    C      : float # capacitance
    Rs     : float # Series resistance connection
    path   : str = '/Users/jjgomezcadenas/Projects/Development/pynextsw/pynext/'

    def __post_init__(self):


        self.fpde = pd.read_csv(self.path+'s13360.csv', delimiter=',')
        self.fdcr = pd.read_csv(self.path+'dcrt2.csv', delimiter=',')
        tf        = np.polyfit(self.fdcr.K.values, np.log(self.fdcr.DCR.values), 2)
        self.TK   = np.poly1d(tf)
        self.F = 1.0 # to achieve the calibrated value o DCR

    @property
    def area(self):
        return self.xsize**2

    def pde(self, lamda : float)->float:
        return np.interp(lamda, self.fpde.WL.values, self.fpde.PDE.values/100)

    def log_dcr(self, t : float)->float: # t in 1/K, K kelvin
        return self.TK(t)

    def dcr_sipm_per_unit_area (self, tC : float)->float:  # tC in Celsius
        tK = tC + 273.15
        ldcr = self.log_dcr(1/tK)
        return self.F* np.exp(ldcr) * hertz / mm2

    def dcr_sipm_per_time (self, tC : float, time: float) ->float:
        return self.area * time * self.dcr_sipm_per_unit_area(tC)

    def __str__(self):
        s =f"""
        sensor ={self.name}, size = {self.xsize/mm} mm, PDE = {self.PDE}
        capacitance = {self.C/pF:.2f} pF;
        """

        return s

    __repr__ = __str__


@dataclass
class TpcEL:
    """
    Defines a EL TPC
    EP = E/P
    dV = drift voltage
    P  = pressure
    d  = EL grid gap
    L  = drift lenght
    Ws = energy to produce a scintillation photon
    Wi = energy to produce a ionization photon
    """
    EP : float = 3.5 * kilovolt / (cm * bar)
    dV : float = 0.5 * kilovolt / cm
    P  : float =  15 * bar
    d  : float =   5 * mm
    L  : float = 120 * cm
    Ws : float =  39.2 * eV
    Wi : float =  21.9 * eV

    def __post_init__(self):
        u = kilovolt / (cm * bar)
        ep = self.EP/u
        self.Vgrid = self.EP * self.d * self.P
        self.Vc    = self.Vgrid + self.dV * self.L
        self.YP    = 140 * ep - 116  # in photons per electron bar^-1 cm^-1
        self.Ng    = self.YP * self.d/cm * self.P/bar #photons/e

    @property
    def e_over_p(self):
        return self.EP

    @property
    def drift_voltage(self):
        return self.dV

    @property
    def grid_voltage(self):
        return self.Vgrid

    @property
    def cathode_voltage(self):
        return self.Vc

    @property
    def pressure(self):
        return self.P

    @property
    def grid_distance(self):
        return self.d

    @property
    def drift_length(self):
        return self.L

    @property
    def optical_gain(self):
        return self.Ng

    def scintillation_photons(self, E):
        return E / self.Ws

    def el_photons(self, E):
        return self.ionization_electrons(E) * self.optical_gain

    def ionization_electrons(self, E):
        return E / self.Wi


    def __str__(self):
        kv_cm_b = kilovolt / (cm * bar)
        kv_cm   = kilovolt / cm
        cm_bar = 1 / (cm * bar)
        s= """
        E/P = %7.2f kV * cm^-1* bar^-1
        dV = drift voltage = %7.2f kV * cm^-1
        P  = pressure = %7.2f bar
        d  = EL grid gap = %7.2f mm
        L  = drift lenght =%7.2f m
        Grid voltage = %7.2f kV
        Cathode voltage = %7.2f kV
        Yield =  %7.2e photons/e

    """%(self.e_over_p / kv_cm_b,
         self.drift_voltage / kv_cm,
         self.pressure / bar,
         self.grid_distance / mm,
         self.drift_length / m,
         self.grid_voltage / kilovolt,
         self.cathode_voltage /kilovolt,
         self.optical_gain)

        s+="""
        Primary scintillation photons per MeV = %7.2e
        Primary ionization electrons per MeV = %7.2e
        EL photons per MeV                   = %7.2e
        """%(self.scintillation_photons(1 * MeV),
             self.ionization_electrons(1 * MeV),
             self.el_photons(1 * MeV))
        s+="""
        Primary scintillation Krypton = %7.2e
        Primary ionization electrons Krypton = %7.2e
        EL photons Krypton                   = %7.2e
        """%(self.scintillation_photons(41.5 * keV),
             self.ionization_electrons(41.5 * keV),
             self.el_photons(41.5 * keV))
        s+="""
        Primary scintillation Qbb = %7.2e
        Primary ionization electrons Qbb = %7.2e
        EL photons Qbb                   = %7.2e
        """%(self.scintillation_photons(2458 * keV),
             self.ionization_electrons(2458 * keV),
             self.el_photons(2458 * keV))

        return s

    __repr__ = __str__

@dataclass
class TpcXe:
    """A cylynder filled with Xenon"""
    cyl  : Cylinder
    name : str = 'Xenon'
    p    : float = 15 * bar


    def __post_init__(self):
        self.ps   = (10,15,20,30, 40)
        self.rhos = (58, 89.9, 124.3, 203.35, (203.35) *(4/3))

    @property
    def pressure(self):
        return self.p

    @property
    def density(self):
        P = self.p/bar
        return np.interp(P, self.ps, self.rhos) * kg / m3

    @property
    def mass(self):
        return self.cyl.volume * self.density

    def __str__(self):
        g_cm3 = g / cm3
        cm2_g = cm2 / g
        icm   = 1 / cm

        s= f"""
        material       = {self.name:s}
        density (rho)  = {self.density/g_cm3:7.2f} g/cm3
        mass           = {self.mass/kg:7.2f}  kg
        Cylinder       = {self.cyl}
    """

        return s

    __repr__ = __str__


@dataclass
class FiberDetector:
    tpcel     : TpcEL
    tpcxe     : TpcXe
    fwls      : FiberWLS  # fiber type
    sipm      : SiPM
    eff_t     : float # transport efficiency of the detector.
    sampling  : float
    adcPerPes : float
    tempC      : float

    def __post_init__(self):
        self.transport = self.eff_t * self.fwls.wls.quantum_efficiency * self.fwls.trapping_efficiency
        self.attenuation = self.fwls.transmittance(self.tpcxe.cyl.length/2)


    @property
    def operating_temperatureC (self):
        return self.tempC

    @property
    def n_fibers (self):
        return int(self.tpcxe.cyl.perimeter / self.fwls.diameter) + 1

    @property
    def efficiency (self):
        return self.transport * self.attenuation

    def __str__(self):
        ren = self.efficiency *  self.sipm.PDE/ self.n_fibers
        detpde = self.efficiency * self.sipm.PDE
        s =f"""
        gas pressure  = {self.tpcxe.pressure/bar:7.2f} bar
        gas density   = {self.tpcxe.density/(g/cm3):7.2f} g/cm3
        gas mass      = {self.tpcxe.mass/kg:7.2f}  kg

        Dimensions(Cylinder)  = {self.tpcxe.cyl}

        Fibers efficiency = {100*self.efficiency:.2f} %
        of which: Transport = {self.transport} % & attenuation = {self.attenuation} %
        SiPM PDE          = {self.sipm.PDE:.2f}
        SiPM size         = {self.sipm.xsize/mm:.2f} mm
        fiber size        = {self.fwls.diameter/mm:.2f} mm
        Sampling S1       = {self.sampling/ns:.2f} ns
        ADC counts 1 PE   = {self.adcPerPes:d}
        number of fibers  = {self.n_fibers:d}

        Primary scintillation Krypton          = {self.tpcel.scintillation_photons(41.5 * keV):.2e}
        EL photons Krypton                     = {self.tpcel.el_photons(41.5 * keV):.2e}
        Primary scintillation Qbb              = {self.tpcel.scintillation_photons(2458 * keV):.2e}
        EL photons Qbb                         = {self.tpcel.el_photons(2458 * keV):.2e}

        Primary scintillation Krypton detected = {self.tpcel.scintillation_photons(41.5 * keV)* detpde:.2e}
        EL photons Krypton detected            = {self.tpcel.el_photons(41.5 * keV)* detpde:.2e}
        Primary scintillation Qbb detected     = {self.tpcel.scintillation_photons(2458 * keV)* detpde:.2e}
        EL photons Qbb detected                = {self.tpcel.el_photons(2458 * keV)* detpde:.2e}

        Primary scintillation Krypton det/fiber = {self.tpcel.scintillation_photons(41.5 * keV) * ren:.2e}
        EL photons Krypton det/fiber            = {self.tpcel.el_photons(41.5 * keV)* ren:.2e}
        Primary scintillation Qbb det/fiber     = {self.tpcel.scintillation_photons(2458 * keV)* ren:.2e}
        EL photons Qbb det/fiber                = {self.tpcel.el_photons(2458 * keV)* ren:.2e}

        Number of DCR photons in the detector for:
         --operating temperature ({self.operating_temperatureC:.2f} C)
         --sampling time of {self.sampling/ns:.2f}
         -- nDCR = {self.sipm.dcr_sipm_per_time(self.operating_temperatureC, self.sampling) * self.n_fibers:.2f}
        """

        return s

    __repr__ = __str__
