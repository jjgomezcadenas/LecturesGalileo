from pynext.system_of_units import *

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
    def __init__(self,
                 EP = 3.5 * kilovolt / (cm * bar),
                 dV = 0.5 * kilovolt / cm,
                 P  =  15 * bar,
                 d  =   5 * mm,
                 L  = 120 * cm,
                 Ws =  24 * eV,
                 Wi =  16 * eV):

        u = kilovolt / (cm * bar)
        ep = EP/u

        self.EP    = EP
        self.dV    = dV
        self.P     =  P
        self.d     =  d
        self.L     =  L
        self.Ws    = Ws
        self.Wi    = Wi
        self.Vgrid = EP * d * P
        self.Vc    = self.Vgrid + dV * L
        self.YP    = 140 * ep - 116  # in photons per electron bar^-1 cm^-1
        self.Ng    = self.YP * self.d/cm * self.P/bar #photons/e

    @property
    def EoverP(self):
        return self.EP

    @property
    def DriftVoltage(self):
        return self.dV

    @property
    def GridVoltage(self):
        return self.Vgrid

    @property
    def CathodeVoltage(self):
        return self.Vc

    @property
    def Pressure(self):
        return self.P

    @property
    def GridDistance(self):
        return self.d

    @property
    def DriftLength(self):
        return self.L

    @property
    def OpticalGain(self):
        return self.Ng

    @property
    def WS(self):
        return self.Wd

    @property
    def WI(self):
        return self.Wi

    def ScintillationPhotons(self, E):
        return E / self.Ws

    def ELPhotons(self, E):
        return self.IonizationElectrons(E) * self.OpticalGain

    def IonizationElectrons(self, E):
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

    """%(self.EoverP / kv_cm_b,
         self.DriftVoltage / kv_cm,
         self.Pressure / bar,
         self.GridDistance / mm,
         self.DriftLength / m,
         self.GridVoltage / kilovolt,
         self.CathodeVoltage /kilovolt,
         self.OpticalGain)

        s+="""
        Primary scintillation photons per MeV = %7.4e
        Primary ionization electrons per MeV = %7.4e
        EL photons per MeV                   = %7.4e
        """%(self.ScintillationPhotons(1 * MeV),
             self.IonizationElectrons(1 * MeV),
             self.ELPhotons(1 * MeV))

        return s

    __repr__ = __str__
