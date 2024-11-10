"""
NEXT100
Defines NEXT100 experiment
"""
from . PhysicalVolume import *
from . Sensors import *
from . TpcEL import *
import pynext.NextData as ND


#Rejection factors

Bi214_RF={0.5:5.6e-7,1.0:1.7e-6,2.0:2.1e-6,3.0:2.5e-6,4.0:2.7e-6}
Tl208_RF={0.5:6.1e-7,1.0:1.7e-6,2.0:4.3e-6,3.0:1.4e-5,4.0:4.5e-5}
Qbb_Xe136 =2462.0*keV
TF= 1./10.
R={0.5:0.03,0.9:0.09}
Qe = e_SI*coulomb
NcPMT=19.
NaPMT=19.


class AngelXe:
    """
    The xenon volume
    R --> fiducial volume radius
    L --> fiducial length
    Rb ---> radius of buffer gas
    Cb --> z of buffer gas in cathode
    Ab --> z of buffer gas in anode
    """
    def __init__(self,R,L,Rb,Cb,Ab,P):
        self.R=R
        self.L=L
        self.gxe = Material('GXe','BI214',P)
        self.P=P*bar

        self.xeShape = Cylinder(R+Rb,L+Cb+Ab)
        self.fidShape = Cylinder(R,L)
        self.xe = PhysicalVolume("XeVolume",self.gxe,self.xeShape)
        self.fid = PhysicalVolume("FiducialVolume",self.gxe,self.fidShape)
        self.Rb = Rb
        self.Cb = Cb
        self.Ab = Ab

        self.LbShape = CylinderShell(R,R+Rb,L+Cb+Ab)
        self.CbShape = FlatPlate(R+Rb,Cb)
        self.AbShape = FlatPlate(R+Rb,Ab)

        self.xeLb = PhysicalVolume("LongBuffer",self.gxe,self.LbShape)
        self.xeCb = PhysicalVolume("CathodeBuffer",self.gxe,self.CbShape)
        self.xeAb = PhysicalVolume("AnodeBuffer",self.gxe,self.AbShape)

    def FiducialRadius(self):
        return self.R

    def XeRadius(self):
        return self.R+self.Rb

    def FiducialLength(self):
        return self.L

    def XeLength(self):
        return self.L+self.Cb+self.Ab

    def FiducialVolume(self):
        return self.fid.Volume()

    def XeVolume(self):
        return self.xe.Volume()

    def FiducialMass(self):
        return self.fid.Mass()

    def XeMass(self):
        return self.xe.Mass()

    def BufferMass(self):
        return self.xeLb.Mass()+self.xeCb.Mass()+self.xeAb.Mass()

    def ReadoutSurface(self):
        return self.xeCb.Surface()

    def __str__(self):
        s= """
        Pressure =%7.2f bar
        Fiducial Radius = %7.2f cm  xe Radius = %7.2f cm
        Fiducial Length = %7.2f cm  Xe Length = %7.2f cm
        """%(self.P/bar,self.FiducialRadius()/cm,self.XeRadius()/cm,
             self.FiducialLength()/cm,self.XeLength()/cm)

        s+="""
        Fiducial Volume = %7.2f m3 Xe Volume = %7.2f m3
        Fiducial mass = %7.2f kg Xe mass = %7.2f kg
        """%(self.FiducialVolume()/m3, self.XeVolume()/m3,
             self.FiducialMass()/kg,self.XeMass()/kg)

        s+= """
        Radial buffer radius = %7.2f cm (E/P ~ 2)
        Cathode buffer thickness = %7.2f cm   (E/P~1)
        Anode buffer thickness = %7.2f cm
        """%(self.Rb/cm, self.Cb/cm,self.Ab/cm)

        s+= """
        Long buffer Volume = %7.2f m3
        Long buffer mass = %7.2f kg
        """%(self.xeLb.Volume()/m3,self.xeLb.Mass()/kg)

        s+= """
        Cathode buffer Volume = %7.2f m3
        Cathode buffer mass = %7.2f kg
    """%(self.xeCb.Volume()/m3,self.xeCb.Mass()/kg)

        s+= """
        Anode buffer Volume = %7.2f m3
        Anode buffer mass = %7.2f kg
    """%(self.xeAb.Volume()/m3,self.xeAb.Mass()/kg)

        s+= """
        Total buffer Mass = %7.2f kg
        Cathode Readout Surface = %7.2f m2
    """%(self.BufferMass()/kg,self.ReadoutSurface()/m2)

        return s

class AngelCan:
    """
    The pressure vessel
    Rin --> inner radius of shell
    Rout --> outer radious of shell
    L ---> length
    tp -->thickness of plate
    The can is made by a cylindrical shell and two flat end plates
    """


    def __init__(self,Rin=ND.PV_IR,Rout=ND.PV_OR,L=ND.PV_Z,Cp=ND.PV_PLATE_TH,Ap=ND.PV_PLATE_TH,
                 cplateType='Flat', aplateType='Flat',
                 FlangeRin=1140*mm/2., FlangeRout=1600*mm/2., FlangeT=90*mm, Mat=ND.PV_MAT):
        self.L=L
        self.Rin=Rin
        self.Rout=Rout
        self.Cp =Cp
        self.Ap =Ap
        self.t= Rout-Rin

        self.FlangeRin = FlangeRin
        self.FlangeRout = FlangeRout
        self.FlangeT = FlangeT

        self.Mat = Material(Mat)

        self.ShellShape = CylinderShell(Rin,Rout,L)
        self.FlangeShape= CylinderShell(FlangeRin,FlangeRout,FlangeT)

        self.Shell =  PhysicalVolume("Shell",self.Mat,self.ShellShape)
        self.Flange =  PhysicalVolume("Flange",self.Mat,self.FlangeShape)

        if cplateType == 'Flat':
            self.cPlateShape = FlatPlate(Rout,Cp)
        else:
            self.cPlateShape = SemiSphereShell(Rout,Cp)
        if aplateType == 'Flat':
            self.aPlateShape = FlatPlate(Rout,Cp)
        else:
            self.aPlateShape = SemiSphereShell(Rout,Cp)

        self.cPlate = PhysicalVolume("cPlate",self.Mat,self.cPlateShape)
        self.aPlate = PhysicalVolume("aPlate",self.Mat,self.aPlateShape)


    def Material(self):
        return self.Mat.Name()

    def InnerRadius(self):
        return self.Rin

    def OuterRadius(self):
        return self.Rout

    def Length(self):
        return self.L

    def CathodePlateThickness(self):
        return self.Cp

    def AnodePlateThickness(self):
        return self.Ap

    def FlangeThickness(self):
        return self.FlangeT


    def Mass(self):
        mass = self.Shell.Mass() + self.cPlate.Mass()+self.aPlate.Mass()
        #mass+=4*self.Flange.Mass()  # no flanges in steel
        return mass

    def Activity(self,FlangeFactor=0.2):
        a= self.Shell.Activity(self.t)/(1./year)
        a+= self.cPlate.Activity(self.Cp)/(1./year)
        a+= self.aPlate.Activity(self.Ap)/(1./year)
        #a+= 4*self.Flange.Activity(self.FlangeT)/(1./year)*FlangeFactor
        return a

    def __str__(self):
        s= """
        Material = %s
        Inner Radius = %7.2f cm
        Outer Radius = %7.2f cm
        Length = %7.2f cm
        Cathode plate thickness = %7.2f cm
        Anode plate thickness = %7.2f cm
        Flange thickness = %7.2f cm
        """%(self.Material(), self.InnerRadius()/cm,self.OuterRadius()/cm,
             self.Length()/cm,
             self.CathodePlateThickness()/cm,
             self.AnodePlateThickness()/cm,
             self.FlangeThickness()/cm)

        s+= """
        Shell Volume = %7.2f m3
        Shell mass = %7.2f kg
        Shell activity = %7.0e c/year
        """%(self.Shell.Volume()/m3,self.Shell.Mass()/kg,
             self.Shell.Activity(self.t)/(1./year))

        s+= """
        Cathode plate Volume = %7.2f m3
        Cathode plate mass = %7.2f kg
        Cathode plate activity = %7.0e c/year
        """%(self.cPlate.Volume()/m3,self.cPlate.Mass()/kg,
             self.cPlate.Activity(self.Cp)/(1./year))

        s+= """
        Anode plate Volume = %7.2f m3
        Anode plate mass = %7.2f kg
        Anode plate activity = %7.0e c/year
        """%(self.cPlate.Volume()/m3,self.aPlate.Mass()/kg,
             self.aPlate.Activity(self.Ap)/(1./year))

        s+= """
        Flange Volume = %7.2f m3
        Flange mass = %7.2f kg
        Flange activity = %7.0e c/year
        """%(self.Flange.Volume()/m3,self.Flange.Mass()/kg,
             self.Flange.Activity(self.FlangeT)/(1./year))


        s+= """
        Can total mass = %7.2f kg
        Can total activity = %7.2e counts/year
        """%(self.Mass()/kg,self.Activity())
        return s



class AngelSiPM:

    def __init__(self, sipm, angelEL, Surface=pi*(53)**2*cm2, pitch=1.1*cm,
                 asic=64, pdm=25):
        """
        Defines a SiPM readout
        """
        self.sipm=sipm
        self.pitch = pitch
        self.el = angelEL
        self.Sp = self.pitch**2
        self.S=Surface
        self.asic=asic
        self.pdm = pdm

    def NumberAsics(self):
        return self.NumberSiPM() / self.asic

    def NumberPDM(self):
        return self.NumberSiPM() / self.pdm

    def QE(self):
        return self.sipm.QE

    def PDE(self):
        return self.sipm.PDE

    def Surface(self):
        return self.S

    def Pitch(self):
        return self.pitch

    def AveragePhotons(self):
        n0 = self.el.dGdx()/(self.pitch)
        R=self.el.GridDistance()
        n=n0/(R**2)*self.S/self.Sp
        return n

    def LinearPackingFactor(self):
        return self.sipm.L/self.pitch

    def NumberSiPM(self):
        return self.LinearPackingFactor()**2*self.S/self.sipm.S

    def AreaSiPM(self):
        return self.sipm.S * self.NumberSiPM()

    def Coverage(self):
        return self.AreaSiPM()/self.Surface()

    def Activity(self):
        a = self.sipm.AU + self.sipm.ATh
        return a*year

    def PesPerPhotonReachingSensorPlane(self):
        return self.Coverage() * self.PDE()



    def __str__(self):
        s= """
        SiPM QE = %7.2f
        SiPM PDE = %7.2f
        SiPM size = %7.2f mm
        SiPM Pitch = %7.2f mm
        Linear Packing factor = %7.2f
        Detection Surface = %7.2g m2
        SiPM Surface = %7.2g m2
        Coverage = %7.2g
        nof SiPM = %7.2g
        nof ASIC = %7.2g
        nof PDM = %7.2g
        """%(self.sipm.QE,
             self.sipm.PDE,
             self.sipm.L/mm,
             self.Pitch()/mm,
             self.LinearPackingFactor(),
             self.Surface()/m2,
             self.AreaSiPM()/m2,
             self.Coverage(),
             self.NumberSiPM(),
             self.NumberAsics(),
             self.NumberPDM()
             )

        return s

class AngelPMT:
    """
    PMT readout
    n1 is the refraction index of the PMT window (quartz)
    n2 is the refraction index of the pressure window (sapphire)
    """
    def __init__(self,pmt,angelXe,coverage=0.35,n1=1.5,n2=1.8):
        self.xe = angelXe
        self.pmt = pmt
        self.coverage = coverage
        self.n1=n1
        self.n2=n2

    def QE(self):
        return self.pmt.QE

    def NumberPMT(self):
        return self.coverage*self.xe.ReadoutSurface()/self.pmt.S

    def CoveragePMT(self):
        return self.coverage

    def WindowT(self):
        R = ((self.n2-self.n1) / (self.n2+ self.n1))**2
        return (1-R)

    def WindowTransmittance(self):
        """
        0.93 ---> pass sapphire window
        0.85 ---> pass quartz window
        0.75 reach photocathode
        """
        return 0.93*0.85*0.75

    def TPBTransmittance(self):
        return 0.5

    def Activity(self):
        a = self.pmt.AU + self.pmt.ATh
        return a*year

    def __str__(self):
        s= """
        PMT name = %s
        PMT QE = %7.2f
        PMT Diameter = %7.2f cm
        PMT Surface = %7.2f cm2
        PMT Window refraction indez = %7.2f
        Optical Window refraction indez = %7.2f
        Optical Transmittance = %7.2f
        nof PMT (coverage=%7.2f) = %7.0f
        """%(self.pmt.name,self.pmt.QE,
             self.pmt.D/cm,self.pmt.S/cm2,
             self.n1,self.n2,self.WindowTransmittance(),
             self.coverage,self.NumberPMT())

        s+="""
        Activity = %7.2e c/year
        """%(self.Activity())
        return s


class AngelEL(TpcEL):
    """
    Defines an Asymmetric EL TPC
    EP = E/P
    dV = drift voltage
    P  = pressure
    d  = EL grid gap
    L  = drift lenght
    """
    def __init__(self,EP=3.0*kilovolt/(cm*bar),
                 dV=0.3*kilovolt/cm,
                 P=15*bar,
                 d=5*mm,
                 L=1.3*m,
                 Ws=76*eV,
                 Wi=24.8*eV,
                 Latt=6000*cm,
                 EPmin=0.6*kilovolt/(cm*bar)):

        TpcEL.__init__(self,EP,dV,P,d,L,Ws,Wi)

        units = kilovolt/(cm*bar)
        ep = EP/units
        epmin = EPmin/units
        self.EPmin=EPmin
        self.YPm = 140*epmin-116.
        self.Ngm = max(0,self.YPm*self.d/cm*self.P/bar) #photons/e
        self.TrkL= 30*cm*P/(10*bar)  # track length of a bb0nu events


    def BufferDistance(self):
        return self.Vc/(self.P*self.EPmin)

    def dEdx(self):
        """
        electrons per unit length
        """
        I = self.IonizationElectrons(Qbb_Xe136)
        return I/self.TrkL

    def dGdx(self):
        """
        photons per unit length
        """
        return self.dEdx()*self.OpticalGain()

    def __str__(self):
        kv_cm_b = kilovolt/(cm*bar)
        kv_cm =kilovolt/cm
        cm_bar = 1./(cm*bar)
        s= """
        E/P = %7.2f kV * cm^-1* bar^-1
        E/P min = %7.2f kV * cm^-1* bar^-1
        dV = drift voltage = %7.2f kV * cm^-1
        P  = pressure = %7.2f bar
        d  = EL grid gap = %7.2f mm
        L  = drift lenght =%7.2f m
        Grid voltage = %7.2f kV
        Cathode voltage = %7.2f kV
        Yield =  %7.2e photons/e
        Yield at EPmin =  %7.2e photons/e
        Buffer length = %7.1f cm
        dE/dx =  %7.2e electrons/cm
        dG/dx =  %7.2e photons/cm
    """%(self.EoverP()/kv_cm_b, self.EPmin/kv_cm_b, self.DriftVoltage()/kv_cm,
         self.Pressure()/bar,
         self.GridDistance()/mm, self.DriftLength()/m, self.GridVoltage()/kilovolt,
         self.CathodeVoltage()/kilovolt,self.OpticalGain(), self.Ngm, self.BufferDistance()/cm,
         self.dEdx()/cm,self.dGdx()/cm)

        s+="""
        Primary scintillation photons per MeV = %7.4e
        Primary ionization electrons per MeV = %7.4e
        """%(self.ScintillationPhotons(1*MeV),
             self.IonizationElectrons(1*MeV))

        return s



class ANGEL:
    """
    Abstraction of ANGEL design
    """
    def __init__(self,angelXe,angelCan,angelPMT,angelEL,angelSiPM):
        """
        We define ANGEL with 4 objects.
        angelPMT for the PMT readout
        angelXe for the xenon gas
        angelCan for the pressure vessel
        angelEL for the electrical and optical properties
        """
        self.PMT = angelPMT
        self.Xe = angelXe
        self.Can = angelCan
        self.EL = angelEL
        self.SiPM = angelSiPM

    def SolidAngle(self):
        R= self.Xe.FiducialRadius()
        L2=self.Can.Length()
        cteta = L2/sqrt(L2**2+R**2)
        return 0.5*(1-cteta)

    def ScintillationSolidAngle(self):
        return 10 * self.SolidAngle()

    def ELSolidAngle(self):
        return self.SolidAngle()

    def CathodePhotoEfficiency(self):
        return (self.PMT.TPBTransmittance() * self.PMT.CoveragePMT() *
                self.PMT.WindowT() * self.PMT.QE())

    def ScintillationReachingCathode(self,E):
        return self.EL.ScintillationPhotons(E) * self.ScintillationSolidAngle()

    def ScintillationPhotoelectrons(self,E):
        return self.ScintillationReachingCathode(E) * self.CathodePhotoEfficiency()

    def ELReachingCathode(self,E):
        return self.EL.ELPhotons(E) * self.ELSolidAngle()

    def ELPhotoelectrons(self,E):
        return self.ELReachingCathode(E) * self.ScintillationSolidAngle()

    def ELPhotoelectronsPerElectron(self):
        #window coated with TPB only 1/2 photon
        # print( """
        # Optical Gain = %7.2e
        # Solid Angle  = %7.2e
        # TPB Transmittance = %7.2e
        # Coverage  = %7.2e
        # Window Transmittance  = %7.2e
        # QE = %7.2e
        # """%(self.EL.OpticalGain(),self.ELSolidAngle(),self.PMT.TPBTransmittance(),
        #      self.PMT.CoveragePMT(),self.PMT.WindowTransmittance(), self.PMT.QE()))

        return self.EL.OpticalGain()*self.ELSolidAngle()*self.CathodePhotoEfficiency()

    def BackgroundRate(self,FWHM=1.0):
        """
        RF is the rejection factor
        TF is the topological rejection factor
        FWHM is the resolution FWHM
        """
        apmt = self.PMT.Activity() # c/year
        acan = self.Can.Activity()
        rf =  0.5*(Bi214_RF[1.0]+Tl208_RF[1.0])
        c_year = (apmt+acan)*rf*TF
        k=(Qbb_Xe136/keV)*FWHM/100.
        c_year_keV = c_year/k
        m=self.Xe.FiducialMass()/kg
        c_year_keV_kg = c_year_keV/m

        return c_year_keV_kg
    def __str__(self):

        s= """
        Optical Gain = %7.2e
        Solid Angle  = %7.2e
        ScintillationSolidAngle = %7.2e
        ELSolidAngle = %7.2e
        Scintillation Photons/MeV = %7.2e
        Scintillation Reaching Cathode/MeV = %7.2e
        Scintillation Photoelectrons/MeV = %7.2e
                   EL Photons/MeV = %7.2e
                   EL Reaching Cathode/MeV = %7.2e
                   EL Photoelectrons/MeV = %7.2e
                   EL pe/electron = %7.2e

        Background rate (c/(keV*kg*year) = %7.2e
        """%(self.EL.OpticalGain(),self.ELSolidAngle(),
             self.ScintillationSolidAngle(),
             self.ELSolidAngle(),
             self.EL.ScintillationPhotons(1*MeV),
             self.ScintillationReachingCathode(1*MeV),
             self.ScintillationPhotoelectrons(1*MeV),
             self.EL.ELPhotons(1*MeV),
             self.ELReachingCathode(1*MeV),
             self.ELPhotoelectrons(1*MeV),
             self.ELPhotoelectronsPerElectron(),
             self.BackgroundRate(FWHM=1.0)
             )

        return s

class NausicaAngel(ANGEL):
    """
    Like ANGEL but with Two SiPM planes
    """
    def __init__(self,xenon,can,EL,energyPlane,trackingPlane):
        super().__init__(xenon,can,energyPlane,EL,trackingPlane)
        self.EP = energyPlane

    def ScintillationPhotoelectrons(self,E):
        return self.ScintillationReachingCathode(E) * self.EP.PesPerPhotonReachingSensorPlane()

    def ELPhotoelectrons(self,E):
        return self.ELReachingCathode(E)* self.EP.PesPerPhotonReachingSensorPlane()

    def ELPhotoelectronsPerElectron(self):
        return self.EL.OpticalGain()* self.ELSolidAngle()* self.EP.PesPerPhotonReachingSensorPlane()


def TitanAngel():

    axe =AngelXe(R=53*cm,L=130*cm,Rb=4*cm,Cb=4*cm,Ab=1*cm,P=15)
    print("AngelXe")
    print(axe)

    print("AngelCan")

    acan = AngelCan(Rin=57*cm,Rout=58*cm,L=145*cm,Cp=0.8*cm,Ap=1.2*cm,
                 cplateType='Flat', aplateType='Flat',
                 FlangeRin=1140*mm/2., FlangeRout=1600*mm/2., FlangeT=9*cm, Mat='Ti')
    print(acan)

    r1141=PMT(name='R1141',gain=1e+6,D=3*2.5*cm, serFWHM=25*ns,filterRT=50*ohm*200*pF,ampliGain=10,
              QE=0.3,AU=3*becquerel*1e-3,ATh=2*becquerel*1e-3)


    print("R1141")
    print(r1141)

    print("AngelPMT")
    apmt =AngelPMT(r1141,axe,coverage=0.20,n1=1.5,n2=1.8)
    print(apmt)

    print("Angel EL")
    ael = AngelEL(EP=3.0*kilovolt/(cm*bar),
                 dV=0.5*kilovolt/cm,
                 P=15*bar,
                 d=5*mm,
                 L=1.3*m,
                 EPmin=1.5*kilovolt/(cm*bar))
    print(ael)


    sipm = SiPM(name='Hamamatsu',L=1.0*mm, QE=0.5,TPB=0.3)
    print("SiPM")
    print(sipm)

    asipm= AngelSiPM(sipm,ael,Surface=pi*(53*cm)**2,pitch=1.0*cm)

    print("ANGEL")
    angel = ANGEL(axe,acan,apmt,ael,asipm)
    print(angel)

def SteelAngel():

    axe =AngelXe()
    print("AngelXe")
    print(axe)

    print("AngelCan")

    acan = AngelCan()
    print(acan)

    r1141=PMT(name='R1141',gain=1e+6,D=3*2.5*cm, serFWHM=25*ns,filterRT=50*ohm*200*pF,ampliGain=10,
              QE=0.3,AU=2*becquerel*1e-3,ATh=2*becquerel*1e-3)


    print("R1141")
    print(r1141)

    print("AngelPMT")
    apmt =AngelPMT(r1141,axe,coverage=0.33,n1=1.5,n2=1.8)
    print(apmt)

    print("Angel EL")
    ael = AngelEL(EP=3.0*kilovolt/(cm*bar),
                 dV=0.3*kilovolt/cm,
                 P=15*bar,
                 d=5*mm,
                 L=1.3*m,
                 EPmin=1.5*kilovolt/(cm*bar))
    print(ael)


    sipm = SiPM(name='Hamamatsu',L=1.0*mm, QE=0.5,TPB=0.5)
    print("SiPM")
    print(sipm)

    asipm= AngelSiPM(sipm,ael,Surface=pi*(53*cm)**2,pitch=1.1*cm)

    print("NEXT100")
    angel = ANGEL(axe,acan,apmt,ael,asipm)
    print(angel)

def NEXT10():

    axe =AngelXe(R=30*cm,L=60*cm,Rb=10*cm,Cb=10*cm,Ab=10*cm,P=15)
    print("AngelXe")
    print(axe)

    print("AngelCan")

    acan = AngelCan(Rin=40*cm,Rout=41*cm,L=80*cm,Cp=0.8*cm,Ap=1.2*cm,
                 cplateType='Flat', aplateType='Flat',
                 FlangeRin=40*cm/2., FlangeRout=41*cm/2., FlangeT=5*cm, Mat='Ti')
    print(acan)

    r1141=PMT('r1141',3*2.5*cm,0.3,3*becquerel*1e-3,2*becquerel*1e-3)

    print("R1141")
    print(r1141)

    print("AngelPMT")
    apmt =AngelPMT(r1141,axe,coverage=0.20,n1=1.5,n2=1.8)
    print(apmt)

    print("Angel EL")
    ael = AngelEL(EP=3.5*kilovolt/(cm*bar),
                 dV=1.0*kilovolt/cm,
                 P=15*bar,
                 d=5*mm,
                 L=60*cm,
                 EPmin=0.5*kilovolt/(cm*bar))
    print(ael)


    sipm = SiPM(name='Hamamatsu',L=0.3*mm, QE=0.5,TPB=0.3)
    print("SiPM")
    print(sipm)


    asipm= AngelSiPM(sipm,ael,Surface=pi*(30)**2*cm**2,pitch=1.0*cm)
    print("ANGEL SiPM")
    print(asipm)

    print("ANGEL")
    angel = ANGEL(axe,acan,apmt,ael,asipm)
    print(angel)

def NEXT1():

    axe =AngelXe(R=16*cm,L=30*cm,Rb=10*cm,Cb=10*cm,Ab=10*cm,P=15)
    print("AngelXe")
    print(axe)

    print("AngelCan")

    acan = AngelCan(Rin=30*cm,Rout=31*cm,L=60*cm,Cp=0.8*cm,Ap=1.2*cm,
                 cplateType='Flat', aplateType='Flat',
                 FlangeRin=40*cm/2., FlangeRout=41*cm/2., FlangeT=5*cm, Mat='Ti')
    print(acan)

    r1141=PMT('r1141',3*2.5*cm,0.3,3*becquerel*1e-3,2*becquerel*1e-3)

    print("R1141")
    print(r1141)

    print("AngelPMT")
    apmt =AngelPMT(r1141,axe,coverage=0.20,n1=1.5,n2=1.8)
    print(apmt)

    print("Angel EL")
    ael = AngelEL(EP=3.5*kilovolt/(cm*bar),
                 dV=1.0*kilovolt/cm,
                 P=15*bar,
                 d=5*mm,
                 L=60*cm,
                 EPmin=0.5*kilovolt/(cm*bar))
    print(ael)


    sipm = SiPM(name='Hamamatsu',L=0.3*mm, QE=0.5,TPB=0.3)
    print("SiPM")
    print(sipm)


    asipm= AngelSiPM(sipm,ael,Surface=pi*(30)**2*cm**2,pitch=1.0*cm)
    print("ANGEL SiPM")
    print(asipm)

    print("ANGEL")
    angel = ANGEL(axe,acan,apmt,ael,asipm)
    print(angel)

def CuAngel():


    print("Copper Can")

    acan = AngelCan(Rin=57*cm,Rout=60*cm,L=145*cm,Cp=5.0*cm,Ap=5.0*cm,
                 cplateType='Flat', aplateType='Flat',
                 FlangeRin=1140*mm/2., FlangeRout=1600*mm/2., FlangeT=14*cm, Mat='Cu')

    print(acan)

def TiTube():
    Lf=11*cm
    Rf=5*cm

    tubeShape = Cylinder(Rf,Lf)
    tubeMat = Material('GXe',15)
    tubePV = PhysicalVolume("TubeVolume",tubeMat,tubeShape)

    print("""
    tube  Radius = %7.2f cm
    tube  Length = %7.2f cm
    """%(Rf/cm,Lf/cm))


    print("""
    tube Volume = %7.3e m3
    tube mass = %7.3e kg
    """%(tubePV.Volume()/m3,tubePV.Mass()/kg))

def NAUSICAA():
    diameter = 250 * cm
    length   = 250 * cm
    l_esipm  = 10 * mm
    l_tsipm  = 1 * mm
    p_esipm  = 11 * mm
    p_tsipm  = 11    * mm

    # diameter = 50 * cm
    # length   = 50 * cm
    # l_esipm  = 10 * mm
    # l_tsipm  = 1 * mm
    # p_esipm  = 10 * mm
    # p_tsipm  = 10 * mm

    axe =AngelXe(R=diameter/2, L=length, Rb=1*cm, Cb=10*cm, Ab=1*cm, P=20)
    print("NAUSICAA")
    print(axe)

    print("NAUSICAA EL")
    ael = AngelEL(EP=2.5*kilovolt/(cm*bar),
                  dV=1.0*kilovolt/cm,
                  P=20*bar,
                  d=5*mm,
                  L=length,
                  EPmin=0.5*kilovolt/(cm*bar))
    print(ael)


    esipm = SiPM(name='Hamamatsu',L=l_esipm, QE=0.5,TPB=0.5)
    print("Energy SiPM")
    print(esipm)

    tsipm = SiPM(name='Hamamatsu',L=l_tsipm, QE=0.5,TPB=0.5)
    print("Tracking SiPM")
    print(tsipm)


    n_esipm= AngelSiPM(esipm, ael, Surface=pi*(diameter/2)**2, pitch=p_esipm)
    print("NAUSICAA Energy SiPM")
    print(n_esipm)

    n_tsipm= AngelSiPM(tsipm, ael, Surface=pi*(diameter/2)**2, pitch=p_tsipm)
    print("NAUSICAA Tracking SiPM")
    print(n_tsipm)




def DEMOPP():

    axe =AngelXe(R=(40/2)*cm,L=52*cm,Rb=1*cm,Cb=1*cm,Ab=1*cm,P=15)
    print("DEMO++")
    print(axe)

    print("Can")

    acan = AngelCan(Rin=64*cm,Rout=65*cm,L=100*cm,Cp=0.8*cm,Ap=1.2*cm,
                 cplateType='Flat', aplateType='Flat',
                 FlangeRin=64*cm/2., FlangeRout=65*cm/2., FlangeT=5*cm, Mat='Fe')
    print(acan)

    r1141=PMT('r1141',3*2.5*cm,0.3,3*becquerel*1e-3,2*becquerel*1e-3)

    print("R1141")
    print(r1141)

    print("PMT")
    apmt =AngelPMT(r1141,axe,coverage=0.20,n1=1.5,n2=1.8)
    print(apmt)

    print("EL")
    ael = AngelEL(EP=2.5*kilovolt/(cm*bar),
                 dV=0.3*kilovolt/cm,
                 P=15*bar,
                 d=5*mm,
                 L=30*cm,
                 EPmin=0.5*kilovolt/(cm*bar))
    print(ael)


    sipm = SiPM(name='Hamamatsu',L=0.3*mm, QE=0.5,TPB=0.3)
    print("SiPM")
    print(sipm)


    asipm= AngelSiPM(sipm,ael,Surface=pi*(30)**2*cm**2,pitch=1.0*cm)
    print( "DEMO++ SiPM")
    print( asipm)

    print( "DEMO++")
    angel = ANGEL(axe,acan,apmt,ael,asipm)
    print( angel)


if __name__ == '__main__':

    print( "Next 100")
    SteelAngel()

    print( "****NEXT1T++****")
    NAUSICAA()


##     print( "TiTube"
##     TiTube()
##     print( "Titan Angel"
##     TitanAngel()
   #CuAngel()
