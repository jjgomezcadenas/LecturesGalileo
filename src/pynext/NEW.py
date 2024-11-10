from . Next100 import *

def NeW2014():

    axe =AngelXe(R=240*mm,L=510*mm,Rb=1*mm,Cb=70*mm,Ab=1*mm,P=15)
    print "NEW 2014: gas"
    print axe

    print "NeW 2014: can"

    acan = AngelCan(640*mm/2.,Rout=664*mm/2.,L=120*mm,Cp=12.5*mm,Ap=12.5*mm,
                    cplateType='Flat', aplateType='Flat',
                    FlangeRin=640*mm/2., FlangeRout=664*mm/2., FlangeT=50*mm, Mat='Fe')

    print acan

    r1141=PMT('r1141',3*2.5*cm,0.3,1*becquerel*1e-3,1*becquerel*1e-3)

    print "R1141"
    print r1141

    print "NEW 2014: PMTs"
    apmt =AngelPMT(r1141,axe,coverage=0.20,n1=1.5,n2=1.5)
    print apmt

    print "Angel EL"
    ael = AngelEL(EP=3.5*kilovolt/(cm*bar),
                 dV=1.0*kilovolt/cm,
                 P=15*bar,
                 d=5*mm,
                 L=510*mm,
                 EPmin=0.5*kilovolt/(cm*bar))
    print ael


    sipm = SiPM(name='SENSL',L=1.0*mm, QE=0.5,TPB=0.5)
    print "SENSL: SiPM"
    print sipm


    asipm= AngelSiPM(sipm,ael,Surface=pi*(240*mm)**2,pitch=1.0*cm)
    print "NEW 2014: SiPM"
    print asipm

    print "NeW 2014"
    new2014 = ANGEL(axe,acan,apmt,ael,asipm)
    print new2014

if __name__ == '__main__':


    print "****NEW 2014++****"
    NeW2014()
