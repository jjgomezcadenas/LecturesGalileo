
#     MATERIALS=['H2O','Fe','Ti','Cu','Pb','LXe','GXe','Cu10','H2OX',
#     'Kevlar','Tensylon','Vectran','Graphene','LSC','LSCTh','LSCU','Poly','Poly2','Peek','PTFE']
#     MU = {'H2O': 0.045 * cm2/g,
#           'H2OX': 0.045 * cm2/g,
#           'Fe' : 0.039 * cm2/g,
#           'Ti' : 0.038 * cm2/g,
#           'Cu' : 0.039 * cm2/g,
#           'Cu10' : 0.039 * cm2/g,
#           'Pb' : 0.044 * cm2/g,
#           'LXe' : 3.8e-2 * cm2/g,
#           'GXe' : 3.8e-2 * cm2/g,
#           'Kevlar' : 0.045  * cm2/g,
#           'Tensylon' : 0.045  * cm2/g,
#           'Vectran' : 0.045  * cm2/g,
#           'Graphene' : 0.045  * cm2/g,
#           'LSC' : 1e-6  * cm2/g,
#           'LSCTh' : 1e-6  * cm2/g,
#           'LSCU' : 1e-6  * cm2/g,
#           'Poly' : 1e-6  * cm2/g,
#           'Poly2' : 1e-6  * cm2/g,
#           'Peek' : 1e-6  * cm2/g,
#           'PTFE' : 1e-6  * cm2/g,
#           }
#
#
#     RHO = {'H2O': 1.00 * g/cm3,
#            'H2OX': 1.00 * g/cm3,
#            'Fe' : 7.87 * g/cm3,
#            'Ti' : 4.54 * g/cm3,
#            'Cu' : 8.96 * g/cm3,
#            'Cu10' : 8.96 * g/cm3,
#            'Pb' : 11.33 * g/cm3,
#            'LXe' :2.98 * g/cm3,
#            'GXe' :5.761 *kg/m3,
#            'Kevlar' : 1440 * kg/m3,
#            'Vectran' : 1.4 * g/cm3,
#            'Tensylon' : 0.97 * g/cm3,
#            'Graphene' : 2.0 * g/cm3,
#            'LSC' : 1.0 * g/cm3,
#            'LSCTh' : 1.0 * g/cm3,
#            'LSCU' : 1.0 * g/cm3,
#            'Poly' : 1.0 * g/cm3,
#            'Poly2' : 1.0 * g/cm3,
#            'Peek' : 1.3 * g/cm3,
#            'PTFE' : 2.0 * g/cm3,
#            }
#
#     AS = {'Fe' : {'BI214':1.9e-3 * becquerel/kg,'TL208':1.e-3/3. * becquerel/kg},
#           'Fe-10': {'BI214':1.e-2 * becquerel/kg,'TL208':1.e-2/3. * becquerel/kg},
#           'Ti': {'BI214':0.93e-3 * becquerel/kg,'TL208':0.22e-3 * becquerel/kg},
#           'Ti-X' : {'BI214':1.e-4 * becquerel/kg,'TL208':1.e-4/3 * becquerel/kg},
#           'Cu10' : {'BI214':12e-6 * becquerel/kg,'TL208':4.e-6/3 * becquerel/kg},
#           'H2O'   : {'BI214':1.e-5 * becquerel/kg,'TL208':1.e-5/3 * becquerel/kg},
#           'H2OX'   : {'BI214':1.e-3 * becquerel/kg,'TL208':1.e-3/3 * becquerel/kg},
#           'LXe'   : {'BI214':1.e-12 * becquerel/kg,'TL208':1.e-12/3 * becquerel/kg},
#           'GXe'   : {'BI214':1.e-12 * becquerel/kg,'TL208':1.e-12/3 * becquerel/kg},
#           'Pb'   : {'BI214':370.0e-6 * becquerel/kg,'TL208':73e-6/3 * becquerel/kg},
#           'Graphene'   : {'BI214':1e-18 * becquerel/kg,'TL208':1.e-18 * becquerel/kg},
#           'Kevlar'   : {'BI214':0.17 * becquerel/kg,'TL208':0.17* becquerel/kg},
#           'Poly'   : {'BI214':2.9e-3 * becquerel/kg,'TL208':1.e-3 * becquerel/kg},
#           'Poly2'   : {'BI214':2.3e-4 * becquerel/kg,'TL208':1.4e-4/3 * becquerel/kg},
#           'Tensylon'   : {'BI214':(0.05/81) * becquerel/kg,'TL208':(0.05/81/3) * becquerel/kg},
#           'Vectran'   : {'BI214':(0.1/81) * becquerel/kg,'TL208':(0.1/81/3) * becquerel/kg},
#           'Peek'   : {'BI214': 36e-3 * becquerel/kg,'TL208':4.3e-3/3 * becquerel/kg},
#           'PTFE'   : {'BI214': 0.025e-3 * becquerel/kg,'TL208':0.031e-3/3 * becquerel/kg},
#           'LSC'   : {'BI214':(3.5e-18/81e-9) * becquerel/kg,'TL208':(5.2e-17/81e-9) * becquerel/kg},
#           }
#
##
# class RMaterial(Material):
#     """
#     RMaterial is a material in which we add properties relevants for its resistance
#     """
#
#     def __init__(self,name,isotope='BI214',S=3200*MPa,p=1):
#
#         Material.__init__(self,name,isotope,p)
#         self.S=S
#
#     def TensileStrength(self):
#         return self.S
#
#
# COPPER=RMaterial('Cu10',isotope='BI214',S=220*MPa)
# TITANIUM=RMaterial('Ti',S=990*MPa)
# STEEL=RMaterial('Fe',S=1860*MPa)
# VECTRAN=RMaterial('Vectran',S=3200*MPa)
# KEVLAR=RMaterial('Kevlar',S=3200*MPa)
# TENSYLON=RMaterial('Tensylon',S=32*MPa)
# LEAD=RMaterial('Pb',S=2*MPa)  #including the effect of creep
# GRAPHENE=RMaterial('Graphene',S=10**5*MPa)
# LXE=RMaterial('LXe',S=0)

#
# class Gas(Material):
#     """
#     Gas is a material in which we add properties relevants for gases
#     vd is the drift velocity
#     mu is the mobility of the electrons
#     mup is the mobility of the positiveions
#     """
#
#     def __init__(self,name,isotope='BI214',p=1,vd=1*mm/1*mus,mue=3e+3*cm2/(volt*second),
#                  mup=4e-3*cm2/(volt*second)):
#
#         Material.__init__(self,name,isotope,p)
#         self.vd=vd
#         self.mue = mue
#         self.mup=mup
#
#     def IonMobility(self):
#         return self.mup
#
#     def ElectronMobility(self):
#         return self.mue
#
#     def DriftVelocity(self):
#         return self.vd
#
# gLXE=Gas('LXe',p=1,vd=5e+5*cm/second,mue=3e+3*cm2/(volt*second),mup=4e-3*cm2/(volt*second))
#
# def vessel(t,R,name):
#
#     vessel = Material(name)
#     S=4*pi*R**2
#     V=(4./3.)*pi*(R+t)**3-(4./3.)*pi*(R)**3
#     M=vessel.rho*V
#     f0=M*vessel.A0
#     N0 = f0*year
#     flux = vessel.A(t)
#
#
#     print("""
#     thickness = %7.2f cm
#     Volume = %7.2f m3
#     Mass = %7.2f kg
#     Activity = %7.2e Bq
#     flux emitted by container (Bq) = %7.1e
#     flux emitted by container (Bq) corrected by shelf-shielding = %7.1e
#     """%(t/cm,V/m3,M/kg,f0/becquerel,N0,flux*S/(1./year)))
#
#
# if __name__ == '__main__':
#     vessel(1*cm,50*cm,'Cu10')
