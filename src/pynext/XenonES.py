from . system_of_units import *
import numpy as np
from math import pi, sqrt, exp, log
import sys
class XenonES:
    """Defines Xenon equation of state.
    https://www.nist.gov/sites/default/files/documents/srd/jpcrd470.pdf

    """
    def __init__(self):
        np.set_printoptions(precision=7)
        np.set_printoptions(suppress=True)
        f = np.array([
             [0.48176287000, -1.17810687000, -0.2405126300, -0.075090970300, -0.21382659500,
              -0.00701074245, -3.56560121],
             [-0.22996818000, 0.61145818500, -0.3576334170, -0.039464195900, 0.85313703300,
              0.03926238030,  8.73657614],
             [0.67678089000, -1.59127084000,  1.3854592100, -0.037603850100, -1.77422585000,
              -0.11095677400, -11.83502506],
             [-0.82687853, 2.18816951999, -2.2662842000,  0.404244668000, 2.11856977,
              0.175796531,  9.56797096],
             [0.59684665, -1.66416249,  1.89132269, -0.462622543, -1.53874755,
              -0.158135317, -4.61464670],
             [-0.26046265, 0.741574506, -0.896367224,  0.24081461, 0.693737872,
              0.07524343650,  1.22917766],
             [0.068249584, -0.192234256, 0.244259366, -0.06689881470, -0.190153915,
              -0.0146753994, -0.13948161],
             [-0.0098938026, 0.0268815303, -0.0356830831,  0.00967590594, 0.0290866982,
              0.,     0.],
             [0.00060938476    , -0.00156464723, 0.00216704074,  -0.000575094683
              , -0.00190789666, 0.,     0.],
             ])
        self.F   = np.asmatrix(f.transpose())
        self.rho = np.array([0., 2.970,0.405,2.325,1.113,0.611,1.485,1.101])
        self.Q   = np.array([-1., 0.,1.,3.,21.,24.83,250.,21.00])
        self.E   = np.array([1.5, 1.94,1.46,1.02,0.98, 1.173,0,0])
        self.Tr = 289.7 # Ref temperature in Kelvin
        self.Rr = 1000. # Ref density in kg/m3
        self.Rm  = 8.31 # Molar gas constant in J/(mol K)
        self.M  = 131.29 # molar xenon mass in g/mol
        self.R  = 1e+3 * self.Rm / self.M # specific gas constant in J((kg k))
        self.T0 = 273.15 # 0 in K
        self.pascal_to_atm = 9.9E-6


        self.MXe = self.M
        self.MHe = 4 # g/mol
        #self.R = (self.Rm / self.M)
        #self.R = (self.Rm / self.M) * J/(K * gram)
        self.RXe = (self.Rm / self.MXe) * J/(K * gram)
        self.RHe = (self.Rm / self.MHe) * J/(K * gram)

        self.k     = 1.38064852 * 1E-23 * J/K
        self.N_A   = 6.02214129E+23 * (1./mol)
        self.ek    = 274 * K # epsilon/k
        self.sigma = 0.3885 * nm

        self.rho_2020 = 124.3 * kg/m3
        self.rho_3020 = 203.35 * kg/m3
        self.rho_1520 = 89.9 * kg/m3
        self.rho_1020 = 58 * kg/m3
        self.Ls       = 0.06    # Lambda*


    def P(self, T, RHO, perfect=False, temp='C'):
        """Computes pressure as a function of density and temperature for xenon
        T in Celsius
        rho in kg/m3
        P in atm
        """
        t = T  / self.Tr
        if temp == 'C':
            t = (T + self.T0) / self.Tr    # relative temperature
        rho = RHO / (kg/m3)
        r = rho / self.Rr  # relative density

        if perfect == False:
            p = r * (t + self.FF_(r, t))  # relative pressure
        else:
            p = r * t

        #return p * self.R * self.Tr * self.Rr * self.pascal_to_atm
        return p * self.R * self.Tr * self.Rr * pascal

    def P2(self, T, RHO, perfect=False, temp='C', gas='Xe'):
        """Pressure a a function of rho and T, using B --secon virial coefficient--"""
        if temp == 'C':
            t = self.T0 + T # kelvin

        b = 0
        if perfect == False:
            b = self.B(t) / (meter3/mol) * RHO /(g/m3)

        if gas == 'Xe':
            return RHO * self.RXe * t * (1  + b / self.MXe)
        elif gas == 'He':
            return RHO * self.RHe * t * (1  + b / self.MHe)
        else:
            print('gas not defined. Exiting')
            sys.exit(0)


    def B(self,T):
        Ts = T/self.ek
        return (2/3) * pi * self.N_A * self.sigma**3 * self.Bs_(Ts)

    def Bs_(self, T):
        return (self.B0s_(T) + self.Ls**2 * self.B1s_(T) + self.Ls**4 * self.B2s_(T) +
                self.Ls**6 * self.B3s_(T)  + self.Ls**3 * self.Bps_(T))

    def B0s_(self, T):
        if T < 1.1:
            return - sqrt(T) * exp(1/T) * (1.18623 + 1.00824 * T + 4.25571 * T**2 -
                                           18.6033 * T**3 + 20.4732 * T**4 - 8.71903 * T**5 +
                                           1.14829 * T**6)

        elif 1.1 < T < 10:
            return - sqrt(T) * exp(1/T) * (0.74685 - 1.0384 * log(T) + 0.31634 * (log(T))**2 -
                                           0.02096 * (log(T))**3 - 0.01498* (log(T))**4)
        else:
            return 0. # can be programmed if needed

    def B1s_(self, T):
        if T < 1.1:
            return T**(-3/2) * exp(1/T) * (0.158192 - 0.037145 * T - 0.0103125 * T**2 +
                                           0.0701852 * T**3 - 0.0328399 * T**4)

        elif 1.1 < T < 10:
            return T**(-3/2)  * exp(1/T) * (0.148 + 0.0241 * log(T) - 0.0123 * (log(T))**2 +
                                            0.0096 * (log(T))**3 - 0.0014* (log(T))**4)
        else:
            return 0. # can be programmed if needed

    def B2s_(self, T):
        if T < 10:
            return T**(-7/2) * exp(1/T) * (0.0152 + 0.0126 * T + 0.0001 * T**2)
        else:
            return 0. # can be programmed if needed

    def B3s_(self, T):
        if T < 1.1:
            return T**(-11/2) * exp(1/T) * (0.001 + 0.006 * T - 0.0197 * T**3 +
                                            0.032 * T**4 - 0.0128 * T**5)
        elif 1.1 < T < 10:
            return T**(-11/2)  * exp(1/T) * (0.0051 - 0.0112 * T - 0.0021 * T**2)
        else:
            return 0. # can be programmed if needed

    def Bps_(self, T):
        return 0

    def psi_(self,i, r):

        def sumij (r):
            sumj = 0
            for k in range(9):
                j = k + 1
                sumj+= j * self.F[i,k] * r**j
            return sumj

        if i < 4 :
            return sumij(r)
        elif i == 4:
            return sumij(self.rho[2])
        elif i == 5:
            return sumij(self.rho[5])
        else:
            print('Error: i = {} out of bounds'.format(i))


    def Xi_(self,i, t):
        if i < 4 :
            return t** (-self.Q[i])
        elif i == 4:
            return t** (-self.Q[4])
        #elif i == 5:
        #    if self.E[4] < t < self.E[5]:
        #        return self.E[2] * t **(-self.Q[6]) + self.E[3] * t**(-self.Q[7])
        #    else:
        #        return 0
        else:
            return 0

    def FF_(self, r, t):
        f = 0
        for i in range(6):
            f += self.psi_(i, r) * self.Xi_(i, t)
        return f


    def __str__(self):
        return '< F = {}\n rho = {}\n Q = {}\n E = {}>'.format(self.F, self.rho, self.Q, self.E)

    __repr__ =     __str__


class XeHe:
    def __init__(self, m, V, tc, xHe, xXe):
        self.R = 8.31441 * J/(K * mol)
        self.T0 = 273.15 # 0 in K
        self.MXe = 131.29 * g/mol
        self.MHe = 4 * g/mol
        self.T = (self.T0 + tc) * K
        self.m = m
        self.V = V
        self.n0Xe = m  / self.MXe  # number of Xe mols for pure xenon.
        self.rho = m / V
        self.P0 = (self.n0Xe * self.R * self.T) / self.V
        self.nXe = self.n0Xe * xXe
        self.nHe = self.n0Xe * xHe
        self.PXe = (self.nXe * self.R * self.T) / self.V
        self.PHe = (self.nHe * self.R * self.T) / self.V
        self.P = self.PXe + self.PHe
        self.xXe = xXe
        self.xHe = xHe

    def tK(self, tc):
        return (273.15 + tc) * K

    def pXe(self, tc):
        return (self.nXe * self.R * self.tK(tc)) / self.V

    def pHe(self, tc):
        return (self.nHe * self.R * self.tK(tc)) / self.V

    def p(self, tc):
        return self.pXe(tc) + self.pHe(tc)

    def pBXe(self, tc):
        return self.pXe(tc) * (1 + self.B(self.tK(tc), self.xHe, self.xXe) * self.nXe/self.V)

    def pBHe(self, tc):
        return self.pHe(tc)  * (1 + self.B(self.tK(tc), self.xHe, self.xXe) * self.nHe/self.V)

    def pB(self, tc):
        return self.pBXe(tc) + self.pBHe(tc)


    def B(self, T, zHe, zXe):
        if zXe == 1:
            return self.BXe(T)
        elif zHe == 1:
             return self.BHe(T)
        else:
            return self.BXeHe(T, zHe, zXe)

    def BXe(self,T):
        b0 = 5.797E-5  * (meter**3/mol)
        L = 1.6940
        ekb = 181.8     #K
        D = exp(ekb/T) - 1
        return b0 * ( 1 - (L**3 -1) * D )

    def BHe(self,T):
        return 11.885370 + 6.6508722E-3 * T - 3.305894E-5 * T**2 + 3.1354694E-8 * T**3

    def BXeHe(self, T, x1=0.15, x2=0.85):
        def B12(T):
            return 2.756941E-4 + (4.617880) / T - 1.741538 / T**3 + 1.364453E+2 / T**3
        return x1**2 * self.BHe(T) + x2**2 * self.BXe(T) + 2 * x1 * x2 * B12(T)


    def __str__(self):
        s = """
            m    (XHe = 1)    = {} kg
            V                 = {} m3
            n0Xe (XHe = 1)    = {} mol
            rho  (XHe = 1)    = {} g/cm3
            P0   (XHe = 1)    = {} atm
            nXe  (XHe = {})   = {} mol
            nHe  (XHe = {})   = {} mol
            PXe               = {} atm
            PHe               = {} atm
            P                 = {} atm
            """.format(self.m / kg, self.V / m3, self.n0Xe/mol, self.rho / (g/cm3), self.P0 / atmosphere,
                       self.nXe/mol, self.xXe, self.nHe/mol, self.xHe,
                       self.PXe/ atmosphere, self.PHe/ atmosphere, self.P/ atmosphere)

        return s

    __repr__ =     __str__
