"""
mathematical functions
"""
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

def attenuation_factor(mu, z):
    att = SelfAtt(mu, z)
    tf, _ = quad(att.f, -np.pi/2 + 0.0001, np.pi/2 -  0.0001)
    return  (tf / (2 * np.pi))
