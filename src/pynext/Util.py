"""
Utilities
"""

from math import *
import scipy.constants as cnt
import system_of_units as units

def wait():
    raw_input("Press a key...")


def drange(start, stop, step):
    r = start
    while r < stop:
        yield r
        r += step


if __name__ == '__main__':
    print('speed of light, c = {} m/s'.format(cnt.c))
    for c in cnt.physical_constants:
        print(c)
