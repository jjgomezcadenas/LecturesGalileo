"""
Shapes
Defines a material shape
"""
from . system_of_units import *
from math import pi
from abc import ABC, abstractmethod

class Shape(ABC):

    """An abstract class representing a geometrical shape.
    In the most general case, the shape is assumed to have a thickness
    (e.g., a sphericall shell or a cylinder shell).

    The interface is defined by 6 methods. To describe them imagine that we want to describe
    a pressure vessel (PV) with gas inside. We can use three shapes.
    A CylinderShell(Rin, Rout, L, ts) to describe the
    PV body; A Disk(R, td) to describe the two end cups. And a Cylinder (Rin, L) to describe
    the volume of gas enclosed by the PV.

    CylinderShell (describes the PV)
    1. inner_volume(): volume enclosed by the CylinderShell
                 = 2 * pi * Rin * L.
    2. shell_volume(): Volume of the shell itself
                 =  2 * pi * (Rout - Rin) * L
    3. inner_surface(): surface in contact with the gas,
                 =  2 * pi * Rin * L
    4. outer_surface(): surface in contact with the air outside the vessel,
                 = 2 * pi * Rout * L
    5. thickness_surface(): surface in contact with the end-cups
                 = pi * (Rout**2 - Rin**2)
    6. thickness() would correspond to the thickness of steel
                 = ts

    Cylinder (describes the volume of gas)
    1. inner_volume()      = 2 * pi * Rin * L
    2. shell_volume()      = 0
    3. inner_surface()     = 2 * pi * Rin * L
    4. outer_surface()     = 2 * pi * Rin * L
    5. thickness_surface() = 0
    6. thickness()         = 0

    Disk (describes the end-cups)
    1. inner_volume()      = pi * R**2 * td
    2. shell_volume()      = pi * R**2 * td
    3. inner_surface()     = pi * R**2
    4. outer_surface()     = pi * R**2
    5. thickness_surface() = 0
    6. thickness()         = 2 * pi * R * td

    Once the interface is defined, the methods volume() and surface()
    as well as the properties V and S come for free.
    """


    @abstractmethod
    def thickness(self):
        pass

    @abstractmethod
    def radius(self):
        pass

    @abstractmethod
    def volume(self):
        pass

    @abstractmethod
    def surface(self):
        pass

    @abstractmethod
    def inner_volume(self):
        pass

    @abstractmethod
    def shell_volume(self):
        pass

    @abstractmethod
    def inner_surface(self):
        pass

    @abstractmethod
    def outer_surface(self):
        pass

    @abstractmethod
    def thickness_surface(self):
        pass

    @property
    def V(self):
        return self.volume()

    @property
    def S(self):
        return self.surface()

    def __str__(self):

        s= """\n
        radius            = {:7.2e} mm
        thickness         = {:7.2e} mm
        inner_volume      = {:7.2e} m3
        shell_volume      = {:7.2e} m3
        inner_surface     = {:7.2e} m2
        outer_surface     = {:7.2e} m2
        thickness_surface = {:7.2e} m2

        volume            = {:7.2e} m3
        surface           = {:7.2e} m2
        """.format(
        self.radius() / mm,
        self.thickness() / mm,
        self.inner_volume() / m3,
        self.shell_volume() / m3,
        self.inner_surface() / m2,
        self.outer_surface() / m2,
        self.thickness_surface() / m2,

        self.volume() / m3,
        self.surface() / m2)

        return s


class Sphere(Shape):

    def __init__(self, R):

        self.R = R

    def inner_volume(self):
        return (4/3) * pi * self.R**3

    def shell_volume(self):
        return 0

    def inner_surface(self):
        return 4 * pi * self.R**2

    def outer_surface(self):
        return 4 * pi * self.R**2

    def thickness_surface(self):
        return 0

    def thickness(self):
        return 0

    def radius(self):
        return self.R

    def volume(self):
        return self.inner_volume()

    def surface(self):
        return self.inner_surface()


    def __str__(self):

        s= """
        Sphere(R = %7.2e m)
        """%(self.R / m)
        s2 = super().__str__()
        return s + s2

    __repr__ = __str__


class SphereShell(Shape):

    def __init__(self, Rin , Rout):

        self.Rin  = Rin
        self.Rout = Rout

    def inner_volume(self):
        return (4/3) * pi * self.Rin**3

    def shell_volume(self):
        return (4/3) * pi * (self.Rout**3 - self.Rin**3)

    def inner_surface(self):
        return 4 * pi * self.Rin**2

    def outer_surface(self):
        return 4 * pi * self.Rout**2

    def thickness_surface(self):
        return 0

    def thickness(self):
        return self.Rout - self.Rin

    def radius(self):
        return self.Rin

    def volume(self):
        return self.shell_volume()

    def surface(self):
        return self.outer_surface()


    def __str__(self):

        s= """
        SphereShell(Rin = %7.2e m , Rout= %7.2e m)
        """%(self.Rin / m, self.Rout / m)
        s2 = super().__str__()
        return s + s2

    __repr__ = __str__


class Cylinder(Shape):

    def __init__(self, R, L):

        self.R = R
        self.L = L

    def inner_volume(self):
        return pi * self.R**2 * self.L

    def shell_volume(self):
        return 0

    def inner_surface(self):
        return 2 * pi * self.R * self.L

    def outer_surface(self):
        return 2 * pi * self.R * self.L

    def thickness_surface(self):
        return 0

    def radius(self):
        return self.R

    def thickness(self):
        return 0

    def volume(self):
        return self.inner_volume()

    def surface(self):
        return self.inner_surface()


    def __str__(self):

        s= """
        Cylinder(R = %7.2e m, L = %7.2e m)
        """%(self.R / m, self.L / m)
        s2 = super().__str__()
        return s + s2

    __repr__ = __str__


class CylinderShell(Shape):

    def __init__(self, Rin, Rout, L):


        self.Rin  = Rin
        self.Rout = Rout
        self.L    = L

    def inner_volume(self):
        return pi * self.Rin**2 * self.L

    def shell_volume(self):
        return pi * (self.Rout**2 - self.Rin**2) * self.L

    def inner_surface(self):
        return 2 * pi * self.Rin * self.L

    def outer_surface(self):
        return 2 * pi * self.Rout * self.L

    def thickness_surface(self):
        return pi * (self.Rout**2 - self.Rin**2)

    def thickness(self):
        return self.Rout - self.Rin

    def radius(self):
        return self.Rin

    def volume(self):
        return self.shell_volume()

    def surface(self):
        return self.inner_surface()

    def __str__(self):

        s= """
        CylinderShell(Rin = %7.2e m, Rout= %7.2e m , L = %7.2e m)
        """%(self.Rin / m, self.Rout / m, self.L / m)
        s2 = super().__str__()
        return s + s2
        return s

    __repr__ = __str__


class Disk(Shape):

    def __init__(self, R, t):
        self.R = R
        self.t = t

    def inner_volume(self):
        return pi * self.R**2 * self.t

    def shell_volume(self):
        return 0

    def inner_surface(self):
        return pi * self.R**2

    def outer_surface(self):
        return pi * self.R**2

    def thickness_surface(self):
        return 2 * pi *  self.R * self.t

    def thickness(self):
        return self.t

    def radius(self):
        return self.Rin

    def volume(self):
        return self.inner_volume()

    def surface(self):
        return self.inner_surface()

    def __str__(self):

        s= """
        Disk(R = %7.2e m, t = %7.2e m)
        """%(self.R / m, self.t / m)
        s2 = super().__str__()
        return s + s2

    __repr__ = __str__


class Brick(Shape):

    def __init__(self, width, heigth, length):

        self.width  = width
        self.heigth = heigth
        self.length = length

    def inner_volume(self):
        return self.width * self.heigth * self.length

    def shell_volume(self):
        return 0

    def inner_surface(self):
        return 2*(self.width * self.heigth + self.width * self.length + self.length * self.heigth)

    def outer_surface(self):
        return 2*(self.width * self.heigth + self.width * self.length + self.length * self.heigth)

    def thickness_surface(self):
        return 0

    def thickness(self):
        return 0

    def radius(self):
        return 0

    def volume(self):
        return self.inner_volume()

    def surface(self):
        return self.inner_surface()

    def __str__(self):

        s= """
        Brick(width= %7.2e m , heigth= %7.2e m , lengthR= %7.2e m)
        """%(self.width / m, self.heigth / m, self.length / m)
        s2 = super().__str__()
        return s + s2

    __repr__ = __str__

FlatPlate = Disk

#
# class Wall:
#     """
#     A wall made of identical bricks
#     """
#     def __init__(self, brick):
#         self.b = brick
#
#     def number_of_bricks(self, w, h, t):
#         nw = w * 1 / self.b.width
#         nh = h * 1 / self.b.heigth
#         nt = t * 1 / self.b.length
#
#         return nw * nh * nt
