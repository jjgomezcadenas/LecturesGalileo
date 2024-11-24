from dataclasses import dataclass, field
import abc
import numpy as np
import pandas as pd
from scipy.linalg import norm
from scipy.special import erfc
from  . system_of_units import *

import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import axes3d


from typing      import Tuple
from typing      import Dict
from typing      import List
from typing      import TypeVar
from typing      import Optional
from enum        import Enum

Number = TypeVar('Number', None, int, float)
Str   = TypeVar('Str', None, str)
Range = TypeVar('Range', None, Tuple[float, float])
Array = TypeVar('Array', List, np.array)
Int = TypeVar('Int', None, int)

class Verbosity(Enum):
    mute     = 0
    concise  = 1
    chat     = 2
    verbose  = 3
    vverbose = 4


def vprint(msg, verbosity, level=Verbosity.mute):
    if verbosity.value <= level.value and level.value >0:
        print(msg)
        

def vpblock(msgs, verbosity, level=Verbosity.mute):
    for msg in msgs:
        vprint(msg, verbosity, level)


@dataclass
class Shape(abc.ABC):
    @property
    def area(self)->float:
        pass
    @property
    def volume(self)->float:
        pass


@dataclass
class Cylinder(Shape):
    r   : float
    zmin: float
    zmax: float
    #p0  : np.array = field(init=False)
    #p1  : np.array = field(init=False)

    def __post_init__(self):
        self.p0 = np.array([0, 0, self.zmin]) # point in one endcup
        self.p1 = np.array([0, 0, self.zmax]) # point in the other

        mag, v, n1, n2 = self.unit_vectors_()
        self.P, self.P2, self.P3 = self.surfaces_(mag,  v, n1, n2)

    def normal_to_barrel(self, P : np.array)->np.array:
        """Normal to the cylinder barrel
        Uses equation of cylynder: F(x,y,z) = x^2 + y^2 - r^2 = 0
        then n = Grad(F)_P /Norm(Grad(F))_P
        n = (2x, 2y, 0)/sqrt(4x^2 + 4y^2) = (x,y,0)/r (P)
        """
        return np.array([P[0], P[1], 0]) / self.r

    def cylinder_equation(self, P : np.array)->float:
        return P[0]**2 + P[1]** 2 - self.r**2

    @property
    def length(self)->float:
        return self.zmax - self.zmin

    @property
    def perimeter(self)->float:
        return 2 * np.pi * self.r

    @property
    def area_barrel(self)->float:
        return 2 * np.pi * self.r * self.length

    @property
    def area_endcap(self)->float:
        return np.pi * self.r **2

    @property
    def area(self)->float:
        return self.area_barrel + 2 * self.area_endcap

    @property
    def volume(self)->float:
        return np.pi * self.r **2 * self.length

    def unit_vectors_(self):
        #vector in direction of axis
        v = self.p1 - self.p0

        #find magnitude of vector
        mag = norm(v)

        #unit vector in direction of axis
        v = v / mag

        # choose (1,0,0) as second axis unless is first axis
        not_v = np.array([1, 0, 0])
        if (v == not_v).all():
            not_v = np.array([0, 1, 0])

        #make vector perpendicular to v and not v
        n1 = np.cross(v, not_v)

        #normalize n1
        n1 /= norm(n1)

        #make unit vector perpendicular to v and n1
        n2 = np.cross(v, n1)

        return mag, v,  n1, n2

    def surfaces_(self, mag, v, n1, n2):
        #surface ranges over t from 0 to length of axis and 0 to 2*pi
        t = np.linspace(0, mag, 2)
        theta = np.linspace(0, 2 * np.pi, 100)
        rsample = np.linspace(0, self.r, 2)

        #use meshgrid to make 2d arrays
        t, theta2 = np.meshgrid(t, theta)

        rsample,theta = np.meshgrid(rsample, theta)

        #generate coordinates for surface
        # "Tube"
        X, Y, Z = [self.p0[i] + v[i] * t + self.r * np.sin(theta2) * n1[i] + self.r * np.cos(theta2) *  n2[i] for i in range(3)]
        # "Bottom"
        X2, Y2, Z2 = [self.p0[i] + rsample[i] * np.sin(theta) * n1[i] + rsample[i] * np.cos(theta) * n2[i] for i in range(3)]
        # "Top"
        X3, Y3, Z3 = [self.p0[i] + v[i]*mag + rsample[i] * np.sin(theta) * n1[i] + rsample[i] * np.cos(theta) * n2[i] for i in range(3)]
        return (X,Y,Z), (X2, Y2, Z2), (X3, Y3, Z3)


@dataclass
class Ray:
    e   : np.array
    d   : np.array
    def ray(self,t):
        return self.e + t * self.d
    @property
    def unit_vector(self):
        v = self.d - self.e
        return v / norm(v)


## Generation of photons and transport

def throw_dice(dice : float)->bool:
    """Throws a random number and compares with value of dice. Returns true if random < dice"""
    cond = False
    if np.random.random() < dice:
        cond = True
    return cond

def sample_spherical(npoints: int, ndim: int=3)->np.array:
    """Generate points distributed in the surface of a unit sphere.
    The points are in a matrix ((x1,x2,x3...xn), (y1,y2,y3...yn), (z1,z2,z3...zn))
    where n is the number of random points

    """
    vec =  np.random.randn(ndim, npoints)
    vec /= np.linalg.norm(vec, axis=0)
    return vec


def vectors_spherical(npoints: int, ndim: int=3)->np.array:
    """Generate vectors distributed in the surface of a unit sphere.
    The vectros are in a matrix ((x1,y1,z1), (x2,y2,z2)... (xn,yn, zn))
    where n is the number of random points

    """
    return sample_spherical(npoints).T


def in_endcaps(c: Cylinder, p : np.array)->bool:
    """Returns True if point in end-caps of cylinder"""
    close = np.isclose(np.array([p[2],p[2]]), np.array([c.zmin, c.zmax]), atol=1e-06)
    return close.any()


def cylinder_intersection_roots(r: Ray, c: Cylinder, eps: float =1e-9)->np.array:
    """Computes intersection roots between a ray and a cylinder"""

    a = r.d[0]**2 + r.d[1]**2
    b = 2 * (r.e[0] * r.d[0] + r.e[1] * r.d[1])
    c = r.e[0]**2 + r.e[1]**2 - c.r**2

    discriminant = b**2 - 4 * a * c
    if discriminant < 0:
        # No real intersection
        return 0
    #print(f"eqr: a = {a}, b = {b}, c = {c}")

    roots =  np.roots([a,b,c])
    proots = [x for x in roots if x>eps]
    if proots:
        return np.min(proots)
    else:
        return 0


def ray_intersection_with_cylinder_end_caps(r: Ray, c: Cylinder, t: float)->np.array:
    """Intersection between a ray and the end-cups of a cylinder"""
    p = r.ray(t)
    if p[2] > c.zmax:
        t = (c.zmax - r.e[2])/r.d[2]
    else:
        t = (c.zmin - r.e[2])/r.d[2]

    return t, r.ray(t)


def ray_intersection_with_cylinder(r: Ray, c:Cylinder)->Tuple[float,np.array]:
    """Intersection between a ray and a cylinder"""
    t = cylinder_intersection_roots(r, c)
    P = r.ray(t)  # proyection to cylinder shell
    z = P[2]
    #if z < c.zmin or z > c.zmax:
    #    t, P = ray_intersection_with_cylinder_end_caps(r, c, t)
    return t, P

## Generate points

def generate_random_points_cylinder_shell(c:Cylinder, num_points:int):
    points = []
    for _ in range(num_points):
        # Generate a point on the cylindrical shell (r = R)
        phi = np.random.uniform(0, 2 * np.pi)
        #z = np.random.uniform(c.zmin, c.zmax)
        z = np.random.uniform(c.zmin, c.zmax)
        x = c.r * np.cos(phi)
        y = c.r * np.sin(phi)
        points.append(np.array((x, y, z)))
    return points


def generate_random_points_cylinder_endcups(c:Cylinder, num_points:int):
    points = []
    for _ in range(num_points):
        # Generate a point on one of the end-cups (z = +/- length/2)
        z = c.zmin if np.random.rand() < 0.5 else c.zmax
        r = np.sqrt(np.random.uniform(0, c.r ** 2))  # Uniform distribution in the circular end-cap
        phi = np.random.uniform(0, 2 * np.pi)
        x = r * np.cos(phi)
        y = r * np.sin(phi)
        points.append(np.array((x, y, z)))
    return points


def generate_and_transport_gammas_lxe_shell(c1 : Cylinder, c2 : Cylinder, 
                                       nphotons: int=10, ndx: int= 100, 
                                       solidAngle="half",
                                       verbose=0, drawRays=True, scale=10)->np.array:
    """
    This function propagates gammas emited by a "background cylinder", c1 (for example a copper vessel)
    to a "fiducial cylinder", c2, that is, the region self-shielded by LXe where signal is sought. 

    The gammas are uniformly  in the shell of c1 and transport them to c2.
    The technique is as follows:
    1. Generate random points in the surface of cylinder c1. This corresponds to 
    gamma emiters in the "background" layer represented by c1.

    2. For each point (gamma emiter), generate random directions. 
       If solidAngle = 'half', count only those photons that point inwards from c1 shell. 
       This case corresponds to situations in which the "background layer" corresponds to 
       a self-shielded volume (like a copper shield), such that only gammas going inwards are counted
       (e.g., the activity of the layer corresponds to gammas going into the detector).
       If solidAngle = 'full' count all photons (this is a normal layer without self-shielding)

    3. shoot a photon with the generated direction (only if inwards), and propagate towards the 
    inner cylinder (e.g., the fiducial region, represented by an inner cylinder). Photons that do not
    reach the inner cylinder have weight zero in the count (they no to contribute to the background).
    Photons that reach the inner cylinder enter the count with exp(-d/latt), where d is the distance
    travelled by the gamma and latt = 8.5 cm is the attenuation length of xenon. 
    Compute a weight LXe, explain. 

    All units assumed to be in mm

    """

    def draw_rays(e, r):
        tt = np.linspace(0, scale, 100)
        xi = e[0] + tt * r[0]
        yi = e[1] + tt * r[1]
        zi = e[2] + tt * r[2]
        ax.plot(xi, yi, zi)
                
    n = int(nphotons)

    if verbose > 0:
        print(f"c1 ={c1}")
        print(f"c2 ={c2}")

    fig = plt.figure(figsize=(6,6))
    ax=plt.subplot(111, projection='3d')
    draw_2cylynder_surface(c1, c2, ax)

    P = generate_random_points_cylinder_shell(c1, n)
    xi,yi,zi = xyz_from_points(P)
    ax.scatter(xi, yi, zi, s=5, c='k', zorder=10)

    WF = []
    DST = []
    latt = 87 # mm, Xenon Latt
    for p in P:
        if verbose >1:
            print(f"Point = {p}, r= {np.sqrt(p[0]**2 + p[1]**2)}, z = {p[2]}")
        W = 0
        icd = 0
        icr = 0
        icz = 0
        ice = 0
        ict = 0

        # generate directions:
        D = vectors_spherical(ndx)

        # loop over directions
        for d in D:
            fce = True  # end-cap fiducial
            fct = True  # end-cap transport

            nx = c1.normal_to_barrel(p) # gamma must go inwards
            if np.dot(nx,d) > 0:
                icd+=1
                continue
                
            r = Ray(p,d)
            t, pt = ray_intersection_with_cylinder(r, c2)
            if t <=0:       # must intersect cylinder
                icr+=1
                continue    # loop away if does not intersect cylinder

            if pt[2] < c2.zmin or pt[2] > c2.zmax: # outside z bounds
                icz+=1
                continue

            # is the end-cap intersection at smaller distance?
            ## This part to be revised

            # t2, pt2 = ray_intersection_with_cylinder_end_caps(r, c2, t)

            # if np.sqrt(pt2[0]**2 + pt2[1]**2) > c2.r: # outside r bounds
            #     ice+=1
            #     fce = False  #not in fiducial 
            
            # if t2 <=0:       # must intersect end-cap
            #     ict+=1
            #     fct = False # transport failed

            if drawRays:
                draw_rays(p, d)
                
            dx = np.linalg.norm(p - pt) # distance traveled to shell

            # if fce and fct: 
            #     dx2 = np.linalg.norm(p - pt2) # distance traveled to end-cap
            #     dx = np.minimum(dx1,dx2)
            # else:
            #     dx = dx1
            
        
            DST.append(dx)
            wx = np.exp(-dx/latt)       # attenuation
            
            if verbose >2:
                print(f"d= {d}, distance = {dx}, wx = {wx}")
                
            W+=wx
        
        if solidAngle=="half":
            den = ndx - icd
        else:
            den = ndx
        WF.append(W/den)
    return WF, DST



def generate_and_transport_gammas_lxe_endcaps(c1 : Cylinder, c2 : Cylinder, 
                                       nphotons: int=10, ndx: int= 100, 
                                       solidAngle="half",
                                       verbose=0, drawRays=True, scale=10)->np.array:
    """
    This function propagates gammas emited by a "background cylinder", c1 (for example a copper vessel)
    to a "fiducial cylinder", c2, that is, the region self-shielded by LXe where signal is sought. 

    The gammas are uniformly  in the end-caps of c1 and transport them to c2.
    The technique is as follows:
    1. Generate random points in the two end-caps of cylinder c1. This corresponds to 
    gamma emiters in the "background" layer represented by c1.

    2. For each point (gamma emiter), generate random directions. 
       If solidAngle = 'half', count only those photons that point inwards from c1 shell. 
       This case corresponds to situations in which the "background layer" corresponds to 
       a self-shielded volume (like a copper shield), such that only gammas going inwards are counted
       (e.g., the activity of the layer corresponds to gammas going into the detector).
       If solidAngle = 'full' count all photons (this is a normal layer without self-shielding)

    3. shoot a photon with the generated direction (only if inwards), and propagate towards the 
    inner cylinder (e.g., the fiducial region, represented by an inner cylinder). Photons that do not
    reach the inner cylinder have weight zero in the count (they no to contribute to the background).
    Photons that reach the inner cylinder enter the count with exp(-d/latt), where d is the distance
    travelled by the gamma and latt = 8.5 cm is the attenuation length of xenon. 
    Compute a weight LXe, explain. 

    All units assumed to be in mm

    """

    def draw_rays(e, r):
        tt = np.linspace(0, scale, 100)
        xi = e[0] + tt * r[0]
        yi = e[1] + tt * r[1]
        zi = e[2] + tt * r[2]
        ax.plot(xi, yi, zi)
                
    n = int(nphotons)

    if verbose > 0:
        print(f"c1 ={c1}")
        print(f"c2 ={c2}")

    
    fig = plt.figure(figsize=(6,6))
    ax=plt.subplot(111, projection='3d')
    draw_2cylynder_surface(c1, c2, ax)

    P = generate_random_points_cylinder_endcups(c1, n)
    xi,yi,zi = xyz_from_points(P)
    ax.scatter(xi, yi, zi, s=5, c='k', zorder=10)

    WF = []
    DST = []
    latt = 87 # mm, Xenon Latt
    for p in P:
        if verbose >1:
            print(f"Point = {p}, r= {np.sqrt(p[0]**2 + p[1]**2)}, z = {p[2]}")
        W = 0
        icd = 0
        icr = 0
        icz = 0
        ice = 0
        ict = 0

        # generate directions:
        D = vectors_spherical(ndx)

        # loop over directions
        for d in D:
            fce = True  # barrel fiducial
            fct = True  # barrel transport

        
            z = p[2]  # either c1.zmin or c1.zmax
            pz = d[2] # must point inwards (-pz in zmax, +pz in zmin)
            if z == c1.zmin:
                zt = c2.zmin 
            elif z == c1.zmax:
                zt = c2.zmax
            else:
                print(f"unexpected: z = {z}")
                continue

            # Calculate the parameter t for the intersection with the end-cap plane
            # if z = c1.zmin zt - z = c2.zmin - c1.zmin > 0 and pz > 0
            # if z = c1.zmmax zt - z = c2.zmax - c1.zmax <0  and pz <0 
            # condition is that t > 0
            t = (zt - z)/ pz 
                
            if t <= 0:
                icd+=1
                #print(f"z = {z}, zt = {zt}, pz = {pz} ")
                continue
    
            
            # 1. Ray intersects with end-caps 
            # Calculate the intersection point
            r = Ray(p,d)
            pt = np.zeros(3)
            pt[0] = p[0] + t * d[0]
            pt[1] = p[1] + t * d[1]
            pt[2] = zt

            # must be within the radius of c2
            if np.sqrt(pt[0]**2 + pt[1]**2) > c2.r:
                icr+=1
                #print(f"after extrap to z: r = {np.sqrt(pt[0]**2 + pt[1]**2)} > {c2.r} ")
                continue

            
            # 2. Ray intersects with cylinder c2
            
            t2, pt2 = ray_intersection_with_cylinder(r, c2)

            if t2 <=0:       # does not intersect cylinder
                fce = False
                
            if pt2[2] < c2.zmin or pt2[2] > c2.zmax: # outside z bounds
                fct = False

            if drawRays:
                draw_rays(p, d)
                
            dx1 = np.linalg.norm(p - pt) # distance traveled to end-cap

            if fce and fct: 
                dx2 = np.linalg.norm(p - pt2) # distance traveled to shell
                dx = np.minimum(dx1,dx2)
            else:
                dx = dx1
            
            DST.append(dx)
            wx = np.exp(-dx/latt)       # attenuation
            
            if verbose >2:
                print(f"d= {d}, distance = {dx}, wx = {wx}")
                
            W+=wx
        
        if solidAngle=="half":
            den = ndx - icd
        else:
            den = ndx
        WF.append(W/den)
    return WF, DST


def generate_and_transport_gammas_gxe_shell(c1 : Cylinder,  
                                       nphotons: int=10, ndx: int= 100, 
                                       solidAngle="half",
                                       verbose=0, drawRays=True, scale=10)->np.array:
    """
    This function propagates gammas emited from the surface of c1 (shell) until they exit 
    the cylinder c1, computing their length and assigning them a weight proportional to it.
    

    The technique is as follows:
    1. Generate random points in the surface of cylinder c1. This corresponds to 
    gamma emiters in the "background" layer represented by c1.

    2. For each point (gamma emiter), generate random directions. 
       If solidAngle = 'half', count only those photons that point inwards from c1 shell. 
       This case corresponds to situations in which the "background layer" corresponds to 
       a self-shielded volume (like a copper shield), such that only gammas going inwards are counted
       (e.g., the activity of the layer corresponds to gammas going into the detector).
       If solidAngle = 'full' count all photons (this is a normal layer without self-shielding)

    3. shoot a photon with the generated direction (only if inwards), and propagate through the cylinder
     until it exists through the shell or the end-caps. Compute the distance d they travel and
     assign them a weight exp(-d/latt), where d is the distance
    travelled by the gamma and latt = 285.2 cm is the attenuation length of xenon at 15 bar. 
    

    All units assumed to be in mm

    """

    def draw_rays(e, r):
        tt = np.linspace(0, scale, 100)
        xi = e[0] + tt * r[0]
        yi = e[1] + tt * r[1]
        zi = e[2] + tt * r[2]
        ax.plot(xi, yi, zi)
                
    n = int(nphotons)

    if verbose > 0:
        print(f"c1 ={c1}")
       

    fig = plt.figure(figsize=(6,6))
    ax=plt.subplot(111, projection='3d')
    draw_cylynder_surface(c1, ax)

    P = generate_random_points_cylinder_shell(c1, n)
    xi,yi,zi = xyz_from_points(P)
    ax.scatter(xi, yi, zi, s=5, c='k', zorder=10)

    WF = []
    DST = []
    latt = 2900 # mm, Xenon Latt
    for p in P:
        if verbose >1:
            print(f"Point = {p}, r= {np.sqrt(p[0]**2 + p[1]**2)}, z = {p[2]}")
        W = 0
        icd = 0
        icr = 0
        icz = 0
        ice = 0
        ict = 0

        # generate directions:
        D = vectors_spherical(ndx)

        # loop over directions
        for d in D:
            fce = False  # end-cap fiducial
            fct = False  # end-cap transport

            nx = c1.normal_to_barrel(p) # gamma must go inwards
            if np.dot(nx,d) > 0:
                icd+=1
                continue
                
            r = Ray(p,d)
            t, pt = ray_intersection_with_cylinder(r, c1)
            if t <=0:       # must intersect cylinder
                icr+=1
                continue    # loop away if does not intersect cylinder

            if pt[2] < c1.zmin or pt[2] > c1.zmax: # outside z bounds

                # does it intercept the end-cap?
                t2, pt2 = ray_intersection_with_cylinder_end_caps(r, c1, t)
                if t2 <=0:       # must intersect cylinder
                    icz+=1
                    continue    # loop away if does not intersect end-cap
                
                if pt2[2] < c1.zmin or pt2[2] > c1.zmax: # outside z bounds
                    ice+=1
                    continue    # loop away if does not intersect end-cap
                fce = True
            else:
                fct = True

            if drawRays:
                draw_rays(p, d)

            if fct:
                dx = np.linalg.norm(p - pt) # distance traveled to end-cap
            elif fce: 
                dx = np.linalg.norm(p - pt2) # distance traveled to shell
            else:
                continue
                
            
            DST.append(dx)
            wx = np.exp(-dx/latt)       # attenuation
            
            if verbose >2:
                print(f"d= {d}, distance = {dx}, wx = {wx}")
                
            W+=wx
        
        if solidAngle=="half":
            den = ndx - icd
        else:
            den = ndx
        WF.append(W/den)
    return WF, DST


## Graphics

def xyz_from_points(P):
    PP = np.array(P)
    xi,yi,zi = PP[:].T
    return xi,yi,zi

def set_fonts(ax, fontsize=20):
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(fontsize)


def draw_cylynder_surface(c: Cylinder, ax, units=mm, alpha=0.2, barrelColor='blue', cupColor='red'):
    ax.plot_surface(c.P[0]/units, c.P[1]/units, c.P[2]/units, color=barrelColor, alpha=alpha)
    ax.plot_surface(c.P2[0]/units, c.P2[1]/units, c.P2[2]/units, color=cupColor, alpha=alpha)
    ax.plot_surface(c.P3[0]/units, c.P3[1]/units, c.P3[2]/units, color=cupColor, alpha=alpha)


def draw_2cylynder_surface(c1: Cylinder, c2: Cylinder, ax, units=mm, alpha=0.2, 
                           barrelColor1='blue', cupColor1='red',
                           barrelColor2='green', cupColor2='orange'):
    
    ax.plot_surface(c1.P[0]/units, c1.P[1]/units, c1.P[2]/units, color=barrelColor1, alpha=alpha)
    ax.plot_surface(c1.P2[0]/units, c1.P2[1]/units, c1.P2[2]/units, color=cupColor1, alpha=alpha)
    ax.plot_surface(c1.P3[0]/units, c1.P3[1]/units, c1.P3[2]/units, color=cupColor1, alpha=alpha)
    ax.plot_surface(c2.P[0]/units, c2.P[1]/units, c2.P[2]/units, color=barrelColor2, alpha=alpha)
    ax.plot_surface(c2.P2[0]/units, c2.P2[1]/units, c2.P2[2]/units, color=cupColor2, alpha=alpha)
    ax.plot_surface(c2.P3[0]/units, c2.P3[1]/units, c2.P3[2]/units, color=cupColor2, alpha=alpha)


def draw_cylinder(c : Cylinder, units=mm, alpha=0.2, barrelColor='blue', cupColor='red',
                  figsize=(16,16), DWORLD=False, WDIM=((-1,1),(-1,1),(-1,1))):


    fig = plt.figure(figsize=figsize)
    ax=plt.subplot(111, projection='3d')
    if DWORLD:
        ax.set_xlim3d(WDIM[0][0], WDIM[0][1])
        ax.set_ylim3d(WDIM[1][0], WDIM[1][1])
        ax.set_zlim3d(WDIM[2][0], WDIM[2][1])
    draw_cylynder_surface(c, ax, units, alpha,  barrelColor, cupColor)
    #ax.plot_surface(c.P[0], c.P[1], c.P[2], color=barrelColor, alpha=alpha)
    #ax.plot_surface(c.P2[0], c.P2[1], c.P2[2], color=cupColor, alpha=alpha)
    #ax.plot_surface(c.P3[0], c.P3[1], c.P3[2], color=cupColor, alpha=alpha)
    plt.show()


def draw_2cylinder(c1: Cylinder, c2: Cylinder, units=mm, alpha=0.2, 
                   barrelColor1='blue', cupColor1='red',
                   barrelColor2='green', cupColor2='orange',
                  figsize=(16,16), DWORLD=False, WDIM=((-1,1),(-1,1),(-1,1))):


    fig = plt.figure(figsize=figsize)
    ax=plt.subplot(111, projection='3d')
    if DWORLD:
        ax.set_xlim3d(WDIM[0][0], WDIM[0][1])
        ax.set_ylim3d(WDIM[1][0], WDIM[1][1])
        ax.set_zlim3d(WDIM[2][0], WDIM[2][1])
    draw_2cylynder_surface(c1, c2, ax, units, alpha,  barrelColor1, cupColor1, barrelColor2, cupColor2)
    #ax.plot_surface(c.P[0], c.P[1], c.P[2], color=barrelColor, alpha=alpha)
    #ax.plot_surface(c.P2[0], c.P2[1], c.P2[2], color=cupColor, alpha=alpha)
    #ax.plot_surface(c.P3[0], c.P3[1], c.P3[2], color=cupColor, alpha=alpha)
    plt.show()


def draw_cylnder_nomal_at_P(P: np.array, c : Cylinder, tscale=1,
                            units=mm, alpha=0.2, barrelColor='blue', cupColor='red',
                            figsize=(16,16)):

    N = c.normal_to_barrel(P)

    def draw_normal(P,N):
        tt = np.linspace(0, tscale, 100)

        xi = P[0] + tt * N[0]
        yi = P[1] + tt * N[1]
        zi = P[2] + tt * N[2]
        ax.plot(xi, yi, zi)

    xi,yi,zi = xyz_from_points(np.array([P]))

    fig = plt.figure(figsize=figsize)
    ax=plt.subplot(111, projection='3d')
    draw_cylynder_surface(c, ax, units, alpha, barrelColor, cupColor)

    draw_normal(P, N)
    ax.scatter(xi, yi, zi, s=25, c='k', zorder=10)
    plt.show()


