# Function to find the intersection of a ray with a cylinder
def intersect_ray_cylinder(p, d, r2, zmin2, zmax2):
    # Extract coordinates of point and direction
    x, y, z = p
    px, py, pz = d

    #print(f"p= {p}, d = {d}")
    
    # Calculate the coefficients of the quadratic equation for intersection
    a = px**2 + py**2
    b = 2 * (x * px + y * py)
    c = x**2 + y**2 - r2**2
    
    # Solve the quadratic equation a*t^2 + b*t + c = 0
    discriminant = b**2 - 4 * a * c

    #print(f" tst: a= {a}, b = {b}, c = {c}, discriminant = {discriminant}")
    
    if discriminant < 0:
        # No real intersection
        return None
    
    # Calculate the two possible solutions for t (parameter along the ray)
    t1 = (-b - np.sqrt(discriminant)) / (2 * a)
    t2 = (-b + np.sqrt(discriminant)) / (2 * a)

    #print(f"tst: t1 = {t1}, t2 = {t2}")
    # We need the closest positive t, since we are interested in the forward direction
    t = min(t for t in [t1, t2] if t > 0)
    
    # Calculate the intersection point
    xi = x + t * px
    yi = y + t * py
    zi = z + t * pz
    
    # Check if the intersection point is within the z bounds of cylinder c2
    if zmin2 <= zi <= zmax2:
        return (xi, yi, zi)
    else:
        # Intersection is outside the z bounds of cylinder c2
        return None

# Example usage
r1 = 10  # Radius of cylinder c1
r2 = 9   # Radius of inner cylinder c2
zmin1 = -10  # Minimum z of cylinder c1
zmax1 = 10   # Maximum z of cylinder c1
zmin2 = -9   # Minimum z of cylinder c2
zmax2 = 9    # Maximum z of cylinder c2

# Point on the surface of cylinder c1
p = (10, 0, 0)  # (x, y, z) such that sqrt(x^2 + y^2) = r1 and zmin1 < z < zmax1
nn = (p[0], p[1], 0)
# Random direction pointing inward

ng = 10 
D = vectors_spherical(ng)

c2 =GCylinder(r2, zmin2, zmax2)
print(f"c2 = {c2}")
       
for vec in D:
    d = vec/np.linalg.norm(vec)

    sgn = np.dot(d, nn)
    print(f"sgn = {sgn}")
    if sgn > 0:
        continue

    intersection = intersect_ray_cylinder(p, d, r2, zmin2, zmax2)
    
    if intersection:
        print(f"Intersection point with cylinder c2: {intersection}")
    else:
        print("No intersection within the bounds of cylinder c2")

    r = Ray(p,d)
    t, Pt = ray_intersection_with_cylinder(r, c2)
    print(f"t1 = {t}, Pt1 = {Pt}")
    print(f"Pt: r = {np.sqrt(Pt[0]**2 + Pt[1]**2)}")
