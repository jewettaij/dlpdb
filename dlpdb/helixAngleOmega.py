from math import sqrt, cos, sin, tan, acos, asin, atan, pi


def length_v(r):
    lsqd = 0.0
    for d in range(0,len(r)):
        lsqd += r[d]*r[d]
    return sqrt(lsqd)


def inner_prod_v(r1,r2):
    result = 0.0
    for d in range(0,len(r1)):
        result += r1[d]*r2[d]
    return result


def cross_prod_v3(a,b):
    c = [0.0,0.0,0.0]
    c[0] = a[1]*b[2] - a[2]*b[1]
    c[1] = a[2]*b[0] - a[0]*b[2]
    c[2] = a[0]*b[1] - a[1]*b[0]
    return c



def CalcOmegaFromThetaPhi(theta, phi):
    """
    This function takes two arguments
    theta (the 3-body bond angle,    1st argument)
    phi   (the 4-body torsion angle, 2nd argument)
    ...and returns the "Omega" rotation angle. 

    Omega is the "exterior-angle" made by the poligonal shadow of the backbone
    trace of a helix, when viewed down it's axis of symmetry ("z-axis")
    and projected onto the plane below (the "xy plane").

    For example the text-diagram below, suppose we are viewing a helix
    down it's symmetry axis, and the 4 consecutive alpha-carbon
    atoms leave the following shadow on the plain below.

       *
        \
     *---------------.*     _
          \         _/      /`
           \      _/       /
            \   _/    <-- Omega
             \ /          angle
              *'

    For an alpha helix (which has roughly 3.6 amino-acids per 360-degree turn,
    Omega ~= 360/3.6 ~= 100.0 degrees

    The Omega angle is related to the 3-body bond angle and 4-body torsion
    angle of these 4 atoms (theta and phi).
    All angles (inputs and outputs) are in units of radians.

    """

    Theta = pi - theta
    Phi = pi - phi

    # Derivation:
    #   This comment does not make sense unless 
    #   you also see the figure that goes with it.
    #
    #tan(phi/2) = deltaY / X
    #        X  = cos(a) sin(Omega)
    #   deltaY  = deltay cos(a)
    #    deltay = sin(a) (1-cos(Omega))
    #   -> deltaY = cos(a)sin(a) (1-cos(Omega))
    # -> tan(phi/2) = sin(a)   (1-cos(Omega)) / sin(Omega)
    #               = sin(a) 2 sin^2(Omega/2) / sin(Omega)
    #               = sin(a) 2 sin^2(Omega/2) / 2 sin(Omega/2) cos(Omega/2)
    #               = sin(a) tan(Omega/2)
    #
    #                 cos(a) = sin(Theta/2) / sin(Omega/2)
    #              -> sin(a) = sqrt(1 - sin^2(Theta/2) / sin^2(Omega/2))
    #
    # -> tan(phi/2) = sqrt(1 - sin^2(Theta/2) / sin^2(Omega/2))  tan(Omega/2)
    #               = sqrt(sin^2(Omega/2) - sin^2(Theta/2)) / cos(Omega/2)
    # ->
    #  tan^2(phi/2) = (sin^2(Omega/2) - Sin^2(Theta/2)) / cos^2(Omega/2)
    # We know phi, and Theta.  Now solve for Omega numerically

    tan2phi2 = tan(0.5*phi)
    tan2phi2 *= tan2phi2

    lower_bound = 0.0
    upper_bound = pi
    tolerance = 0.000001
    while ((upper_bound - lower_bound) > tolerance):
        Omega = lower_bound + 0.5*(upper_bound - lower_bound)
        cos2Ome2 =  cos(0.5*Omega)
        cos2Ome2 *= cos2Ome2
        sin2Ome2 =  sin(0.5*Omega)
        sin2Ome2 *= sin2Ome2
        sin2The2 =  sin(0.5*Theta)
        sin2The2 *= sin2The2

        if (sin2Ome2 - sin2The2) / cos2Ome2 > tan2phi2:
            upper_bound = Omega
        else:
            lower_bound = Omega

    # Deal with possible negative values of Omega (left-handed helices):
    if phi < 0.0:
        Omega = -Omega

    return Omega





def CalcOmega(r0, r1, r2, r3):
    r10 = [0.0, 0.0, 0.0]
    r21 = [0.0, 0.0, 0.0]
    r32 = [0.0, 0.0, 0.0]
    for d in range(0,3):
        r10[d] = r1[d] - r0[d]
        r21[d] = r2[d] - r1[d]
        r32[d] = r3[d] - r2[d]
    l10 = length_v(r10)
    l21 = length_v(r21)
    l32 = length_v(r32)

    n012 = cross_prod_v3(r10, r21)
    n123 = cross_prod_v3(r21, r32)


    # The torsion-angle or 4-body angle is named "angle0124"
    cos_phi = inner_prod_v(n012, n123) /(length_v(n012)*length_v(n123))

    # There is a problem whenever 4 consecutive atoms are coplanar:
    #
    #            *---*
    #                |      (all 4 atoms are coplanar, and phi = 0)
    #            *---*
    #
    # In this example, the torsion angle phi is well defined and =0.
    # The problem is that, due to floating point roundoff
    # "cos_phi" can sometimes slightly exceed 1.
    # This causes a NAN when you calculate acos(cos_phi).

    if (cos_phi > 1.0):
        cos_phi = 1.0
    elif (cos_phi < -1.0):
        cos_phi = -1.0

    phi = acos(cos_phi)

    # This formula does not distinguish positive and negative phi.
    #
    # Negative torsion angles:
    #
    # Check if  the position of atom i+3 is above the phi=0 plane
    # (in the region of positive phi), or below the phi=0 plane.
    # It is above the phi=0 plane if the bond from atom i+2 to i+3
    # points in the same direction as the negative-phi-tangent-vector 
    # for atom i (not i+3)  (...which points in the n012 direction)
    if inner_prod_v(n012, r32) < 0.0:
        phi = -phi

    # The two bond-angles or 3-body angles are named "angle012" and "angle123"
    angle012 = acos( -inner_prod_v(r10, r21) / (l10 * l21) )
    angle123 = acos( -inner_prod_v(r21, r32) / (l21 * l32) )
    # (The negative sign above comes from the fact that we are really 
    #  interested in the angle between r21 and r01 (which is -r10).)

    # Convert these angles to degrees, and print it out
    #sys.stdout.write(str(angle012*180.0/pi)+' ')
    #sys.stdout.write(str(angle123*180.0/pi)+' ')

    # Let theta = the average of the two 3-body angles
    theta = 0.5 * (angle012 + angle123) 

    # Omega (usually 360/3.6 ~= 100 degrees) is the helix rotation angle.
    Omega = CalcOmega(theta, phi)

    return Omega
