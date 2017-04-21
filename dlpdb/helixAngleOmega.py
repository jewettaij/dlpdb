
from math import sqrt, cos, sin, tan, acos, asin, atan, pi

def CalcOmega(theta, phi):
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


