import sys
from math import sqrt, cos, sin, tan, acos, asin, atan, pi, floor

def DotProd(va, vb):
    assert(len(va) == len(vb))
    result = 0.0
    for i in range(0, len(va)):
        result += va[i]*vb[i]
    return result

def LengthVect(v):
    return math.sqrt(DotProd(v, v))

def AddVect(va, vb):
    assert(len(va) == len(vb))
    result = [0.0 for i in range(0, len(va))]
    for i in range(0, len(va)):
        result[i] = va[i] + vb[i]
    return result

def SubtractVect(va, vb):
    assert(len(va) == len(vb))
    result = [0.0 for i in range(0, len(va))]
    for i in range(0, len(va)):
        result[i] = va[i] - vb[i]
    return result

def ScaleVect(c, v):
    result = [0.0 for i in range(0, len(v))]
    for i in range(0, len(v)):
        result[i] = c*v[i]
    return result

def CrossProd3(a,b):
    c = [0.0,0.0,0.0]
    c[0] = a[1]*b[2] - a[2]*b[1]
    c[1] = a[2]*b[0] - a[0]*b[2]
    c[2] = a[0]*b[1] - a[1]*b[0]
    return c


def PrintVect(v):
    for d in range(0, len(v)):
        sys.stdout.write(str(v[d]))
        if d+1 < len(v):
            sys.stdout.write(' ')
    sys.stdout.write('\n')


def ClosestLinePoints(ra0, rb0, va, vb):
    """
    Return the closest pair of points in 3D along two lines of infinite length
    pointing in directions va,vb and passing through points ra0,rb0.
    """
    assert(len(ra0) == len(rb0) == len(va) == len(vb))
    # Derivation:
    # I use "bra" "ket" notation to denote the vectors for each point.
    # The two points are denoted |ra> and |rb>.
    # |ra> is a point on a line along the direction |va> passing through |ra0>.
    # |rb> is a point on a line along the direction |vb> passing through |rb0>.
    # The formula for the closest pair of points, |ra> and |rb>, is:
    # |ra> = |ra0> + ta |va>
    # |rb> = |rb0> + tb |vb>
    # Where:
    # ta = (vb^2 <va| - <va|vb> <vb|) (|rb0> - |ra0>) / (va^2 vb^2 - <va|vb>^2)
    # tb = (va^2 <vb| - <va|vb> <va|) (|ra0> - |rb0>) / (va^2 vb^2 - <va|vb>^2)
    # Notation used:
    # 1) The notation <V|W> denotes the dot product between vectors |V> and |W>
    # 2) Similarly, notation like: (<va| + <vb|) (|ra0> - |rb>) denotes the
    #   dot product between vector sum (<va| + <vb|) and the sum (|ra0> - |rb0>)
    #   (and also equals <va|ra0> + <vb|ra0> - <va|rb0> - <vb|rb0>)
    # 3) va^2 is shorthand for the squared length of vector va (=<va|va>)
    rab = SubtractVect(ra0, rb0)  # = |ra> - |rb>
    va2 = DotProd(va,va)
    vb2 = DotProd(vb,vb)
    va_vb = DotProd(va, vb)
    descr = va2*vb2 - va_vb*va_vb
    if descr == 0.0:
        # Then the lines are parallel.
        # In that case, try to pick a point along each line
        # which lies half-way in between ra0 and rb0
        scaleb = 0.5*DotProd(rab, vb) / vb2
        deltab = ScaleVect(scaleb, vb)
        rb = AddVect(rb0, deltab)
        ra = SubtractVect(ra0, deltab)
        return (ra, rb)
    ta = -DotProd(SubtractVect(ScaleVect(vb2,va),ScaleVect(va_vb,vb)),rab)/descr
    tb =  DotProd(SubtractVect(ScaleVect(va2,vb),ScaleVect(va_vb,va)),rab)/descr
    ra = AddVect(ra0, ScaleVect(ta, va))
    rb = AddVect(rb0, ScaleVect(tb, vb))
    return (ra, rb)
