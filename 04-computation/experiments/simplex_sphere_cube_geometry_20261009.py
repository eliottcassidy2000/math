"""Three exact norm sections of one marked positive ray.

Geometry transports authenticated source coordinates, not ROOT proofs.
Sphere coordinates are stored exactly as integer numerators / sqrt(radicand).
"""
from dataclasses import dataclass, replace
from fractions import Fraction as F
from functools import reduce
from itertools import product
from math import gcd, isqrt, lcm

CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(message)


def integer(n, least=0):
    if type(n) is not int or n < least:
        raise ValueError('exact integer in declared domain required')
    return n


def rational(n):
    if type(n) not in (int,F):
        raise ValueError('exact rational required')
    return F(n)


def vector(v):
    if type(v) is not tuple or len(v) < 2:
        raise ValueError('tuple vector of dimension at least two required')
    return tuple(rational(x) for x in v)


def primitive(v):
    v = vector(v)
    if any(x < 0 for x in v) or not any(v):
        raise ValueError('nonzero nonnegative ray required')
    denominator = lcm(*(x.denominator for x in v))
    nums = tuple(int(x*denominator) for x in v)
    common = reduce(gcd,nums)
    return tuple(x//common for x in nums)


@dataclass(frozen=True)
class Section:
    kind: str
    numerator: tuple
    radicand: int


def encode(v, kind):
    """Actual point is numerator/sqrt(radicand), in a named norm section."""
    if type(kind) is not str or kind not in ('simplex','sphere','cube'):
        raise ValueError('declared section required')
    nums = primitive(v)
    r = {'simplex':sum(nums)**2,
         'sphere':sum(x*x for x in nums),
         'cube':max(nums)**2}[kind]
    return Section(kind,nums,r)


def audit(point):
    if type(point) is not Section or type(point.numerator) is not tuple:
        raise ValueError('typed exact Section required')
    for x in point.numerator:
        integer(x)
    integer(point.radicand,1)
    if point != encode(point.numerator,point.kind):
        raise ValueError('noncanonical or forged section')
    return point


def transition(point,kind):
    audit(point)
    return encode(point.numerator,kind)


def rational_coordinates(point):
    audit(point)
    root = isqrt(point.radicand)
    if root*root != point.radicand:
        raise ValueError('these exact spherical coordinates are irrational')
    return tuple(F(x,root) for x in point.numerator)


def decode_monitor(point):
    audit(point)
    if len(point.numerator) != 3 or point.numerator[-1] == 0:
        raise ValueError('marked last-coordinate affine monitor required')
    x,n,z = point.numerator
    x,n = F(x,z),F(n,z)
    if x.denominator != 1 or n.denominator != 1 or x <= 0 or n <= 0:
        raise ValueError('monitor is not a positive integral source pair')
    return int(x),int(n)


def marker_square(point):
    """Exact squared last coordinate; the chart convention marks that coordinate."""
    audit(point)
    return F(point.numerator[-1]**2,point.radicand)


def height_from_floor(eta):
    """Uniform integer height consequence, conditional on a marker-square floor."""
    eta = rational(eta)
    if not 0 < eta <= 1:
        raise ValueError('positive exact squared-marker floor at most one required')
    return isqrt(eta.denominator//eta.numerator)


def marked_height_bound(point,eta):
    """Check the floor at this one monitor; this is not a whole-orbit assertion."""
    bound = height_from_floor(eta)
    x,n = decode_monitor(point)
    if marker_square(point) < eta:
        raise ValueError('the supplied marker floor does not hold at this point')
    if max(x,n) > bound:
        raise ArithmeticError('norm-section height bound failed')
    return x,n,bound


def sign_guard(point,coefficients):
    """Sign of a homogeneous affine inequality in the marked source chart."""
    audit(point)
    coefficients = vector(coefficients)
    if len(coefficients) != len(point.numerator):
        raise ValueError('guard and point dimensions differ')
    v = sum(a*b for a,b in zip(coefficients,point.numerator))
    return (v > 0)-(v < 0)


def boundary_labels(point):
    audit(point)
    zeros = tuple(i for i,x in enumerate(point.numerator) if x == 0)
    if point.kind == 'cube':
        maximum = max(point.numerator)
        return zeros,tuple(i for i,x in enumerate(point.numerator) if x == maximum)
    return zeros,()


def matrix_word(word):
    if type(word) is not tuple:
        raise ValueError('tuple valuation word required')
    P,Q,B = 1,1,0
    for a in word:
        integer(a,1)
        P,Q,B = 3*P,(1<<a)*Q,3*B+Q
    return ((P,0,B),(0,Q,0),(0,0,Q))


def act(matrix,point):
    audit(point)
    if (type(matrix) is not tuple or len(matrix) != len(point.numerator)
        or any(type(row) is not tuple or len(row) != len(point.numerator) for row in matrix)):
        raise ValueError('exact square matrix in the same dimension required')
    matrix = tuple(tuple(rational(a) for a in row) for row in matrix)
    image = tuple(sum(a*b for a,b in zip(row,point.numerator)) for row in matrix)
    return encode(image,point.kind)


def guarded_step(point,valuation):
    integer(valuation,1)
    x,n = decode_monitor(point)
    if x % 2 == 0 or n % 2 == 0 or x == 1:
        raise ValueError('nonterminal odd current and odd original source required')
    y = 3*x+1
    if (y & -y).bit_length()-1 != valuation:
        raise ValueError('valuation guard fails at the original marked current')
    result = act(matrix_word((valuation,)),point)
    if decode_monitor(result) != (y>>valuation,n):
        raise ArithmeticError('source marker was not retained')
    return result


def guarded_word(point,word):
    matrix_word(word)  # validate before any interpretation
    initial = audit(point)
    for a in word:
        point = guarded_step(point,a)
    if point != act(matrix_word(word),initial):
        raise ArithmeticError('normalized composition is not associative')
    return point


def dot(a,b):
    return sum(x*y for x,y in zip(a,b))


def simplex(d):
    integer(d,1)
    m = d+1
    vertices = tuple(tuple(F(int(i==j))-F(1,m) for j in range(m)) for i in range(m))
    normals = tuple(tuple(-m*x for x in v) for v in vertices)
    return vertices,normals


def cuboid(widths):
    widths = vector(widths)
    if any(a <= 0 for a in widths):
        raise ValueError('positive widths required')
    vertices = tuple(tuple(a*s for a,s in zip(widths,signs))
                     for signs in product((-1,1),repeat=len(widths)))
    polar = tuple(tuple(F(sign,widths[i]) if j==i else F(0)
                        for j in range(len(widths)))
                  for i in range(len(widths)) for sign in (-1,1))
    return vertices,polar


def rejected(thunk):
    try:
        thunk()
    except (ValueError,TypeError):
        return True
    return False


def main():
    global CHECKS
    CHECKS = 0
    kinds = ('simplex','sphere','cube')
    section_cases = 0
    for v in product(range(5),repeat=3):
        if not any(v):
            continue
        for kind in kinds:
            q = encode(v,kind)
            need(encode(tuple(7*x for x in v),kind) == q,'positive scaling is the same section point')
            for other in kinds:
                need(transition(transition(q,other),kind) == q,'section transition is lossless')
            need(boundary_labels(q)[0] == tuple(i for i,x in enumerate(v) if x==0),'all zero-coordinate face labels retained')
            section_cases += 1
    for x in range(1,64,2):
        for n in range(1,64,2):
            for kind in kinds:
                q = encode((x,n,1),kind)
                need(decode_monitor(q) == (x,n),'original affine chart recovered')
                need(sign_guard(q,(1,-1,0)) == (x>n)-(x<n),'original-source payment wall retained')
                if x != 1:
                    z = 3*x+1
                    a = (z&-z).bit_length()-1
                    result = guarded_step(q,a)
                    need(decode_monitor(result) == (z>>a,n),'actual branch transported')
                    need(rejected(lambda:guarded_step(q,a+1)),'wrong native branch rejected')

    for kind in kinds:
        # Removing the unit marker aliases distinct ordinary sources.
        need(encode((3,5),kind) == encode((9,15),kind),'unmarked direction loses integer source')
        need(encode((3,5,1),kind) != encode((9,15,1),kind),'marked source repair')
        a = guarded_word(encode((151,151,1),kind),(1,1,10))
        b = guarded_word(encode((75,75,1),kind),(1,2,8))
        need(decode_monitor(a) == (1,151) and decode_monitor(b) == (1,75),'supplied factory ROOT words transported')
        need(a != b,'ROOT does not erase the immutable original source')
        tie = encode((1,1,1),kind)
        if kind == 'cube':
            need(boundary_labels(tie)[1] == (0,1,2),'cube corner needs all incident facet labels')
    need(rejected(lambda:rational_coordinates(encode((3,5,1),'sphere'))),'sqrt35 is not a rational lattice map')
    need(rational_coordinates(encode((3,4),'sphere')) == (F(3,5),F(4,5)),'rational spherical control')
    need((3*3+1)&-(3*3+1) == 2 and (3*9+1)&-(3*9+1) == 4,'same unmarked ray changes actual valuation')

    height_cases = 0
    for x in range(1,32):
        for n in range(1,32):
            M = max(x,n)
            points = tuple(encode((x,n,1),kind) for kind in kinds)
            squares = tuple(marker_square(q) for q in points)
            need(squares == (F(1,(x+n+1)**2),F(1,x*x+n*n+1),F(1,M*M)),
                 'exact marker formulas')
            need(F(1,(2*M+1)**2) <= squares[0] <= squares[1] <= squares[2],
                 'quantitative comparison of gauges')
            need(F(1,2*M*M+1) <= squares[1],'spherical lower comparison')
            for q,square in zip(points,squares):
                result = marked_height_bound(q,square)
                need(result[:2] == (x,n) and result[2] >= M,'independent floor implies finite height')
                need(rejected(lambda q=q,square=square:marked_height_bound(q,square+F(1,2))),
                     'unsupported higher floor rejected')
                height_cases += 1
    for kind in kinds:
        previous = F(1)
        for j in range(1,33):
            q = encode(((1<<j)+1,3,1),kind)
            s = marker_square(q)
            need(0 < s < previous,'positive interior marker can tend to zero')
            need(s <= F(1,((1<<j)+1)**2),'explicit escape bound, not an orbit claim')
            previous = s
        boundary = encode((1,0,0),kind)
        need(marker_square(boundary) == 0 and rejected(lambda:decode_monitor(boundary)),
             'marker-zero compactification point is not an ordinary integer monitor')
    need(height_from_floor(F(1,101)) == 10,'non-square reciprocal rounded exactly')
    need(marked_height_bound(encode((3,5,1),'cube'),F(1,25)) == (3,5,5),'sharp common height bound')

    polar_cases = 0
    for d in range(1,7):
        vertices,normals = simplex(d)
        for i,v in enumerate(vertices):
            for j,a in enumerate(normals):
                need(dot(v,a) == (-d if i==j else 1),'regular simplex reflected polar incidence')
        for values in product(range(-2,3),repeat=d):
            x = tuple(F(a,7) for a in values)
            x = x+(-sum(x),)
            need(all(dot(x,a)<=1 for a in normals) == all(t>=-F(1,d+1) for t in x),'simplex facet membership')
            polar_cases += 1
        if d >= 2:
            widths = tuple(range(1,d+1))
            cube,dual = cuboid(widths)
            need(len(cube) == 2**d and len(dual) == 2*d,'cube and crosspolytope vertex counts')
            for a in dual:
                need(max(dot(a,v) for v in cube) == 1,'cuboid dual support')
            # A generic y is in the polar iff its weighted l1 norm is <=1.
            for y in product((F(-1,3),F(0),F(1,3)),repeat=d):
                need(max(dot(y,v) for v in cube) == sum(a*abs(t) for a,t in zip(widths,y)),'exact weighted polar norm')
    # Square polar is a rotated/rescaled square: this rational map avoids sqrt2.
    square,_ = cuboid((1,1))
    rotated = {((x-y)/2,(x+y)/2) for x,y in square}
    need(rotated == {(F(1),F(0)),(F(-1),F(0)),(F(0),F(1)),(F(0),F(-1))},'dimension-two exception')
    need(2**3 != 2*3,'cube is not combinatorially selfdual in dimension three')
    # Inverse-transpose transport preserves pairings, but not the standard lattice.
    need(dot((2,3),(F(1,2),1)) == dot((1,3),(1,1)),'contragredient pairing')
    need(F(1,2).denominator != 1,'real duality needs an integer-lattice sidecar')

    hostile = (
        lambda:encode((True,1,1),'sphere'),
        lambda:encode((1.,1,1),'simplex'),
        lambda:encode((0,0,0),'cube'),
        lambda:encode((-1,1,1),'sphere'),
        lambda:audit(replace(encode((1,1,1),'sphere'),radicand=True)),
        lambda:audit(replace(encode((1,1,1),'cube'),numerator=(True,1,1))),
        lambda:decode_monitor(encode((1,2,0),'simplex')),
        lambda:decode_monitor(encode((1,2,3),'simplex')),
        lambda:guarded_step(encode((3,5,1),'sphere'),True),
        lambda:guarded_word(encode((1,3,1),'cube'),(2,)),
        lambda:guarded_word(encode((3,5,1),'simplex'),(1.0,)),
        lambda:sign_guard(encode((3,5,1),'sphere'),(1,-1)),
        lambda:height_from_floor(True),
        lambda:height_from_floor(0.1),
        lambda:height_from_floor(F(0)),
        lambda:height_from_floor(F(-1,2)),
        lambda:height_from_floor(F(2)),
    )
    for thunk in hostile:
        need(rejected(thunk),'malformed chart or guard rejected')
    print('Three gauges: positive simplex / spherical octant / cubical patches, with exact marked transition.')
    print(f'Boundary section controls: {section_cases}; polar simplex grid controls: {polar_cases}.')
    print('Marked (current,source,1) retains every tested native valuation and payment sign.')
    print('Unmarked (3,5)=(9,15) as a direction, but actual current valuations differ: 1 versus 2.')
    print('Factory 151 and75 retain separate source markers after their supplied ROOT words.')
    print('Polar checks: regular simplex reflected scaling; cube versus crosspolytope; square exception.')
    print(f'Marker-height controls: {height_cases}; squared floor eta gives max(current,source) <= isqrt(floor(1/eta)).')
    print('Interior markers can approach zero; a pointwise chart membership is not a uniform orbit-height bound.')
    print('Hostiles: boundary ties, irrational sphere coordinates, nonunimodular dual transport, malformed packets.')
    print(f'ALL {CHECKS} CHECKS PASSED')
    print('Scope: exact representation transport and classical polarity; no new ROOT coverage or uniform orbit bound.')


if __name__ == '__main__':
    main()
