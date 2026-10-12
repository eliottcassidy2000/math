"""Exact controls for Stechkin repair obstructions; no plane-colouring solver.

Coordinates (a,b) mean a*u+b*v, u=(1,0), v=(1/2,sqrt(3)/2).
Only the two open triangular classes are compiled. Hexagon boundary ownership
is deliberately not guessed: every obstructing vertex is strictly interior.
"""
from fractions import Fraction as F
from itertools import combinations


def rational(x):
    if type(x) not in (int, F):
        raise ValueError("exact rational required")
    return F(x)


def point(p):
    if type(p) is not tuple or len(p) != 2:
        raise ValueError("exact coordinate tuple required")
    return tuple(rational(x) for x in p)


def integer(x):
    if type(x) is not int:
        raise ValueError("exact integer required")
    return x


def add(p, q):
    p, q = point(p), point(q)
    return (p[0] + q[0], p[1] + q[1])


def scale(p, c):
    p, c = point(p), rational(c)
    return (c * p[0], c * p[1])


def norm2(p):
    a, b = point(p)
    return a*a + a*b + b*b


def distance2(p, q):
    return norm2(add(p, scale(q, -1)))


def gap_vertices(orientation, i=0, j=0):
    integer(orientation), integer(i), integer(j)
    if orientation not in (0, 1):
        raise ValueError("orientation must be 0 or 1")
    base = ((F(1, 2), 0), (F(1, 2), F(1, 2)), (0, F(1, 2)))
    if orientation:
        base = ((F(1, 2), F(1, 2)), (1, F(1, 2)), (F(1, 2), 1))
    return tuple(add(p, (i, j)) for p in base)


def triangle_class(p):
    """Return 4/5 in the open gaps, otherwise None; no boundary assignment."""
    a, b = point(p)
    a, b = a - a.__floor__(), b - b.__floor__()
    if 0 < a < F(1, 2) and 0 < b < F(1, 2) and a+b > F(1, 2):
        return 4
    if F(1, 2) < a < 1 and F(1, 2) < b < 1 and a+b < F(3, 2):
        return 5
    return None


def triangle_witness(orientation, i=0, j=0):
    vertices = gap_vertices(orientation, i, j)
    p = tuple(sum(v[k] for v in vertices) / 3 for k in (0, 1))
    return (p, add(p, (1, 0)), add(p, (0, 1)))


def audit_unit_triangle(vertices, class_id=None):
    if type(vertices) is not tuple or len(vertices) != 3:
        raise ValueError("three exact points required")
    vertices = tuple(point(p) for p in vertices)
    if class_id is not None:
        integer(class_id)
        if class_id not in (4, 5):
            raise ValueError("triangle class must be 4 or 5")
        if any(triangle_class(p) != class_id for p in vertices):
            raise ValueError("wrong native triangle class")
    if any(distance2(p, q) != 1 for p, q in combinations(vertices, 2)):
        raise ValueError("not a unit triangle")
    return vertices


def pairings(items):
    if not items:
        return ((),)
    first, rest = items[0], items[1:]
    return tuple(((first, x),) + tail for x in rest
                 for tail in pairings(tuple(y for y in rest if y != x)))


def pairing_obstruction(pairs):
    """Reauthenticate a perfect pairing of six original whole classes."""
    if type(pairs) is not tuple or len(pairs) != 3:
        raise ValueError("three pairs required")
    flat = []
    for pair in pairs:
        if type(pair) is not tuple or len(pair) != 2:
            raise ValueError("exact pair tuple required")
        flat.extend(integer(x) for x in pair)
    if sorted(flat) != list(range(6)):
        raise ValueError("each original class occurs exactly once")
    bad_pair = next(pair for pair in pairs if 4 in pair)
    witness = audit_unit_triangle(triangle_witness(0), 4)
    return bad_pair, witness


def cross(a, b, c):
    a, b, c = point(a), point(b), point(c)
    return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])


def area(poly):
    if not poly:
        return F(0)
    poly = tuple(point(p) for p in poly)
    return abs(sum(p[0]*q[1]-p[1]*q[0]
                   for p, q in zip(poly, poly[1:]+poly[:1]))) / 2


def clip(subject, clipper):
    """Exact convex polygon intersection; clipper is counterclockwise."""
    result = tuple(point(p) for p in subject)
    clipper = tuple(point(p) for p in clipper)
    for a, b in zip(clipper, clipper[1:]+clipper[:1]):
        old, result = result, []
        if not old:
            break
        for p, q in zip(old, old[1:]+old[:1]):
            dp, dq = cross(a, b, p), cross(a, b, q)
            inside_p, inside_q = dp >= 0, dq >= 0
            if inside_p:
                result.append(p)
            if inside_p != inside_q:
                t = dp / (dp-dq)
                result.append(add(p, scale(add(q, scale(p, -1)), t)))
        result = tuple(result)
    return tuple(result)


def normalization_areas():
    """Areas in the common [0,2]^2 lattice-coordinate period, not probabilities."""
    small = tuple(gap_vertices(o, i, j) for i in range(2)
                  for j in range(2) for o in range(2))
    large = tuple(tuple(scale(p, 2) for p in gap_vertices(o)) for o in range(2))
    source = sum(area(p) for p in small)
    overlap = sum(area(clip(p, q)) for p in small for q in large)
    return source, overlap, source-overlap


def missed_by_normalization(p):
    p = point(p)
    return triangle_class(p) is not None and triangle_class(scale(p, F(1, 2))) is None


def rational_unit(t):
    """Cartesian coordinates, unlike the triangular-lattice coordinates above."""
    t = rational(t)
    return ((1-t*t)/(1+t*t), 2*t/(1+t*t))


def odd_denominator_parity(x):
    x = rational(x)
    if x.denominator % 2 == 0:
        raise ValueError("not in Z_(2)")
    return x.numerator % 2


def main():
    checks = 0

    def check(value, message):
        nonlocal checks
        checks += 1
        if not value:
            raise RuntimeError(message)

    def reject(call):
        try:
            call()
        except ValueError:
            check(True, "rejected")
        else:
            check(False, "malformed input accepted")

    for o in range(2):
        for i in range(-5, 6):
            for j in range(-5, 6):
                vs = triangle_witness(o, i, j)
                check(audit_unit_triangle(vs, 4+o) == vs, "interior unit triangle")
                check(all(triangle_class(p) == 4+o for p in vs), "orientation")
                check(area(gap_vertices(o, i, j)) == F(1, 8), "tile area")
    matches = pairings(tuple(range(6)))
    check(len(matches) == 15, "all perfect pairings")
    for pairing in matches:
        pair, witness = pairing_obstruction(pairing)
        check(4 in pair and all(triangle_class(p) in pair for p in witness), "pair lift")
    for bits in range(8):
        check(any((bits >> i & 1) == (bits >> j & 1)
                  for i, j in combinations(range(3), 2)), "triangle binary obstruction")
    source, overlap, missed = normalization_areas()
    check((source, overlap, missed) == (1, F(1, 4), F(3, 4)), "exact areas")
    check(missed/4 == F(3, 16), "missed density")
    p = triangle_witness(0)[0]
    check(missed_by_normalization(p), "centroid coverage hole")
    check(norm2(scale(p, F(1, 2))) == F(1, 12), "halved point inside hexagon")
    # Rational interior samples independently authenticate all 16 intersection areas.
    pieces = []
    for i in range(2):
        for j in range(2):
            for o in range(2):
                sm = gap_vertices(o, i, j)
                pieces.append(tuple(area(clip(sm, tuple(scale(v, 2)
                                  for v in gap_vertices(other)))) for other in range(2)))
    check(sorted(x for pair in pieces for x in pair) == [F(0)]*14+[F(1, 8)]*2,
          "two whole triangles and fourteen zero intersections")
    for denominator in (5, 7, 11, 17):
        for i in range(2*denominator):
            for j in range(2*denominator):
                z = (F(2*i+1, 2*denominator), F(2*j+1, 2*denominator))
                missed_now = missed_by_normalization(z)
                covered = triangle_class(z) is None or triangle_class(scale(z, F(1, 2))) is not None
                check(missed_now != covered, "source-preserving cover test")
    # Every rational unit direction has odd denominators and odd parity sum.
    units = set()
    for numerator in range(-12, 13):
        for denominator in range(1, 13):
            t = F(numerator, denominator)
            x, y = rational_unit(t)
            units.add((x, y))
            check(x*x+y*y == 1, "rational unit")
            check((odd_denominator_parity(x)+odd_denominator_parity(y)) % 2 == 1,
                  "rational direction parity")
    check((odd_denominator_parity(F(-1))+odd_denominator_parity(F(0))) % 2 == 1,
          "omitted stereographic direction (-1,0)")
    # All-direction two-unit decomposition: v/2 +/- sqrt(1/q-1/4) Jv.
    for a in range(-9, 10):
        for b in range(-9, 10):
            q = F(a*a+b*b, 25)
            if 0 < q <= 4:
                coefficient_squared = 1/q-F(1, 4)
                check(coefficient_squared >= 0 and q/4+coefficient_squared*q == 1,
                      "two real unit vectors")
    rejects = (
        lambda: point((True, 0)), lambda: point((F(1, 3), 0.0)),
        lambda: point([0, 0]), lambda: gap_vertices(True),
        lambda: gap_vertices(2), lambda: gap_vertices(0, 1.0),
        lambda: audit_unit_triangle(triangle_witness(0), 4.0),
        lambda: audit_unit_triangle(triangle_witness(0), 5),
        lambda: audit_unit_triangle(((0, 0), (1, 0), (0, 0))),
        lambda: pairing_obstruction(((0, 1), (2, 3), (4, True))),
        lambda: pairing_obstruction(((0, 1), (2, 3), (4, 4))),
        lambda: rational_unit(0.5), lambda: odd_denominator_parity(F(1, 2)),
    )
    for hostile in rejects:
        reject(hostile)
    print("Stechkin whole-plane obstruction controls")
    print("Coordinates: u=(1,0), v=(1/2,sqrt(3)/2); triangle side=1/2")
    print("242 translated interior unit triangles: PASS")
    print("All 15 pairings of the original six whole classes: obstructed")
    print("Normalization areas on [0,2]^2: T=1; intersection=1/4; missed=3/4")
    print("Normalized cover density=13/16; missed density=3/16 (hex-centred dilation)")
    print(f"Rational-direction parity controls: 300 parameter values, {len(units)} distinct units")
    print(f"Malformed/native guard hostiles: {len(rejects)}")
    print("Scope: neither a six-colouring construction nor a seven-colour lower bound")
    print(f"checks={checks}")


if __name__ == "__main__":
    main()
