"""Finite exact witness interfaces; no CH, HN, JC, or D19 solution claimed."""
from dataclasses import dataclass, replace
from fractions import Fraction as F
from itertools import product, permutations
from math import gcd, isqrt
import json

import plane_algebraic_witness_20261011 as plane


def need(ok, message):
    if not ok:
        raise ValueError(message)


def integer(x):
    need(type(x) is int, "exact integer required")
    return x


def rational(x):
    need(type(x) in (int, F), "exact rational required")
    return F(x)


def corner_norms(x, y, width, height):
    for z in (x, y, width, height):
        integer(z)
    need(width > 0 and height > 0, "positive side lengths")
    return (x*x+y*y, (x-width)**2+y*y,
            x*x+(y-height)**2, (x-width)**2+(y-height)**2)


@dataclass(frozen=True)
class RectangleReceipt:
    x: int
    y: int
    width: int
    height: int
    distances: tuple


def verify_rectangle(packet, *, inside=False):
    need(type(packet) is RectangleReceipt and type(inside) is bool, "typed rectangle")
    norms = corner_norms(packet.x, packet.y, packet.width, packet.height)
    need(type(packet.distances) is tuple and len(packet.distances) == 4,
         "four marked distances")
    need(all(type(a) is int and a >= 0 for a in packet.distances),
         "nonnegative exact distances")
    need(tuple(a*a for a in packet.distances) == norms, "all four exact equations")
    if inside:
        need(0 < packet.x < packet.width and 0 < packet.y < packet.height,
             "strict inside guard")
    return norms


def verify_square(packet, *, inside=False, primitive=True):
    need(type(primitive) is bool, "boolean normalization flag")
    norms = verify_rectangle(packet, inside=inside)
    need(packet.width == packet.height, "same source-owned side in both directions")
    if primitive:
        need(gcd(gcd(packet.x, packet.y), packet.width) == 1, "primitive source")
    return norms


def primitive_local_cell(x, y, side, modulus=24):
    for z in (x, y, side, modulus):
        integer(z)
    need(modulus > 0, "positive modulus")
    residues = {i*i % modulus for i in range(modulus)}
    norms = (x*x+y*y, (x-side)**2+y*y,
             x*x+(y-side)**2, (x-side)**2+(y-side)**2)
    return gcd(gcd(x, y), gcd(side, modulus)) == 1 and all(
        z % modulus in residues for z in norms)


def mod24_description(x, y, side):
    for z in (x, y, side):
        integer(z)
    parity = (x % 4 == 0 and y % 2 == 1) or (y % 4 == 0 and x % 2 == 1)
    ternary = (x % 3 == 0) != (y % 3 == 0)
    return side % 12 == 0 and parity and ternary


@dataclass(frozen=True)
class Polynomial:
    n: int
    terms: tuple


def poly(n, terms):
    integer(n)
    need(n >= 1, "positive dimension")
    out = {}
    for exponent, coefficient in terms:
        need(type(exponent) is tuple and len(exponent) == n and
             all(type(e) is int and e >= 0 for e in exponent), "exact multiindex")
        c = rational(coefficient)
        out[exponent] = out.get(exponent, F(0)) + c
    return Polynomial(n, tuple((a, c) for a, c in sorted(out.items()) if c))


def checked(p):
    need(type(p) is Polynomial and type(p.terms) is tuple, "polynomial packet")
    need(poly(p.n, p.terms) == p, "canonical nonzero sorted polynomial terms")
    return p


def constant(n, c):
    return poly(n, [((0,)*n, c)])


def variable(n, j):
    integer(j)
    need(0 <= j < n, "marked variable")
    return poly(n, [(tuple(int(i == j) for i in range(n)), 1)])


def add(p, q):
    checked(p); checked(q)
    need(p.n == q.n, "same source ring")
    return poly(p.n, p.terms + q.terms)


def scale(p, c):
    checked(p); c = rational(c)
    return poly(p.n, [(a, c*b) for a, b in p.terms])


def mul(p, q):
    checked(p); checked(q)
    need(p.n == q.n, "same source ring")
    return poly(p.n, [(tuple(x+y for x, y in zip(a, b)), c*d)
                      for a, c in p.terms for b, d in q.terms])


def power(p, k):
    checked(p); integer(k)
    need(k >= 0, "nonnegative power")
    out = constant(p.n, 1)
    for _ in range(k):
        out = mul(out, p)
    return out


def derivative(p, j):
    checked(p); integer(j)
    need(0 <= j < p.n, "derivative variable")
    return poly(p.n, [(tuple(e-int(i == j) for i, e in enumerate(a)), c*a[j])
                      for a, c in p.terms if a[j]])


def evaluated(p, point):
    checked(p)
    need(type(point) is tuple and len(point) == p.n, "marked rational point")
    point = tuple(rational(v) for v in point)
    ans = F(0)
    for a, c in p.terms:
        term = c
        for v, e in zip(point, a):
            term *= v**e
        ans += term
    return ans


def map_checked(mapping):
    need(type(mapping) is tuple and len(mapping) >= 1, "polynomial self-map")
    for p in mapping:
        checked(p)
        need(p.n == len(mapping), "self-map source dimension")
    return mapping


def compose(p, mapping):
    checked(p); map_checked(mapping)
    need(p.n == len(mapping), "matching source variables")
    out = constant(p.n, 0)
    for a, c in p.terms:
        term = constant(p.n, c)
        for q, e in zip(mapping, a):
            term = mul(term, power(q, e))
        out = add(out, term)
    return out


def jacobian(mapping):
    map_checked(mapping)
    n = len(mapping)
    out = constant(n, 0)
    for order in permutations(range(n)):
        sign = (-1)**sum(order[i] > order[j] for i in range(n) for j in range(i+1,n))
        term = constant(n, sign)
        for i, j in enumerate(order):
            term = mul(term, derivative(mapping[i], j))
        out = add(out, term)
    return out


@dataclass(frozen=True)
class CollisionReceipt:
    characteristic: int
    mapping: tuple
    left: tuple
    right: tuple
    coordinate: int
    inverse_difference: F


def verify_collision(packet):
    need(type(packet) is CollisionReceipt, "collision receipt")
    integer(packet.characteristic)
    need(packet.characteristic == 0, "this verifier is characteristic zero only")
    mapping = map_checked(packet.mapping)
    need(jacobian(mapping) == constant(len(mapping), 1), "normalized Keller identity")
    j = integer(packet.coordinate)
    need(0 <= j < len(mapping), "marked differing coordinate")
    left = tuple(evaluated(p, packet.left) for p in mapping)
    right = tuple(evaluated(p, packet.right) for p in mapping)
    need(left == right, "same image")
    u = rational(packet.inverse_difference)
    need(u*(rational(packet.left[j])-rational(packet.right[j])) == 1,
         "nonzero-coordinate witness")
    return left


def verify_inverse(mapping, inverse):
    map_checked(mapping); map_checked(inverse)
    need(len(mapping) == len(inverse), "same inverse source dimension")
    n = len(mapping)
    need(jacobian(mapping) == constant(n, 1), "normalized Keller map")
    identity = tuple(variable(n, j) for j in range(n))
    need(tuple(compose(p, inverse) for p in mapping) == identity, "right inverse")
    need(tuple(compose(p, mapping) for p in inverse) == identity, "left inverse")
    return True


def main():
    checks = 0
    def check(ok, label):
        nonlocal checks
        need(ok, label); checks += 1
    def rejects(fn):
        nonlocal checks
        try:
            fn()
        except (ValueError, TypeError):
            checks += 1
            return
        raise ValueError("hostile accepted")

    cells = []
    for x, y, side in product(range(24), repeat=3):
        actual = primitive_local_cell(x, y, side)
        check(actual == mod24_description(x,y,side), "complete primitive local ledger")
        if actual:
            cells.append((x,y,side))
            check(all(a % 24 == 1 for a in (x*x+y*y,(x-side)**2+y*y,
                       x*x+(y-side)**2,(x-side)**2+(y-side)**2)), "all distances are units mod6")
    check(len(cells) == 128, "mod24 exact finite universe")
    hostile = (3,4,12)
    norms = corner_norms(*hostile, hostile[2])
    check(primitive_local_cell(*hostile) and norms == (25,97,73,145), "local is not global")
    check(tuple(isqrt(v)**2 == v for v in norms) == (True,False,False,False), "three missing rational distances")
    centre_norm = F(1,2)
    check(centre_norm == F(1,2)**2 + F(1,2)**2, "algebraic centre norm")
    check(isqrt(centre_norm.denominator)**2 != centre_norm.denominator,
          "algebraic distance need not be rational")
    rectangle = RectangleReceipt(3,4,6,8,(5,5,5,5))
    check(verify_rectangle(rectangle,inside=True) == (25,)*4, "honest rectangle positive")
    rejects(lambda: verify_square(rectangle,inside=True))
    rejects(lambda: verify_square(RectangleReceipt(3,4,12,12,(5,10,9,12))))
    rejects(lambda: verify_rectangle(replace(rectangle,x=3.0)))
    rejects(lambda: verify_rectangle(replace(rectangle,distances=(True,5,5,5))))
    rejects(lambda: primitive_local_cell(False,4,12))
    rejects(lambda: primitive_local_cell(3,4,12,0))

    X,Y = variable(2,0),variable(2,1)
    inverse_cases = 0
    for a,b in product(range(-2,3),repeat=2):
        for d,e in product(range(2,4),repeat=2):
            v = add(Y,scale(power(X,d),b))
            mapping = (add(X,scale(power(v,e),a)),v)
            u = add(X,scale(power(Y,e),-a))
            inverse = (u,add(Y,scale(power(u,d),-b)))
            check(verify_inverse(mapping,inverse), "two exact shear compositions")
            inverse_cases += 1
    same = CollisionReceipt(0,(X,Y),(0,0),(0,0),0,F(1))
    rejects(lambda: verify_collision(same))
    fake = CollisionReceipt(0,(power(X,2),Y),(-1,0),(1,0),0,F(-1,2))
    rejects(lambda: verify_collision(fake))
    pointwise = CollisionReceipt(0,(add(X,power(X,2)),Y),(0,0),(-1,0),0,F(1))
    check(evaluated(jacobian(pointwise.mapping),(0,0)) == 1, "pointwise determinant hostile")
    rejects(lambda: verify_collision(pointwise))
    rejects(lambda: verify_inverse((X,Y),(Y,X)))
    rejects(lambda: verify_collision(replace(same,characteristic=False)))
    rejects(lambda: poly(2,[((True,0),1)]))
    rejects(lambda: evaluated(X,(1.0,2)))
    prime_hostiles=[]
    for p in (2,3,5,7,11):
        t = variable(1,0)
        f = add(t,scale(power(t,p),-1))
        df = derivative(f,0)
        check(all((c.numerator % p == (1 if a == (0,) else 0))
                  for a,c in df.terms), "characteristic-p derivative1")
        check(evaluated(f,(0,)) == evaluated(f,(1,)) == 0, "characteristic-p collision")
        check(df != constant(1,1), "same integer map is not char0 Keller")
        packet = CollisionReceipt(p,(f,),(0,),(1,),0,F(-1))
        rejects(lambda packet=packet: verify_collision(packet))
        prime_hostiles.append(p)

    points = plane.moser_points()
    keep, spindle = plane.extract_algebraic((64,0,-28,0,1),5,6,points,4)
    check(keep == tuple(range(7)) and plane.verify_algebraic(spindle) == 139,
          "inherited exact finite algebraic witness")
    check(len(spindle.critical.edges) == 11, "actual induced spindle geometry")
    rejects(lambda: plane.verify_algebraic(replace(spindle,points=points[:-1]+(points[0],))))
    rejects(lambda: plane.verify_algebraic(replace(spindle,critical=replace(spindle.critical,colours=6))))

    return {"status":"FINITE-EXACT; no finite HN6, JC collision or D19 square witness supplied",
            "checks":checks,"mod24":{"residue_triples":24**3,"accepted":len(cells),
                "side_counts":{str(l):sum(c[2]==l for c in cells) for l in (0,12)},
                "local_not_global":{"source":[3,4,12],"norms":list(norms)}},
            "rectangle_control":{"sides":[6,8],"point":[3,4],"distances":[5]*4},
            "positive_inverse_maps":inverse_cases,"characteristic_hostiles":prime_hostiles,
            "algebraic_graph_control":{"vertices":7,"edges":11,"chromatic":4,
                                         "refutation_nodes":139},
            "set_theory_scope":"absoluteness is a mathematical theorem, not a finite computation"}


if __name__ == "__main__":
    print(json.dumps(main(),indent=2,sort_keys=True))
