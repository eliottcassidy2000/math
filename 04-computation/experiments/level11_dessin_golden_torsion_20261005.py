"""Exact level-11 dessin, marked Paley chart, and golden E[3] clock.

Standard library only. No files are written. Every check survives python -O.
The note separates cited modular identifications from elementary proofs.
"""
from collections import deque
from itertools import product


def require(value, message):
    if not value:
        raise AssertionError(message)


def matmul(a, b, modulus=None):
    c = (a[0]*b[0]+a[1]*b[2], a[0]*b[1]+a[1]*b[3],
         a[2]*b[0]+a[3]*b[2], a[2]*b[1]+a[3]*b[3])
    return tuple(x % modulus for x in c) if modulus else c


I = (1, 0, 0, 1)


def matpow(a, k, modulus=None):
    out = I
    while k:
        if k & 1:
            out = matmul(out, a, modulus)
        a = matmul(a, a, modulus)
        k //= 2
    return out


def matrix_order(a, modulus, limit):
    b = I
    for k in range(1, limit+1):
        b = matmul(b, a, modulus)
        if b == I:
            return k
    raise ValueError("order exceeds declared bound")


def canonical(a, p):
    a = tuple(x % p for x in a)
    return min(a, tuple(-x % p for x in a))


def mobius(a, z, p):
    # The integer p denotes infinity; finite coordinates are 0,...,p-1.
    if z == p:
        return p if a[2] == 0 else a[0]*pow(a[2], -1, p) % p
    bottom = (a[2]*z+a[3]) % p
    return p if bottom == 0 else (a[0]*z+a[1])*pow(bottom, -1, p) % p


def cycles(perm):
    unseen = set(range(len(perm)))
    out = []
    while unseen:
        x = min(unseen)
        cycle = []
        while x in unseen:
            unseen.remove(x)
            cycle.append(x)
            x = perm[x]
        require(x == cycle[0], "input is not a permutation")
        out.append(tuple(cycle))
    return out


def dessin_controls():
    p = 11
    s = canonical((0, -1, 1, 0), p)
    t = canonical((1, 1, 0, 1), p)
    r = canonical(matmul(s, t, p), p)
    generators = (s, t)
    group = {canonical(I, p)}
    pending = deque(group)
    while pending:
        g = pending.popleft()
        for h in generators:
            v = canonical(matmul(g, h, p), p)
            if v not in group:
                group.add(v)
                pending.append(v)
    all_matrices = {canonical(a, p) for a in product(range(p), repeat=4)
                    if (a[0]*a[3]-a[1]*a[2]) % p == 1}
    require(group == all_matrices and len(group) == 660, "PSL2 closure")
    perms = {g: tuple(mobius(g, x, p) for x in range(p+1)) for g in group}
    require(len(set(perms.values())) == 660, "projective action faithful")
    passport = tuple(tuple(sorted(map(len, cycles(perms[g])))) for g in (s, r, t))
    require(passport == ((2,)*6, (3,)*4, (1, 11)), "dessin passport")
    require(6+4+2-12 == 0, "genus one")
    pair_orbit = {(v[p], v[0]) for v in perms.values()}
    require(len(pair_orbit) == 132, "two-transitivity")
    ordered = sorted(group)
    index = {g: i for i, g in enumerate(ordered)}
    regular_cycles = []
    for h in (s, r, t):
        reg = tuple(index[canonical(matmul(g, h, p), p)] for g in ordered)
        regular_cycles.append(len(cycles(reg)))
    require(regular_cycles == [330, 220, 60], "regular dessin passport")
    require(2-sum(regular_cycles)+660 == 52, "regular cover genus26")
    borel = [g for g in ordered if mobius(g, p, p) == p]
    require(len(borel) == 55, "marked stabilizer size")
    squares = {x*x % p for x in range(1, p)}
    arcs = {(x, y) for x in range(p) for y in range(p)
            if (y-x) % p in squares}
    require(len(arcs) == 55 and not any((b, a) in arcs for a, b in arcs), "tournament")
    for g in borel:
        require({(mobius(g, a, p), mobius(g, b, p)) for a, b in arcs} == arcs,
                "Borel preserves marked Paley tournament")
    require((1, 2) in arcs and (mobius(s, 1, p), mobius(s, 2, p)) == (10, 5)
            and (10, 5) not in arcs, "unmarked inversion hostile")
    print("DESSIN PSL2(F11)=660; P1=12; passport=2^6,3^4,11*1; genus=1")
    print("REGULAR COVER: degree660; cycles330,220,60; genus26; quotient degree55")
    print("MARKED PALEY: Borel55; 55 arcs; 3025 arc-image checks; inversion reverses1->2")
    print("DESSIN GENERATOR CYCLES (11 denotes infinity):", passport,
          tuple(cycles(perms[g]) for g in (s, r, t)))


# Exact binary field F_256 in the polynomial basis modulo X^8+X^4+X^3+X+1.
def fm(a, b):
    out = 0
    while b:
        if b & 1:
            out ^= a
        a <<= 1
        if a & 256:
            a ^= 0x11B
        b >>= 1
    return out


def fp(a, n):
    out = 1
    while n:
        if n & 1:
            out = fm(out, a)
        a = fm(a, a)
        n //= 2
    return out


def curve(point):
    if point is None:
        return True
    x, y = point
    return fm(y, y) ^ y == fm(fm(x, x), x) ^ fm(x, x)


def add(a, b):
    if a is None:
        return b
    if b is None:
        return a
    x, y = a
    u, v = b
    if x == u and y != v:
        return None
    slope = fm(x, x) if a == b else fm(y ^ v, fp(x ^ u, 254))
    xx = fm(slope, slope) ^ 1 ^ x ^ u
    yy = fm(slope, xx ^ x) ^ y ^ 1
    return xx, yy


def scalar(k, point):
    out = None
    for _ in range(k):
        out = add(out, point)
    return out


def frobenius(point):
    return None if point is None else (fm(point[0], point[0]), fm(point[1], point[1]))


def torsion_controls():
    require(all(fp(a, 255) == 1 for a in range(1, 256)), "binary modulus is a field")
    points = {(x, y) for x in range(256) for y in range(256) if curve((x, y))}
    require(len(points) == 224, "independent F256 point count")
    torsion = {p for p in points if add(add(p, p), p) is None}
    quartic = {x for x in range(256) if fp(x, 4) ^ x ^ 1 == 0}
    require(len(torsion) == 8 and len(quartic) == 4, "3-torsion sizes")
    require(torsion == {p for p in points if p[0] in quartic}, "tangent division criterion")
    require(all(fp(x, 16) == x for x in quartic), "x coordinates in F16")
    require(not any(fp(y, 16) == y for _, y in torsion), "full points not in F16")
    p = min(torsion)
    q = frobenius(p)
    register = {(a, b): add(scalar(a, p), scalar(b, q)) for a, b in product(range(3), repeat=2)}
    require(set(register.values()) == torsion | {None}, "basis P,Frob(P)")
    for (a, b), point in register.items():
        require(frobenius(point) == register[(b, (a+b) % 3)], "golden intertwiner")
    for a, b in product(register, repeat=2):
        require(add(register[a], register[b]) == register[((a[0]+b[0]) % 3, (a[1]+b[1]) % 3)],
                "additive basis map")
    orbit = []
    z = p
    while z not in orbit:
        orbit.append(z)
        z = frobenius(z)
    require(len(orbit) == 8 and z == p and set(orbit) == torsion, "all nonzero torsion one orbit")
    for z in torsion:
        y = z
        for _ in range(4):
            y = frobenius(y)
        require(y == (z[0], z[1] ^ 1), "fourth Frobenius is negation")
    print("TORSION: exhaustive65536 affine pairs overF256; E(F256)=225; E[3]=9")
    print("TORSION: P,Frob(P) basis=", p, q, "; 81 additive controls; golden matrix on all9 points")
    print("TORSION ORBIT:", orbit)
    print("TORSION COORDINATES: x^4+x+1 irreducible overF2; four x inF16; eight points requireF256")


def clock_controls():
    a = (-2, -2, 1, 0)
    m = (0, 1, 1, 1)
    swap = (0, 1, 1, 0)
    require(matpow(a, 4) == (-4, 0, 0, -4), "integral local Frobenius relation")
    require(matmul(matmul(swap, m, 3), swap, 3) == tuple(x % 3 for x in a), "mod3 conjugacy")
    require(matrix_order(a, 3, 8) == matrix_order(m, 3, 8) == 8, "eight-clock")
    require(matpow(m, 4, 3) == (2, 0, 0, 2), "projective four-clock")
    require(matrix_order(a, 9, 24) == matrix_order(m, 9, 24) == 24, "same order at mod9")
    require((a[0]+a[3]) % 9 != (m[0]+m[3]) % 9, "no mod9 similarity: trace obstruction")
    require((a[0]*a[3]-a[1]*a[2]) % 9 != (m[0]*m[3]-m[1]*m[2]) % 9,
            "no mod9 similarity: determinant obstruction")
    # Every binary quartic with constant1 is checked, independently of field embedding.
    def remainder(poly, divisor):
        while poly.bit_length() >= divisor.bit_length():
            poly ^= divisor << (poly.bit_length()-divisor.bit_length())
        return poly
    require(all(remainder(0b10011, f) for f in (0b10, 0b11, 0b111)), "quartic irreducibility")
    print("CLOCK: A^4=-4I; A mod3 is swap-conjugate to goldenM; both order8")
    print("HOSTILE mod9: both order24, but traces7/1 and determinants2/8 prevent conjugacy")


def main():
    dessin_controls()
    clock_controls()
    torsion_controls()
    print("PASS: exact finite controls; modular identification is separately cited; no Collatz guard transferred")


if __name__ == "__main__":
    main()
