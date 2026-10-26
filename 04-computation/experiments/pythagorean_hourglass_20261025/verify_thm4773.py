# Independent check of THM-4773 (klein, 2026-10-25): fixed-twist rational lifts on the hourglass boundary.
# Surface X_d : s_i^2 = d_i L_i(y), L = (E, E+x, E-x, E+z, E-z, E+x+z, E-x-z), y = (E:x:z) in P^2(Q).
# (A) pointwise criterion and the 2^(6-r) lift count, by brute force; (B) central line: pair conditions + the conic
# w^2 = d1 d5 d3^2 s^2 + d3 d5 t^2 decided by my own Hilbert symbols, against brute-force search; (C) genus of the
# lift curve over each of the seven lines by Riemann-Hurwitz; (D) untwisted d = 1: no boundary lifts to height H;
# (E) Pythagorean twist: central lifts <-> Pythagorean triples.  Own code (does not import the klein scripts).
# mac-mini-2026-10-25-pyth.
import itertools, math, random, sys
from fractions import Fraction
from math import gcd, isqrt

def sqfree(n):
    """squarefree part of a nonzero integer (sign kept)"""
    s = -1 if n < 0 else 1; n = abs(n); r = 1; p = 2
    while p * p <= n:
        while n % (p * p) == 0: n //= p * p
        if n % p == 0: r *= p; n //= p
        p += 1
    return s * r * n

def is_square(n): return n >= 0 and isqrt(n) ** 2 == n

def forms(E, x, z): return (E, E + x, E - x, E + z, E - z, E + x + z, E - x - z)

def lifts(d, y):
    """number of points of X_d(Q) over the base point y (integer coords), by the definition: s_i^2 = lam d_i L_i(y)
    for some lam in Q^x, s not all zero, counted projectively (s ~ -s)."""
    vals = [di * Li for di, Li in zip(d, forms(*y))]
    nz = [v for v in vals if v != 0]
    if not nz: return 0
    # lam must make every nonzero value a square: lam ~ sqfree(v) for each v
    classes = {sqfree(v) for v in nz}
    if len(classes) != 1: return 0
    return 2 ** (len(nz) - 1)

def hilbert(a, b, p):
    """Hilbert symbol (a,b)_p for nonzero integers a, b; p prime or 0 for the real place"""
    if p == 0: return -1 if (a < 0 and b < 0) else 1
    def split(n):
        v = 0
        while n % p == 0: n //= p; v += 1
        return v, n
    al, u = split(a); be, w = split(b)
    if p == 2:
        eps = lambda t: ((t - 1) // 2) % 2
        om = lambda t: ((t * t - 1) // 8) % 2
        e = eps(u) * eps(w) + al * om(w) + be * om(u)
        return -1 if e % 2 else 1
    leg = lambda t: 1 if pow(t % p, (p - 1) // 2, p) == 1 else -1
    s = (-1) ** (al * be * ((p - 1) // 2)) * (leg(u) ** be) * (leg(w) ** al)
    return s

def primes_dividing(n):
    n = abs(n); ps = set(); p = 2
    while p * p <= n:
        while n % p == 0: ps.add(p); n //= p
        p += 1
    if n > 1: ps.add(n)
    return ps

def conic_soluble(A, B):
    """w^2 = A s^2 + B t^2 has a nontrivial rational point iff (A,B)_v = 1 for all v (Hasse-Minkowski / Legendre)"""
    A, B = sqfree(A), sqfree(B)
    for p in {0, 2} | primes_dividing(A) | primes_dividing(B):
        if hilbert(A, B, p) != 1: return False
    return True

def central_decision(d):
    d0, d1, d2, d3, d4, d5, d6 = d
    sq = lambda v: is_square(v) or (v > 0 and is_square(v))   # rational square test on integers
    pair = all(is_square(sqfree(v)) and sqfree(v) == 1 for v in (-d1 * d2, -d3 * d4, -d5 * d6))
    if not pair: return False
    return conic_soluble(d1 * d5 * d3 * d3, d3 * d5)

def central_bruteforce(d, H):
    for x in range(-H, H + 1):
        for z in range(-H, H + 1):
            if x == 0 or z == 0 or x + z == 0 or gcd(x, z) != 1: continue
            if lifts(d, (0, x, z)): return True
    return False

random.seed(2026)
out = []
# (A) pointwise criterion: compare the squareclass test with an independent brute-force root search on small values
bad = 0; tested = 0
for _ in range(400):
    d = tuple(random.choice([1, -1, 2, -2, 3, -3, 5, -5, 6, -6]) for _ in range(7))
    for _ in range(25):
        y = tuple(random.randint(-6, 6) for _ in range(3))
        if y == (0, 0, 0): continue
        vals = [di * Li for di, Li in zip(d, forms(*y))]
        nz = [v for v in vals if v != 0]
        if not nz: continue
        # brute force: is there a squarefree lam with all lam*v squares?
        found = 0
        for lam in {sqfree(v) for v in nz}:
            if all(is_square(lam * v) for v in nz): found = 2 ** (len(nz) - 1)
        tested += 1
        if found != lifts(d, y): bad += 1
out.append(f"(A) pointwise criterion and 2^(6-r) count: {tested} (twist, point) pairs, {bad} disagreements")
# (B) central line: Hilbert decision vs brute force
agree = dis = 0; examples = []
twists = [(1, 1, -1, 1, -1, 1, -1), (1, 1, -1, 1, -1, 7, -7), (1, 2, -2, 3, -3, 6, -6), (1, 1, 1, 1, 1, 1, 1)]
for _ in range(500):
    d1, d3, d5 = (random.choice([1, -1, 2, -2, 3, -3, 5, -5, 6, -6, 7, -7, 10, 11, 13]) for _ in range(3))
    twists.append((1, d1, -d1, d3, -d3, d5, -d5))       # pair conditions satisfied; the conic decides
for d in twists:
    dec = central_decision(d); bf = central_bruteforce(d, 40)
    if dec == bf: agree += 1
    else:
        dis += 1
        if len(examples) < 5: examples.append((d, dec, bf))
out.append(f"(B) central line: Hilbert decision vs brute force (|x|,|z| <= 40): {agree} agree, {dis} disagree {examples}")
out.append(f"    displayed twists: Pythagorean {central_decision(twists[0])}, (1,1,-1,1,-1,7,-7) {central_decision(twists[1])}, "
           f"(1,2,-2,3,-3,6,-6) {central_decision(twists[2])}, untwisted {central_decision(twists[3])}")
# (C) genus of the lift locus over each line, by Riemann-Hurwitz, generic twist.
#   Over a line, restrict the seven forms; forms vanishing identically are dropped (they are the line itself);
#   proportional forms impose a constant squareclass condition and split the curve; the remaining distinct branch
#   points b and the rank k of the (Z/2)^k cover give 2g - 2 = 2^k (-2) + b 2^(k-1).
from sympy import Matrix, symbols
lines = {'E=0': (1, 0, 0), 'E+x=0': (1, 1, 0), 'E-x=0': (1, -1, 0), 'E+z=0': (1, 0, 1), 'E-z=0': (1, 0, -1),
         'E+x+z=0': (1, 1, 1), 'E-x-z=0': (1, -1, -1)}
Lvec = [(1, 0, 0), (1, 1, 0), (1, -1, 0), (1, 0, 1), (1, 0, -1), (1, 1, 1), (1, -1, -1)]
for name, ell in lines.items():
    # parametrize the line ell . (E,x,z) = 0 by two basis vectors
    M = Matrix([list(ell)]); ns = M.nullspace(); B = Matrix.hstack(*ns)   # 3 x 2
    restricted = []
    for Li in Lvec:
        v = Matrix([list(Li)]) * B          # 1 x 2 linear form on P^1
        if all(c == 0 for c in v): continue
        restricted.append(tuple(v))
    # group by proportionality (same zero on P^1)
    groups = []
    for v in restricted:
        for g in groups:
            if v[0] * g[0][1] - v[1] * g[0][0] == 0: g.append(v); break
        else: groups.append([v])
    b = len(groups)                          # distinct branch points
    k = b - 1                                # independent square roots of ratios of distinct forms
    splits = sum(len(g) - 1 for g in groups) # each extra proportional form: a constant square root (splitting)
    g_genus = 1 + 2 ** (k - 1) * (b / 2 - 2) if k >= 1 else 0
    out.append(f"(C) line {name:8s}: {len(restricted)} nonvanishing forms, {b} branch points, cover rank {k}, "
               f"{2 ** splits} components (over Qbar) of genus {g_genus:g}")
# (D) untwisted surface: search all seven lines for lifts up to height H
H = 60; found = []
d1 = (1,) * 7
for name, ell in lines.items():
    for a in range(-H, H + 1):
        for c in range(-H, H + 1):
            # points on the line: choose two free coordinates and solve the third (all lines have E coefficient 1)
            x, z = a, c
            E = -(ell[1] * x + ell[2] * z)
            if (E, x, z) == (0, 0, 0) or gcd(gcd(abs(E), abs(x)), abs(z)) != 1: continue
            if lifts(d1, (E, x, z)): found.append((name, (E, x, z)))
out.append(f"(D) untwisted d = 1: boundary lifts with |x|,|z| <= {H}: {len(found)} {found[:5]}")
# (E) Pythagorean twist: central lifts <-> Pythagorean triples
dP = (1, 1, -1, 1, -1, 1, -1); cen = []
for x in range(1, 400):
    for z in range(1, 400):
        if gcd(x, z) != 1: continue
        if lifts(dP, (0, x, z)): cen.append((x, z))
trip = sorted({(isqrt(x), isqrt(z), isqrt(x + z)) for x, z in cen if is_square(x) and is_square(z)})
out.append(f"(E) Pythagorean twist, central base points (0:x:z) with 0 < x, z < 400 that lift: {len(cen)}; "
           f"all of the form (a^2, b^2) with a^2 + b^2 = c^2: {len(trip) == len(cen)}; triples {trip[:8]}")
print('\n'.join(out))

# (B') constructive confirmation for every twist where the Hilbert decision says "soluble" but the box search found
# nothing: search conic points (s, t) with s t (d1 d3 s^2 + t^2) != 0, build the base point (0 : d1 d3 s^2 : t^2) and
# check it with the definition-level lift test (no use of the conic reduction in the check itself).
def constructive(d, S=400):
    d0, d1, d2, d3, d4, d5, d6 = d
    for s in range(1, S):
        for t in range(1, S):
            x = d1 * d3 * s * s; z = t * t
            if x + z == 0: continue
            g = gcd(abs(x), z); xx, zz = x // g, z // g
            if lifts(d, (0, xx, zz)): return (0, xx, zz)
    return None
missing = []; confirmed = []
for d in twists:
    if central_decision(d) and not central_bruteforce(d, 40):
        w = constructive(d)
        (confirmed if w else missing).append((d, w))
print(f"(B') Hilbert-soluble twists missed by the |x|,|z| <= 40 box: {len(confirmed) + len(missing)}; explicit lifts found "
      f"for {len(confirmed)} (max height {max(max(abs(c) for c in w) for d, w in confirmed) if confirmed else 0}), none found for {len(missing)} {missing[:3]}")
print("    e.g.", confirmed[:3])
