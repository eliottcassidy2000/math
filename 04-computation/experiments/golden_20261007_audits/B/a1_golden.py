#!/usr/bin/env python3
"""Independent audit of THM-4592 (golden reading of standard-Collatz periodic points).

Own conventions (deliberately different code path from golden_cycle_dictionary.py):
  * C(x) = 3x+1 (x odd), x/2 (x even) on Z_(2) cap Q.
  * periodic point from word via  C^n(x) = (3^k x + B)/2^e,  B <- 3B + 2^e on odd, e <- e+1 on even.
  * kappa_n(w) = sum_j w_j phi^(n-1-j) in Z[phi] = Z^2 (coords of 1, phi), reduced to a canonical
    representative modulo the ideal (phi^n - 1) via a Hermite normal form computed here.
  * equivariance is checked against the word of C(x) recomputed from the point itself.
"""
from fractions import Fraction
from itertools import product
from math import gcd

def egcd(a, b):
    if b == 0:
        return (abs(a), (1 if a >= 0 else -1), 0)
    g, x, y = egcd(b, a % b)
    return (g, y, x - (a // b) * y)

def fibpair(n):
    # phi^n = F_{n-1} + F_n phi  (n >= 1), phi^0 = 1
    if n == 0:
        return (1, 0)
    a, b = 0, 1  # F0, F1
    for _ in range(n - 1):
        a, b = b, a + b
    return (a, b)  # (F_{n-1}, F_n)

def mulphi(x):
    # (m + k phi) * phi = k + (m + k) phi
    m, k = x
    return (k, m + k)

def mul(x, y):
    a, b = x; c, d = y
    return (a * c + b * d, a * d + b * c + b * d)

class Quot:
    """Z[phi]/(g) with g = a + b phi, via HNF of lattice spanned by g, g*phi."""
    def __init__(self, g):
        a, b = g
        v1 = (a, b); v2 = mulphi(g)          # (b, a+b)
        G, u, v = egcd(v1[0], v2[0])          # u*v1[0] + v*v2[0] = G
        r1 = (u * v1[0] + v * v2[0], u * v1[1] + v * v2[1])
        assert r1[0] == G
        if G != 0:
            r2 = ((v2[0] // G) * v1[0] - (v1[0] // G) * v2[0], (v2[0] // G) * v1[1] - (v1[0] // G) * v2[1])
        else:
            raise ValueError
        assert r2[0] == 0
        d2 = abs(r2[1])
        self.d1, self.d2 = G, d2
        self.x = r1[1] % d2
        self.N = G * d2
    def red(self, z):
        m, k = z
        q = m // self.d1
        m -= q * self.d1; k -= q * self.x
        return (m, k % self.d2)

def golden_words(n):
    out = []
    for w in product((0, 1), repeat=n):
        if all(not (w[i] and w[(i + 1) % n]) for i in range(n)):
            out.append(w)
    return out

def point(w):
    B, e, k = 0, 0, 0
    for b in w:
        if b:
            B = 3 * B + 2 ** e; k += 1
        else:
            e += 1
    return Fraction(B, 2 ** e - 3 ** k)

def C(x):
    assert x.denominator % 2 == 1
    return 3 * x + 1 if x.numerator % 2 else x / 2

def word_of(x, n):
    w = []
    for _ in range(n):
        w.append(x.numerator % 2); x = C(x)
    return tuple(w), x

def kappa(w):
    n = len(w); s = (0, 0)
    for j, b in enumerate(w):
        if b:
            p = fibpair(n - 1 - j); s = (s[0] + p[0], s[1] + p[1])
    return s

def lucas(n):
    a, b = 2, 1
    for _ in range(n):
        a, b = b, a + b
    return a

print("=== Statement 1/2: equivariance and bijectivity, n = 1..22 ===")
for n in range(1, 23):
    g = fibpair(n); g = (g[0] - 1, g[1])
    Q = Quot(g)
    W = golden_words(n)
    assert len(W) == lucas(n), (n, len(W), lucas(n))
    assert Q.N == lucas(n) - 1 - (-1) ** n or n == 1, (n, Q.N)
    cls = {}
    equi_ok = True
    for w in W:
        x = point(w)
        ww, back = word_of(x, n)
        assert ww == w and back == x, (w, x)
        c = Q.red(kappa(w))
        cls.setdefault(c, []).append(x)
        # equivariance using the point C(x): its word recomputed from the orbit
        y = C(x)
        wy, _ = word_of(y, n)
        if Q.red(kappa(wy)) != Q.red(mulphi(c)):
            equi_ok = False
    fib = sorted((sorted(v) for v in cls.values() if len(v) > 1))
    nonzero_bad = [v for v in fib if v != [Fraction(-2), Fraction(-1), Fraction(0)]]
    print(f"n={n:2d} L_n={lucas(n):6d} |R_n|={Q.N:6d} classes={len(cls):6d} "
          f"surj={len(cls)==Q.N} collapsed={[list(map(str,v)) for v in fib]} equivariant={equi_ok} other_collisions={len(nonzero_bad)}")

print("\n=== Statement 3: n = 5 labels in F_11 (phi -> 4) ===")
H = {1, 3, 4, 5, 9}
negH = {(-h) % 11 for h in H}
g5 = fibpair(5); Q5 = Quot((g5[0] - 1, g5[1]))
print("R_5 HNF d1,d2 =", Q5.d1, Q5.d2, " phi^5-1 =", (g5[0] - 1, g5[1]),
      " phi^3*(4-phi) =", mul(fibpair(3), (4, -1)))
lab = {}
for w in golden_words(5):
    x = point(w); c = kappa(w)
    l = (c[0] + 4 * c[1]) % 11
    lab[x] = l
for start in (Fraction(1, 13), Fraction(-5)):
    orb = []; x = start
    for _ in range(5):
        orb.append(x); x = C(x)
    print("cycle (C-order):", [str(t) for t in orb], "labels:", [lab[t] for t in orb],
          "label set == H:", {lab[t] for t in orb} == H, " == -H:", {lab[t] for t in orb} == negH)
print("table order 1/13,2/13,4/13,8/13,16/13 ->", [lab[Fraction(a, 13)] for a in (1, 2, 4, 8, 16)])
print("table order -5,-14,-7,-20,-10 ->", [lab[Fraction(a)] for a in (-5, -14, -7, -20, -10)])
print("label(0) =", lab[Fraction(0)], "; 4 has order", next(k for k in range(1, 11) if pow(4, k, 11) == 1))
print("C acts as x4 on labels:", all(lab[C(x)] == 4 * lab[x] % 11 for x in lab))
print("4^4 mod 11 =", pow(4, 4, 11), "; <3> =", sorted({pow(3, i, 11) for i in range(10)}),
      "; <4> =", sorted({pow(4, i, 11) for i in range(10)}), "; <5> =", sorted({pow(5, i, 11) for i in range(10)}))

print("\n=== Statement 4: G_5 ===")
def gen_group(gens, p):
    idt = tuple(range(p)); G = {idt}; fr = [idt]
    while fr:
        new = []
        for s in fr:
            for g in gens:
                t = tuple(g[s[i]] for i in range(p))
                if t not in G:
                    G.add(t); new.append(t)
        fr = new
    return G
p = 11
S = tuple((i + 1) % p for i in range(p))
G_phi = gen_group([S, tuple(4 * i % p for i in range(p))], p)
G_3 = gen_group([S, tuple(3 * i % p for i in range(p))], p)
Aff_H = {tuple((a * i + b) % p for i in range(p)) for a in H for b in range(p)}
print("|<x+1,4x>| =", len(G_phi), " |<x+1,3x>| =", len(G_3), " equal:", G_phi == G_3, " == {ax+b: a in H}:", G_phi == Aff_H)
pairs = [frozenset((i, j)) for i in range(p) for j in range(i + 1, p)]
orb = {frozenset((g[0], g[1])) for g in G_phi}
stab = [g for g in G_phi if frozenset((g[0], g[1])) == frozenset((0, 1))]
print("orbit of {0,1}:", len(orb), "of", len(pairs), " stabilizer size:", len(stab))
arcs = {(g[0], g[1]) for g in G_phi}
paley = {(x, y) for x in range(p) for y in range(p) if (y - x) % p in H}
print("orbit of arc (0,1) == Paley QR_11 arcs:", arcs == paley, " |arcs| =", len(arcs))
# Borel of PSL(2,11): matrices [[a,b],[0,d]], ad=1, acting x -> (a x + b)/d = a^2 x + a b, modulo +-1
borel = set()
for a in range(1, p):
    d = pow(a, -1, p)
    for b in range(p):
        borel.add(tuple((a * x + b) * d % p for x in range(p)))
print("Borel(PSL(2,11)) on F_11 == G_5:", borel == G_phi, " |Borel| =", len(borel))

print("\n=== Statement 5: n = 10 ===")
g10 = fibpair(10); g10 = (g10[0] - 1, g10[1])
Q10 = Quot(g10)
print("phi^10-1 =", g10, " = 11*phi^5 =", tuple(11 * t for t in fibpair(5)), " HNF d1,d2,x =", Q10.d1, Q10.d2, Q10.x)
def ab(c):
    return ((c[0] + 4 * c[1]) % 11, (c[0] + 8 * c[1]) % 11)
W10 = golden_words(10)
pts = {}
for w in W10:
    x = point(w); pts[x] = ab(kappa(w))
def period(x):
    y = C(x); t = 1
    while y != x:
        y = C(y); t += 1
    return t
prim = [x for x in pts if period(x) == 10]
print("#points period|10 =", len(pts), " #primitive period 10 =", len(prim))
cells_prim = [pts[x] for x in prim]
print("primitive points -> distinct cells:", len(set(cells_prim)) == 110, " all beta != 0:", all(c[1] != 0 for c in cells_prim))
beta0 = {x: pts[x] for x in pts if pts[x][1] == 0}
print("points in beta=0 cells:", sorted(str(x) for x in beta0), "count", len(beta0))
# C^5 directly
ok = True
for x in prim:
    y = x
    for _ in range(5):
        y = C(y)
    al, be = pts[x]
    if pts[y] != (al, (-be) % 11):
        ok = False
print("C^5 acts as (a,b)->(a,-b) on the 110 primitive points (direct iteration):", ok)
# period-5 points: cell (2*label, 0)?
twist = all(pts[x] == ((2 * lab[x]) % 11, 0) for x in lab)
print("period-5 point with R_5-label l sits in R_10 cell (2l, 0):", twist)
print("  e.g. 1/13 (label 3, in H) -> cell", pts[Fraction(1, 13)], "; paper projection (a,b)->a gives base cell",
      pts[Fraction(1, 13)][0], "in H?", pts[Fraction(1, 13)][0] in H)
# Psi_- bijection on the 55 C^5-orbits
orbits = set()
for x in prim:
    y = x
    for _ in range(5):
        y = C(y)
    orbits.add(frozenset((x, y)))
imgs = set()
for o in orbits:
    x = next(iter(o)); al, be = pts[x]; bi = pow(be, -1, 11)
    imgs.add(frozenset(((al + bi) % 11, (al - bi) % 11)))
print("#C^5-orbits =", len(orbits), " Psi_- images distinct pairs =", len(imgs))

print("\n=== Remark: R_18 and the -17 cycle ===")
g18 = fibpair(18); g18 = (g18[0] - 1, g18[1])
print("phi^18-1 =", g18, " 76*phi^9 =", tuple(76 * t for t in fibpair(9)), " L_9 =", lucas(9))
w, back = word_of(Fraction(-17), 18)
print("-17: C^18(-17) == -17:", back == -17, " minimal period:", period(Fraction(-17)), " odd steps:", sum(w), " even:", 18 - sum(w))
print("19 splits in Z[phi]: x^2-x-1 roots mod 19:", [r for r in range(19) if (r * r - r - 1) % 19 == 0])
Q18 = Quot(g18); print("|R_18| =", Q18.N, "= 76^2:", Q18.N == 76 ** 2, " HNF:", Q18.d1, Q18.d2)
