#!/usr/bin/env python3
"""Rank-drop locus of the arc ideal in the one-parameter families of HYP-9298 (audit2).
q=10: f1 = x y^2 - r x^4, f2 = y^3 - x^3 y   (r = beta/gamma; VA value 1/12)
q=11: f1 = y^3 - r x^3 y, f2 = x^2 y^2 - x^5  (r = d/e;      VA value 15/2)
At the first weight W0 where the generic HS falls below the product (12 resp. 14), the product needs a rank drop, so
the set of r with HS_r = product through any weight >= W0 is contained in the roots of
g(r) = gcd over blocks/projections of det(R M_block(r) C).  We compute g mod P, factor it, and test each F_P-root."""
import sys, random
from fractions import Fraction
import flint
from arc_audit import Arc, a2_product

P = 2147483647
q = int(sys.argv[1]); Wtest = int(sys.argv[2]) if len(sys.argv) > 2 else 22
W0 = {10: 12, 11: 14}[q]
def rels(r):
    if q == 10:
        return [{(1, 2): 1, (4, 0): (-r) % P}, {(0, 3): 1, (3, 1): -1}]
    return [{(0, 3): 1, (3, 1): (-r) % P}, {(2, 2): 1, (5, 0): -1}]

def block_matrix(A, W, L):
    cols = A.monos(W, L); ci = {m: i for i, m in enumerate(cols)}
    rows = []
    for (wf, F) in A.Fs:
        for n in range(0, L + 1):
            Wm, Lm = W - wf - n, L - n
            if Wm < 0 or n >= len(F): continue
            for m in A.monos(Wm, Lm):
                row = [0] * len(cols)
                for mo, co in F[n].items():
                    k = ci[tuple(sorted(mo + m))]; row[k] = (row[k] + co) % P
                rows.append(row)
    return rows, cols

random.seed(2026)
rgen = random.randrange(2, P - 1)
Agen = Arc(rels(rgen), W0, P)
polys = []
for L in range(0, W0 + 1):
    rows, cols = block_matrix(Agen, W0, L)
    if not rows or not cols: continue
    rho = flint.nmod_mat(rows, P).rank()
    if rho == 0: continue
    # evaluate det(R M(r) C) at rho+2 points for two random projections
    for trial in range(2):
        R = flint.nmod_mat([[random.randrange(P) for _ in range(len(rows))] for _ in range(rho)], P)
        C = flint.nmod_mat([[random.randrange(P) for _ in range(rho)] for _ in range(len(cols))], P)
        xs, ys = [], []
        for t in range(rho + 2):
            rv = random.randrange(P)
            At = Arc(rels(rv), W0, P)
            Mt, _ = block_matrix(At, W0, L)
            d = (R * flint.nmod_mat(Mt, P) * C).det()
            xs.append(rv); ys.append(int(d))
        # Lagrange interpolation via flint
        poly = flint.nmod_poly([0], P)
        for i in range(len(xs)):
            num = flint.nmod_poly([1], P); den = 1
            for j in range(len(xs)):
                if j != i:
                    num *= flint.nmod_poly([(-xs[j]) % P, 1], P); den = den * (xs[i] - xs[j]) % P
            poly += num * (ys[i] * pow(den, P - 2, P) % P)
        polys.append((L, trial, poly))
# per block gcd of the two projections
byL = {}
for L, t, p in polys:
    byL.setdefault(L, []).append(p)
allroots = set()
for L, ps in sorted(byL.items()):
    g = ps[0].gcd(ps[1]) if len(ps) > 1 else ps[0]
    if g.degree() <= 0:
        continue
    fac = g.factor()
    desc = []
    for f, e in fac[1]:
        if f.degree() == 1:
            root = (-int(f.coeffs()[0])) * pow(int(f.coeffs()[1]), P - 2, P) % P
            allroots.add(root)
            desc.append(('root', root, e))
        else:
            desc.append(('irreducible factor of degree', f.degree(), e))
    print(f'q={q} W={W0} block L={L}: rank-drop polynomial degree {g.degree()}: {desc}')
def ratrec(a):
    import math
    r0, r1, s0, s1 = P, a % P, 0, 1
    B = math.isqrt(P // 2)
    while r1 > B:
        qq = r0 // r1; r0, r1 = r1, r0 - qq * r1; s0, s1 = s1, s0 - qq * s1
    return Fraction(r1, s1)
ref = a2_product(q, Wtest)
for root in sorted(allroots):
    A = Arc(rels(root), Wtest, P)
    hs = A.hs_list(Wtest)
    diff = [w for w in range(Wtest + 1) if hs[w] != ref[w]]
    print(f'  root r = {ratrec(root)} (mod P {root}):', 'HS = product through %d' % Wtest if not diff else
          'differs first at %d (HS %d vs product %d)' % (diff[0], hs[diff[0]], ref[diff[0]]))
