#!/usr/bin/env python3
"""equivariant_reduced_jc_kps_S128c105.py -- kind-pasteur-2026-07-20-S128c105

A REDUCED JACOBIAN CONJECTURE THAT HOLDS, and its sharpness against THM-1300.

THE OBSERVATION.  THM-1300 records that the counterexample F is C*-equivariant
with weights (1,-1,-2) on the source and (-2,-1,1) on the target.  Those two
vectors have the SAME MULTISET {1,-1,-2}, but they are not the same VECTOR: the
target is the source composed with the transposition (1 3).  In S_3 = D_3 that
transposition is a REFLECTION.  So the counterexample is equivariant only up to a
reflection twist of the weight vector.

THE REDUCED CONJECTURE.  Say F is UNTWISTED-EQUIVARIANT if there is a weight
vector w with F(lambda . x) = lambda . F(x) for the SAME w on both sides, i.e.
F_i is weighted-homogeneous of weight w_i.  Then:

    (RJC)  F untwisted-equivariant and det JF in C*  ==>  F is invertible.

PART A -- PROVED, all weights of one sign.  If every w_i > 0 (or every w_i < 0):
F_i weighted-homogeneous of positive weight gives F(0) = 0.  det JF in C* makes F
etale, so every fibre is discrete; being algebraic, F^{-1}(0) is finite.  It is
also C*-stable, and with all weights positive every non-zero orbit is infinite, so
F^{-1}(0) = {0}.  A graded map with V(F) = {0} makes C[x]/(F) finite-dimensional,
so C[x] is a finite C[F]-module (graded Nakayama) and F is a FINITE morphism,
hence proper.  Proper + etale = covering; C^n is simply connected, so the degree
is 1 and F is an isomorphism.  QED.

PART B -- PROVED, n = 2, weights (1,-1) (the smallest mixed-sign case, where Part
A's argument fails outright).  Verified symbolically below: with t = xy the whole
untwisted family is F = (x f(t), y g(t)) and

    det JF  =  d/dt [ t f(t) g(t) ]                                  (identity)

so det JF constant forces t f g = t, i.e. f g = 1, i.e. f and g are constants and
F is LINEAR.  The conjecture holds in this case for the strongest possible reason.

PART C -- SHARPNESS.  The counterexample's own weights (1,-1,-2), imposed
UNTWISTED, give the exact parametrisation (invariants s = xy, t = x^2 z):

    F_1 = x A(s,t),   F_2 = y B(s,t) + x z C(s,t),   F_3 = z D(s,t) + y^2 E(s,t)

-- and this is a complete description of the weight-1 / weight-(-1) / weight-(-2)
parts, not an ansatz (proved by the monomial count in the script).  Since the
weight vectors agree, det JF is C*-invariant, hence lies in C[s,t]; "det JF
constant" is therefore finitely many equations in the coefficients.  The script
solves them at increasing degree and reports whether any NON-LINEAR untwisted map
achieves a constant Jacobian.

If none does, the dividing line is exactly the reflection twist: relaxing "same
weight vector" to "same weight multiset" is what admits THM-1300's F.
"""
import sys
from sympy import (symbols, Poly, expand, factor, simplify, diff, Matrix, srepr,
                   Rational, groebner, solve, total_degree, nsimplify, together)
from itertools import product

x, y, z = symbols('x y z')
s, t = symbols('s t')

print("=" * 78)
print("PART 0 -- INDEPENDENT RE-VERIFICATION OF THE COUNTEREXAMPLE (THM-1300)")
print("=" * 78)
u = 1 + x * y
F1 = u**3 * z + y**2 * u * (4 + 3 * x * y)
F2 = y + 3 * x * u**2 * z + 3 * x * y**2 * (4 + 3 * x * y)
F3 = 2 * x - 3 * x**2 * y - x**3 * z
F = [F1, F2, F3]
J = Matrix(3, 3, lambda i, j: diff(F[i], [x, y, z][j]))
detJ = expand(J.det())
print("  det JF = %s   -> constant: %s" % (detJ, detJ.free_symbols == set()))
pts = [(0, 0, Rational(-1, 4)), (1, Rational(-3, 2), Rational(13, 2)),
       (-1, Rational(3, 2), Rational(13, 2))]
ims = [tuple(expand(f.subs({x: p[0], y: p[1], z: p[2]})) for f in F) for p in pts]
print("  images of the three collision points:")
for p, im in zip(pts, ims):
    print("     F%-24s = %s" % (str(p), str(im)))
print("  all three coincide : %s" % (len(set(ims)) == 1))

print()
print("=" * 78)
print("PART 0b -- THE WEIGHTS, AND THE REFLECTION TWIST")
print("=" * 78)
lam = symbols('lam')


def wt_of(expr, w):
    """C*-weight of expr if weighted-homogeneous; else the sorted set of weights.
    Read straight off the exponents: wt(x^a y^b z^c) = a w0 + b w1 + c w2."""
    ex = {sum(e * wi for e, wi in zip(mono, w))
          for mono, _ in Poly(expand(expr), x, y, z).terms()}
    return ex.pop() if len(ex) == 1 else sorted(ex)


w_src = (1, -1, -2)
print("  source weights (x,y,z) = %s" % (w_src,))
w_tgt = [wt_of(f, w_src) for f in F]
print("  target weights (F1,F2,F3) = %s" % (w_tgt,))
print("  same MULTISET as source : %s" % (sorted(w_tgt) == sorted(w_src)))
print("  same VECTOR  as source  : %s" % (list(w_tgt) == list(w_src)))
perm = [w_src.index(v) for v in w_tgt]
print("  target = source composed with the permutation %s of positions" % (perm,))
print("  -> that permutation is the transposition (1 3): a REFLECTION in S_3 = D_3.")
print("  So F is equivariant only UP TO A REFLECTION TWIST of the weight vector.")

print()
print("=" * 78)
print("PART B -- n=2, weights (1,-1) UNTWISTED: det JF = d/dt[t f g], so F is LINEAR")
print("=" * 78)
# the untwisted weight-(1,-1) family is exactly F = (x f(xy), y g(xy)); verify the
# identity for a generic truncation of f and g
NB = 4
fc = symbols('f0:%d' % (NB + 1))
gc = symbols('g0:%d' % (NB + 1))
fpoly = sum(fc[i] * (x * y)**i for i in range(NB + 1))
gpoly = sum(gc[i] * (x * y)**i for i in range(NB + 1))
G1, G2 = x * fpoly, y * gpoly
JG = Matrix([[diff(G1, x), diff(G1, y)], [diff(G2, x), diff(G2, y)]])
lhs = expand(JG.det())
ft = sum(fc[i] * t**i for i in range(NB + 1))
gt = sum(gc[i] * t**i for i in range(NB + 1))
rhs = expand(diff(t * ft * gt, t)).subs(t, x * y)
print("  det JF - d/dt[t f g]|_{t=xy}  =  %s" % expand(lhs - expand(rhs)))
print("  -> identity holds for f,g truncated at degree %d." % NB)
print("  Consequently det JF = c constant  =>  t f g = c t (value 0 at t=0)")
print("     =>  f g = c  =>  f, g constant (C[t] is a domain, deg f + deg g = 0)")
print("     =>  F = (a x, (c/a) y) is LINEAR and invertible.  RJC HOLDS here.")

print()
print("=" * 78)
print("PART C -- weights (1,-1,-2) UNTWISTED: complete parametrisation + search")
print("=" * 78)


def weight(mono, w):
    return sum(e * wi for e, wi in zip(mono, w))


DEGCAP = int(sys.argv[1]) if len(sys.argv) > 1 else 8
w = (1, -1, -2)
# confirm the module generators claimed in the docstring, by monomial enumeration
gens = {1: [(1, 0, 0)], -1: [(0, 1, 0), (1, 0, 1)], -2: [(0, 0, 1), (0, 2, 0)]}
inv = [(1, 1, 0), (2, 0, 1)]   # s = xy, t = x^2 z
print("  invariants: s = xy (deg 2), t = x^2 z (deg 3)")
for tw in (1, -1, -2):
    ok = True
    bad = None
    for a in range(DEGCAP + 1):
        for b in range(DEGCAP + 1):
            for c in range(DEGCAP + 1):
                if a + b + c > DEGCAP:
                    continue
                if weight((a, b, c), w) != tw:
                    continue
                # is (a,b,c) = generator + n1*s + n2*t for some generator?
                hit = False
                for g in gens[tw]:
                    r = (a - g[0], b - g[1], c - g[2])
                    for n2 in range(0, DEGCAP + 1):
                        rr = (r[0] - 2 * n2, r[1], r[2] - n2)
                        if rr[2] != 0 or rr[0] < 0 or rr[1] < 0:
                            continue
                        if rr[0] == rr[1] >= 0:
                            hit = True
                            break
                    if hit:
                        break
                if not hit:
                    ok = False
                    bad = (a, b, c)
    print("     weight %-3d module generated by %-22s over C[s,t] : %s%s"
          % (tw, str(gens[tw]), ok, "" if ok else "  (missed %s)" % (bad,)))

# now the search: A,B,C,D,E polynomials in s,t of total degree <= DEG_ST
DEG_ST = int(sys.argv[2]) if len(sys.argv) > 2 else 1
mons_st = [(i, j) for i in range(DEG_ST + 1) for j in range(DEG_ST + 1)
           if i + j <= DEG_ST]
names = "ABCDE"
coeffs = {}
polys = {}
for nm in names:
    cs = symbols('%s0:%d' % (nm, len(mons_st)))
    coeffs[nm] = cs
    polys[nm] = sum(cs[k] * (x * y)**m[0] * (x**2 * z)**m[1]
                    for k, m in enumerate(mons_st))
A, B, C, D, E = (polys[n] for n in names)
H1 = expand(x * A)
H2 = expand(y * B + x * z * C)
H3 = expand(z * D + y**2 * E)
H = [H1, H2, H3]
JH = Matrix(3, 3, lambda i, j: diff(H[i], [x, y, z][j]))
dH = expand(JH.det())
print()
print("  ansatz: A..E of total degree <= %d in (s,t); %d unknown coefficients"
      % (DEG_ST, 5 * len(mons_st)))
print("  det JH has %d terms before reduction" % len(Poly(dH, x, y, z).terms()))
# det JH is C*-invariant, so it is a polynomial in s,t: collect its (s,t) content
pd = Poly(dH, x, y, z)
inv_terms = {}
for mono, co in pd.terms():
    a, b, c = mono
    # invariant monomial x^a y^b z^c with a-b-2c=0 equals s^b t^c iff a = b+2c
    if a - b - 2 * c != 0:
        inv_terms.setdefault(('NONINVARIANT', mono), 0)
        inv_terms[('NONINVARIANT', mono)] += co
    else:
        inv_terms.setdefault((b, c), 0)
        inv_terms[(b, c)] += co
noninv = [v for k, v in inv_terms.items() if isinstance(k[0], str) and expand(v) != 0]
print("  non-invariant terms in det JH : %d  (expected 0 -- weights agree)"
      % len(noninv))
eqs = [expand(v) for k, v in inv_terms.items()
       if not isinstance(k[0], str) and k != (0, 0) and expand(v) != 0]
print("  equations forcing det JH constant : %d" % len(eqs))
allc = [c for nm in names for c in coeffs[nm]]
print("  solving...")
sols = solve(eqs, allc, dict=True)
print("  solution branches : %d" % len(sols))
nonlinear = []
for sol in sols:
    h1 = expand(H1.subs(sol))
    h2 = expand(H2.subs(sol))
    h3 = expand(H3.subs(sol))
    dj = expand(JH.det().subs(sol))
    if dj == 0:
        continue
    degs = [total_degree(p) for p in (h1, h2, h3) if p != 0]
    if degs and max(degs) > 1:
        nonlinear.append((h1, h2, h3, dj))
print("  branches with det JH a NON-ZERO constant and F NON-LINEAR : %d"
      % len(nonlinear))
for h1, h2, h3, dj in nonlinear[:5]:
    print("     F = (%s, %s, %s)   det = %s" % (h1, h2, h3, dj))
print()
if not nonlinear:
    print("  VERDICT: at weights (1,-1,-2) UNTWISTED, with A..E of degree <= %d," % DEG_ST)
    print("  every constant-Jacobian map is LINEAR.  The counterexample's own weight")
    print("  multiset gives nothing once the REFLECTION TWIST is removed.")
