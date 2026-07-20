#!/usr/bin/env python3
"""
death-star-2026-07-19-S59n (HYP-8080) -- the equivariant family computation.

Ansatz (fully general for weights (1,-1,-2), E linear in s):
  F1 = z A(t) + y^2 B(t),  F2 = y C(t) + x z D(t),  F3 = x (E0(t) - E1(t) s)
Work in Q[x,y,z]; det J computed exactly; det is automatically weight-0, i.e.
a polynomial in t = xy, s = x^2 z.  We:
  (1) regression: given functions => det == -2;
  (2) express det symbolically over unknown coefficients (multivariate dict
      engine, variables = (x,y,z) + unknown coeffs) and extract the constraint
      system det - const == 0;
  (3) tangent space at the known solution: nullity of the linearization
      (the infinitesimal deformation space WITHIN the ansatz);
  (4) attempt exact continuation: solve the constraints with v (the cube root
      of A) deformed, to produce NEW exact counterexample candidates.
"""
from fractions import Fraction as Fr
import itertools, sys

NV = None  # set later: number of variables = 3 + n_unknowns

def pmul(a, b):
    r = {}
    for ka, ca in a.items():
        for kb, cb in b.items():
            k = tuple(p+q for p, q in zip(ka, kb))
            v = r.get(k, 0) + ca*cb
            if v: r[k] = v
            elif k in r: del r[k]
    return r
def padd(*ps):
    r = {}
    for p in ps:
        for k, c in p.items():
            v = r.get(k, 0) + c
            if v: r[k] = v
            elif k in r: del r[k]
    return r
def psc(p, s): return {k: c*s for k, c in p.items() if c*s != 0}
def var(i):
    k = [0]*NV; k[i] = 1
    return {tuple(k): 1}
def pdiff(p, i):
    r = {}
    for k, c in p.items():
        if k[i] > 0:
            k2 = list(k); k2[i] -= 1
            r[tuple(k2)] = c*k[i]
    return r

def build_F(coeffs_ABCDE):
    """coeffs: dict name -> list of engine-polys (coefficients of t^j).
       Returns F1,F2,F3 as engine polys in x,y,z (+unknown vars)."""
    A, B, C, D, E0, E1 = (coeffs_ABCDE[n] for n in "A B C D E0 E1".split())
    x, y, z = var(0), var(1), var(2)
    t = pmul(x, y)
    s = pmul(pmul(x, x), z)
    def poly_t(cl):
        out = {}
        tp = {tuple([0]*NV): 1}
        for c in cl:
            out = padd(out, pmul(c, tp))
            tp = pmul(tp, t)
        return out
    F1 = padd(pmul(z, poly_t(A)), pmul(pmul(y, y), poly_t(B)))
    F2 = padd(pmul(y, poly_t(C)), pmul(pmul(x, z), poly_t(D)))
    F3 = pmul(x, padd(poly_t(E0), psc(pmul(s, poly_t(E1)), -1)))
    return F1, F2, F3

def detJ(F1, F2, F3):
    J = [[pdiff(F, i) for i in range(3)] for F in (F1, F2, F3)]
    return padd(
        pmul(J[0][0], padd(pmul(J[1][1], J[2][2]), psc(pmul(J[1][2], J[2][1]), -1))),
        psc(pmul(J[0][1], padd(pmul(J[1][0], J[2][2]), psc(pmul(J[1][2], J[2][0]), -1))), -1),
        pmul(J[0][2], padd(pmul(J[1][0], J[2][1]), psc(pmul(J[1][1], J[2][0]), -1))))

# ---------- (1) regression with the given functions ----------
NV = 3
K = lambda c: {tuple([0]*NV): c} if c else {}
given = {
    "A":  [K(1), K(3), K(3), K(1)],          # (1+t)^3
    "B":  [K(4), K(7), K(3)],                # (1+t)(4+3t)
    "C":  [K(1), K(12), K(9)],
    "D":  [K(3), K(6), K(3)],                # 3(1+t)^2
    "E0": [K(2), K(-3)],
    "E1": [K(1)],
}
F1, F2, F3 = build_F(given)
d = detJ(F1, F2, F3)
print("(1) regression: det J of six-function form =", d, "== -2:", d == {(0,0,0): -2})

# ---------- (2)+(3) symbolic: unknowns as extra variables ----------
# unknown layout: A: deg<=4 (5), B: deg<=3 (4), C: deg<=3 (4), D: deg<=3 (4),
#                 E0: deg<=2 (3), E1: deg<=1 (2)   => 22 unknowns
layout = [("A",5), ("B",4), ("C",4), ("D",4), ("E0",3), ("E1",2)]
NU = sum(n for _, n in layout)
NV = 3 + NU
names = []
coeffs = {}
i = 3
for nm, n in layout:
    coeffs[nm] = [var(i+j) for j in range(n)]
    names += [f"{nm}{j}" for j in range(n)]
    i += n
F1, F2, F3 = build_F(coeffs)
D_sym = detJ(F1, F2, F3)   # polynomial in x,y,z and the 22 unknowns (cubic in unknowns)
print("(2) symbolic det: ", len(D_sym), "terms")

# known solution vector
sol = {}
for nm, n in layout:
    src = given[nm] + [{}]*10
    for j in range(n):
        sol[f"{nm}{j}"] = src[j].get((0,0,0), 0) if src[j] else 0
sol_vec = [sol[nm] for nm in names]

def eval_unknowns(p, vals):
    """substitute unknown vars (indices 3..) with numbers; return dict over (x,y,z)."""
    out = {}
    for k, c in p.items():
        m = c
        for j in range(NU):
            e = k[3+j]
            if e:
                m *= vals[j]**e
        if m:
            kk = (k[0], k[1], k[2])
            out[kk] = out.get(kk, 0) + m
    return {k: c for k, c in out.items() if c}

chk = eval_unknowns(D_sym, sol_vec)
print("(3a) det at known solution:", chk == {(0,0,0): -2})

# linearization: d(det)/d(unknown_j) at sol, for each j: a poly in (x,y,z).
# constraint: det must be CONSTANT => all nonconstant (x,y,z)-coefficients of
# the perturbed det vanish to first order.
rows = {}   # (xyz-monomial) -> row over unknowns
for j in range(NU):
    dp = pdiff(D_sym, 3+j)
    dpv = eval_unknowns(dp, sol_vec)
    for k, c in dpv.items():
        if k == (0,0,0):
            continue   # constant part free (det may drift to another constant)
        rows.setdefault(k, [0]*NU)
        rows[k][j] = c
M = list(rows.values())
print(f"(3b) linearized constraint matrix: {len(M)} x {NU}")

def rank_frac(Mat):
    Mat = [ [Fr(c) for c in row] for row in Mat ]
    r = 0; cols = len(Mat[0]) if Mat else 0
    for c in range(cols):
        piv = None
        for i in range(r, len(Mat)):
            if Mat[i][c] != 0: piv = i; break
        if piv is None: continue
        Mat[r], Mat[piv] = Mat[piv], Mat[r]
        pv = Mat[r][c]
        Mat[r] = [x / pv for x in Mat[r]]
        for i in range(len(Mat)):
            if i != r and Mat[i][c] != 0:
                f = Mat[i][c]
                Mat[i] = [a - f*b for a, b in zip(Mat[i], Mat[r])]
        r += 1
        if r == len(Mat): break
    return r

rk = rank_frac(M)
print(f"(3c) rank = {rk}, NULLITY (infinitesimal deformations within ansatz) = {NU - rk}")
print("     unknown order:", names)
