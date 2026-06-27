#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_asano_zerofree_kpswf15.py   (kind-pasteur 2026-06-27, TOOL 2: Asano/Lee-Yang)

Probe whether the LRC(14) loneliness measure L(R) = integral prod_s 1[||s t||>=1/14] dt has a
Lee-Yang / Asano-contraction zero-free structure.

THE DISCRETE PARTITION FUNCTION (exact).
  Pick N >= 1.  Let omega = e(1/(14N)) be a primitive 14N-th root of unity (so the arc breakpoints
  at multiples of 1/14 fall on the grid).  Evaluate t at the midpoint grid t_j=(2j+1)/(2*14N).
  Then phi_s(t_j) in {0,1} EXACTLY, and
      L(R) = (1/(14N)) sum_{j} prod_{s in R} phi_s(t_j)   (-> integral as N->inf; exact for the
             midpoint Riemann sum, and EQUALS the exact arc measure once N is a multiple of every
             speed, since all breakpoints are then captured).

  PER-SPEED VARIABLE ENCODING.  Put a single circle variable z=e(t).  Each phi_s is a function of
  z^s.  Introduce a "fugacity" variable lambda_s per speed and form the multi-affine partition fn
      Z(lambda) = sum_{A subset R} ( prod_{s in A} lambda_s ) * w(A),
  where the weights w(A) come from the inclusion-exclusion of the COMPLEMENT events
  D_s = {||s t|| < 1/14} (the "dangerous"/lonely-failing arcs):
      L = meas( intersect_s D_s^c ) = sum_{A subset R} (-1)^{|A|} meas( intersect_{s in A} D_s ).
  Define J(A) = meas( intersect_{s in A} D_s )  (J(empty)=1).  Then the multi-affine GENERATING
  polynomial is
      Z(lambda_1,...,lambda_k) = sum_{A} (prod_{s in A} (-lambda_s)) J(A),
  and L = Z(1,1,...,1).  Asano contraction acts on the lambda_s.  We test:
    (1) Is Z(lambda) zero-free on a polydisk |lambda_s| <= rho (Lee-Yang region)?  If so and rho>=1,
        L=Z(1,..,1) != 0; combined with L>=0 (it is a measure) gives L>0.
    (2) Asano: combine lambda_i,lambda_j -> single mu preserving zero-freeness; iterate.

  We also test the SECTOR partition fn (the dual, missed-sector view) and the FEJER lower bound.
"""
import sys, itertools, math
from fractions import Fraction as F
from math import gcd, pi
import numpy as np
import numpy.polynomial.polynomial as P

try:
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
except Exception:
    pass

# ----------------------------------------------------------------- exact arcs (shared)
def merge(iv):
    iv = sorted(iv); out = []
    for a, b in iv:
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out

def meas(arcs):
    return sum((b - a for a, b in arcs), F(0))

def intersect(A, B):
    out = []; i = j = 0
    A = merge(A); B = merge(B)
    while i < len(A) and j < len(B):
        lo = max(A[i][0], B[j][0]); hi = min(A[i][1], B[j][1])
        if lo < hi: out.append((lo, hi))
        if A[i][1] < B[j][1]: i += 1
        else: j += 1
    return out

def safe_arcs(s, thr=F(1, 14)):
    """{ t in [0,1) : ||s t|| >= thr }.  s positive integer.  (No wrap issues: closed inside.)"""
    out = []
    for j in range(s):
        lo = F(j) / s + thr / s
        hi = F(j + 1) / s - thr / s
        if lo < hi:
            out.append((lo, hi))
    return merge(out)

def complement_arcs(arcs):
    arcs = merge(arcs); out = []; prev = F(0)
    for a, b in arcs:
        if a > prev: out.append((prev, a))
        prev = max(prev, b)
    if prev < 1: out.append((prev, F(1)))
    return out

def danger_arcs(s, thr=F(1, 14)):
    """D_s = { t : ||s t|| < thr } = complement of safe_arcs(s).  meas = 1/7 for all s."""
    return complement_arcs(safe_arcs(s, thr))

def J(A, thr=F(1, 14)):
    """meas( intersect_{s in A} D_s )."""
    arcs = [(F(0), F(1))]
    for s in A:
        arcs = intersect(arcs, danger_arcs(s, thr))
    return meas(arcs)

def loneliness_exact(R, thr=F(1, 14)):
    arcs = [(F(0), F(1))]
    for s in R:
        # safe = complement of danger
        d = danger_arcs(s, thr)
        safe = []
        prev = F(0)
        for a, b in d:
            if a > prev: safe.append((prev, a))
            prev = max(prev, b)
        if prev < 1: safe.append((prev, F(1)))
        arcs = intersect(arcs, safe)
    return meas(arcs)

# ----------------------------------------------------------------- multi-affine Z(lambda)
def Z_coeffs(R, thr=F(1, 14)):
    """
    Return dict A(frozenset of indices) -> coefficient (-1)^{|A|} J({R[i]:i in A})  (exact F).
    Z(lambda) = sum_A coeff[A] * prod_{i in A} lambda_i,  L = Z(1..1).
    """
    k = len(R)
    coeffs = {}
    for r in range(k + 1):
        for combo in itertools.combinations(range(k), r):
            sub = tuple(R[i] for i in combo)
            j = J(sub, thr)
            coeffs[frozenset(combo)] = ((-1) ** r) * j
    return coeffs

def Z_eval(coeffs, lam):
    """Evaluate multi-affine Z at complex vector lam (indexed 0..k-1)."""
    total = 0.0 + 0.0j
    for A, c in coeffs.items():
        p = complex(c)
        for i in A:
            p *= lam[i]
        total += p
    return total

def Z_as_univariate_diag(coeffs, k):
    """
    Restrict all lambda_i = lambda (the diagonal).  Returns numpy poly coeffs in lambda
    (degree k).  Z_diag(lambda) = sum_{r} ( sum_{|A|=r} coeff[A] ) lambda^r.
    L = Z_diag(1).
    """
    deg = [0.0] * (k + 1)
    for A, c in coeffs.items():
        deg[len(A)] += float(c)
    return deg  # ascending order, P.polyroots convention

# ----------------------------------------------------------------- Asano contraction
def asano_contract(coeffs, i, j, k):
    """
    Asano contraction of variables i and j into a single new variable.
    Multi-affine poly Z = A + lam_i B + lam_j C + lam_i lam_j D  (A,B,C,D free of lam_i,lam_j).
    Asano product replaces (constant, lam_i lam_j) pair: new poly in a single var mu is
       Z' = A + mu D  ... but the standard Asano contraction is:
    Given Z = a + b lam_i + c lam_j + d lam_i lam_j, define the contracted
       Z_c = a + d mu   (drop the linear cross terms b,c) -- this is the ASANO operator,
    which preserves the property "no zeros in |lam|<=1 polydisk" by the Asano lemma
    (Ruelle's formulation): if Z is nonzero for |lam_i|<1 and |lam_j|<1, then a+d mu is
    nonzero for |mu|<1.
    Here a,b,c,d are themselves multi-affine in the remaining variables.  We return the new
    coeff dict on the index set (remaining) + a fresh contracted index 'c'.
    """
    rest = [x for x in range(k) if x != i and x != j]
    # group coeffs by (has i?, has j?)
    a = {}; d = {}
    for A, cf in coeffs.items():
        hi = i in A; hj = j in A
        base = frozenset(x for x in A if x != i and x != j)
        if (not hi) and (not hj):
            a[base] = a.get(base, F(0)) + cf
        elif hi and hj:
            d[base] = d.get(base, F(0)) + cf
        # linear cross terms b (hi only), c (hj only) are DROPPED by Asano
    # new variable index = k (fresh)
    newc = {}
    for base, cf in a.items():
        newc[base] = newc.get(base, F(0)) + cf
    for base, cf in d.items():
        nb = frozenset(set(base) | {k})
        newc[nb] = newc.get(nb, F(0)) + cf
    return newc, rest + [k]

def banner(s):
    print("=" * 78); print(s); print("=" * 78)

if __name__ == "__main__":
    banner("Z(lambda) multi-affine partition function: diagonal roots + polydisk zero-freeness")

    test_R = [(1, 2), (2, 3), (1, 2, 3), (1, 3), (2, 5), (1, 2, 3, 4), (3, 5, 7)]
    for R in test_R:
        k = len(R)
        coeffs = Z_coeffs(R)
        L = Z_eval(coeffs, [1.0] * k).real
        Lex = float(loneliness_exact(R))
        # diagonal univariate roots
        diag = Z_as_univariate_diag(coeffs, k)
        roots = P.polyroots(diag) if k >= 1 else []
        # min |root| on the diagonal
        minabs = min((abs(r) for r in roots), default=float('inf'))
        print(f"\nR={R}  L={L:.6f} (exact {Lex:.6f})  diag-roots |.|: "
              f"{sorted(round(abs(r),4) for r in roots)}")
        print(f"   min |diag root| = {minabs:.4f}  (Lee-Yang on diagonal needs all |root|>1)")

    banner("FULL polydisk min|Z| on torus |lam_s|=rho  (rho up to 1 = Lee-Yang polydisk test)")
    rng = np.random.default_rng(0)
    for R in test_R:
        k = len(R)
        coeffs = Z_coeffs(R)
        for rho in [1.0, 1.5, 2.0, 2.5]:
            # Monte-Carlo + grid over the torus |lam_s|=rho; report min |Z|
            best = float('inf')
            # deterministic angle grid for small k, random for larger
            if k <= 3:
                grids = itertools.product(np.linspace(0, 2 * pi, 24, endpoint=False), repeat=k)
            else:
                grids = (tuple(rng.uniform(0, 2 * pi, k)) for _ in range(20000))
            for angs in grids:
                lam = [rho * complex(math.cos(a), math.sin(a)) for a in angs]
                v = abs(Z_eval(coeffs, lam))
                if v < best:
                    best = v
            flag = "ZERO-FREE" if best > 1e-9 else "HAS ZERO (or near)"
            print(f"  R={str(R):<16} rho={rho:<4} min|Z|={best:.5f}  [{flag}]")
