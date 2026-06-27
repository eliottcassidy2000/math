#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_asano_contraction_proof_kpswf15.py   (kind-pasteur 2026-06-27, TOOL 2 Asano/Lee-Yang)

THE ASANO-CONTRACTION POSITIVITY ROUTE for the LRC(14) loneliness floor.

PARTITION FUNCTION (independent-set / polymer form).
  Danger events  D_s = { t : ||s t|| < 1/14 },   meas(D_s) = 1/7  for EVERY speed s.
  Safe loneliness  L(R) = meas( intersect_s D_s^c ) = integral prod_s (1 - 1_{D_s}(t)) dt.
  Introduce fugacities lambda_s.  Define the GRAND PARTITION FUNCTION
      Xi(lambda) = integral_0^1 prod_s ( 1 - lambda_s 1_{D_s}(t) ) dt
                 = sum_{A subset R} prod_{s in A}(-lambda_s) * J(A),   J(A)=meas(cap_{s in A} D_s).
  This is MULTI-AFFINE in lambda, and  L(R) = Xi(1,...,1).

ASANO LEMMA (Ruelle 1971; Asano 1970).
  If a multi-affine polynomial P = a + b x + c y + d x y has NO zero in {|x|<r}x{|y|<s}, then the
  ASANO CONTRACTION  P# = a + d (xy->u) := a + d u  has no zero in {|u| < r s}.
  (Proof: a+d u = 0 => u = -a/d; the polydisk-zero-freeness of P forces |a/d| >= r s.)

SINGLE-SPEED FACTOR.  Xi_s(lambda_s) = 1 - lambda_s/7   (since J({s})=meas(D_s)=1/7).
  Zero only at lambda_s = 7.  So EACH factor (1 - lambda_s 1_{D_s}) -- as a polynomial in lambda_s
  at a GENERIC FIXED t with 1_{D_s}(t) in {0,1} -- is zero-free for |lambda_s| < 1 (if 1_{D_s}=1,
  root at 1) UNION the whole plane (if 1_{D_s}=0).  The AVERAGE Xi is what we contract.

THE PROOF STRATEGY (Asano-merge the speeds one at a time):
  We DO NOT have a product of single-variable polynomials (the integral couples them).  Instead
  we use the ASANO MERGE on the integrand viewed as a multilinear form, exactly as in the
  Lee-Yang / Ising literature: write Xi as the diagonal (all lambda_s = lambda) and bound the
  root via successive Asano contractions starting from the pair structure.  We TEST:

  (Test 1) Single-speed root = 7 (|root|>1).  [TRIVIAL zero-free, the seed.]
  (Test 2) For each pair (s,s'), is the 2-var restriction of Xi (other lambdas = 0)
           P_2 = 1 - lambda_s/7 - lambda_s'/7 + lambda_s lambda_s' J(s,s')
           zero-free on |lambda|<r x |lambda|<r for r up to ?  -> gives the Asano radius.
  (Test 3) The Asano-contracted univariate polynomial and its root radius, iterated, vs the TRUE L.
  (Test 4) The KEY RIGOROUS BOUND: Xi(lambda,...,lambda) has all roots OUTSIDE |lambda|<=R0 where
           R0 is certified by the GROENEVELD / penrose-tree / Dobrushin condition
           sum_{s} meas(D_s) <= ... .  We compute the certified R0 and compare to the empirical
           min|root|.
"""
import sys, itertools, math
from fractions import Fraction as F
from math import gcd, pi
import numpy as np
import numpy.polynomial.polynomial as Poly

try:
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
except Exception:
    pass

# reuse exact arc machinery
sys.path.insert(0, "04-computation")
from lrc14_asano_zerofree_kpswf15 import (merge, meas, intersect, safe_arcs,
                                          complement_arcs, danger_arcs, J, loneliness_exact)

# ----------------------------------------------------------------- 2-var pair restriction
def pair_poly(s1, s2, thr=F(1, 14)):
    """P_2(x,y) = 1 - x*meas(D_s1) - y*meas(D_s2) + x*y*J(s1,s2).  Coeffs (a,b,c,d) exact."""
    a = F(1)
    b = -J((s1,), thr)
    c = -J((s2,), thr)
    d = J((s1, s2), thr)
    return a, b, c, d

def pair_zerofree_radius(a, b, c, d, samples=400):
    """
    Empirically find the largest r such that a + b x + c y + d x y != 0 on |x|<=r,|y|<=r.
    Multi-affine: min over torus |x|=|y|=r of |P|.  Binary search r.
    """
    def minmod(r):
        best = float('inf')
        for ia in range(samples):
            ang1 = 2 * pi * ia / samples
            x = r * complex(math.cos(ang1), math.sin(ang1))
            # for fixed x, P is affine in y: P = (a+b x) + (c + d x) y; min over |y|=r is
            # | |a+b x| - r|c+d x| |  (affine-in-y minimum on circle)
            A = complex(a) + complex(b) * x
            B = complex(c) + complex(d) * x
            m = abs(abs(A) - r * abs(B))
            if m < best:
                best = m
        return best
    lo, hi = 0.0, 50.0
    for _ in range(60):
        mid = (lo + hi) / 2
        if minmod(mid) > 1e-9:
            lo = mid
        else:
            hi = mid
    return lo

# ----------------------------------------------------------------- diagonal partition fn roots
def Xi_diag_coeffs(R, thr=F(1, 14)):
    """ascending coeffs of Xi(lambda,...,lambda) = sum_r (sum_{|A|=r} (-1)^r J(A)) lambda^r."""
    k = len(R)
    deg = [F(0)] * (k + 1)
    for r in range(k + 1):
        for combo in itertools.combinations(range(k), r):
            sub = tuple(R[i] for i in combo)
            deg[r] += ((-1) ** r) * J(sub, thr)
    return deg

def banner(s):
    print("=" * 78); print(s); print("=" * 78)

if __name__ == "__main__":
    banner("Test 1+2: single-speed root=7; pairwise Asano radius r (P_2 zero-free on |.|<r)")
    pairs = [(1, 2), (2, 3), (1, 3), (2, 4), (1, 6), (2, 5), (3, 5), (5, 7), (11, 13)]
    print(f"{'pair':<10}{'meas D':>9}{'J(pair)':>12}{'pair Asano r':>16}{'product 7*7=49?':>18}")
    for s1, s2 in pairs:
        a, b, c, d = pair_poly(s1, s2)
        r = pair_zerofree_radius(a, b, c, d)
        # Asano contraction of P_2 -> P# = a + d u, root u = -a/d = -1/J; |root|=1/J
        contracted_root = float(F(1) / d) if d != 0 else float('inf')
        print(f"({s1},{s2})".ljust(10) + f"{float(-b):>9.5f}{float(d):>12.5f}{r:>16.4f}"
              f"{contracted_root:>18.4f}")

    banner("Test 3+4: diagonal Xi roots vs certified radius; L=Xi(1)")
    test_R = [(1, 2), (1, 2, 3), (1, 2, 3, 4), (1, 2, 3, 4, 5), (1, 2, 3, 4, 5, 6),
              (2, 3, 5, 7, 11, 13), (3, 5, 7), (84, 154)]
    print(f"{'R':<26}{'k':>3}{'sum meas D':>12}{'L=Xi(1)':>12}{'min|root|':>12}{'1/(k/7)':>10}")
    for R in test_R:
        k = len(R)
        coeffs = [float(c) for c in Xi_diag_coeffs(R)]
        roots = Poly.polyroots(coeffs)
        minabs = min((abs(r) for r in roots), default=float('inf'))
        L = float(loneliness_exact(R))
        summeas = k * (1.0 / 7.0)
        # naive Dobrushin-style radius guess 1/(sum meas) -- compare
        print(f"{str(R):<26}{k:>3}{summeas:>12.5f}{L:>12.5f}{minabs:>12.5f}{(7.0/k):>10.4f}")
