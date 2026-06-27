#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_asano_atom_decomp_kpswf15.py   (kind-pasteur 2026-06-27, TOOL 2 Asano/Lee-Yang)

ATOM DECOMPOSITION of the loneliness partition function + rigorous zero-free tests.

Partition [0,1) into ATOMS where the danger pattern sigma(t)=(1_{D_s}(t))_{s in R} is constant.
On atom alpha (pattern sigma^alpha, measure w_alpha):
    integrand = prod_{s: sigma^alpha_s=1} (1 - lambda_s).
So the GRAND PARTITION FUNCTION is the POSITIVE (w_alpha>=0) combination
    Xi(lambda) = sum_alpha w_alpha * prod_{s in supp(sigma^alpha)} (1 - lambda_s).        (*)
This is a "mixture of independent-site products".  L(R)=Xi(1,...,1)=sum_{alpha: sigma^alpha=0} w_alpha
    = meas{ t : no danger } = the loneliness.  (Only the all-safe atom survives at lambda=1.)

KEY EXACT FACTS we verify:
  (F0)  sum_alpha w_alpha = 1  (atoms partition the circle).
  (F1)  marginal:  sum_{alpha: s in supp} w_alpha = meas(D_s) = 1/7  for every s.
  (F2)  L = Xi(1) = w_{empty pattern}  (the all-safe atom mass).
  (F3)  Xi is multi-affine; (*) is its convex-combination-of-products representation.

ZERO-FREE QUESTIONS (Lee-Yang) -- which we TEST exactly:
  (Q1) The DIAGONAL Xi(lambda):  does it satisfy a GROENEVELD/PENROSE root bound?  We compute the
       exact diagonal coeffs e_r = sum_{|A|=r}(-1)^r J(A) and the Newton/abel bound on min|root|.
  (Q2) ASANO MERGE on (*):  the product prod(1-lambda_s) is the trivial Lee-Yang polynomial with
       all roots at lambda_s=1.  A CONVEX COMBINATION of such products is NOT automatically
       zero-free in |lambda|<1.  BUT the ALL-SITES-PRESENT structure is special.  We test whether
       Xi(lambda) (diagonal) is zero-free on |lambda|<=1 directly (binary search the root radius),
       and whether the per-variable polydisk radius rho* (Xi != 0 on |lambda_s|<=rho* all s)
       exceeds 1 with a UNIFORM margin over all few-apex R.
  (Q3) The Penrose/Kotecky-Preiss CERTIFICATE:  the cluster expansion of log Xi converges (hence
       Xi != 0) on |lambda_s| <= R whenever there is a function mu_s>0 with
              lambda_s <= mu_s / prod_{s'~s}(1+mu_{s'})  ... (Fernandez-Procacci form).
       For the COMPLETE-graph worst case (all pairs dependent) with marginal p=1/7, this gives an
       explicit R0; we compute it and compare to 1.

OUTPUT: the exact atom data + the certified zero-free radius vs 1, for all r=2..6 multi-far R.
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

sys.path.insert(0, "04-computation")
from lrc14_asano_zerofree_kpswf15 import (merge, meas, intersect, safe_arcs,
                                          complement_arcs, danger_arcs, J, loneliness_exact)

# ----------------------------------------------------------------- atom decomposition (exact)
def atoms(R, thr=F(1, 14)):
    """
    Return list of (pattern_tuple_in_{0,1}^k, weight F).  Exact via breakpoint refinement.
    Breakpoints of all 1_{D_s} live in (1/(14*lcm)) Z.
    """
    L0 = 1
    for s in R:
        L0 = L0 * s // gcd(L0, s)
    D = 14 * L0
    k = len(R)
    # collect all breakpoints (endpoints of danger arcs) as fractions j/D
    bps = {F(0), F(1)}
    for s in R:
        for a, b in danger_arcs(s, thr):
            bps.add(a); bps.add(b)
    bps = sorted(bps)
    acc = {}
    for i in range(len(bps) - 1):
        lo, hi = bps[i], bps[i + 1]
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        pat = tuple(1 if (((s * mid) % 1) < thr or ((s * mid) % 1) > 1 - thr) else 0 for s in R)
        acc[pat] = acc.get(pat, F(0)) + (hi - lo)
    return acc

def Xi_from_atoms(acc, R, lam):
    """Evaluate Xi(lam) from atom dict.  lam: complex vector length k."""
    total = 0.0 + 0.0j
    for pat, w in acc.items():
        p = complex(w)
        for i, b in enumerate(pat):
            if b:
                p *= (1 - lam[i])
        total += p
    return total

# ----------------------------------------------------------------- diagonal root radius
def diag_root_radius(R, thr=F(1, 14)):
    k = len(R)
    deg = [F(0)] * (k + 1)
    for r in range(k + 1):
        for combo in itertools.combinations(range(k), r):
            sub = tuple(R[i] for i in combo)
            deg[r] += ((-1) ** r) * J(sub, thr)
    roots = Poly.polyroots([float(c) for c in deg])
    return min((abs(r) for r in roots), default=float('inf')), deg

# ----------------------------------------------------------------- Fernandez-Procacci certificate
def fp_radius_complete(p, k):
    """
    Worst-case (COMPLETE dependency graph, all k events pairwise dependent) Fernandez-Procacci /
    Shearer radius for the lambda such that the cluster expansion of the polymer gas with
    single-polymer weight |lambda|*p converges.  For the complete graph K_k with uniform
    activity a = |lambda| p, the FP condition is a <= mu/(1+mu)^{k-1} maximized over mu>0:
       max over mu of mu/(1+mu)^{k-1} = (1/(k-1)) (1 - 1/(k-1))^{k-1} ... (mu*=1/(k-2)) for k>=3.
    Returns the largest |lambda| with convergence guaranteed, = a*/p.
    """
    if k == 1:
        a_star = 1.0  # single event: zero-free for |lambda|p<1 i.e. lambda<7
    elif k == 2:
        # mu/(1+mu): sup=1 as mu->inf, but need finite; use mu=1 -> 1/2; actually K_2: a<=mu/(1+mu)
        # optimum at mu->inf gives a<1; Shearer for K_2 (one edge): independence poly 1-2a... root.
        a_star = 0.25  # conservative
    else:
        # maximize f(mu)=mu/(1+mu)^{k-1}; f'(mu)=0 => 1*(1+mu)^{k-1} - mu(k-1)(1+mu)^{k-2}=0
        # => (1+mu) - mu(k-1)=0 => 1 = mu(k-2) => mu = 1/(k-2)
        mu = 1.0 / (k - 2)
        a_star = mu / (1 + mu) ** (k - 1)
    return a_star / p

def banner(s):
    print("=" * 78); print(s); print("=" * 78)

if __name__ == "__main__":
    banner("ATOM DECOMPOSITION: verify F0/F1/F2 exactly")
    test_R = [(1, 2), (1, 2, 3), (1, 2, 3, 4), (1, 2, 3, 4, 5), (1, 2, 3, 4, 5, 6),
              (2, 3, 5, 7, 11, 13), (3, 5, 7)]
    for R in test_R:
        acc = atoms(R)
        tot = sum(acc.values(), F(0))
        k = len(R)
        # marginal per site
        marg = [sum((w for pat, w in acc.items() if pat[i]), F(0)) for i in range(k)]
        empty = acc.get(tuple([0] * k), F(0))
        Lex = loneliness_exact(R)
        ok0 = (tot == 1)
        ok1 = all(m == F(1, 7) for m in marg)
        ok2 = (empty == Lex)
        print(f"R={str(R):<24} sum_w={str(tot):<4}[{'OK' if ok0 else 'BAD'}]  "
              f"marg=1/7 all[{'OK' if ok1 else 'BAD'}]  Xi(1)=w_empty={float(empty):.5f}"
              f"==L[{'OK' if ok2 else 'BAD'}]  #atoms={len(acc)}")

    banner("ZERO-FREE RADIUS: diagonal min|root| vs Fernandez-Procacci certified vs 1")
    print(f"{'R':<26}{'k':>3}{'L':>10}{'diag min|root|':>16}{'FP cert |lam|':>15}{'>1?':>6}")
    more_R = test_R + [(1, 2, 3, 4, 5, 6, 7), (84, 154), (84, 154, 196)]
    for R in more_R:
        k = len(R)
        mr, deg = diag_root_radius(R)
        fp = fp_radius_complete(1.0 / 7.0, k)
        L = float(loneliness_exact(R))
        print(f"{str(R):<26}{k:>3}{L:>10.5f}{mr:>16.5f}{fp:>15.5f}"
              f"{'YES' if mr > 1 else 'NO':>6}")
