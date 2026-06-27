#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_asano_joint_floor_kpswf15.py   (kind-pasteur 2026-06-27, TOOL 2 Asano/Lee-Yang)

THE JOINT FEW-APEX FLOOR as a Lee-Yang / Asano partition function.

Covering set S = R u 14Q, |S|=13.  R = 14-free big part (13-r speeds), Q = apexes (r=|Q| in 2..6).
Joint loneliness = meas(R-safe cap Q-lonely) = meas( cap_{s in R u 14Q} D_s^c ) = P(M_full=0),
where M_full(t) = #{ s in R u 14Q : ||s t|| < 1/14 } and EVERY danger event has meas(D_s)=1/7.

THE DANGER-COUNT PGF (the partition function).
   G(z) = E[ z^{M_full} ] = sum_j P(M_full=j) z^j,    L = meas(both) = G(0) = P(M_full=0).
   With separate fugacities for the two blocks:
   Xi(lambda, mu) = integral_0^1 prod_{s in R}(1 - lambda 1_{D_s}) * prod_{m in Q}(1 - mu 1_{D_{14m}}) dt.
   meas(R-safe)   = Xi(1, 0)         (Q-fugacity off: Q-block contributes 1)
   meas(Q-lonely) = Xi(0, 1)
   meas(both)     = Xi(1, 1).
   R' = Xi(1,1) / ( Xi(1,0) Xi(0,1) ).

LEE-YANG / ASANO QUESTIONS.
  (1) ZERO-FREE bidisk:  for which (rho_lam, rho_mu) is Xi(lambda,mu) != 0 on |lambda|<=rho_lam,
      |mu|<=rho_mu?  If (1,1) is strictly inside the zero-free region and Xi>=0 on the real box,
      then meas(both)=Xi(1,1)>0.
  (2) ASANO MERGE of the R-block and Q-block: collapse each block's fugacity to one variable, then
      Asano-contract the two blocks.  Track the certified radius.
  (3) The HALF-PLANE (Hurwitz) form: G(z) has roots z=zeta; meas(both)=G(0)!=0 iff 0 is not a root,
      which holds iff P(M=0)>0.  We measure the CLEARANCE of the PGF roots from 0 and from the
      Lee-Yang circle |1-zeta|=1 (the image of the lambda-unit-circle).

We TEST on the documented few-apex covering sets and synthetic r=2..6 variants, report R' and the
zero-free clearances, and check whether a UNIFORM clearance c>0 holds.
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
                                          complement_arcs, danger_arcs, loneliness_exact)

# ----------------------------------------------------------------- joint atom decomposition
def joint_atoms(R, Q, thr=F(1, 14)):
    """
    Atoms of [0,1) by the pattern (a_s for s in R ; b_m for m in Q) of danger membership.
    Returns dict pattern -> weight.  pattern = (rbits..., qbits...).
    """
    speedsR = list(R)
    speedsQ = [14 * m for m in Q]
    speeds = speedsR + speedsQ
    L0 = 1
    for s in speeds:
        L0 = L0 * s // gcd(L0, s)
    D = 14 * L0
    bps = {F(0), F(1)}
    for s in speeds:
        for a, b in danger_arcs(s, thr):
            bps.add(a); bps.add(b)
    bps = sorted(bps)
    acc = {}
    for i in range(len(bps) - 1):
        lo, hi = bps[i], bps[i + 1]
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        pat = tuple(1 if (((s * mid) % 1) < thr or ((s * mid) % 1) > 1 - thr) else 0 for s in speeds)
        acc[pat] = acc.get(pat, F(0)) + (hi - lo)
    return acc, len(speedsR), len(speedsQ)

def Xi_joint(acc, nR, nQ, lam, mu):
    total = 0.0 + 0.0j
    for pat, w in acc.items():
        p = complex(w)
        for i in range(nR):
            if pat[i]:
                p *= (1 - lam)
        for j in range(nQ):
            if pat[nR + j]:
                p *= (1 - mu)
        total += p
    return total

def joint_quantities(R, Q, thr=F(1, 14)):
    acc, nR, nQ = joint_atoms(R, Q, thr)
    measR = Xi_joint(acc, nR, nQ, 1.0, 0.0).real
    measQ = Xi_joint(acc, nR, nQ, 0.0, 1.0).real
    both = Xi_joint(acc, nR, nQ, 1.0, 1.0).real
    Rp = both / (measR * measQ) if measR > 1e-15 and measQ > 1e-15 else None
    return acc, nR, nQ, measR, measQ, both, Rp

def bidisk_minmod(acc, nR, nQ, rho_l, rho_m, n=60):
    """min |Xi(lambda,mu)| over the bidisk torus |lambda|=rho_l, |mu|=rho_m."""
    best = float('inf')
    for ia in range(n):
        al = 2 * pi * ia / n
        lam = rho_l * complex(math.cos(al), math.sin(al))
        for ib in range(n):
            am = 2 * pi * ib / n
            mu = rho_m * complex(math.cos(am), math.sin(am))
            v = abs(Xi_joint(acc, nR, nQ, lam, mu))
            if v < best:
                best = v
    return best

def pgf_roots(R, Q, thr=F(1, 14)):
    """Danger-count PGF G(z)=E[z^M]; return roots and clearances."""
    acc, nR, nQ = joint_atoms(R, Q, thr)
    dist = {}
    for pat, w in acc.items():
        j = sum(pat)
        dist[j] = dist.get(j, F(0)) + w
    k = nR + nQ
    coeffs = [float(dist.get(j, F(0))) for j in range(k + 1)]
    roots = Poly.polyroots(coeffs)
    return roots, dist

def banner(s):
    print("=" * 78); print(s); print("=" * 78)

if __name__ == "__main__":
    banner("Joint few-apex floor R' and danger-count PGF clearance (documented + synthetic)")
    cases = [
        ("r=1 doc", tuple(range(1, 12)) + (13,), (6,)),
        ("r=2 doc", tuple(range(1, 11)) + (13,), (6, 11)),
        ("r=2 alt", tuple(range(1, 11)) + (12,), (6, 13)),
        ("r=3", tuple(range(1, 10)) + (11,), (6, 12, 13)),
        ("r=4", tuple(range(1, 9)) + (10,), (6, 11, 12, 13)),
        ("r=5", tuple(range(1, 8)) + (9,), (6, 10, 11, 12, 13)),
        ("r=6", tuple(range(1, 7)) + (8,), (6, 9, 10, 11, 12, 13)),
    ]
    print(f"{'tag':<9}{'|R|':>4}{'r':>3}{'measR':>10}{'measQ':>10}{'both':>10}{'R prime':>10}"
          f"{'min|root|':>11}{'min|1-z|':>10}")
    rows = []
    for tag, R, Q in cases:
        acc, nR, nQ, mR, mQ, both, Rp = joint_quantities(R, Q)
        roots, dist = pgf_roots(R, Q)
        minabs = min((abs(r) for r in roots), default=float('inf'))
        min1mz = min((abs(1 - r) for r in roots), default=float('inf'))
        Rps = f"{Rp:.5f}" if Rp is not None else "n/a"
        print(f"{tag:<9}{nR:>4}{nQ:>3}{mR:>10.5f}{mQ:>10.5f}{both:>10.5f}{Rps:>10}"
              f"{minabs:>11.4f}{min1mz:>10.4f}")
        rows.append((tag, Rp, min1mz))

    banner("Bidisk zero-free test: is (lambda,mu)=(1,1) strictly inside the zero-free region?")
    for tag, R, Q in cases[:5]:
        acc, nR, nQ = joint_atoms(R, Q)
        print(f"  {tag}: min|Xi| on bidisk torus")
        for (rl, rm) in [(1.0, 1.0), (1.1, 1.1), (1.0, 1.5), (1.5, 1.0)]:
            mm = bidisk_minmod(acc, nR, nQ, rl, rm, n=48)
            print(f"     |lam|={rl},|mu|={rm}: min|Xi|={mm:.6f}  [{'zero-free' if mm>1e-7 else 'ZERO'}]")
