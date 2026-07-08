#!/usr/bin/env python3
r"""
lrc_k10_union_wide_kps_S75.py  (kind-pasteur-2026-07-07-S75, HYP-5247)

THE CLEAN k=10 SPLIT (sidesteps the teeth AND the average form):
  - diam <= D0:  intersection ledger (kps-S60), conditional, covers compact families.
  - diam >  D0:  UNION BOUND rho*(P,E) >= meas(G_P) + mu(E) - 1, using the UNCONDITIONAL
                 mu(E) = meas{maxgap(E) > 1/7}.  Closes iff mu(E) >= 1 - meas(G_P) + m_P.

The union bound uses UNCONDITIONAL mu -- NO G_P teeth, NO resonance-window bookkeeping.  The
only question: is inf over WIDE primitive 10-clusters of mu(E) >= 1 - min_{|P|=3} meas(G_P) +
m_P?  Wide families have HIGH mu (spread), so this should hold with room; the compact families
(low mu, e.g. AP_10 mu=0.388) are exactly the ledger's job.  This file finds:
  (1) min_{|P|=3} meas(G_P) exact -> the union requirement mu >= mu_req;
  (2) min mu over primitive 10-clusters as a function of the diameter lower bound D0
      (structured/adversarial search) -> the crossover D0* where min mu(diam>D0) >= mu_req;
  (3) compare D0* to the ledger reach (kps-S60: k=10 diam <= 17).  If D0* <= 17 (ledger reach),
      k=10 is CLOSED: ledger below, union bound above, contiguous.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations
import numpy as np
import random

M_P = F(14249, 252252)

def gp_measure_exact(P):
    bad = []
    for p in P:
        w = F(1, 14 * p)
        for j in range(p + 1):
            lo, hi = F(j, p) - w, F(j, p) + w
            bad.append((max(lo, F(0)), min(hi, F(1))))
    bad.sort()
    tot = F(0); cl = ch = None
    for lo, hi in bad:
        if cl is None: cl, ch = lo, hi
        elif lo <= ch: ch = max(ch, hi)
        else: tot += ch - cl; cl, ch = lo, hi
    if cl is not None: tot += ch - cl
    return 1 - tot

# (1) min meas(G_P) over |P|=3
best = F(2); bestP = None
for P in combinations(range(1, 14), 3):
    m = gp_measure_exact(P)
    if m < best: best, bestP = m, P
mu_req = 1 - best + M_P
print("=" * 88)
print("k=10 CLEAN SPLIT: intersection ledger (diam<=17) + union bound (diam>17)")
print("=" * 88)
print(f"  min_|P|=3 meas(G_P) = {best} = {float(best):.4f} at P={bestP}")
print(f"  => union bound needs UNCONDITIONAL mu(E) >= 1 - meas(G_P) + m_P = {float(mu_req):.4f}")
print(f"     (worst over P; for larger meas(G_P) the requirement is easier)")

# (2) min mu over primitive 10-clusters vs diameter lower bound
def mu17(E, n=120000):
    Ev = np.array(E, dtype=np.float64)
    xs = (np.arange(n) + 0.5) / n
    ph = np.mod(np.outer(xs, Ev), 1.0); ph.sort(axis=1)
    g = np.diff(ph, axis=1); wrap = ph[:, 0] + 1.0 - ph[:, -1]
    mg = np.maximum(g.max(axis=1), wrap)
    return float(np.mean(mg > 1.0 / 7.0))

def primitive(E):
    g = 0
    for a in E[1:]: g = gcd(g, a - E[0])
    return g == 1

def min_mu_wide(D0, rng, tries=2500):
    """minimize mu over primitive 10-clusters with diam >= D0 (structured + random + descent)."""
    best = (1.0, None)
    # structured seeds: dilated-AP-like (lowest mu among wide), 2-AP+bump, AP+far
    seeds = []
    for d in range(2, D0 // 2 + 3):
        base = list(range(0, 10 * d, d))[:9]
        for bump in (D0, D0 + 3, 10 * d - d + 1):
            E = sorted(set(base + [bump]))[:10]
            if len(E) == 10 and E[-1] - E[0] >= D0 and primitive([e - E[0] for e in E]):
                seeds.append([e - E[0] for e in E])
    for _ in range(400):
        k = rng.randint(8, 9)
        core = sorted(random.Random(rng.random()).sample(range(0, 14), k))
        bump = rng.randint(D0, D0 + 25)
        E = sorted(set(core + [bump]))
        if len(E) == 10 and E[-1] - E[0] >= D0:
            E = [e - E[0] for e in E]
            if primitive(E): seeds.append(E)
    for s in seeds:
        v = mu17(s)
        if v < best[0]: best = (v, s)
    # descent
    cur = best
    for _ in range(tries):
        if cur[1] is None: break
        E = list(cur[1]); idx = rng.randrange(1, 10)
        E[idx] += rng.choice([-2, -1, 1, 2])
        E = sorted(set(E))
        if len(E) != 10 or E[0] < 0 or E[-1] - E[0] < D0: continue
        E = [e - E[0] for e in E]
        if not primitive(E): continue
        v = mu17(E)
        if v < cur[0]: cur = (v, E)
    return cur

rng = random.Random(75)
print(f"\n  min mu over primitive 10-clusters, by diameter lower bound D0:")
print(f"  {'D0':>4} {'min mu (wide)':>14} {'>= mu_req?':>10}  {'argmin E'}")
crossover = None
for D0 in (11, 13, 15, 17, 18, 20, 24, 30, 40):
    best = (1.0, None)
    for s in range(5):
        v = min_mu_wide(D0, random.Random(75 + s))
        if v[0] < best[0]: best = v
    ok = best[0] >= float(mu_req)
    if ok and crossover is None: crossover = D0
    print(f"  {D0:>4} {best[0]:>14.4f} {str(ok):>10}  {best[1]}")
print()
if crossover is not None:
    print(f"  => union bound (unconditional mu) closes wide families for diam >= {crossover}.")
    if crossover <= 18:
        print(f"     Ledger covers diam <= 17 (kps-S60); {crossover} <= 18 => CONTIGUOUS => k=10 CLOSED")
        print(f"     via [ledger diam<=17] + [union bound diam>=18], NO average form, NO teeth.")
    else:
        print(f"     GAP: ledger reaches 17, union needs diam >= {crossover}; band [18,{crossover-1}]")
        print(f"     needs the average form / D_q-window / exhaustion.")
else:
    print("  => union bound does not reach mu_req at tested D0; wide families still dip below.")
