#!/usr/bin/env python3
"""
klein-2026-07-07-S157 -- THE VOLTAGE-LIFT READING: lift-with-defect classification of the
mean sidecar (does inf E[maxgap] over c-lift families dip below T* = 1/7 + (6/7)m_P?).

FRAME (the tiling analogy, voltage-lift face): the tournament SC-blowup / double round-robin
lifts a base tournament to 2n vertices by sheet doubling (CLAUDE.md, THM-378). The LRC analog:
a c-LIFT of a base set A is c*A (orbit of c*A at x = orbit of A at cx: a c-sheeted cover in
time), and a LIFT-WITH-DEFECT is F = c*A u B (small defect set B breaking sheets). monad's
E[maxgap] record 2*{1..11} u {11,13} IS a defective 2-lift of AP_11; death-star's prim-sat
2*{1..12} u {13} likewise. QUESTION: over the systematic lift-defect class, what is the exact
min of E[maxgap], and does it stay above T* = 56291/294294 ~ 0.191275 (the honest mean bar,
monad HYP-4787)? If some family dips below T*, the mean sidecar at the m_P bar DIES; if the
class bottoms out above, the record mechanism is bounded and the sidecar stands.

METHOD: numeric screen (numpy grid) over F = c*{1..a} u B, c in {2,3}, |F| = 13,
B subset of [1..Bmax] \ c*{1..a} with |B| = 13 - a; then EXACT rational E[maxgap]
(death-star-style corrected integrator: breakpoints at all m/d, d <= kdenom = max pairwise
difference, with per-cell linear-order + sub-breakpoints for gap-max changes) for the top
candidates. Also mu_{1/7} for winners (tail-vs-mean divergence check).
"""
import numpy as np
from fractions import Fraction as F
from itertools import combinations
from math import gcd
from functools import reduce

TSTAR = F(1,7) + F(6,7)*F(14249,252252)

def Emaxgap_exact(E):
    E = list(E)
    kden = max(max(abs(j) for j in E), max(abs(a-b) for a in E for b in E))
    bps = {F(0), F(1)}
    for d in range(1, kden+1):
        for m in range(1, d):
            bps.add(F(m, d))
    bps = sorted(bps)
    total = F(0)
    for a, b in zip(bps, bps[1:]):
        mid = (a+b)/2
        fl = {j: (j*mid).__floor__() for j in E}
        order = sorted(E, key=lambda j: (j*mid - fl[j]))
        n = len(order); gaps = []
        for s_ in range(n):
            j1 = order[s_]; j2 = order[(s_+1) % n]
            if s_ < n-1:
                c_ = F(j2-j1); b0 = F(-(fl[j2]-fl[j1]))
            else:
                c_ = F(order[0]-order[-1]); b0 = F(-(fl[order[0]]-fl[order[-1]])+1)
            gaps.append((c_, b0))
        subbp = {a, b}
        for i in range(n):
            ci, bi = gaps[i]
            for jx in range(i+1, n):
                cj, bj = gaps[jx]
                if ci != cj:
                    xc = (bj-bi)/(ci-cj)
                    if a < xc < b: subbp.add(xc)
        for u, v in zip(sorted(subbp), sorted(subbp)[1:]):
            m2 = (u+v)/2
            cb, bb = max(gaps, key=lambda cb_: cb_[0]*m2+cb_[1])
            total += (cb*u+bb + cb*v+bb)/2*(v-u)
    return total

def Emg_numeric(E, x):
    A = np.asarray(E, float)
    P = np.sort((A[None, :] * x[:, None]) % 1.0, axis=1)
    G = np.diff(P, axis=1); wrap = P[:, :1] + 1.0 - P[:, -1:]
    mg = np.maximum(G.max(axis=1), wrap[:, 0])
    return float(mg.mean()), float((mg > 1/7).mean())

def primitive(E):
    g = reduce(gcd, E)
    return tuple(sorted(e//g for e in E))

if __name__ == "__main__":
    x = (np.arange(50021) + 0.5)/50021
    print(f"T* = {TSTAR} = {float(TSTAR):.6f}   (known record 12907/65520 = {float(F(12907,65520)):.6f})")

    print("\n=== screen: F = c*{1..a} u B, c in {2,3}, |F|=13, exact-confirm the leaders ===")
    cands = []
    seen = set()
    for c in (2, 3):
        for a in range(13 - 4, 13):           # a = 9..12 base length, nb = 13-a defects 1..4
            nb = 13 - a
            base = [c*i for i in range(1, a+1)]
            pool = [v for v in range(1, c*a + c + 2) if v not in base]
            for B in combinations(pool, nb):
                Fam = sorted(base + list(B))
                pf = primitive(Fam)
                if pf in seen: continue
                seen.add(pf)
                cands.append(pf)
    print(f"  candidates (primitive-deduped): {len(cands)}")
    vals = []
    for E in cands:
        emg, mu = Emg_numeric(list(E), x)
        vals.append((emg, mu, E))
    vals.sort()
    print(f"  numeric leaders (E[mg] ascending):")
    for emg, mu, E in vals[:10]:
        print(f"    E[mg]~{emg:.6f}  mu~{mu:.4f}  {E}")
    print("\n  exact confirmation of the top 6:")
    best = (F(10), None)
    for emg, mu, E in vals[:6]:
        ex = Emaxgap_exact(list(E))
        flag = "BELOW T*  <<<< MEAN SIDECAR DIES" if ex < TSTAR else f"above T* by {float(ex-TSTAR):+.6f}"
        if ex < best[0]: best = (ex, E)
        print(f"    {str(E):>52}: EXACT = {ex} = {float(ex):.6f}  [{flag}]")
    print(f"\n  CLASS MIN (exact, this sweep): {best[0]} = {float(best[0]):.6f} at {best[1]}")
    print(f"  record check: is monad's 2*{{1..11}}u{{11,13}} in/beaten by the class? "
          f"{'BEATEN' if best[0] < F(12907,65520) else 'record stands within class'}")

    print("\n=== defective 2-lift mechanism check: winners' tail mu_{1/7} (tail-vs-mean divergence) ===")
    for emg, mu, E in vals[:4]:
        print(f"    {str(E):>52}: mean~{emg:.5f}  tail mu~{mu:.4f}  (bar m_P=0.0565: {'OK' if mu>=0.0565 else 'FAIL'})")
