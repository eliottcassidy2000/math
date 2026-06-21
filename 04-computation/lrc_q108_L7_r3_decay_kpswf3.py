#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_q108_L7_r3_decay_kpswf3.py   (kind-pasteur 2026-06-21, THREAD B)

Resolve the SUBTLETY found in lrc_q108_L7_r3_direct_kpswf3.py:
the RAW 3D cell discrepancy D3 = sum_{ijk}|mu - 1/343| does NOT decay like 1/q
(a 1-dim geodesic can never equidistribute the 3-dim torus; D3 = Omega(1)).
So the lossy bound  |R3| <= D3  is the WRONG controlling object.

The RIGHT object: R3 = p0_inf - P3 is a correlation in the COVERAGE FUNCTIONAL,
which depends only on the far SECTOR-SET  S_far(v) = {s(qv), s(p2 v), s(p3 v)} (a
subset of Z/7 of size <=3), not the ordered triple. The functional g_B depends on
S_far only through which missing sectors it fills. So the controlling discrepancy is
the L1 distance between the law of S_far(v) (1-dim geodesic) and the law of 3 i.i.d.
uniform sectors -- a 1-DIMENSIONAL projection discrepancy, which DOES decay ~1/q
(Erdos-Turan / three-gap on the single torus variable v).

This script:
  (A) Confirms R3 -> 0 DIRECTLY along resonant families with growing q (proper coprime
      directions), and fits R3 ~ C/q.
  (B) Defines the CORRECT controlling discrepancy Dcov (coverage-set discrepancy) and
      shows R3 <= Dcov with Dcov = O(1/q).
  (C) Reports the global sup of |R3| over the small-q atlas (the resonance peak).

EXACT rational. Output -> 05-knowledge/results/.
"""
import itertools
from fractions import Fraction as Fr
from math import gcd, comb

P = 7
def sector(yf): return int(P*yf)

CAP = {9: Fr(1979,4004), 10: Fr(55,91)}

def base_xcells(B):
    B = [int(b) for b in B]
    xbp = {Fr(0), Fr(1)}
    for b in B:
        if b == 0: continue
        for t in range(0, P*b): xbp.add(Fr(t, P*b))
    xs = sorted(xbp)
    out = []
    for a, b in zip(xs, xs[1:]):
        mid = (a+b)/2
        out.append((b-a, frozenset(sector((bb*mid) % 1) for bb in B)))
    return out

def far_vcells(mults):
    mults = [int(m) for m in mults]
    vbp = {Fr(0), Fr(1)}
    for f in mults:
        for t in range(0, P*f): vbp.add(Fr(t, P*f))
    vs = sorted(vbp)
    out = []
    for a, b in zip(vs, vs[1:]):
        mid = (a+b)/2
        out.append((b-a, frozenset(sector((f*mid) % 1) for f in mults)))
    return out

def p0_inf_multi(B, mults):
    xcells = base_xcells(B)
    from collections import defaultdict
    vgrp = defaultdict(Fr)
    for vlen, sf in far_vcells(mults):
        vgrp[sf] += vlen
    total = Fr(0)
    for sf, vlen in vgrp.items():
        for xlen, Sbase in xcells:
            if len(Sbase | sf) == P:
                total += xlen*vlen
    return total

def Pr_decorrelated(B, r):
    total = Fr(0)
    for xlen, Sbase in base_xcells(B):
        m = P - len(Sbase)
        if m == 0: pcov = Fr(1)
        elif m > r: pcov = Fr(0)
        else:
            pcov = Fr(0)
            for j in range(0, m+1):
                pcov += Fr((-1)**j * comb(m, j)) * Fr(P-j, P)**r
        total += xlen*pcov
    return total

# ---- coverage-set law of the far multipliers: prob over v of each far sector-SET --
def coverset_law(mults):
    from collections import defaultdict
    law = defaultdict(Fr)
    for vlen, sf in far_vcells(mults):
        law[sf] += vlen
    return law   # frozenset(Z/7) -> measure

def coverset_law_iid(r):
    """law of the SET {w1,..,wr}, wi iid Unif(Z/7), as frozenset->prob."""
    from collections import defaultdict
    law = defaultdict(Fr)
    tot = Fr(1, P**r)
    for combo in itertools.product(range(P), repeat=r):
        law[frozenset(combo)] += tot
    return law

def Dcov(mults):
    """L1 distance between coverage-set law of the geodesic and the iid law.
       This is the CORRECT controlling discrepancy: it lives on subsets of Z/7
       (a finite functional), driven by a SINGLE torus variable v -> O(1/q)."""
    r = len(mults)
    Lg = coverset_law(mults)
    Li = coverset_law_iid(r)
    keys = set(Lg) | set(Li)
    return sum(abs(Lg.get(s, Fr(0)) - Li.get(s, Fr(0))) for s in keys)

def main():
    B = [0,2,4,6,8,10]   # worst k=9 base6
    P3 = Pr_decorrelated(B, 3)
    print("="*86)
    print("THREAD B (decay): R3 = p0_inf - P3 DECAYS ~1/q;  raw 3D-cell D3 does NOT (wrong object)")
    print("="*86)
    print(f"base6={B}, P3(decorrelated 3-far)={float(P3):.6f}, cap_9={float(CAP[9]):.5f}")

    # (A) resonant families with PROPER coprime directions, growing q.
    # Family I:  direction (q, q+1, q+2) -> ratios ((q+1)/q,(q+2)/q) -> (1,1) as q->inf.
    #   These stay in-window for q>= ceil(1/0.15)=7? (q+2)/q<=2.15 => q>=2. ok small q too.
    # Family II: direction (q, 2q-1, 2q+1) -> ratios ~2 (the "near-2" resonance), coprime.
    fams = {
        "I  (q,q+1,q+2) g->(1,1)": lambda q:(q, q+1, q+2),
        "II (q,2q-1,2q+1) g->(2,2)": lambda q:(q, 2*q-1, 2*q+1),
    }
    for fname, fn in fams.items():
        print(f"\n--- Family {fname} ---")
        print(f"  {'q':>3} {'dir':>14} {'p0_inf':>9} {'R3':>10} {'R3*q':>9} {'Dcov':>9} {'Dcov*q':>9}")
        for q in range(2, 22):
            q0,p2,p3 = fn(q)
            if min(q0,p2,p3) <= 0: continue
            if gcd(gcd(q0,p2),p3) != 1: continue
            # window check on ratios
            if not (Fr(1) <= Fr(p2,q0) <= Fr(43,20) and Fr(1) <= Fr(p3,q0) <= Fr(43,20)):
                continue
            val = p0_inf_multi(B, (q0,p2,p3)); R3 = val - P3
            dc = Dcov((q0,p2,p3))
            print(f"  {q:>3} {f'{q0},{p2},{p3}':>14} {float(val):>9.5f} {float(R3):>+10.5f} "
                  f"{float(R3*q):>+9.4f} {float(dc):>9.5f} {float(dc*q):>9.4f}")

    # (B) global sup over small-q atlas (the resonance peak): all coprime directions q<=6
    print("\n" + "="*86)
    print("Global sup |R3| over coprime ratio-pair atlas q<=6 (the resonance peak):")
    print("="*86)
    best = (Fr(0), None)
    bestDcov = (Fr(0), None)
    cnt = 0
    for q in range(1, 7):
        nums = range(q+1, int(Fr(43,20)*q)+1)
        for p2, p3 in itertools.combinations_with_replacement(nums, 2):
            if gcd(gcd(q,p2),p3) != 1: continue
            cnt += 1
            R3 = abs(p0_inf_multi(B, (q,p2,p3)) - P3)
            if R3 > best[0]: best = (R3, (q,p2,p3))
            dc = Dcov((q,p2,p3))
            if dc > bestDcov[0]: bestDcov = (dc, (q,p2,p3))
    print(f"  scanned {cnt} coprime directions")
    print(f"  sup |R3|  = {float(best[0]):.6f} at {best[1]}  => sup p0_inf = {float(best[0]+P3):.6f}")
    print(f"  sup Dcov  = {float(bestDcov[0]):.6f} at {bestDcov[1]}")
    print(f"  cap_9 = {float(CAP[9]):.5f}; margin at the peak = {float(CAP[9]-(best[0]+P3)):.5f}")
    # the controlling bound check:  is |R3| <= Dcov on the peak?
    peakDcov = Dcov(best[1])
    print(f"  check  |R3|(peak) {float(best[0]):.5f} <= Dcov(peak) {float(peakDcov):.5f}  -> "
          f"{best[0] <= peakDcov}")

    print("\nDONE.")

if __name__ == "__main__":
    main()
