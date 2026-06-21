#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_krawtchouk_convex_order_macmini.py  (mac-mini 2026-06-21, THREAD D, lead ii decisive)

DECISIVE TEST of the Krawtchouk-moment governing-scalar lead (HYP-2737 candidate).

Finding so far: E[K_4(N)] has Kendall concordance 0.862 with measS7 (>> additive
energy 0.691); E[K_2(N)] has 0.743 and K_2 is EXACTLY CONVEX (2nd diff = 4), so
E[K_2(N)] is a SCHUR-CONVEX functional => a convex-order/majorization argument can
prove its maximizer WITHOUT a finite atlas.

DECISIVE QUESTIONS:
 (Q1) Is consec [1..k] the UNIQUE maximizer of the convex moment E[K_2(N)] over
      primitive k-sets?  If YES, AND we can show measS7 is monotone in E[K_2(N)] (or
      that the LP dual = nonneg combo of K_j with K_2 dominant), we get a non-finite
      route.  If NO (some set has larger E[K_2(N)] but smaller measS7), K_2 alone is
      insufficient and we need the full nonneg Krawtchouk combo (= the Delsarte LP).
 (Q2) The Delsarte dual at k=8 is c=(10,0,0,1,0,0,10)/scale (HYP-2726d): the bound is
      10*q0 <= 10 K_0 + 1 K_3 + 10 K_6 ... wait that's the DUAL of the q-vector. Let me
      directly form the moment L_y = sum_j c_j E[K_j(N)] with the VERIFIED nonneg c_j
      and check: is consec the maximizer of L_y? (THM-534 says p0 <= L_y per-E; the
      crux is consec maximizes L_y.)  If consec maximizes EACH nonneg-weighted K_j
      moment that appears, done.
 (Q3) Recover the equal-AE pair (consec vs odd-AP[1,3..15]) at max=15 and confirm
      which Krawtchouk moment separates them in the right direction.

This is the concrete variational target for closing consec-max non-finitely.
"""
import itertools
from fractions import Fraction as Fr
from math import gcd, comb
from collections import defaultdict

P = 7
def sector(yf): return int(P * yf)

def breakpoints(E):
    bp = {Fr(0), Fr(1)}
    for e in E:
        if e == 0: continue
        for t in range(0, P * e):
            bp.add(Fr(t, P * e))
    return sorted(bp)

def empty_count_law(E):
    E = [int(e) for e in E if int(e) != 0]
    xs = breakpoints(E); law = defaultdict(lambda: Fr(0))
    for a, b in zip(xs, xs[1:]):
        mid = (a + b) / 2
        N = P - len(set(sector((e * mid) % 1) for e in E))
        law[N] += (b - a)
    return law

def krawtchouk(j, x, m=6):
    return sum((-1)**i * comb(x, i) * comb(m - x, j - i) for i in range(j + 1))

def Kmoment(E, j):
    law = empty_count_law(E)
    return sum(krawtchouk(j, N) * m for N, m in law.items())

def measS7(E):
    law = empty_count_law(E)
    return law.get(0, Fr(0))  # N=0 means 0 empty => all 7 hit

def main():
    print("#"*80)
    print("# KRAWTCHOUK CONVEX-ORDER decisive test (THREAD D)")
    print("#"*80)

    consec = tuple(range(1, 9))

    # population: primitive 8-sets min=1, max<=15 (recovers odd-AP[1,3..15])
    print("\nbuilding primitive 8-sets, min=1, max<=15 ...")
    rows = []
    for combo in itertools.combinations(range(2, 16), 7):
        E = (1,) + combo
        g = 0
        for e in E: g = gcd(g, e)
        if g != 1: continue
        rows.append(E)
    print(f"population = {len(rows)}")

    # precompute moments + p0
    data = {}
    for E in rows:
        law = empty_count_law(E)
        p0 = law.get(0, Fr(0))
        data[E] = {
            "p0": p0,
            "K2": sum(krawtchouk(2, N) * m for N, m in law.items()),
            "K4": sum(krawtchouk(4, N) * m for N, m in law.items()),
            "K6": sum(krawtchouk(6, N) * m for N, m in law.items()),
            "K0": sum(krawtchouk(0, N) * m for N, m in law.items()),
            "K3": sum(krawtchouk(3, N) * m for N, m in law.items()),
        }

    # (Q1) is consec the maximizer of E[K_2(N)] (convex moment)?
    print("\n=== Q1: maximizer of the CONVEX moment E[K_2(N)] ===")
    bestK2 = max(rows, key=lambda E: data[E]["K2"])
    print(f"  argmax E[K2(N)] = {bestK2}  K2={float(data[bestK2]['K2']):.5f}  p0={float(data[bestK2]['p0']):.5f}")
    print(f"  consec K2 = {float(data[consec]['K2']):.5f}  p0={float(data[consec]['p0']):.5f}")
    print(f"  consec maximizes E[K2(N)]? {bestK2 == consec}")
    # top 4 by K2
    for E in sorted(rows, key=lambda E:-data[E]["K2"])[:4]:
        print(f"    K2={float(data[E]['K2']):.5f} p0={float(data[E]['p0']):.5f} E={E}")

    # (Q2) the VERIFIED k=8 Delsarte dual: 10*q0 <= 10 K0 + 1 K3 + 10 K6  (HYP-2726d gK8)
    #   so the moment majorant is L = (10 K0 + K3 + 10 K6)/10. Is consec the maximizer of L?
    print("\n=== Q2: maximizer of the k=8 Delsarte moment majorant L = (10 K0 + K3 + 10 K6)/10 ===")
    def Lmaj(E):
        d = data[E]
        return (10*d["K0"] + d["K3"] + 10*d["K6"]) / 10
    bestL = max(rows, key=lambda E: Lmaj(E))
    print(f"  argmax L = {bestL}  L={float(Lmaj(bestL)):.5f}")
    print(f"  consec L = {float(Lmaj(consec)):.5f}")
    print(f"  consec maximizes L (the Delsarte majorant)? {bestL == consec}")
    for E in sorted(rows, key=lambda E:-Lmaj(E))[:4]:
        print(f"    L={float(Lmaj(E)):.5f} p0={float(data[E]['p0']):.5f} E={E}")
    # sanity: L >= p0 per E? (Delsarte majorant)
    viol = sum(1 for E in rows if Lmaj(E) < data[E]["p0"] - Fr(1,10**12))
    print(f"  Delsarte majorant L >= p0 violations: {viol} (should be 0)")

    # (Q3) the equal-AE pair separation
    print("\n=== Q3: equal-AE pair consec[1..8] vs odd-AP[1,3..15] ===")
    odap = (1,3,5,7,9,11,13,15)
    for name, E in [("consec", consec), ("oddAP", odap)]:
        d = data[E]
        print(f"  {name:7s} p0={float(d['p0']):.5f} K2={float(d['K2']):.4f} "
              f"K4={float(d['K4']):.4f} K6={float(d['K6']):.4f}")
    dc, do = data[consec], data[odap]
    print(f"  separation: K2 favors consec? {dc['K2']>do['K2']}  "
          f"K4 favors consec? {dc['K4']>do['K4']}  K6 favors consec? {dc['K6']>do['K6']}")

    # (Q4) GLOBAL: does consec maximize p0? and is the FULL nonneg-Krawtchouk combo
    #      (Delsarte LP majorant) ALSO maximized at consec?  (the clean claim)
    print("\n=== Q4: global p0 maximizer & is it consec? ===")
    bestp0 = max(rows, key=lambda E: data[E]["p0"])
    print(f"  global argmax p0 = {bestp0} p0={float(data[bestp0]['p0']):.5f}  "
          f"{'== consec' if bestp0==consec else 'NOT consec!'}")

    print("\nDONE.")

if __name__ == "__main__":
    main()
