#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_q108_threadC_sharp_bound_kpswf4.py  (kind-pasteur 2026-06-21, THREAD C)

THE SHARP 1D BOUND.  We have reduced D_{p,q} = (1/(pq)) * S, where
    S = S(p,q) = sum_{j=0}^{6} | r_j - pq/7 |,
and r = row 0 = histogram of the pq points  x_{k,t} = 7 frac(pk/q) + (2t+1)/(2q) (mod 7)
into the 7 unit bins.  Target: D <= 12/(7q)  <=>  S <= 12 p/7  <=>  7 S <= 12 p.
Since 7 r_j and pq/7*7 = pq are involved, let e_j = 7 r_j - pq (integer, sum_j e_j = 0).
Then 7 S = sum_j |7 r_j - pq| = sum_j |e_j|.  So the SHARP target is

        SHARP:   sum_j |e_j|  <=  12 p,     where e_j = 7 r_j - pq,  sum e_j = 0.

We also observed sup(D p) = 20/7 i.e. sum_j|e_j| <= 20 q is the OTHER sharp face.
This script:
  (1) confirms 7S = sum|e_j| and checks the two sharp faces 12p and 20q over a big atlas;
  (2) prints the EXACT achievers of each face;
  (3) decomposes S over the q arcs to expose WHY 12 appears (each arc's contribution).
"""
from fractions import Fraction as Fr
from math import gcd

P = 7

def row0(p, q):
    r = [0] * P
    for k in range(q):
        Ak = (7 * Fr(p * k, q)) % 7
        for t in range(p):
            x = (Ak + Fr(2 * t + 1, 2 * q)) % 7
            r[int(x)] += 1
    return r

def evec(p, q):
    r = row0(p, q)
    return [7 * rj - p * q for rj in r]   # integers, sum 0

def main():
    print("THREAD C: the sharp 1D bound  sum_j |e_j| <= 12 p   (e_j = 7 r_j - pq)")
    print("=" * 72)
    QMAX = 200
    viol_p = []   # sum|e| > 12 p
    viol_q = []   # sum|e| > 20 q
    achieve_p = []  # sum|e| == 12p
    achieve_q = []  # sum|e| == 20q
    best_pface = (0, None)  # max sum|e|/p
    best_qface = (0, None)
    n = 0
    for q in range(1, QMAX + 1):
        for p in range(q + 1, int(Fr(43, 20) * q) + 1):
            if gcd(p, q) != 1 or not (Fr(1) < Fr(p, q) <= Fr(43, 20)):
                continue
            e = evec(p, q)
            S7 = sum(abs(x) for x in e)
            n += 1
            if S7 > 12 * p:
                viol_p.append((p, q, S7))
            if S7 > 20 * q:
                viol_q.append((p, q, S7))
            if S7 == 12 * p:
                achieve_p.append((p, q))
            if S7 == 20 * q:
                achieve_q.append((p, q))
            rp = Fr(S7, p)
            if rp > best_pface[0]:
                best_pface = (rp, (p, q, S7))
            rq = Fr(S7, q)
            if rq > best_qface[0]:
                best_qface = (rq, (p, q, S7))
    print(f"checked {n} window ratios, q<={QMAX}")
    print(f"  violations of sum|e| <= 12p : {len(viol_p)}  {viol_p[:5]}")
    print(f"  violations of sum|e| <= 20q : {len(viol_q)}  {viol_q[:5]}")
    print(f"  max sum|e|/p = {best_pface[0]} ({float(best_pface[0]):.5f}) at {best_pface[1]}  (target sup=12)")
    print(f"  max sum|e|/q = {best_qface[0]} ({float(best_qface[0]):.5f}) at {best_qface[1]}  (target sup=20)")
    print(f"  ratios achieving sum|e|=12p (face P): {achieve_p[:12]}{'...' if len(achieve_p)>12 else ''} (count {len(achieve_p)})")
    print(f"  ratios achieving sum|e|=20q (face Q): {achieve_q[:12]}{'...' if len(achieve_q)>12 else ''} (count {len(achieve_q)})")

    # ---- arc-decomposition of S to expose the constant 12 --------------------
    print("\n" + "=" * 72)
    print("PER-ARC view: each arc contributes to r; e = sum over q arcs of (arc hist*7 - p)")
    print("Show e and the e-vector for the face-P achievers.")
    print("=" * 72)
    for (p, q) in achieve_p[:6]:
        e = evec(p, q)
        print(f"  p/q={p}/{q}: e={e} sum|e|={sum(abs(x) for x in e)} =12p={12*p}? {sum(abs(x) for x in e)==12*p}")

    # ---- the extreme ratio 3/2 in detail ------------------------------------
    print("\nExtreme ratio p/q=3/2 (the global sup D*q=12/7):")
    r = row0(3, 2); e = evec(3, 2)
    print(f"  r={r}  pq/7={float(Fr(6,7)):.4f}  e=7r-6={e}  sum|e|={sum(abs(x) for x in e)}  12p=36  20q=40")
    print(f"  D = sum|e|/(7pq) = {Fr(sum(abs(x) for x in e), 7*6)} = {float(Fr(sum(abs(x) for x in e),42)):.5f}")

if __name__ == "__main__":
    main()
