#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_q108_threadC_arc_contrib_kpswf4.py  (kind-pasteur 2026-06-21, THREAD C)

WHERE DOES 12 COME FROM?  Decompose e = 7r - pq into per-arc contributions.

Each arc k contributes a histogram h^k_j (j=0..6) of its p points to the 7 bins,
sum_j h^k_j = p.  Define the per-arc deviation  d^k_j = 7 h^k_j - p  (sum_j d^k_j = 0).
Then r_j = sum_k h^k_j, and e_j = 7 r_j - pq = sum_k (7 h^k_j - p) = sum_k d^k_j.

CLAIM A (per-arc total variation): for EVERY arc,  sum_j |d^k_j| <= 12 ... ? test.
If sum_j |d^k_j| <= 12 for each of the q arcs, then by triangle inequality on e_j=sum_k d^k_j,
   sum_j |e_j| <= sum_k sum_j |d^k_j| <= 12 q.
That gives the 20q? no, 12q -> D <= 12/(7q)... wait. Let's see: sum|e| <= 12q gives
D = sum|e|/(7pq) <= 12q/(7pq) = 12/(7p). Hmm that is the OTHER scaling (per-p).

Let's just MEASURE: per-arc sum_j |d^k_j| (call it tv_k). Its max and the achiever.
Also measure the BAND view differently. The point: each arc is a length-(p/q) sub-interval
of the circle sampled by p equally spaced points (spacing 1/q). A length-L interval sampled
by p equal points lands in floor or ceil of (its bin overlaps). Quantify tv per arc.
"""
from fractions import Fraction as Fr
from math import gcd

P = 7

def arc_hists(p, q):
    """list over k of the per-arc bin histogram h^k (length 7), sum = p each."""
    out = []
    for k in range(q):
        base = 14 * ((p * k) % q)
        h = [0] * P
        for t in range(p):
            val = base + (2 * t + 1)
            j = (val // (2 * q)) % P
            h[j] += 1
        out.append(h)
    return out

def main():
    print("THREAD C: per-arc deviation total variation tv_k = sum_j |7 h^k_j - p|")
    print("=" * 72)
    max_tv = (0, 1, None)
    # also collect distribution of tv_k values
    from collections import Counter
    cnt = Counter()
    for q in range(1, 120):
        for p in range(q + 1, int(Fr(43, 20) * q) + 1):
            if gcd(p, q) != 1 or not (Fr(1) < Fr(p, q) <= Fr(43, 20)):
                continue
            for h in arc_hists(p, q):
                tv = sum(abs(7 * x - p) for x in h)
                # normalize tv/p to find sup
                if tv * max_tv[1] > max_tv[0] * p:
                    max_tv = (tv, p, (p, q, h, tv))
                cnt[Fr(tv, p)] += 1
    print(f"  sup over all arcs of tv_k / p = {max_tv[0]}/{max_tv[1]} = {max_tv[0]/max_tv[1]:.5f}")
    print(f"     at p/q={max_tv[2][0]}/{max_tv[2][1]} arc hist={max_tv[2][2]} tv={max_tv[2][3]}")
    print(f"  Note: bound sum_j|e_j| <= sum_k tv_k <= q * sup(tv) is TOO LOSSY (cancellation).")

    # Better: the # of bins each arc touches, and that consecutive arcs' deviations CANCEL.
    # Let's look at the e-vector vs the sum of |d^k|.
    print("\nCancellation factor: sum_j|e_j|  vs  sum_k sum_j|d^k_j|  (ratio):")
    for (p, q) in [(3,2),(5,3),(9,5),(5,4),(11,10),(43,20),(40,19),(23,13)]:
        if gcd(p,q)!=1: continue
        H = arc_hists(p,q)
        e = [7*sum(H[k][j] for k in range(q)) - p*q for j in range(P)]
        Se = sum(abs(x) for x in e)
        Sd = sum(sum(abs(7*H[k][j]-p) for j in range(P)) for k in range(q))
        print(f"  p/q={p}/{q:<3d} sum|e|={Se}  sum_k tv_k={Sd}  ratio={Se/Sd:.4f}  "
              f"sum|e|/p={Se/p:.4f}")

    # The REAL structure: each arc is a contiguous interval; its deviation d^k is supported on
    # ~3 consecutive bins. Print the support pattern.
    print("\nSupport pattern of per-arc deviation d^k (which bins nonzero), p/q=9/5:")
    p,q=9,5
    for k,h in enumerate(arc_hists(p,q)):
        d=[7*x-p for x in h]
        supp=[j for j in range(P) if h[j]>0]
        print(f"  arc {k}: h={h} d={d} support bins={supp} (contiguous mod 7)")

if __name__ == "__main__":
    main()
