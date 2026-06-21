#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_q108_threadC_row0_structure_kpswf4.py  (kind-pasteur 2026-06-21, THREAD C)

Goal: understand the integer row-0 vector r_j = c[0][j] well enough to prove
    sum_j |r_j - pq/7| <= 12p/7    (equivalently D <= 12/(7q)).

ROW-0 STRUCTURE (band decomposition).
Row 0 = { v : sector(qv)=0 } = disjoint union of q bands B_k = [k/q, k/q+1/(7q)), k=0..q-1.
On the 7pq grid, band B_k consists of exactly p consecutive subintervals (length 1/(7q) = p/(7pq)).
So r_j = (# subintervals in row-0 bands whose p-sector is j).
Total subintervals in row 0 = q*p = pq = rowsum. Good.

On band B_k, parameter v=k/q + w/(7q), w in [0,1). Then
    7 p v = 7pk/q + (p/q) w.
The p-sector j = floor(7 p v) mod 7 = floor(7pk/q + (p/q)w) mod 7.
As w runs over the p subintervals (w = (2t+1)/(2p), t=0..p-1, midpoints), the integer part
floor(7pk/q + (p/q)w) takes p consecutive(ish) integer values mod 7.

KEY: define the band STARTING phase. Let A_k = 7pk/q mod 7 = 7p k/q  reduced mod 7
(a value in [0,7)). The p subintervals of band k have 7pv values
    A_k + (p/q)*w,  w in [0,1)  ->  sweep [A_k, A_k + p/q).
So band k contributes a CONTIGUOUS ARC (mod 7) of length p/q to the "7pv mod 7" circle,
sampled at p equally spaced midpoints. Its j-histogram = how a length-(p/q) arc starting at
phase A_k, sampled by p equal points, distributes over the 7 unit-sectors.

This is a 1D *interval-on-circle* occupancy: r = sum over q arcs, each arc length p/q,
starting phases A_k = frac7(7 p k/q) = 7*frac(pk/q), k=0..q-1.
Because gcd(p,q)=1, {pk mod q} = {0..q-1}, so {frac(pk/q)} = {0,1/q,...,(q-1)/q} -- the
q starting phases A_k are the q points {7*m/q : m=0..q-1} = {0, 7/q, 14/q, ...} EQUALLY
SPACED on [0,7) with gap 7/q.

So ROW 0 IS: q equally spaced arcs (gap 7/q on the length-7 circle), EACH of length p/q,
each sampled by p equal subinterval-midpoints, accumulated into 7 unit bins.

We verify this exact reformulation and print the arc model, to set up the sharp bound.
"""
from fractions import Fraction as Fr
from math import gcd

P = 7

def cmatrix_row0(p, q):
    """row 0 of the integer cell matrix on the 7pq grid."""
    N = 7 * p * q
    r = [0] * P
    # row 0: subintervals k with sector(qv)=0, v=(2k+1)/(2N)
    for k in range(N):
        num_q = q * (2 * k + 1)
        i = (P * (num_q % (2 * N))) // (2 * N)
        if i != 0:
            continue
        num_p = p * (2 * k + 1)
        j = (P * (num_p % (2 * N))) // (2 * N)
        r[j] += 1
    return r

def row0_arc_model(p, q):
    """Independent computation via the arc model:
       q arcs, starting phases A_k = 7*frac(p k / q) (equally spaced, gap 7/q),
       each arc of length p/q on the length-7 circle, sampled at p midpoints
       at offsets (2t+1)/(2p) * (p/q) = (2t+1)/(2q), t=0..p-1, added to A_k.
       Bin = floor( (A_k + offset) mod 7 ).
    """
    r = [0] * P
    for k in range(q):
        Ak = (7 * Fr(p * k, q)) % 7   # in [0,7)
        for t in range(p):
            # midpoint of t-th subinterval: 7pv = A_k + (p/q)*((2t+1)/(2p))
            val = Ak + Fr(p, q) * Fr(2 * t + 1, 2 * p)   # = A_k + (2t+1)/(2q)
            j = int(val % 7)
            r[j] += 1
    return r

def main():
    print("THREAD C: row-0 ARC MODEL")
    print("=" * 70)
    print("Row 0 = q arcs, equally spaced starts (gap 7/q), each length p/q,")
    print("sampled by p midpoints. Verify arc model == direct grid row 0.")
    print()
    ok = True
    for q in range(1, 30):
        for p in range(q + 1, int(Fr(43, 20) * q) + 1):
            if gcd(p, q) != 1 or not (Fr(1) < Fr(p, q) <= Fr(43, 20)):
                continue
            if row0_arc_model(p, q) != cmatrix_row0(p, q):
                ok = False
                print("MISMATCH", p, q, row0_arc_model(p, q), cmatrix_row0(p, q))
    print("arc model == grid row 0:", "YES" if ok else "NO")

    # Insight: the arc has length p/q. Write p = a*q + b, 0<b<q (since 1<p/q<=43/20<3,
    # a in {1,2}). p/q = a + b/q. So each arc covers EXACTLY a full unit-bins worth
    # of length plus a fractional b/q. The sampling by p points...
    print("\nArc length p/q = a + b/q  decomposition (a=floor(p/q) in {1,2}):")
    for (p, q) in [(3,2),(2,1),(4,3),(5,3),(5,4),(9,5),(11,10),(17,10)]:
        a, b = divmod(p, q)
        print(f"  p/q={p}/{q}: a={a} b={b}  arclen={float(Fr(p,q)):.4f}")

    # The contribution of one arc of length p/q (sampled by p midpoints) to the 7 bins:
    # since the arc is length < 3 and bins are length 1, it touches at most ceil(p/q)+1 bins.
    print("\nPer-arc bin histograms (showing arcs touch few bins):")
    for (p, q) in [(3,2),(5,3),(9,5)]:
        print(f" p/q={p}/{q}:")
        for k in range(q):
            r = [0]*P
            Ak = (7*Fr(p*k,q))%7
            for t in range(p):
                val = Ak + Fr(p,q)*Fr(2*t+1,2*p)
                r[int(val%7)] += 1
            touched = sum(1 for x in r if x>0)
            print(f"   arc k={k}: start phase={float(Ak):.3f} hist={r} touches {touched} bins")

if __name__ == "__main__":
    main()
