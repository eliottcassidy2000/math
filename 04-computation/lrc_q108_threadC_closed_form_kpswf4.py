#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_q108_threadC_closed_form_kpswf4.py  (kind-pasteur 2026-06-21, THREAD C)

Seek a CLOSED FORM for r_j (row 0) and prove sum_j|7 r_j - pq| <= 12p combinatorially.

Reformulation as point counting.  The pq points are
    P = { 7 frac(pk/q) + (2t+1)/(2q) mod 7 : 0<=k<q, 0<=t<p }.
r_j = #(P cap [j, j+1)).

Equivalent lattice description.  Scale by 2q: the points 2q*x are
    y_{k,t} = 14 (pk mod q) + (2t+1)   (an ODD integer), taken mod 14q.
Bin j  <=>  2q j <= (y mod 14q) < 2q(j+1)  <=>  floor( (y mod 14q) / (2q) ) = j.
So r_j counts odd integers y in the multiset { 14(pk mod q) + (2t+1) } landing in
the length-2q window [2qj, 2q(j+1)) mod 14q.

Since {pk mod q : k} = {0,..,q-1} (gcd=1), and (2t+1) ranges over odds 1,3,..,2p-1,
    P-scaled = { 14 a + b : a in Z/q, b in {1,3,...,2p-1} }  mod 14q.
For fixed a, the b-values 14a+1, 14a+3, ..., 14a+(2p-1) are p consecutive odd integers
=> an arithmetic progression step 2, length p, starting at 14a+1.  This is the arc.

Now r_j = sum_{a=0}^{q-1} #{ t in [0,p) : 14a + 2t+1 in [2qj, 2q(j+1)) (mod 14q) }.
The inner count = # of the p-term AP (step 2, from 14a+1) inside a length-2q window.
A length-2q window contains floor(2q/2)=q or q+/-... of a step-2 AP -- about q points,
but capped at p (only p terms exist) and shifted by the window position.

GOAL of this script: get the EXACT inner count formula and the global r_j, then the bound.
We brute force to discover the pattern, focusing on the two regimes a=floor(p/q) in {1,2}.
"""
from fractions import Fraction as Fr
from math import gcd

P = 7

def r_direct(p, q):
    r = [0] * P
    for k in range(q):
        base = 14 * ((p * k) % q)
        for t in range(p):
            j = ((base + 2 * t + 1) // (2 * q)) % P
            r[j] += 1
    return r

def r_via_a(p, q):
    """same thing but summing over a = pk mod q directly (a runs over all residues)."""
    r = [0] * P
    for a in range(q):
        base = 14 * a
        for t in range(p):
            j = ((base + 2 * t + 1) // (2 * q)) % P
            r[j] += 1
    return r

def main():
    print("THREAD C: closed-form hunt for r_j")
    print("=" * 72)
    # confirm r_via_a == r_direct (reindexing a=pk mod q is a bijection)
    ok = all(r_direct(p, q) == r_via_a(p, q)
             for q in range(1, 40) for p in range(q + 1, int(Fr(43, 20) * q) + 1)
             if gcd(p, q) == 1 and Fr(1) < Fr(p, q) <= Fr(43, 20))
    print("r_direct == r_via_a (a=pk mod q reindex):", "YES" if ok else "NO")

    # Now r_via_a no longer depends on p EXCEPT through the count p and the modulus interplay!
    # Actually r_via_a(p,q) sums over ALL a in Z/q, b in odds<2p, of bin(14a+b).
    # This is purely a function of (p, q) lattice geometry. Let's see the e-vector and its
    # decomposition: e_j = 7 r_j - pq.
    print("\nThe map a -> floor((14a + 2t+1)/(2q)) mod 7 -- structure of contributions.")
    print("Count, for each residue class, how many (a,t) land in bin j.")
    # Key: 14a + 2t + 1 ranges over an interval of integers as (a,t) vary.
    # full set of values {14a+2t+1 : 0<=a<q, 0<=t<p} -- these are NOT all distinct mod 14q.
    # Let's see the bound directly and the per-(p,q) e-vector with sum|e|.
    print("\np/q : e=7r-pq : sum|e| : 12p : 20q : sum|e|<=min?")
    rows = []
    for q in range(1, 50):
        for p in range(q + 1, int(Fr(43, 20) * q) + 1):
            if gcd(p, q) != 1 or not (Fr(1) < Fr(p, q) <= Fr(43, 20)):
                continue
            r = r_direct(p, q)
            e = [7 * x - p * q for x in r]
            S = sum(abs(x) for x in e)
            rows.append((p, q, e, S, 12 * p, 20 * q))
    # print extremes (largest sum|e|/p)
    rows.sort(key=lambda t: Fr(t[3], t[0]), reverse=True)
    for p, q, e, S, f12, f20 in rows[:12]:
        print(f"  {p}/{q:<3d} e={e} S={S} 12p={f12} 20q={f20} "
              f"S/p={S/p:.4f} ok={S<=min(f12,f20)}")

    # observe: which e-vectors appear? all entries multiples of something? gcd of e?
    print("\ngcd of e-vector entries (suggests a quantum):")
    from math import gcd as g
    for p, q, e, S, _, _ in rows[:8]:
        gg = 0
        for x in e:
            gg = g(gg, abs(x))
        print(f"  {p}/{q}: e={e} gcd={gg}")

if __name__ == "__main__":
    main()
