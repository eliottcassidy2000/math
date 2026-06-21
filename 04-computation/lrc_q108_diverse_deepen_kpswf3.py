#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_q108_diverse_deepen_kpswf3.py   (kind-pasteur 2026-06-21, THREAD D follow-up)

Deepen the angles that showed signal in lrc_q108_diverse_angles_kpswf3.py:

 (A) SCALING INVARIANCE: measS7(c*E) = measS7(E)?  ([1..8] == [2,4..16] exactly.)
     -> if true, the cover bound only depends on E up to GLOBAL scaling => reduces the
        search space and is a genuine structural lemma. PROVE-check exactly.

 (B) BONFERRONI / inclusion-exclusion structure.  B2 upper bound was useless (>1) but
     the TRUE p0 has an EXACT inclusion-exclusion p0 = sum_{j} (-1)^j S_j where S_j =
     sum over j-subsets T of (measure all sectors in T empty).  Test which truncation
     ALTERNATES and brackets, and whether a SMARTER even/odd Bonferroni on the
     COMPLEMENT (miss) lattice gives a usable cap.  Also: the miss-lattice is graded by
     |T|; compute the full S_j vector and the "miss polynomial".

 (C) TOURNAMENT-SIDE collision count: is p0 monotone / bounded by E[#collisions]?
     The decorrelated coupon depends only on k; the EXCESS corr correlates with the
     element structure. Test corr vs E[#collisions] and vs the GCD-structure.

 (D) The complement/miss generating function M(z) = sum_j S_j z^j and p0 = (-1 eval).
     Is M(z) a nice product (independence across sectors)?  QR(7) test for anomaly.
"""
import itertools
from fractions import Fraction as Fr
from math import gcd, comb

P = 7
def sector(yf): return int(P * yf)

def breakpoints(E):
    bp = {Fr(0), Fr(1)}
    for e in E:
        if e == 0: continue
        for t in range(0, P * e): bp.add(Fr(t, P * e))
    return sorted(bp)

def measS7(E):
    E = [int(e) for e in E if int(e) != 0]
    xs = breakpoints(E); total = Fr(0)
    for a, b in zip(xs, xs[1:]):
        mid = (a + b) / 2
        if len(set(sector((e*mid)%1) for e in E)) == P: total += (b-a)
    return total

def miss_prob(E, missing):
    E = [int(e) for e in E if int(e) != 0]
    xs = breakpoints(E); miss = set(missing); total = Fr(0)
    for a, b in zip(xs, xs[1:]):
        mid = (a + b) / 2
        secs = set(sector((e*mid)%1) for e in E)
        if not (secs & miss): total += (b-a)
    return total

def Svec(E):
    """S_j = sum over j-subsets T of Z/7 of measure( all sectors in T empty ).
       p0 = sum_j (-1)^j S_j  (inclusion-exclusion on 'sector empty')."""
    S = [Fr(0)] * (P + 1)
    S[0] = Fr(1)
    for j in range(1, P + 1):
        for T in itertools.combinations(range(P), j):
            S[j] += miss_prob(E, set(T))
    return S

# ----------------------------------------------------------------------------
def main():
    print("#" * 80)
    print("# THREAD D follow-up: scaling, Bonferroni lattice, collision, miss-GF")
    print("#" * 80)
    CAP = {8: Fr(2243,5880), 9: Fr(1979,4004), 10: Fr(55,91)}

    # ---- (A) SCALING INVARIANCE ----
    print("\n=== (A) SCALING INVARIANCE  measS7(c*E) == measS7(E) ? ===")
    tests = [
        [1,2,3,4,5],
        [1,3,4,5,9],         # QR(7) set itself shifted
        [2,3,5,7,11],
        [1,2,4,8,16],
    ]
    for E in tests:
        base = measS7(E)
        ok = True
        scaled = {}
        for c in [2,3,5,7]:
            v = measS7([c*e for e in E])
            scaled[c] = v
            if v != base: ok = False
        # also test c=7 (a multiple of P) which may DIFFER
        v7 = measS7([7*e for e in E])
        print(f"  E={E}: base={float(base):.5f}  c2={float(scaled[2]):.5f}"
              f" c3={float(scaled[3]):.5f} c5={float(scaled[5]):.5f}"
              f"  scale-inv(c coprime 7)={'YES' if ok else 'NO'}")
        if v7 != base:
            print(f"      NOTE c=7: measS7(7E)={float(v7):.5f} != base  (7|c breaks it, as expected)")
    print("  => If YES universally for gcd(c,7)=1: measS7 depends only on E mod global")
    print("     scaling by a 7-coprime unit. Big search-space reduction.")

    # ---- (B)+(D) BONFERRONI / MISS GENERATING FUNCTION ----
    print("\n=== (B/D) MISS-LATTICE  S_j  and bracketing Bonferroni ===")
    bases = {
        "k=8 [1..8]":   [1,2,3,4,5,6,7,8],
        "k=9 [1..9]":   [1,2,3,4,5,6,7,8,9],
        "k=9 2clust":   [1,2,3,4, 50,51,52,53,54],
        "k=8 QR-fl":    [1,2,3,4,5,6,7,9],
    }
    for name, E in bases.items():
        S = Svec(E)
        # inclusion-exclusion partial sums P_t = sum_{j<=t} (-1)^j S_j
        partial = []
        acc = Fr(0)
        for j in range(P+1):
            acc += (-1)**j * S[j]
            partial.append(acc)
        p0 = partial[-1]
        cap = CAP.get(len(E))
        print(f"\n  {name}: p0={float(p0):.5f}  cap={float(cap):.4f}" if cap else f"\n  {name}: p0={float(p0):.5f}")
        print(f"    S_j  = [{', '.join(f'{float(s):.4f}' for s in S)}]")
        print(f"    IE partial sums (Bonferroni bracketing):")
        print(f"      {[f'{float(p):.4f}' for p in partial]}")
        # Which partial sums are valid UPPER bounds (even truncation) and how tight?
        # Bonferroni: P_{2t} >= p0 >= P_{2t+1}. find tightest even >= cap-violating?
        evens = [partial[t] for t in range(0,P+1,2)]
        odds  = [partial[t] for t in range(1,P+1,2)]
        tightest_upper = min(evens)
        print(f"    tightest even-trunc UPPER on p0 = {float(tightest_upper):.4f}"
              f"  ({'BEATS cap' if cap and tightest_upper<=cap else 'no good'})")

    # ---- (C) collision vs corr ----
    print("\n=== (C) TOURNAMENT-SIDE: corr vs E[#collisions], and gcd structure ===")
    def collisions_and_p0(E):
        E = [int(e) for e in E if int(e)!=0]
        xs = breakpoints(E); Eco=Fr(0); p0=Fr(0)
        for a,b in zip(xs,xs[1:]):
            mid=(a+b)/2
            secs=[sector((e*mid)%1) for e in E]
            cnt={}
            for s in secs: cnt[s]=cnt.get(s,0)+1
            Eco += (b-a)*sum(v*(v-1)//2 for v in cnt.values())
            if len(set(secs))==P: p0+=(b-a)
        return p0,Eco
    def coupon(k):
        return sum((-1)**j*comb(P,j)*Fr((P-j)**k,P**k) for j in range(P+1))
    fam = {
        "[1..8]":[1,2,3,4,5,6,7,8],
        "[1..9]":[1,2,3,4,5,6,7,8,9],
        "[2,4..16]":[2,4,6,8,10,12,14,16],
        "primes[2,3,5,7,11,13,17,19]":[2,3,5,7,11,13,17,19],
        "2clust k9":[1,2,3,4,50,51,52,53,54],
        "sep k9 [1,2,3,40,41,42,500,501,502]":[1,2,3,40,41,42,500,501,502],
    }
    print(f"  {'family':38s}{'k':>3}{'p0':>9}{'coupon':>9}{'corr':>9}{'E[coll]':>9}")
    for name,E in fam.items():
        p0,Eco=collisions_and_p0(E); k=len(E); c=coupon(k)
        print(f"  {name:38s}{k:>3}{float(p0):>9.4f}{float(c):>9.4f}{float(p0-c):>+9.4f}{float(Eco):>9.4f}")
    print("  hypothesis check: does corr DECREASE as clusters separate (-> coupon)?")

    print("\nDONE.")

if __name__ == "__main__":
    main()
