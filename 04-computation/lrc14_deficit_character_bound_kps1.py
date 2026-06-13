#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LRC(14): the covering DEFICIT as a character sum, the (6/7)^13 main term, and the
empirical RESOURCE BOUND (the ladder-height ceiling that would finish t-0124).
kind-pasteur-2026-06-13-S1.  Builds on lrc14_dilated_band_covering_kps1.py.

REFRAME (validated in the companion script).  A speed set S (13 speeds, with a
multiple of 14) is LOOSE at shell q iff the 13 dilated bands  A_v = v^{-1} B_q
(B_q = {r: ||r/q|| <= floor(q/14)/q}, a centered interval of 2 floor(q/14)+1 ~ q/7
residues) do NOT cover Z/q.  Define the COVERING DEFICIT

    D(q,S) = #{ a in Z/q : v a mod q not in B_q for all v in S }     (full group)
    D*(q,S)= same but a ranges over units (Z/q)*  (the strict-witness count)

so D(q,S) > 0  <=>  shell-q witness exists  <=>  S loose via shell q.

WHY A WITNESS MUST EVENTUALLY EXIST (the heuristic the character sum makes rigorous):
each band has relative size beta_q = (2 floor(q/14)+1)/q -> 1/7.  If the 13 rotated
intervals were "independent/equidistributed", D(q,S) ~ q * prod_v (1 - beta_q) ~
q (6/7)^13 ~ 0.135 q  >  0.  The cover can be COMPLETE (D=0) only if the intervals
are ANTI-equidistributed -- a correlation = an additive relation among the v's.
The deviation D - q(6/7)^13 is a multilinear character sum; for prime q it is bounded
by incomplete-character (Polya-Vinogradov / Weil) estimates ~ sqrt(q) polylog, which
the main term 0.135 q beats for q > q0.  Bounding q0 uniformly over the 13 speeds
(and excluding the B'-dominant escape) is exactly the resource bound f(13) of t-0124.

THIS SCRIPT (all exact integer arithmetic):
 (A) main-term tracking: D(q,S) vs q(6/7)^13 for the hard families across bands.
 (B) the DEVIATION delta(q,S) = D(q,S) - q*prod(1-beta_q); is it O(sqrt q)?  fit.
 (C) the RESOURCE BOUND: min over a large primitive-non-dominant config search of
     D(q,S) at each q -- the first q where the per-q minimum deficit is > 0 and
     STAYS > 0 is the candidate finite-check ceiling K.  (Honest: bounded entries.)
 (D) the additive-relation signature: at a near-cover (small deficit), measure the
     pairwise-overlap excess sum_{v<v'} |A_v cap A_v'| vs its independent value
     C(13,2) beta^2 q -- the correlation that the cover spends.
"""

import itertools, time, random
from math import gcd, floor
from fractions import Fraction


def band_indicator(q, n=14):
    h = q // n
    B = [False] * q
    for r in range(q):
        if min(r, q - r) <= h:
            B[r] = True
    return B, (2 * h + 1)


def deficit(S, q, units_only=False, n=14):
    B, _ = band_indicator(q, n)
    cnt = 0
    for a in range(q):
        if units_only and gcd(a, q) != 1:
            continue
        if all(not B[(v * a) % q] for v in S):
            cnt += 1
    return cnt


def main_term(q, n=14):
    h = q // n
    beta = Fraction(2 * h + 1, q)
    return q * (1 - beta) ** 13, float(beta)


# ----------------------------------------------------------- config families

def five_evaders():
    base = [7 * k for k in range(1, 13)]
    return {r: sorted(base + [r]) for r in (611, 702, 793, 962, 1053)}


def is_primitive_mult14(S):
    from functools import reduce
    return reduce(gcd, S) == 1 and any(v % 14 == 0 for v in S)


def is_Bprime_dominant(S, n=14):
    """One runner > (n-1) * (all others) => loose by THM-398 Cor B2; exclude these."""
    s = sorted(S)
    return s[-1] > (n - 1) * s[-2]


# ----------------------------------------------------------------- the lab

def part_A():
    print("=== A. deficit D(q,S) vs main term q(6/7)^13 (~0.1348 q) — five evaders ===", flush=True)
    print("   (band-1 ceiling q=27; evaders fully cover band-1, leak at band-2 q=40/41)", flush=True)
    for r, S in five_evaders().items():
        row = []
        for q in (23, 27, 28, 40, 41, 55, 69, 83, 97):
            D = deficit(S, q)
            mt, beta = main_term(q)
            row.append(f"q={q}:D={D}(mt~{mt:.1f})")
        print(f"   r={r}: " + "  ".join(row), flush=True)


def part_B():
    print("\n=== B. deviation delta(q)=D-q*prod(1-beta) vs sqrt(q) — evader r=611 over primes ===", flush=True)
    S = five_evaders()[611]
    import math
    print("      q   D     mainterm   delta    delta/sqrt(q)", flush=True)
    for q in (29, 41, 43, 53, 67, 71, 83, 97, 113, 127, 151, 181, 211):
        D = deficit(S, q)
        mt, _ = main_term(q)
        delta = D - float(mt)
        print(f"     {q:4d} {D:4d}   {float(mt):7.2f}   {delta:7.2f}   {delta/math.sqrt(q):7.3f}", flush=True)
    print("      (bounded delta/sqrt(q) => character-sum/Polya-Vinogradov regime; "
          "main term linear beats it => D>0 forced for q>q0)", flush=True)


def part_C(Hmax=120):
    print("\n=== C. RESOURCE BOUND: per-q MINIMUM deficit over primitive non-dominant configs ===", flush=True)
    print("   (search 13-speed sets w/ a multiple of 14, entries<=80, excluding B'-dominant;", flush=True)
    print("    the q where min-deficit first becomes >0 and stays >0 = candidate finite ceiling K)", flush=True)
    rng = random.Random(20260613)
    configs = []
    tries = 0
    while len(configs) < 3000 and tries < 400000:
        tries += 1
        S = sorted(rng.sample(range(1, 81), 13))
        if is_primitive_mult14(S) and not is_Bprime_dominant(S):
            configs.append(S)
    print(f"   sampled {len(configs)} configs.", flush=True)
    # include the hard families explicitly
    for r, S in five_evaders().items():
        if not is_Bprime_dominant(S):
            configs.append(S)
    print("      q : min-deficit (over configs)  #configs-with-deficit-0-at-this-q", flush=True)
    last_zero_q = 0
    for q in range(14, Hmax + 1):
        mind = None; zeros = 0
        for S in configs:
            D = deficit(S, q)
            if D == 0:
                zeros += 1
            if mind is None or D < mind:
                mind = D
        if mind == 0:
            last_zero_q = q
        if q <= 50 or mind == 0 or q % 10 == 0:
            print(f"     {q:4d} : {mind:4d}            {zeros}", flush=True)
    print(f"   --> largest q with ANY config still fully covered (deficit 0) = {last_zero_q}", flush=True)
    print(f"       (above this, EVERY sampled config has a witness at this single q; "
          f"the union over q<=K of witnesses covers all configs at much smaller K)", flush=True)
    # the real finite-check quantity: max over configs of ladder height
    heights = []
    for S in configs:
        h = None
        for q in range(2, Hmax + 1):
            if deficit(S, q) > 0:
                h = q; break
        heights.append(h if h else Hmax + 1)
    heights.sort()
    print(f"   ladder heights: min={heights[0]} median={heights[len(heights)//2]} "
          f"max={heights[-1]} (max = empirical resource bound f(13) on this box)", flush=True)


def part_D():
    print("\n=== D. additive-relation signature at near-cover (evaders, critical shell) ===", flush=True)
    for r, S in five_evaders().items():
        q = 41 if r in (611, 702) else 40
        B, blen = band_indicator(q)
        beta = blen / q
        # pairwise overlaps
        Asets = []
        for v in S:
            Asets.append(set(a for a in range(q) if B[(v * a) % q]))
        ov = 0
        for i in range(13):
            for j in range(i + 1, 13):
                ov += len(Asets[i] & Asets[j])
        indep = (13 * 12 // 2) * beta * beta * q
        D = deficit(S, q)
        print(f"   r={r} q={q}: deficit={D}, pairwise-overlap sum={ov}, "
              f"independent~{indep:.1f}, excess={ov-indep:+.1f} "
              f"({'OVER' if ov>indep else 'UNDER'}-correlated => cover {'tighter' if ov<indep else 'looser'})", flush=True)


def main():
    t0 = time.time()
    part_A()
    part_B()
    part_C()
    part_D()
    print(f"\n=== TOTAL {time.time()-t0:.1f}s ===", flush=True)


if __name__ == "__main__":
    main()
