#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LRC(14): the COMPRESSIBILITY REDUCTION — m runners sharing a common factor d
reduce the config to LRC-on-the-quotient-block + dodging the (13-m) strangers.
kind-pasteur-2026-06-13-S2.  Generalizes codex's HYP-2465 (which is the d=7,
core=7*{1..12} special case) to ARBITRARY common-factor blocks.

THE LEMMA TO TEST/STATE.
If m of the 13 runners are divisible by d (write that block as d*A, |A|=m, A a set
of m distinct positive integers; the other 13-m are "strangers"), then by dilation
M(d*A as a standalone) = M(A) >= 1/(m+1) (LRC(m+1), proven for m+1<=13). At the
witness t*/d for A, the m block-runners all exceed 1/(m+1) >= 1/14 (for m<=13). The
config is loose iff the (13-m) strangers can be simultaneously dodged in the safe
window of the block — a Criterion-B' / dominance dodge whose difficulty grows with
the number of strangers (13-m). So HIGH compressibility m (FEW strangers) => the
dodge is easy => loose with an EARLY witness; the reduction degrades as m drops.

EXPERIMENTS (exact band criterion):
 (1) THRESHOLD: for configs with exactly m runners sharing some factor d, the
     ladder height h(S) as a function of (m, #strangers=13-m). Does high m force
     low h (early witness)?  Where is the threshold?
 (2) The block-witness: for high-m configs, is the witness shell related to the
     block-quotient's own structure (the A-set), confirming the reduction?
 (3) Compare to codex HYP-2465 (d=7 core, m>=9 => Q27 witness): does the general
     d-block version also give early witnesses at m>=9?
"""

import random, time
from math import gcd
from functools import reduce
from collections import defaultdict

BANDS = {}
def band(q, n=14):
    if q in BANDS: return BANDS[q]
    h = q // 14
    B = [min(r, q - r) <= h for r in range(q)]
    BANDS[q] = B
    return B

def deficit(S, q):
    B = band(q); c = 0
    for a in range(1, q):
        if gcd(a, q) != 1: continue
        if all(not B[(v * a) % q] for v in S): c += 1
    return c

def ladder_height(S, Hmax=130):
    for q in range(2, Hmax + 1):
        if deficit(S, q) > 0: return q
    return Hmax + 1

def primitive_mult14(S):
    return reduce(gcd, S) == 1 and any(v % 14 == 0 for v in S)
def is_dominant(S, n=14):
    s = sorted(S); return s[-1] > 13 * s[-2]

def max_block(S):
    """(m, d): largest m runners sharing a common factor d>1."""
    best = (0, 1)
    for d in range(2, max(S) + 1):
        m = sum(1 for v in S if v % d == 0)
        if m > best[0]:
            best = (m, d)
    return best


def build_block_config(d, A, strangers):
    """d*A union strangers, as a sorted set."""
    return sorted(set([d * a for a in A] + list(strangers)))


def main():
    t0 = time.time()
    rng = random.Random(13)
    print("=== (1) THRESHOLD: ladder height h vs compressibility m (=runners sharing a factor) ===", flush=True)
    print("    constructed family: d*A (|A|=m, A subset {1..14}) + (13-m) random strangers, primitive, non-dom", flush=True)
    bym = defaultdict(list)
    tries = 0
    while sum(len(v) for v in bym.values()) < 2400 and tries < 300000:
        tries += 1
        m = rng.randint(4, 13)
        d = rng.choice([2, 3, 5, 6, 7, 10, 14])
        A = rng.sample(range(1, 15), m)
        block = [d * a for a in A]
        nst = 13 - m
        strangers = []
        guard = 0
        while len(strangers) < nst and guard < 200:
            guard += 1
            s = rng.randint(1, 200)
            if s not in block and s not in strangers and s % d != 0:
                strangers.append(s)
        if len(strangers) != nst:
            if nst == 0:
                strangers = []
            else:
                continue
        S = sorted(set(block + strangers))
        if len(S) != 13: continue
        if not primitive_mult14(S) or is_dominant(S): continue
        mm, dd = max_block(S)
        h = ladder_height(S)
        bym[mm].append(h)
    print("      m (max block) | n | h: min  median  max  | %with h<=27 (band-1)", flush=True)
    for m in sorted(bym):
        hs = sorted(bym[m])
        pct = 100.0 * sum(1 for x in hs if x <= 27) / len(hs)
        print(f"        {m:2d}           |{len(hs):4d}| {hs[0]:3d}   {hs[len(hs)//2]:3d}    {hs[-1]:3d}  |  {pct:4.0f}%", flush=True)

    print("\n=== (2) high-m configs: is the witness early & block-driven? (m=12, vary stranger) ===", flush=True)
    # d*{1..12} + stranger, d in {2,7}: the evader-style family generalized
    for d in (2, 3, 7):
        hs = []
        for r in range(1, 400):
            block = [d * a for a in range(1, 13)]
            if r % d == 0 or r in block: continue
            S = sorted(block + [r])
            if len(S) != 13 or not primitive_mult14(S) or is_dominant(S): continue
            hs.append(ladder_height(S))
        if hs:
            hs.sort()
            print(f"   d={d}, d*{{1..12}}+stranger: {len(hs)} valid strangers, "
                  f"h: min={hs[0]} median={hs[len(hs)//2]} max={hs[-1]} "
                  f"(all loose; max h = the family's hardest)", flush=True)

    print("\n=== (3) general-d analog of HYP-2465 (m>=9 => early witness?) ===", flush=True)
    big = [h for m in bym for h in bym[m] if m >= 9]
    small = [h for m in bym for h in bym[m] if m <= 6]
    if big and small:
        big.sort(); small.sort()
        print(f"   m>=9 ({len(big)} configs): median h={big[len(big)//2]}, max={big[-1]}, "
              f"%band-1(<=27)={100*sum(1 for x in big if x<=27)/len(big):.0f}%", flush=True)
        print(f"   m<=6 ({len(small)} configs): median h={small[len(small)//2]}, max={small[-1]}, "
              f"%band-1(<=27)={100*sum(1 for x in small if x<=27)/len(small):.0f}%", flush=True)
        print("   (codex HYP-2465: d=7 core m>=9 => Q27 witness; here the GENERAL-d m>=9 trend)", flush=True)
    print(f"\n=== TOTAL {time.time()-t0:.1f}s ===", flush=True)


if __name__ == "__main__":
    main()
