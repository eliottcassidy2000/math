#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LRC(14): is the BLOCKING HEIGHT (the resource) governed by the MULTIPLICATIVE
RESONANCE of the speed set?  The structure-vs-witness dichotomy.
kind-pasteur-2026-06-13-S2.  Builds on THM-497 (counting can't prove C'(14); hard
configs OVER-CORRELATE their dilated bands).

THE THEORY.  Band criterion (THM-492): a/q is a strict witness iff every v has
va mod q not in B_q = +-{0..floor(q/14)}.  Covering deficit D(q,S) = #escaping
units; D>0 <=> shell-q witness.  THM-497: the cover-completeness (D=0) is forced
NOT by cardinality but by ADDITIVE ALIGNMENT of the dilated bands {v^{-1}B_q}.

THE NEW CLAIM TO TEST.  Two bands A_v, A_v' over-overlap exactly when the RATIO
v/v' mod q is RESONANT (close to a small-denominator rational), since
|A_v cap A_v'| counts a with both va, v'a in [-k,k] -- large when v'a is pinned by
va, i.e. when v'/v is near m/m' small.  So:
  blocking height h(S)  ~  how RESONANT the speed set's ratio-multiset is.
Dichotomy: GENERIC speed sets (few resonant ratios) have D>0 at low q (early
witness); only STRUCTURED sets (many resonant ratios, e.g. c*{1..12}) climb.
If clean, this reduces C'(14) to the high-resonance/structured class -- exactly
where the elementary witnesses (dominance dodge, divisor clocks, low-q clocks) live.

EXPERIMENTS (exact arithmetic; validated band criterion):
 (1) LINK: across configs, correlate h(S) with structural measures of the speed
     set: (a) over-correlation excess at the top covering shell; (b) resonance
     count R(S,q) = #{i<j : min over small m of ||(v_i - m v_j) mod q|| is small};
     (c) the common-factor compressibility (evaders = 7*small).
 (2) DICHOTOMY: bin configs by a q-free genericity score; show low-genericity =>
     low h; only high-structure climbs.
 (3) The evaders / known hard configs vs random: are they resonance-extreme?
"""

import random, time
from math import gcd
from functools import reduce
from itertools import combinations

BANDS = {}
def band(q, n=14):
    if q in BANDS: return BANDS[q]
    h = q // n
    B = [min(r, q - r) <= h for r in range(q)]
    BANDS[q] = B
    return B

def deficit(S, q):
    B = band(q); c = 0
    for a in range(1, q):
        if gcd(a, q) != 1: continue
        if all(not B[(v * a) % q] for v in S): c += 1
    return c

def ladder_height(S, Hmax=120):
    for q in range(2, Hmax + 1):
        if deficit(S, q) > 0: return q
    return Hmax + 1

def primitive_mult14(S):
    return reduce(gcd, S) == 1 and any(v % 14 == 0 for v in S)
def is_dominant(S, n=14):
    s = sorted(S); return s[-1] > (n - 1) * s[-2]

# ---- structural measures of the speed set ----

def overcorrelation(S, q):
    """sum over pairs of (overlap - independent-expectation) at shell q."""
    B = band(q); k = q // 14
    indep = (2 * k + 1) ** 2 / q
    Asets = []
    for v in S:
        Asets.append(set(a for a in range(q) if B[(v * a) % q]))
    exc = 0.0
    for i, j in combinations(range(len(S)), 2):
        exc += len(Asets[i] & Asets[j]) - indep
    return exc

def resonance_count(S, q, mmax=6):
    """# pairs (i,j) whose ratio is near a small rational: exists small m,m'<=mmax
    with m*v_i ≡ +- m'*v_j mod q (a small-denominator resonance)."""
    cnt = 0
    for i, j in combinations(range(len(S)), 2):
        vi, vj = S[i] % q, S[j] % q
        hit = False
        for m in range(1, mmax + 1):
            for mp in range(1, mmax + 1):
                if (m * vi - mp * vj) % q == 0 or (m * vi + mp * vj) % q == 0:
                    hit = True; break
            if hit: break
        if hit: cnt += 1
    return cnt

def compressibility(S):
    """largest subset sharing a common factor d>1, scaled down = small integers.
    measure = (max over d>1 of #{v in S: d|v}) ; evaders 7*{1..12} -> 12."""
    best = 0
    for d in range(2, max(S) + 1):
        c = sum(1 for v in S if v % d == 0)
        if c > best: best = c
    return best

def additive_energy_mod(S, q):
    """E = #{(a,b,c,d) in S^4 : a+b ≡ c+d mod q}/normalizer; report raw count
    of additive quadruples of the RESIDUES (structure detector)."""
    res = [v % q for v in S]
    from collections import Counter
    sums = Counter((res[i] + res[j]) % q for i in range(len(res)) for j in range(len(res)))
    return sum(c * c for c in sums.values())


def main():
    t0 = time.time()
    rng = random.Random(20260613)
    print("=== (1)+(2) LINK + DICHOTOMY: blocking height vs structure ===", flush=True)
    rows = []
    tries = 0
    while len(rows) < 1500 and tries < 200000:
        tries += 1
        S = sorted(rng.sample(range(1, 100), 13))
        if not primitive_mult14(S) or is_dominant(S): continue
        h = ladder_height(S)
        comp = compressibility(S)
        rows.append((h, comp, S))
    rows.sort()
    # bin by compressibility
    from collections import defaultdict
    byc = defaultdict(list)
    for h, comp, S in rows:
        byc[comp].append(h)
    print(f"   {len(rows)} primitive non-dominant configs (entries<100).", flush=True)
    print("   compressibility c (max #runners sharing a factor) -> ladder-height stats:", flush=True)
    print("      c | n | h: min  median  max", flush=True)
    for c in sorted(byc):
        hs = sorted(byc[c])
        print(f"     {c:2d} | {len(hs):4d} | {hs[0]:3d}   {hs[len(hs)//2]:3d}    {hs[-1]:3d}", flush=True)
    # the high climbers: are they high-compressibility?
    rows_byh = sorted(rows, reverse=True)
    print("\n   TOP-10 climbers (highest h): their compressibility:", flush=True)
    for h, comp, S in rows_byh[:10]:
        print(f"      h={h:3d}  comp={comp:2d}  gcd-reduced-spread={max(S)//reduce(gcd,[v for v in S if comp and v%_dominant_factor(S)==0] or S)}  S={S}", flush=True)

    print("\n=== (3) known hard configs (evaders) vs random — resonance-extreme? ===", flush=True)
    evaders = {r: sorted([7*k for k in range(1,13)] + [r]) for r in (611,702,793,962,1053)}
    for r, S in evaders.items():
        h = ladder_height(S)
        comp = compressibility(S)
        q = h if h <= 120 else 41
        oc = overcorrelation(S, min(q, 50))
        print(f"   evader r={r}: h={h} comp={comp} (=12 means 12 runners share factor 7) "
              f"overcorr@q={min(q,50)}={oc:.1f}", flush=True)
    # random comparison
    samp = [S for _,_,S in rows[:200]]
    avgcomp = sum(compressibility(S) for S in samp)/len(samp)
    print(f"   random configs avg compressibility = {avgcomp:.2f} (evaders have 12)", flush=True)
    print(f"\n=== TOTAL {time.time()-t0:.1f}s ===", flush=True)

def _dominant_factor(S):
    best=1; bestc=0
    for d in range(2, max(S)+1):
        c=sum(1 for v in S if v%d==0)
        if c>bestc: bestc=c; best=d
    return best

if __name__ == "__main__":
    main()
