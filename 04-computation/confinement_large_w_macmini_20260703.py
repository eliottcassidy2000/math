#!/usr/bin/env python3
"""
CONFINEMENT stress test: large-speed primitive tight q*=28 families? (mac-mini-2026-07-03-S33)
The anti-correlation obstruction (S32) bounds only min(w1,w2); one tightener could be LARGE. Prior searches
capped at speed ~120. Test HARD: fix even part E=2U (U binds at t=1/14, loose globally), brute-force odd
tighteners w1,w2 up to LARGE speed, check tightness (M=1/14 exact). A hit => confinement FALSE (big).
Also |F|=3 (10 even + 3 odd). numpy float pre-filter + exact confirm.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations
import numpy as np

def gcd_all(xs): return reduce(gcd, xs)
_G = np.arange(1, 16800)/16800.0   # resolves 1/16800 << 1/28, fine enough for large w up to ~800
def approxM(sp):
    v = np.array(sp, dtype=np.float64); ph = np.outer(v, _G) % 1.0
    return np.minimum(ph, 1.0-ph).min(axis=0).max()
def nd(x):
    x = x % 1
    return x if x <= 1-x else 1-x
def M_exact(sp):
    vs = sorted(set(sp)); Q = set()
    for i in range(len(vs)):
        Q.add(2*vs[i])
        for j in range(i+1, len(vs)): Q.add(vs[i]+vs[j]); Q.add(vs[j]-vs[i])
    best = F(0)
    for q in Q:
        if q < 2: continue
        for a in range(1, q):
            m = min(nd(v*F(a,q)) for v in sp)
            if m > best: best = m
    return best
target = F(1,14); tf = 1.0/14
def is_tight(S):
    if len(set(S))!=13 or gcd_all(S)!=1: return False
    if abs(approxM(S)-tf) > 1.5e-4: return False
    return M_exact(S)==target

if __name__ == "__main__":
    import sys
    def out(*a): print(*a); sys.stdout.flush()
    # even parts E=2U with U an 11-set binding at 1/14 (contains a runner ≡ ±1 mod 14) and loose
    Us = [list(range(1,12)), list(range(1,11))+[13], list(range(1,10))+[12,13],
          [1,2,3,4,5,6,7,8,9,10,25], [1,3,5,7,9,11,13,2,4,6,8]]
    hits = []
    WMAX = 799  # large odd tighteners
    out(f"(|F|=2) E=2U (11 even) + w1,w2 odd up to {WMAX}. Checking tightness...")
    for U in Us:
        U = sorted(set(U))[:11]
        if len(U)!=11: continue
        E = [2*u for u in U]
        odds = [w for w in range(1, WMAX+1, 2) if w not in E]
        cnt = 0; found_U = 0
        # float-prefilter: precompute per-w min-dist profile is heavy; just brute pairs with early approxM
        for i in range(len(odds)):
            w1 = odds[i]
            for j in range(i+1, len(odds)):
                w2 = odds[j]
                S = E + [w1, w2]
                if gcd_all(S) != 1: continue
                cnt += 1
                if is_tight(S):
                    hits.append(sorted(S)); found_U += 1
                    out(f"   TIGHT q*=28 primitive! {sorted(S)}  (E=2*{U})")
        out(f"   U={U}: tested {cnt} (w1,w2) pairs up to {WMAX}; tight found: {found_U}")
    out(f"\n(|F|=3) E=2U (10 even) + 3 odd tighteners up to 200 (sampled)")
    import random
    rng = random.Random(33)
    for U in [list(range(1,11)), list(range(1,10))+[13]]:
        E = [2*u for u in U]; cnt=0; f=0
        for _ in range(200000):
            w = sorted(rng.sample([x for x in range(1,201,2) if x not in E], 3))
            S = E + w
            if gcd_all(S)!=1: continue
            cnt+=1
            if is_tight(S): hits.append(sorted(S)); f+=1; out(f"   TIGHT |F|=3! {sorted(S)}")
        out(f"   U={U}: tested {cnt}; tight found: {f}")
    out("\n" + "="*72)
    out(f"TOTAL primitive tight q*=28 families found (large-w): {len(hits)}")
    if hits:
        out("  !!! CONFINEMENT FALSE — primitive tight q*=28 exists (override THM-612 confinement conjecture):")
        for S in hits[:10]: out(f"     {S}")
    else:
        out(f"  => NONE up to w={WMAX} (|F|=2) and sampled |F|=3. Confinement robust to large tighteners.")
        out("     The anti-correlation obstruction holds empirically even at large speed => not a small-speed artifact.")
