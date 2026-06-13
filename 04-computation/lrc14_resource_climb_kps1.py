#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LRC(14): how HIGH can 13 runners climb the band ladder? — the resource-bound climb.
kind-pasteur-2026-06-13-S1.  (Refutes the band-2 ceiling f(13)=41: found a primitive
non-dominant config blocking all q<=41 with first witness q=43.  So: is there ANY
finite ceiling for NON-DOMINANT configs, or does the ladder climb forever unless a
runner goes B'-dominant?  This is t-0124 / HYP-2438.)

THE EXPERIMENT.  For a target ceiling K in {27,41,55,69,83,97,...}, greedily build a
primitive multiple-of-14 13-set that BLOCKS every shell q<=K (no strict witness
<=K, i.e. ladder height > K).  Measure:
  - the max K actually achievable (does it plateau? = a ceiling),
  - the DOMINANCE RATIO rho = max/second-max of the achieved blocker (does it grow
    with K? => climbing forces dominance => B' eventually fires => ladder u B' closes,
    the HYP-2438 mechanism made quantitative),
  - the ENTRY GROWTH (largest speed needed to block up to K).
Honest: greedy is a heuristic lower bound on achievable height; a NON-finding of high
climbers is weak evidence, a finding of forced-dominance-growth is the positive signal.

All exact; band criterion (THM-492): a/q witness iff all v: va mod q not in
B_q = +-{0..floor(q/14)}.
"""

import random, time
from math import gcd
from functools import reduce

BANDS = {}
def band(q, n=14):
    if q in BANDS: return BANDS[q]
    h = q // n
    B = [min(r, q - r) <= h for r in range(q)]
    BANDS[q] = B
    return B

def escapers(S, Q):
    """# escaping witness-numerators summed over q=2..Q (0 => blocks all shells<=Q)."""
    tot = 0
    for q in range(2, Q + 1):
        B = band(q)
        for a in range(1, q):
            if gcd(a, q) != 1: continue
            if all(not B[(v * a) % q] for v in S): tot += 1
    return tot

def first_witness(S, Hmax=400):
    for q in range(2, Hmax + 1):
        B = band(q)
        for a in range(1, q):
            if gcd(a, q) != 1: continue
            if all(not B[(v * a) % q] for v in S):
                return q
    return None

def primitive_mult14(S):
    return reduce(gcd, S) == 1 and any(v % 14 == 0 for v in S)

def dominance(S):
    s = sorted(S); return s[-1] / s[-2]

def is_dominant(S, n=14):
    s = sorted(S); return s[-1] > (n - 1) * s[-2]


def greedy_block_to(K, emax, restarts, seed):
    """Best (fewest residual escapers over q<=K) primitive 13-set; returns the
    config achieving escapers==0 with SMALLEST dominance ratio, else best-effort."""
    rng = random.Random(seed)
    pool = [x for x in range(1, emax + 1)]
    best_blocker = None  # (dominance, S) among full blockers
    best_effort = None   # (escapers, S)
    for _ in range(restarts):
        S = [14 * rng.randint(1, max(1, emax // 14))]
        while len(S) < 13:
            cands = rng.sample(pool, min(50, len(pool)))
            scored = []
            for c in cands:
                if c in S: continue
                scored.append((escapers(S + [c], K), c))
            if not scored: break
            scored.sort()
            # tie-break toward smaller entries (avoid forcing dominance)
            m = scored[0][0]
            choices = [c for e, c in scored if e == m]
            S.append(min(choices))
        S = sorted(set(S))
        if len(S) != 13 or not primitive_mult14(S):
            continue
        e = escapers(S, K)
        if best_effort is None or e < best_effort[0]:
            best_effort = (e, S)
        if e == 0:
            rho = dominance(S)
            if best_blocker is None or rho < best_blocker[0]:
                best_blocker = (rho, S)
    return best_blocker, best_effort


def main():
    t0 = time.time()
    print("=== RESOURCE CLIMB: how high can a primitive 13-set block the band ladder? ===", flush=True)
    print("   K = target ceiling (block all shells q<=K).  band-k ceiling = (k+1)*14-1.", flush=True)
    print("   For each K: can we block all q<=K? with what dominance ratio rho=max/2nd?", flush=True)
    print("   rho > 13 => B'-DOMINANT (THM-398 Cor B2 fires => loose anyway).", flush=True)
    print(flush=True)
    print("    K  band  blocked?  min-dominance-rho  non-dom-blocker?  maxentry  true-height", flush=True)
    results = []
    for K in (27, 41, 55, 69, 83):
        bandk = (K + 1) // 14
        # scale emax with K (need bigger entries to block higher)
        emax = max(120, K * 18)
        blk, eff = greedy_block_to(K, emax=emax, restarts=240, seed=100 + K)
        if blk is not None:
            rho, S = blk
            nondom = rho <= 13
            h = first_witness(S)
            results.append((K, S, rho, nondom, h))
            print(f"   {K:3d}  ~{bandk:2d}   YES       rho={rho:6.2f}          "
                  f"{'YES' if nondom else 'NO(B-dom)'}        {max(S):5d}     {h}", flush=True)
        else:
            e, S = eff
            print(f"   {K:3d}  ~{bandk:2d}   no(res={e})   --                  --             "
                  f"{max(S):5d}     best-effort", flush=True)
            results.append((K, S, None, None, None))

    print("\n=== verdict ===", flush=True)
    nd = [(K, rho, h, max(S)) for (K, S, rho, nd_, h) in results if rho is not None and nd_]
    if nd:
        print("   NON-DOMINANT blockers found at K =", [x[0] for x in nd], flush=True)
        print("   dominance rho vs K (does rho -> 13 as K grows? => forced dominance):", flush=True)
        for K, rho, h, me in nd:
            print(f"      K={K}: min non-dom rho achievable = {rho:.2f}, true height = {h}, max entry {me}", flush=True)
        khi = max(x[0] for x in nd)
        print(f"   HIGHEST K with a NON-DOMINANT full blocker = {khi} "
              f"(band-{(khi+1)//14}); ladder height there > {khi}.", flush=True)
        print("   => the band ladder has NO finite ceiling for non-dominant configs at least up", flush=True)
        print(f"      to band-{(khi+1)//14}; consistent with HYP-2438 (closure needs B'(any)).", flush=True)
    print(f"\n=== TOTAL {time.time()-t0:.1f}s ===", flush=True)


if __name__ == "__main__":
    main()
