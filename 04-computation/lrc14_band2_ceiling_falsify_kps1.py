#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LRC(14): the BAND-2 CEILING conjecture and an adversarial attempt to FALSIFY it.
kind-pasteur-2026-06-13-S1.

CONJECTURE (sharpened resource bound, the t-0124 finite-check target).
  Every primitive, non-B'-dominant speed set S of 13 integers containing a multiple
  of 14 has a strict loneliness witness t = a/q at some shell q <= 41 = 3n-1 (the
  BAND-2 CEILING).  [Empirically max ladder height = 41 over the companion search;
  claudebox-S7's five evaders leak at exactly q in {40,41}; codex's eight-core at
  {31,33}.]  If true, C'(14) -- hence LRC(14) -- is a FINITE residue check mod
  lcm(1..41), since the band criterion at shells <= 41 depends only on {v_i mod q}.

FALSIFICATION = find a primitive non-dominant multiple-of-14 config whose ladder
height EXCEEDS 41 (no witness at any q <= 41 => it must block ALL of band-1 and
band-2).  Such a config would push the ceiling to band-3 (4n-1 = 55) and refute
f(13)=41.  We attack hard:
  (1) GREEDY COVER-BLOCKER: build S speed-by-speed to keep Z/q covered for as many
      q<=41 as possible (each new speed chosen to maximally reduce the # of escaping
      witnesses summed over q=2..41).
  (2) STRUCTURED FAMILIES: 7*{subset} u {strangers}, scaled APs, 13|r evader-style,
      mod-27/mod-41 aligned strangers.
  (3) RANDOM with rejection: many primitive non-dominant configs, record the max
      ladder height and any that reach >41.
A height >41 found => conjecture REFUTED (report the config). None found across a
hard search => evidence the ceiling is 41 (NOT a proof; honest).

Uses the validated band criterion (THM-492): a/q witness iff all v: va mod q not in
B_q = +-{0..floor(q/14)}.  Witness-at-q iff some unit escapes all bands.
"""

import itertools, time, random
from math import gcd
from functools import reduce


def band(q, n=14):
    h = q // n
    return [min(r, q - r) <= h for r in range(q)]


BANDS = {q: band(q) for q in range(2, 60)}


def has_witness(S, q):
    B = BANDS[q]
    for a in range(1, q):
        if gcd(a, q) != 1:
            continue
        if all(not B[(v * a) % q] for v in S):
            return True
    return False


def ladder_height(S, Hmax=58):
    for q in range(2, Hmax + 1):
        if has_witness(S, q):
            return q
    return None  # no witness up to Hmax => climbs past Hmax (FALSIFIES if Hmax>=41)


def primitive_mult14(S):
    return reduce(gcd, S) == 1 and any(v % 14 == 0 for v in S)


def dominant(S, n=14):
    s = sorted(S)
    return s[-1] > (n - 1) * s[-2]


def escapers_upto(S, Q=41):
    """total # escaping witness-numerators over q=2..Q (0 = blocks all shells <=Q)."""
    tot = 0
    for q in range(2, Q + 1):
        B = BANDS[q]
        for a in range(1, q):
            if gcd(a, q) != 1:
                continue
            if all(not B[(v * a) % q] for v in S):
                tot += 1
    return tot


# --------------------------------------------- (1) greedy cover-blocker

def greedy_blocker(Q=41, pool=range(1, 600), seed=0, restarts=200):
    """Try to assemble 13 speeds (incl. a multiple of 14) that BLOCK every shell
    q<=Q (escapers_upto == 0 => ladder height > Q => FALSIFIES)."""
    rng = random.Random(seed)
    best = None
    pool = list(pool)
    for _ in range(restarts):
        # seed with a multiple of 14
        S = [14 * rng.randint(1, 6)]
        while len(S) < 13:
            # candidate batch; pick the one minimizing remaining escapers
            cands = rng.sample(pool, min(40, len(pool)))
            scored = []
            for c in cands:
                if c in S:
                    continue
                T = S + [c]
                scored.append((escapers_upto(T, Q), c))
            if not scored:
                break
            scored.sort()
            S.append(scored[0][1])
        S = sorted(set(S))
        if len(S) != 13:
            continue
        e = escapers_upto(S, Q)
        if best is None or e < best[0]:
            best = (e, S)
        if e == 0 and primitive_mult14(S) and not dominant(S):
            return ("FALSIFIED", S, e)
    return ("best", best[1], best[0]) if best else ("none", None, None)


# --------------------------------------------- (2) structured families

def structured_search(Q=41):
    hits = []
    base = [7 * k for k in range(1, 13)]
    # evader-style: 7*{1..12} u {r}, r climbing, 13|r, mod-27/41 aligned
    for r in range(28, 4000, 1):
        S = sorted(base + [r])
        if len(set(S)) != 13:
            continue
        if not primitive_mult14(S) or dominant(S):
            continue
        h = ladder_height(S)
        if h is None or h > 41:
            hits.append((r, S, h))
    # scaled APs c*{1..13} with a multiple of 14
    for c in range(1, 60):
        S = sorted(c * k for k in range(1, 14))
        if not primitive_mult14(S) or dominant(S):
            continue
        h = ladder_height(S)
        if h is None or h > 41:
            hits.append((('AP', c), S, h))
    return hits


# --------------------------------------------- (3) random rejection

def random_search(trials=20000, emax=400, seed=7):
    rng = random.Random(seed)
    maxh = 0; argmax = None; over = []
    done = 0
    for _ in range(trials):
        S = sorted(rng.sample(range(1, emax + 1), 13))
        if not primitive_mult14(S) or dominant(S):
            continue
        done += 1
        h = ladder_height(S)
        if h is None:
            over.append(S)            # climbs past 58 -- very strong falsifier
        elif h > maxh:
            maxh, argmax = h, S
    return done, maxh, argmax, over


def main():
    t0 = time.time()
    print("=== BAND-2 CEILING conjecture: f(13) = 41 = 3n-1.  Try to FALSIFY. ===", flush=True)

    print("\n--- (2) structured families (7*{1..12}u{r}, scaled APs) ---", flush=True)
    hits = structured_search()
    over = [h for h in hits if (h[2] is None or h[2] > 41)]
    print(f"   configs with ladder height > 41 (FALSIFIERS): {len(over)}", flush=True)
    for tag, S, h in over[:10]:
        print(f"      {tag}: height={h}  S={S}", flush=True)
    # report the max height among the family
    allh = [(h[2] if h[2] else 99, h[0]) for h in hits]
    if allh:
        allh.sort()
        print(f"      max structured ladder height = {allh[-1][0]} at {allh[-1][1]}", flush=True)

    print("\n--- (3) random primitive non-dominant configs (entries<=400) ---", flush=True)
    done, maxh, argmax, rover = random_search(trials=30000)
    print(f"   {done} valid configs; MAX ladder height = {maxh}", flush=True)
    print(f"   configs climbing past shell 58: {len(rover)} (any => strong falsifier)", flush=True)
    if maxh > 41:
        print(f"   *** FALSIFIER (height {maxh} > 41): {argmax} ***", flush=True)
    else:
        print(f"   max height {maxh} <= 41: consistent with the band-2 ceiling. argmax={argmax}", flush=True)

    print("\n--- (1) greedy cover-blocker (adversarial: block ALL shells q<=41) ---", flush=True)
    status, S, e = greedy_blocker(restarts=120)
    print(f"   {status}: best residual escapers over q<=41 = {e}", flush=True)
    if status == "FALSIFIED":
        print(f"   *** FALSIFIER: blocks all shells <=41: {S} ***", flush=True)
    else:
        print(f"   best blocker config (escapers {e} > 0 => still has a witness <=41): {S}", flush=True)
        if S:
            print(f"      its ladder height = {ladder_height(S)}", flush=True)

    print(f"\n=== TOTAL {time.time()-t0:.1f}s ===", flush=True)


if __name__ == "__main__":
    main()
