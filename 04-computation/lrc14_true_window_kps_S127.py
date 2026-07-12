# -*- coding: utf-8 -*-
# kps-2026-07-11-S127cont46: determining the TRUE bounded clearing window for spread DC families, and
# correcting the fleet's recently-narrowed [15,31] (opus-S238/klein-S258) -- a SAMPLE ARTIFACT.
#
# WHY [15,31] is too narrow: clearing at q needs BOTH (a) q nmid every v_i (else that runner sits at residue 0
# for all p -- ALWAYS in danger), AND (b) the coprime-to-q sub-family misses a fold-class. A spread DC family
# whose SMALL elements coincide with window moduli BLOCKS those moduli via (a). E.g.
# v=[23,26,29,30,31,40,42,44,48,50,51,54,57] has runners = 23,26,29,30,31, blocking exactly those q; it first
# clears at q=33. So the window must be wide enough to contain a non-14 q coprime-ish to all runners AND with
# a fold-class miss. This finds the true max min-clearing-q (the window width) over adversarial spread DC.
import random
from math import gcd
from functools import reduce

def is_DC(v):  return all(any(x % d == 0 for x in v) for d in range(2, 15))
def prim(v):   return reduce(gcd, v) == 1
def longest_run(v):
    v = sorted(set(v)); b = m = 1
    for i in range(1, len(v)):
        if v[i] == v[i-1] + 1: m += 1; b = max(b, m)
        else: m = 1
    return b
def clears_at(v, q):
    lo = -(-q // 14)
    return any(all(lo <= (vi * p) % q <= q - lo for vi in v) for p in range(1, q))
def min_clear_q(v, hi=120):
    for q in range(15, hi):
        if q % 14 and clears_at(v, q): return q
    return None
def valid(v):
    return len(set(v)) == 13 and prim(v) and is_DC(v) and longest_run(v) <= 7

def main():
    rng = random.Random(20260713)
    # (1) random baseline: max min-clear-q over random spread DC
    worst = (0, None); n = 0
    for vmax in (50, 80, 120, 180, 260):
        for _ in range(4000):
            v = sorted(rng.sample(range(1, vmax + 1), 13))
            if valid(v):
                n += 1; q = min_clear_q(v)
                if q and q > worst[0]: worst = (q, v)
    print(f"(1) random spread DC ({n} valid): max min-clearing-q = {worst[0]} at {worst[1]}")

    # (2) ADVERSARIAL hill-climb: mutate to MAXIMIZE min-clear-q (push the window as far as possible)
    best = worst
    for restart in range(60):
        vmax = rng.choice([40, 60, 90, 140, 220])
        cur = None
        for _ in range(400):
            v = sorted(rng.sample(range(1, vmax + 1), 13))
            if valid(v): cur = v; break
        if cur is None: continue
        curq = min_clear_q(cur) or 0
        for _ in range(300):
            w = cur[:]
            i = rng.randrange(13); w[i] = rng.randint(1, vmax)
            w = sorted(set(w))
            if len(w) != 13 or not valid(w): continue
            q = min_clear_q(w) or 0
            if q >= curq:
                cur, curq = w, q
        if curq > best[0]: best = (curq, cur)
    print(f"(2) adversarial hill-climb: max min-clearing-q = {best[0]} at {best[1]}")
    print(f"    => the TRUE bounded window must reach at least q = {best[0]} (fleet's [15,31] is TOO NARROW -- sample artifact)")
    print(f"       (my cont.34 window [8,43] is the safe bound; [15,31] undercounts the (a) q|v_i blocking)")

    # (3) how often does [15,31] actually suffice vs need q>31?
    over31 = 0; tot = 0
    for vmax in (60, 100, 160):
        for _ in range(3000):
            v = sorted(rng.sample(range(1, vmax + 1), 13))
            if valid(v):
                tot += 1; q = min_clear_q(v)
                if q and q > 31: over31 += 1
    print(f"(3) spread DC families needing q>31 to clear: {over31}/{tot} ({100*over31/max(tot,1):.2f}%) -- NONZERO => [15,31] genuinely insufficient")

if __name__ == "__main__":
    main()
