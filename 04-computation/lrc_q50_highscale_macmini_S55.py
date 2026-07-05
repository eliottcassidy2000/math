#!/usr/bin/env python3
"""
mac-mini-2026-07-05-S55 -- HYP-4119: Q50 at UNBOUNDED heights (adversarial sampling).

The census verifies Q50 (profile => witness at q <= 50) exhaustively to height 48.
THIS: random profile-passing families at heights ~10^4..10^7, built by CRT so the
pinning holds by construction -- the heights no census will ever reach.  If Q50
holds there too, the conjecture is height-independent evidence-wise (as the
periodicity reduction predicts: the witness condition only sees residues).

Construction per sample:
  - draw 12 residues mod P = 5354228880 (= lcm(2..25) / 5 to keep ints small?
    NO: use L = lcm(2..25) directly, python bigints fine)
    ... simpler: draw the 12 elements as x_i = r_i + L * t_i with r_i a
    PROFILE-PASSING template drawn from the census-like generator:
      * start from a random census-style low family (height <= 48) that passes
        the full profile (we re-generate via rejection from random 12-subsets
        of [1,48] -- cheap since we know the pass rate ~ 5e5/2e9 ~ too low!
      => INSTEAD: lift ACTUAL census survivors: take a survivor W0 (re-derived
        by a mini-census at B=36), and lift each element x -> x + L * t_i with
        random t_i in [0, T]; residues mod every q <= 25 UNCHANGED => pinning,
        covering (q <= 12) still hold; spread improves; 24-compression and
        pair-38 rechecked; primitivity rechecked.
  - test: witness at q <= 50?  (Q50); also record first-q.
"""
from math import gcd, lcm
from functools import reduce
from itertools import combinations
import random, sys, time

T0 = time.time()
def log(m=""):
    print(m, flush=True)

L = lcm(*range(2, 26))
random.seed(5555)

def dq(x, q):
    x %= q
    return min(x, q - x)

def profile_pass(W):
    Ws = sorted(W)
    for m in range(2, 13):
        if not any(v % m == 0 for v in W):
            return False
    if 2 * Ws[-1] <= 23 * Ws[0]:
        return False
    if Ws[-1] > 24 * Ws[-2]:
        return False
    if Ws[-1] + Ws[-2] < 38:
        return False
    if reduce(gcd, W) != 1:
        return False
    for q in range(2, 26):
        if any(v % q == 0 for v in W):
            continue
        for b in range(1, q // 2 + 1):
            if gcd(b, q) != 1:
                continue
            if not any(v % q in (b, q - b) for v in W):
                return False
    return True

def first_witness_q(W, qmax=50):
    for q in range(8, qmax + 1):
        if any(v % q == 0 for v in W):
            continue
        for a in range(1, q // 2 + 1):
            if gcd(a, q) != 1:
                continue
            if all(25 * dq(a * v, q) >= 2 * q for v in W):
                return q
    return 0

# -- step 1: harvest ~500 low profile survivors by direct search in [1,36]
log("harvesting low-profile seeds in [1,36]...")
seeds = []
vals = list(range(1, 37))
tries = 0
while len(seeds) < 500 and tries < 3_000_000:
    tries += 1
    W = random.sample(vals, 12)
    if max(W) >= 25 and profile_pass(W):
        seeds.append(sorted(W))
log(f"seeds: {len(seeds)} (from {tries} tries)  [{time.time()-T0:.0f}s]")

# -- step 2: lift seeds to high scales, keep profile, test Q50
stats = dict(lifted=0, profile_ok=0, q50_ok=0, q50_fail=0)
firstq = {}
fails = []
for trial in range(20000):
    W0 = random.choice(seeds)
    # lift a random subset of elements by L * t (residues mod <=25 unchanged)
    W = []
    for x in W0:
        t = random.choice([0, 0, 1, 2, 5, 17, 100, 4321]) if random.random() < 0.7 else 0
        W.append(x + L * t)
    W = sorted(set(W))
    if len(W) < 12:
        continue
    stats['lifted'] += 1
    if not profile_pass(W):
        continue                    # 24-compression/pair38 can break; skip
    stats['profile_ok'] += 1
    q = first_witness_q(W)
    if q:
        stats['q50_ok'] += 1
        firstq[q] = firstq.get(q, 0) + 1
    else:
        stats['q50_fail'] += 1
        fails.append(list(W))
        if len(fails) <= 5:
            log(f"  !! Q50 FAILURE: {W}")
log(f"\nhigh-scale Q50 test: {stats}  [{time.time()-T0:.0f}s]")
log(f"first-q histogram: {dict(sorted(firstq.items()))}")
log(f"heights reached: up to ~{4321 * L}")
log("VERDICT: " + ("Q50 HOLDS on all high-scale profile samples -- height-independent, "
                   "as the periodicity reduction predicts (witness sees only residues)."
                   if stats['q50_fail'] == 0 else
                   f"Q50 FAILS on {stats['q50_fail']} samples -- STRUCTURAL NEWS, see lists."))
