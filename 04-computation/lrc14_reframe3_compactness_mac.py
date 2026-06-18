#!/usr/bin/env python3
"""
REFRAME 3 (mac-mini): compactness / finite reduction for LRC(14), gap side.

Goal: reduce "M(S) >= 1/14 for ALL covering 13-sets" (infinite family) to a
FINITE check, by bounding the essential speeds of an M-minimizer.

Core mechanism (DECOUPLING LEMMA, established empirically + provable easy half):
  Fix a 12-subset core C of S. Let m = M(C) and let
      P(C) = { t in [0,1) : g(C,t) = m }   (the plateau of optimal times for C).
  Adding a 13th speed w gives M(C u {w}) = max_t min(g(C,t), ||w t||) <= M(C).
  EASY HALF (provable, 0 counterexamples): if SOME t in P(C) has ||w t|| >= m,
  then M(C u {w}) = m (w does not lower M at all).
  Whether ||w t|| >= m for a plateau point t = p/Q depends ONLY on (w mod Q),
  where Q = lcm of plateau denominators. So the *large-w* behaviour of M(C u {w})
  is periodic in w mod Q -- a huge speed cannot beat the smallest representative
  of its residue class. THIS is what bounds the minimizer's speeds.

This script:
  (1) verifies the easy-half implication on random cores;
  (2) measures, for the canonical covering core {1..11,13}, that M(C u {84m})
      is INCREASING in m (so the minimizer wants the SMALLEST covering big elt);
  (3) the main reframe-3 experiment: search covering 13-sets, sorted by M,
      and check whether the inf is attained at BOUNDED max-speed (and stabilizes
      at 7/89). Report the explicit speed bound B observed and whether the
      finite check is feasible.
"""
from fractions import Fraction as F
import random
from math import gcd
from functools import reduce

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r

def g(S, t):
    return min(nrm(v * t) for v in S)

def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2 * k + 1, 2 * v) <= F(1, 2):
            C.add(F(2 * k + 1, 2 * v)); k += 1
    for i in range(len(S)):
        for j in range(i + 1, len(S)):
            for d in (S[i] + S[j], S[j] - S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2):
                        C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C

def M(S):
    b = F(0)
    for t in cand(S):
        v = g(S, t)
        if v > b: b = v
    return b

def plateau(C):
    m = M(C)
    return [t for t in cand(C) if g(C, t) >= m], m

def is_covering(S):
    for q in range(2, 15):
        if not any(v % q == 0 for v in S):
            return False
    return True

def primitive(S):
    return reduce(gcd, S) == 1


def part1_verify_easy_half():
    random.seed(7)
    fails = 0; tested = 0
    for _ in range(120):
        C = sorted(random.sample(range(1, 40), 12))
        pl, m = plateau(C)
        for w in random.sample(range(1, 300), 12):
            if w in C:
                continue
            tested += 1
            survives = any(g(C, t) >= m and nrm(w * t) >= m for t in pl)
            if survives and M(C + [w]) < m:
                fails += 1
    print(f"[1] easy-half implication: tested={tested}, failures={fails}")


def part2_canonical_family():
    print("[2] canonical covering family {1..11,13,84m}:")
    for m in range(1, 8):
        S = [1,2,3,4,5,6,7,8,9,10,11,13, 84 * m]
        print(f"    m={m}: bigelt={84*m} M={M(S)}={float(M(S)):.5f} covering={is_covering(S)}")


def part3_search(ntrials=20000, wmax=2000):
    random.seed(3)
    best = F(1); results = []
    seen = set()
    for _ in range(ntrials):
        base = sorted(random.sample(range(1, 30), 12))
        w = random.randint(1, wmax)
        S = tuple(sorted(set(base + [w])))
        if len(S) < 13 or S in seen:
            continue
        seen.add(S)
        if is_covering(S) and primitive(S):
            m = M(list(S))
            results.append((m, max(S), S))
            if m < best: best = m
    results.sort()
    print(f"[3] random search: {len(results)} primitive covering 13-sets")
    for m, mx, S in results[:10]:
        print(f"    M={m}={float(m):.5f} maxspeed={mx} S={list(S)}")
    print(f"    global-min M = {best} = {float(best):.6f} (7/89={float(F(7,89)):.6f})")
    return results


if __name__ == "__main__":
    part1_verify_easy_half()
    part2_canonical_family()
    part3_search()
