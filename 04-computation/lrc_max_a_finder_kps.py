#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LRC second-point structure: find the LARGEST a such that a gcd-1 k-set S achieves
max-min M(S) = a/(a(k+1)-1), the a-th Stern-Brocot mediant above the floor 1/(k+1).
kind-pasteur-2026-06-19-S9.

g(k) = sigma_2(k) - 1/(k+1) = 1/((a_max(k+1)-1)(k+1)) = Theta(1/(a_max * k^2)).
LIVE QUESTION: does a_max(k) grow with k (gap dips below Theta(1/k^2)) or stay bounded?

Two engines:
 (1) RESIDUE CONSTRUCTION: at t*=a/q (q=a(k+1)-1), a speed v is "safe at level a"
     iff (v*a mod q) in [a, q-a].  Take the k smallest gcd-1 reps of safe residues
     and check M(S) == a/q exactly (global max-min, not just the bound at t*).
 (2) LOCAL DESCENT: start from doubled apex, greedily apply moves (double a speed,
     +-1 a speed, replace) that strictly decrease M, toward the floor.
"""
import sys, itertools, random
from fractions import Fraction
from functools import reduce
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout, 'reconfigure') else None

def frac_norm(x):
    r = x - (x.numerator // x.denominator)
    return min(r, 1 - r)

def maxmin(S):
    S = sorted(set(S)); cand = set()
    for i in range(len(S)):
        for j in range(i, len(S)):
            for d in (S[i] + S[j], S[i] - S[j]):
                if d == 0: continue
                d = abs(d)
                for m in range(1, d): cand.add(Fraction(m, d))
    best = Fraction(0); bestt = Fraction(0)
    for t in cand:
        mv = min(frac_norm(v * t) for v in S)
        if mv > best: best, bestt = mv, t
    return best, bestt

def primitive(S): return reduce(gcd, S) == 1

def construct_level_a(k, a, max_extra=200):
    """Build a candidate k-set whose speeds are smallest reps of safe residues at level a.
       Returns S (sorted) or None.  q=a(k+1)-1, safe residue band [a, q-a]."""
    q = a * (k + 1) - 1
    if q <= 2 * a: return None
    # speed v is safe iff (a*v mod q) in [a, q-a]
    speeds = []
    v = 1
    while len(speeds) < k and v <= q + max_extra:
        r = (a * v) % q
        if a <= r <= q - a:
            speeds.append(v)
        v += 1
    if len(speeds) < k: return None
    # take the k smallest safe speeds; ensure speed 1 reachable / gcd 1
    S = tuple(sorted(speeds[:k]))
    return S

def best_construction(k, a):
    """Try several safe-speed selections at level a; return (M, S) with M minimal that
       actually realizes a/q, else the min M found."""
    q = a * (k + 1) - 1
    target = Fraction(a, q)
    floor = Fraction(1, k + 1)
    if q <= 2 * a: return None
    # gather many safe speeds, then pick k of them in a few ways
    safe = [v for v in range(1, q + 200) if a <= (a * v) % q <= q - a]
    results = []
    # selection 1: k smallest
    for sel in [safe[:k]]:
        if len(sel) == k and primitive(sel):
            M, t = maxmin(tuple(sorted(sel)))
            results.append((M, tuple(sorted(sel))))
    # selection 2: smallest including both boundary-residue speeds (bind at t*)
    lo = [v for v in safe if (a * v) % q == a]
    hi = [v for v in safe if (a * v) % q == q - a]
    if lo and hi:
        forced = {lo[0], hi[0]}
        rest = [v for v in safe if v not in forced]
        sel = sorted(forced) + rest[:k - len(forced)]
        sel = sorted(set(sel))[:k]
        if len(sel) == k and primitive(sel):
            M, t = maxmin(tuple(sorted(sel)))
            results.append((M, tuple(sorted(sel))))
    if not results: return None
    results.sort()
    return results[0], target, floor

def local_descent(k, iters=4000, restarts=40):
    """Random local search for the config minimizing M(S) (the true sigma_2)."""
    floor = Fraction(1, k + 1)
    best_M = Fraction(1); best_S = None
    rng = random.Random(12345 + k)
    for _ in range(restarts):
        # seed: doubled apex, or random
        if rng.random() < 0.5:
            S = list(range(1, k)) + [2 * k]
        else:
            S = sorted(rng.sample(range(1, 4 * k + 4), k))
            if 1 not in S: S[0] = 1; S = sorted(set(S))
            while len(S) < k:
                c = rng.randint(1, 4*k+4)
                if c not in S: S.append(c)
            S = sorted(S)
        if not primitive(S): continue
        M, _ = maxmin(tuple(S))
        cur_M, cur_S = M, tuple(S)
        for _ in range(iters // restarts):
            i = rng.randrange(k)
            S2 = list(cur_S)
            move = rng.randrange(4)
            if move == 0: S2[i] = S2[i] + 1
            elif move == 1: S2[i] = max(1, S2[i] - 1)
            elif move == 2: S2[i] = S2[i] * 2
            else: S2[i] = rng.randint(1, 4 * k + 4)
            S2 = sorted(set(S2))
            if len(S2) != k or not primitive(S2): continue
            M2, _ = maxmin(tuple(S2))
            if M2 < cur_M or (M2 == cur_M and rng.random() < 0.3):
                cur_M, cur_S = M2, tuple(S2)
        if cur_M < best_M and cur_M > floor:
            best_M, best_S = cur_M, cur_S
    return best_M, best_S

def express_as_mediant(M, k):
    """If M = a/(a(k+1)-1) return a, else None. Also return excess e=p(k+1)-q."""
    p, q = M.numerator, M.denominator
    e = p * (k + 1) - q
    a = None
    if q == p * (k + 1) - 1:   # e == p  -> wait: a/(a(k+1)-1): p=a, q=a(k+1)-1 => e = a(k+1)-(a(k+1)-1)=1
        pass
    if e == 1:
        a = p
    return a, e

if __name__ == "__main__":
    print("=== construct level-a families {smallest safe reps mod q=a(k+1)-1} ===")
    print("    checking which a are realizable as exact max-min a/(a(k+1)-1)\n")
    for k in [5, 6, 7, 8, 9, 10, 11, 12, 13]:
        floor = Fraction(1, k + 1)
        achievable = []
        for a in range(2, k + 2):
            res = best_construction(k, a)
            if res is None: continue
            (Mbest, Sbest), target, fl = res
            hit = (Mbest == target)
            if hit:
                achievable.append((a, Sbest))
        if achievable:
            amax, Smax = achievable[-1]
            q = amax * (k + 1) - 1
            g = Fraction(amax, q) - floor
            print(f"  k={k:2d}: achievable a = {[a for a,_ in achievable]}  | a_max={amax}  "
                  f"sigma_2<= {amax}/{q}  g={g}  g*k^2={float(g)*k*k:.4f}  witness={Smax}")
        else:
            print(f"  k={k:2d}: no level-a construction realized exactly (search wider)")

    print("\n=== local descent: true minimal M(S) found (lower sigma_2 => larger a) ===")
    for k in [5, 6, 7, 8, 9, 10, 11, 12, 13]:
        M, S = local_descent(k)
        floor = Fraction(1, k + 1)
        a, e = express_as_mediant(M, k)
        g = M - floor
        tag = f"= {a}/(a(k+1)-1)" if a else f"excess e={e}"
        print(f"  k={k:2d}: min M={M} ({float(M):.5f})  floor={floor}  g={g}  g*k^2={float(g)*k*k:.4f}  [{tag}]  S={S}")
