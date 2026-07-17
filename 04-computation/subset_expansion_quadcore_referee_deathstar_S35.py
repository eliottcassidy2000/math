#!/usr/bin/env python3
"""death-star-2026-07-16-S35 (HYP-7172): referee for LRCB5SubsetExpansion.lean and
LRCBlockSplitLift.lean.

SUBSET EXPANSION: B5 = sum_(k<=5) (-1)^k sum_(|T|=k) N_T with N_T = joint failure
counts; the deviation identity B5 = (q-1)*2052/16807 + signed deviation sum EXACTLY
(rational arithmetic). Referee: direct numeric check on random (v, q).

QUAD CORE: the triple-block mass fee at eps = 3/(2 w[j+3]):
  sum_{u in triple} [2d/7 + 3/(7u) + eps*(2du + 3)] < 2d - eps,  d = (13-j)/(14(j+1)B).
Referee: lifted-dichotomy exhaustiveness (singles | pair | triple-block | QuadDenseCore)
+ planted block families close + the fee's honest measure on random dense cores."""
import random
from fractions import Fraction as Fr
from itertools import combinations

def in_band(vi, q, p):
    r = (vi * p) % q
    return q <= 14 * r <= 13 * q

def B5_direct(v, q):
    tot = 0
    for p in range(1, q):
        c = sum(1 for vi in v if not in_band(vi, q, p))
        for d in range(6):
            if d <= c:
                # C(c, d)
                from math import comb
                tot += (-1) ** d * comb(c, d)
    return tot

def B5_subset(v, q):
    from math import comb
    tot = 0
    for k in range(6):
        for T in combinations(range(13), k):
            NT = sum(1 for p in range(1, q)
                     if all(not in_band(v[i], q, p) for i in T))
            tot += (-1) ** k * NT
    return tot

def deviation_identity(v, q):
    from math import comb
    dev_sum = Fr(0)
    for k in range(6):
        for T in combinations(range(13), k):
            NT = sum(1 for p in range(1, q)
                     if all(not in_band(v[i], q, p) for i in T))
            dev_sum += (-1) ** k * (Fr(NT) - Fr(q - 1, 7 ** k))
    return Fr(q - 1) * Fr(2052, 16807) + dev_sum

def referee_subset(trials=40, seed=35):
    rnd = random.Random(seed)
    ok = True
    for _ in range(trials):
        v = [rnd.randint(1, 400) for _ in range(13)]
        q = rnd.randint(20, 90)
        b_direct = B5_direct(v, q)
        b_subset = B5_subset(v, q)
        b_ident = deviation_identity(v, q)
        if b_direct != b_subset or Fr(b_direct) != b_ident:
            ok = False
            print(f"  FAIL v={v} q={q}: direct {b_direct} subset {b_subset} ident {b_ident}")
    print(f"subset-expansion referee ({trials} random (v,q)): {'PASS' if ok else 'FAIL'}"
          " (moment=subset=equilibrium+deviation, exact rationals)")

def dense_core_certificate(w):
    bad = [j for j in range(12) if w[j + 1] < 3 * w[j]]
    if not bad:
        return None
    js = max(bad)
    if not all(3 * w[k] <= w[k + 1] for k in range(js + 1, 12)):
        return None
    if not all(2 * (12 - k) * w[k + 1] < 21 * (k + 2) * w[k] for k in range(js, 12)):
        return None
    return js

def singles_split_works(w, k):
    for jj in range(k, 12):
        if w[jj + 1] < 3 * w[jj]:
            return False
    return k == 0 or 21 * (k + 1) * w[k - 1] <= 2 * (13 - k) * w[k]

def pair_split_works(w, j):
    if not w[j + 1] < 3 * w[j]:
        return False
    if j >= 1 and not 13 * (j + 1) * w[j - 1] <= (13 - j) * w[j]:
        return False
    if j + 2 <= 12:
        if not 21 * w[j + 1] <= w[j + 2]:
            return False
        for k in range(j + 2, 12):
            if w[k + 1] < 3 * w[k]:
                return False
    return True

def triple_block_works(w, j):
    """the S35 fee: j <= 9, eps = 3/(2 w[j+3]), B = w[j-1] (1 at j=0)."""
    if j > 9:
        return False
    B = 1 if j == 0 else w[j - 1]
    d = Fr(13 - j, 14 * (j + 1) * B)
    eps = Fr(3, 2 * w[j + 3])
    fee = sum(2 * d / 7 + Fr(3, 7 * u) + eps * (2 * d * u + 3)
              for u in (w[j], w[j + 1], w[j + 2]))
    if not fee < 2 * d - eps:
        return False
    for k in range(j + 3, 12):  # tail ratios (pair indices >= j+3 > j: from cert)
        if w[k + 1] < 3 * w[k]:
            return False
    return True

def referee_quadcore(trials=150000, seed=135):
    rnd = random.Random(seed)
    ok = True
    n_s = n_p = n_b = n_q = 0
    for _ in range(trials):
        style = rnd.random()
        if style < 0.3:
            w = sorted(rnd.randint(1, 10**5) for _ in range(13))
        else:
            # dense-core-leaning families
            j = rnd.randint(0, 11)
            w = sorted(rnd.randint(1, 50) for _ in range(j))
            base = (w[-1] if w else 1) * rnd.randint(1, 3) + rnd.randint(0, 5)
            w.append(base)
            w.append(base + max(1, int(base * rnd.uniform(0.1, 1.8))))
            if w[-1] >= 3 * w[-2]:
                continue
            good = True
            for k in range(j + 1, 12):
                lo, hi = 3 * w[-1], (21 * (k + 2) * w[-1]) // (2 * (12 - k))
                if lo > hi - 1:
                    good = False
                    break
                w.append(rnd.randint(lo, max(lo, hi - 1)))
            if not good or len(w) != 13 or len(set(w)) != 13:
                continue
        if any(singles_split_works(w, k) for k in range(13)):
            n_s += 1
            continue
        js = dense_core_certificate(w)
        if js is None:
            ok = False
            print(f"  FAIL no split, no cert: {w}")
            continue
        if pair_split_works(w, js):
            n_p += 1
        elif triple_block_works(w, js):
            n_b += 1
        else:
            n_q += 1  # QuadDenseCore (fee fail or j >= 10)
    tot = n_s + n_p + n_b + n_q
    print(f"quad-core referee: {'PASS' if ok else 'FAIL'} "
          f"[singles {n_s} | pair {n_p} | triple-block {n_b} | quad-core {n_q}] of {tot}")

if __name__ == "__main__":
    referee_subset()
    referee_quadcore()
