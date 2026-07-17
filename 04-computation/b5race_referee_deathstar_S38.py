#!/usr/bin/env python3
"""death-star-2026-07-17-S38 (HYP-7182): referee for LRCB5Race.lean.

(1) counting lemma: doubling pairs in an injective positive 13-family <= 12
    (pairs = m - chains; top-determined; min never a top);
(2) trap confinement: on chain-dense families every doubling top sits <= j+1;
(3) the scoreboard inequality: B5 >= (q-1)*2052/16807 - 78(q-1)/49 - 715(q-1)/2401
    - tail3 - tail5 with tail3/tail5 the ACTUAL odd sums (so the inequality is tight
    against the proved floors) -- verified exactly on random (v, q)."""
import random
from fractions import Fraction as Fr
from itertools import combinations
from math import gcd

def in_band(vi, q, p):
    r = (vi * p) % q
    return q <= 14 * r <= 13 * q

def NT(v, q, T):
    return sum(1 for p in range(1, q) if all(not in_band(v[i], q, p) for i in T))

def B5(v, q):
    from math import comb
    tot = 0
    for p in range(1, q):
        c = sum(1 for vi in v if not in_band(vi, q, p))
        for d in range(6):
            tot += (-1) ** d * comb(c, d)
    return tot

def referee_counting(trials=30000, seed=38):
    rnd = random.Random(seed)
    ok = True
    for _ in range(trials):
        S = rnd.sample(range(1, 3000), 13)
        pairs = sum(1 for x in S for y in S if y == 2 * x)
        if pairs > 12:
            ok = False
    print(f"counting lemma (<=12): {'PASS' if ok else 'FAIL'}")

def referee_confinement(trials=4000, seed=138):
    rnd = random.Random(seed)
    ok = True
    checked = 0
    for _ in range(trials):
        # chain-dense: dense pair at j, ratios >= 3 above
        j = rnd.randint(0, 11)
        w = sorted(rnd.sample(range(1, 60), j)) if j else []
        base = (w[-1] if w else 1) * rnd.randint(1, 3) + rnd.randint(0, 4)
        w.append(base)
        w.append(base + max(1, int(base * rnd.uniform(0.2, 1.9))))
        if w[-1] >= 3 * w[-2]:
            continue
        good = True
        for k in range(j + 1, 12):
            lo = 3 * w[-1]
            w.append(lo + rnd.randint(0, lo))
            if False:
                good = False
        if len(w) != 13 or len(set(w)) != 13:
            continue
        checked += 1
        for si in range(13):
            for ti in range(13):
                if w[ti] == 2 * w[si] and ti > j + 1:
                    ok = False
                    print(f"  FAIL confinement: j={j} t={ti} w={w}")
    print(f"trap confinement (top <= j+1): {'PASS' if ok else 'FAIL'} ({checked} families)")

def referee_scoreboard(trials=25, seed=238):
    rnd = random.Random(seed)
    ok = True
    for _ in range(trials):
        q = rnd.randint(30, 90)
        v = [rnd.randint(1, 10**5) for _ in range(13)]
        if any(gcd(x, q) != 1 for x in v):
            continue
        b5 = B5(v, q)
        E = Fr(q - 1) * Fr(2052, 16807)
        tail3 = sum(Fr(NT(v, q, T)) - Fr(q - 1, 7**3) for T in combinations(range(13), 3))
        tail5 = sum(Fr(NT(v, q, T)) - Fr(q - 1, 7**5) for T in combinations(range(13), 5))
        floor = E - 78 * Fr(q - 1, 49) - 715 * Fr(q - 1, 2401) - tail3 - tail5
        if not floor <= b5:
            ok = False
            print(f"  FAIL scoreboard: q={q} b5={b5} floor={float(floor):.3f}")
    print(f"scoreboard inequality: {'PASS' if ok else 'FAIL'}")

if __name__ == "__main__":
    referee_counting()
    referee_confinement()
    referee_scoreboard()
