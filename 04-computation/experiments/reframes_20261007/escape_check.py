#!/usr/bin/env python3
"""Finite checks for the px+1 non-coalescence theorem (p odd >= 5).
Chain (THM-4581 table with 3 -> p): sigma = e mod 2 (2-adic), beta a fair coin:
  (0,0): (k, e/2); (0,1): (k, (p e + 1 - p^k)/2); (1,0): (k+1, (p e + 1)/2); (1,1): (k-1, (e - p^(k-1))/2).
(A) Growth lemma: from (0, e), e integer != 0, the path 'xp/2 at every run step, at k = 0 depart with beta = 1
    (to k = -1), return toward 0 at the first flip' reaches (0, e*) with |e*| >= (p/4)|e| - 1, never visiting (0,0).
    Checked exactly for all 1 <= |e| <= 2000 and p in {5,7,9,11,13,15,17,19,21,23,25,27,29,31}.
(B) Small starts: for p = 5 (threshold |e| > 4) and p = 7 (threshold |e| > 4/3), every integer start 1 <= |e| <= 4
    reaches some (0, e') with |e'| above the threshold along a coin sequence avoiding (0,0) (BFS, depth <= 24).
(C) Departure multiplier: at k = 0 with e odd both coins give |f'| in [|e|/2 - 1/(2p), |e|/2 + 1/(2p)]; every other
    step multiplies |f| by 1/2 or p/2 (one coin value each) up to an additive error <= 1/2 + 1/(2p). Checked on random states."""
from fractions import Fraction as Fr
import random
def par(e): return (e.numerator * pow(e.denominator, -1, 2)) % 2
def step(p, k, e, b):
    s = par(e)
    if s == 0:
        return (k, e / 2) if b == 0 else (k, (p * e + 1 - Fr(p) ** k) / 2)
    return (k + 1, (p * e + 1) / 2) if b == 0 else (k - 1, (e - Fr(p) ** (k - 1)) / 2)
def f_of(p, k, e): return e / Fr(p) ** max(k, 0)
def growth_path(p, e):
    k, e = 0, Fr(e); steps = 0
    # runs at k = 0 with the xp/2 outcome (beta = 1: both odd -> p e/2) until e odd
    while par(e) == 0:
        k, e = step(p, k, e, 1); steps += 1
        assert k == 0 and e != 0
    k, e = step(p, k, e, 1); steps += 1       # departure to k = -1
    assert k == -1
    while True:
        if par(e) == 0:
            k, e = step(p, k, e, 1); steps += 1   # run at k = -1, both odd: xp/2
            assert k == -1
        else:
            k, e = step(p, k, e, 0); steps += 1   # flip toward 0 (u odd): xp/2
            assert k == 0
            return e, steps
        if steps > 10000: return None, steps
okA = True
for p in range(5, 33, 2):
    for e0 in list(range(1, 2001)) + list(range(-2000, 0)):
        es, st = growth_path(p, e0)
        if es is None or es == 0 or abs(es) < Fr(p, 4) * abs(e0) - 1:
            okA = False; print("FAIL A", p, e0, es)
print(f"(A) growth lemma |e*| >= (p/4)|e| - 1 with e* != 0, all odd p in [5,31], 1 <= |e| <= 2000: {'PASS' if okA else 'FAIL'}")
def bfs(p, e0, thresh, depth):
    frontier = {(0, Fr(e0)): ()}
    seen = set(frontier)
    for dep in range(depth):
        new = {}
        for (k, e), path in frontier.items():
            for b in (0, 1):
                k2, e2 = step(p, k, e, b)
                if k2 == 0 and e2 == 0: continue
                if k2 == 0 and abs(e2) > thresh: return path + (b,), e2
                if (k2, e2) not in seen and abs(f_of(p, k2, e2)) < 10**6:
                    seen.add((k2, e2)); new[(k2, e2)] = path + (b,)
        frontier = new
    return None, None
okB = True
for p, thresh in ((5, 4), (7, Fr(4, 3))):
    for e0 in (1, 2, 3, 4, -1, -2, -3, -4):
        path, e2 = bfs(p, e0, thresh, 24)
        print(f"(B) p={p} start (0,{e0:2d}): {'reaches (0,%s) via coins %s' % (e2, ''.join(map(str, path))) if path else 'NO PATH FOUND'}")
        okB &= path is not None
print(f"(B) small starts: {'PASS' if okB else 'FAIL'}")
rnd = random.Random(3); okC = True
for _ in range(20000):
    p = rnd.choice(range(5, 33, 2)); k = rnd.randint(-6, 6)
    e = Fr(rnd.randint(-10**6, 10**6), p ** max(0, -k))
    f = abs(f_of(p, k, e))
    outs = [abs(f_of(p, *step(p, k, e, b))) for b in (0, 1)]
    if k == 0 and par(e) == 1:
        good = all(abs(o - f / 2) <= Fr(1, 2 * p) for o in outs)
    else:
        lo, hi = sorted(outs)
        good = abs(lo - f / 2) <= Fr(1, 2) + Fr(1, 2 * p) and abs(hi - Fr(p, 2) * f) <= Fr(1, 2) + Fr(1, 2 * p)
    okC &= good
print(f"(C) step multipliers (20000 random states, |k| <= 6): {'PASS' if okC else 'FAIL'}")
