#!/usr/bin/env python3
"""
boxeph-S14: CHAIN-COARSENED liveness for hB5.
Runner-level degree-1 Bonferroni fails always (13/7 > 1 over-covering,
forcing klein's quintic B5). COARSEN to chains: C'(p) = # maximal doubling
chains with >= 1 member in the danger band at p. Then
    LM(q) >= (q-1) - Sum_p C'(p)     [degree-1, chain level]
and the over-covering ratio is ~ Sum_i (1 - mu_{L_i}) (my S11 budget),
< 1 for W0 >= 11 ALWAYS, and per-instance often far below the partition
worst case. TEST: on covering sets with doubling structure, at pair-sum
rulers q, compare runner-level degree-1 (always negative), chain-level
degree-1, and true LM(q). Exact integers.
"""
import random
from math import gcd

def is_covering(S):
    return all(any(v % q == 0 for v in S) for q in range(2, 15))

def chains(S):
    Sset = set(S)
    out = []
    for v in sorted(S):
        if v % 2 == 0 and v // 2 in Sset:
            continue
        c = [v]
        while c[-1] * 2 in Sset:
            c.append(c[-1] * 2)
        out.append(c)
    return out

def danger(v, p, q):
    r = (v * p) % q
    return 14 * min(r, q - r) < q

def analyze(S, q):
    ch = chains(S)
    lm = 0
    sumC = 0
    sumCp = 0
    for p in range(1, q):
        C = sum(1 for v in S if danger(v, p, q))
        Cp = sum(1 for c in ch if any(danger(v, p, q) for v in c))
        if C == 0:
            lm += 1
        sumC += C
        sumCp += Cp
    d1_runner = (q - 1) - sumC
    d1_chain = (q - 1) - sumCp
    return lm, d1_runner, d1_chain, len(ch)

random.seed(14)
print("chain-coarsened degree-1 liveness at pair-sum rulers (exact)")
print("W0  m  q      LM(q)  d1_runner  d1_chain")
rows = 0
pos_chain = 0
for trial in range(60000):
    if rows >= 14:
        break
    a = random.randint(3, 40)
    # build a doubling-rich covering set: two or three chains + fill
    S = set()
    for c0, ln in [(a, random.randint(3, 6)), (random.randint(3, 50), random.randint(2, 5))]:
        x = c0
        for _ in range(ln):
            S.add(x); x *= 2
    pool = list(range(2, 300)); random.shuffle(pool)
    for v in pool:
        if len(S) >= 13: break
        S.add(v)
    S = sorted(S)
    if len(S) != 13 or not is_covering(S): continue
    g = 0
    for v in S: g = gcd(g, v)
    if g != 1: continue
    ch = chains(S)
    W0 = 13 - len(ch)
    if W0 < 5: continue
    # a pair-sum ruler from the largest two speeds
    q = S[-1] + S[-2]
    lm, d1r, d1c, m = analyze(S, q)
    rows += 1
    if d1c > 0: pos_chain += 1
    print(f"{W0:2d} {m:2d} {q:5d} {lm:6d} {d1r:9d} {d1c:9d}  "
          f"{'CHAIN-D1 POSITIVE' if d1c > 0 else ''}")
print(f"\nchain-degree-1 positive on {pos_chain}/{rows} instances "
      f"(runner-degree-1 positive on 0 by over-covering)")
