#!/usr/bin/env python3
"""
HYP-5862 (boxeph-S13): hB5 via the DOUBLING-PAIR-SUM ruler q = 3a.
A doubling pair (a, 2a) has pair-sum 3a; the 3a-ruler hosts the chain
anchor p = a (p/q = 1/3), where every chain element a*2^e clears at 1/3.
TEST: on covering 13-sets containing a doubling pair, is the liveness
count LM(3a) = #{p in [1,q-1] : all residues in the safe band} > 0?
If yes on the adversarial bank, hB5 splits: [doubling-rich: 3a-ruler]
+ [doubling-poor: small E3-diagonal, klein's signed box]. Exact integers.
"""
import random
from math import gcd

def is_covering(S):
    return all(any(v % q == 0 for v in S) for q in range(2, 15))

def LM(S, q):
    c = 0
    for p in range(1, q):
        if all(14 * min((v * p) % q, q - (v * p) % q) >= q for v in S):
            c += 1
    return c

def anchor_live(S, a):
    q = 3 * a
    p = a
    return all(14 * min((v * p) % q, q - (v * p) % q) >= q for v in S)

random.seed(13)
print("hB5 via the 3a-ruler: exact LM(3a) on covering sets with doubling pairs")
worst = None
tested = 0
fails = 0
for trial in range(30000):
    a = random.randint(2, 60)
    base = [a, 2 * a]
    pool = list(range(2, 200))
    random.shuffle(pool)
    S = set(base)
    for v in pool:
        if len(S) >= 13:
            break
        S.add(v)
    S = sorted(S)
    if len(S) != 13 or not is_covering(S):
        continue
    g = 0
    for v in S:
        g = gcd(g, v)
    if g != 1:
        continue
    tested += 1
    lm = LM(S, 3 * a)
    if worst is None or lm < worst[0]:
        worst = (lm, S, a)
    if lm == 0:
        fails += 1
        if fails <= 3:
            print(f"  DEAD 3a-ruler: a={a} S={S}")
print(f"tested {tested} covering+primitive sets with doubling pair")
print(f"dead-ruler count: {fails}")
if worst:
    print(f"worst LM(3a) = {worst[0]} at a={worst[2]}, S={worst[1]}")
    lm, S, a = worst
    print(f"  anchor p=a live at worst case: {anchor_live(S, a)}")
