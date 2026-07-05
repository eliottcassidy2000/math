#!/usr/bin/env python3
"""
mac-mini-2026-07-05-S55 -- HYP-4119 item 3: the l >= 7 BOUNDED STRATUM sweep.

opus-S81's gap descent (MISTAKE-105 fix) dodges SPREAD tops (consecutive ratio
>= 26/11) with bottom entry >= ~134, for ANY top count.  The complementary
bounded piece: l >= 7 lifts with ALL lifted values in [14, 133].

STRUCTURE (from opus-S80's kernel lemmas): l >= 7 lifted coords meet the
unique-multiple positions {7..12}; a lifted r in {7..12} needs r | k_r
(height forcing), so its value is 14r >= 98; for r in {10,11,12}: 14r >= 140
> 133 -- IMPOSSIBLE in the window.  Hence bounded-stratum patterns have
C n {10,11,12} = EMPTY, so C is a subset of {1..9} with |C| >= 7:
C(9,7)+C(9,8)+C(9,9) = 36+9+1 = 46 patterns; forced coords {7,8,9} n C at
k_r = r (values 98, 112, 126); free coords C n {1..6} at k <= 10 (value <= 133).

Sweep: sieve (F1) + primitivity + witness scan at 2/25 + exact-M escalation.
Zero sub-2/25 expected => the bounded stratum is CLEAN and opus's descent
covers the rest (their ratio-cluster leg for high clustered tops).
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations, product
import sys, time

sys.path.insert(0, '04-computation')
from lonely_profile import profile

T0 = time.time()
def log(m=""):
    print(m, flush=True)

BETA = F(2, 25)
AP = list(range(1, 13))

def M_exact(S):
    for cap in (11, 8, 6, 4, 3, 2):
        p = profile(sorted(S), F(1, cap))
        m = p.M()
        if m is not None:
            return m
    return None

def dq(x, q):
    x %= q
    return min(x, q - x)

def sieve_ok(W):
    return all(any(v % m == 0 for v in W) for m in range(2, 13))

QLIB = [q for q in range(26, 61)] + [13 * u for u in range(2, 15)] + [25, 50] + [q for q in range(8, 26)]

def scan25(W):
    for q in QLIB:
        if any(v % q == 0 for v in W):
            continue
        for a in range(1, q // 2 + 1):
            if gcd(a, q) != 1:
                continue
            if all(25 * dq(a * v, q) >= 2 * q for v in W):
                return q, a
    return None

stats = dict(total=0, sieved=0, nonprim=0, scanned=0, exact=0, sub=0)
worst = None
patterns = []
for l in (7, 8, 9):
    for C in combinations(range(1, 10), l):
        patterns.append(C)
log(f"l>=7 bounded stratum: {len(patterns)} patterns (C subset {{1..9}}, |C| in 7..9)")
for C in patterns:
    forced = [r for r in C if r >= 7]          # k_r = r, value 14r
    free = [c for c in C if c <= 6]            # k in 1..10, value c+13k <= 133
    base = [v for v in AP if v not in C]
    kranges = [range(1, min(10, (133 - c) // 13) + 1) for c in free]
    for kv in product(*kranges):
        W = sorted(base + [14 * r for r in forced] + [c + 13 * k for c, k in zip(free, kv)])
        stats['total'] += 1
        if len(set(W)) < 12:
            continue
        if not sieve_ok(W):
            stats['sieved'] += 1
            continue
        if reduce(gcd, W) != 1:
            stats['nonprim'] += 1
            continue
        hit = scan25(W)
        if hit:
            stats['scanned'] += 1
            continue
        stats['exact'] += 1
        M = M_exact(W)
        if worst is None or M < worst[0]:
            worst = (M, tuple(W))
        if M < BETA:
            stats['sub'] += 1
            log(f"  << SUB-2/25: M={M} W={list(W)}")
log(f"\nl>=7 bounded stratum sweep: {stats}  [{time.time()-T0:.0f}s]")
if worst:
    log(f"worst exact case: M = {worst[0]} at {list(worst[1])}")
log("VERDICT: " + ("CLEAN -- every l>=7 lift with all lifted values <= 133 has margin >= 2/25;"
                   " with opus-S81's descent (spread tops dodged at any count, any scale),"
                   " the l>=7 stratum's remaining open piece is exactly their high ratio-cluster leg."
                   if stats['sub'] == 0 else "SUB-2/25 FOUND -- see above"))
