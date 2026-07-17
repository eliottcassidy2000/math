#!/usr/bin/env python3
"""Recon for HYP-7233 (death-star-S47): the locked-chain joint count.

(a) PAIR IFF: for v2 = M*v1, M <= 13: (both fail at p) <=> exists w with
    14*M*|v1*p - w*q| < q.  Also chain version: members u*M_i, M_i <= 13:
    all fail <=> 14*M_max*|u*p - w*q| < q.
(b) EXACT COUNT at gcd(u,q)=1: #{p in (0,q) : narrow-fail} = 2*floor((q-1)/(14*M)).
(c) LOCKED PAIR DEVIATION: D = N_pair - (q-1)/49 law across M = 2..13.
Sharpness probe: does the iff break at M = 14?
"""
from math import gcd
import random
random.seed(47)

def fails(v, q, p):
    r = (v * p) % q
    return 14 * r < q or 14 * r > 13 * q

def narrow_fails(u, q, p, M):
    r = (u * p) % q
    return 14 * M * r < q or 14 * M * r > (14 * M - 1) * q

iff_checked = iff_fails = 0
for _ in range(300000):
    M = random.randint(2, 13)
    u = random.randint(1, 300)
    q = random.randint(2, 800)
    p = random.randint(1, q - 1)
    lhs = fails(u, q, p) and fails(M * u, q, p)
    rhs = narrow_fails(u, q, p, M)
    iff_checked += 1
    if lhs != rhs:
        iff_fails += 1
        if iff_fails <= 3: print("IFF FAIL", M, u, q, p, lhs, rhs)
print(f"(a) pair iff M<=13: {iff_checked} checked, fails={iff_fails}",
      "PASS" if iff_fails == 0 else "FAIL")

# chain version: u, u*M2, u*M3 with 1 < M2 < M3 <= 13
chain_checked = chain_fails = 0
for _ in range(150000):
    M3 = random.randint(3, 13); M2 = random.randint(2, M3 - 1)
    u = random.randint(1, 200); q = random.randint(2, 600); p = random.randint(1, q - 1)
    lhs = fails(u, q, p) and fails(M2 * u, q, p) and fails(M3 * u, q, p)
    rhs = narrow_fails(u, q, p, M3)
    chain_checked += 1
    if lhs != rhs:
        chain_fails += 1
        if chain_fails <= 3: print("CHAIN FAIL", M2, M3, u, q, p, lhs, rhs)
print(f"(a') chain iff (3 members): {chain_checked} checked, fails={chain_fails}",
      "PASS" if chain_fails == 0 else "FAIL")

# M = 14 sharpness of the iff
break14 = []
for q in range(2, 400):
    for u in range(1, 40):
        if len(break14) >= 3: break
        for p in range(1, q):
            lhs = fails(u, q, p) and fails(14 * u, q, p)
            if lhs != narrow_fails(u, q, p, 14):
                break14.append((u, q, p)); break
print(f"(a'') iff at M=14: breaks found {len(break14)} {break14[:2]}")

# (b) exact count at coprime moduli
cnt_checked = cnt_fails = 0
for _ in range(4000):
    M = random.randint(2, 13); u = random.randint(1, 300); q = random.randint(2, 700)
    if gcd(u, q) != 1: continue
    cnt_checked += 1
    N = sum(1 for p in range(1, q) if narrow_fails(u, q, p, M))
    pred = 2 * ((q - 1) // (14 * M))
    if N != pred:
        cnt_fails += 1
        if cnt_fails <= 3: print("COUNT FAIL", M, u, q, N, pred)
print(f"(b) exact count 2*floor((q-1)/(14M)) at gcd=1: {cnt_checked} checked, fails={cnt_fails}",
      "PASS" if cnt_fails == 0 else "FAIL")

# (c) the deviation law sign pattern
print("(c) locked-pair deviation D(M) = 2*floor((q-1)/(14M)) - (q-1)/49 at q=9800:")
q = 9800
for M in range(2, 14):
    D = 2 * ((q - 1) // (14 * M)) - (q - 1) / 49
    print(f"    M={M:2d}: D = {D:9.1f}   ({'+' if D>0 else '-' if D<0 else '0'})")
