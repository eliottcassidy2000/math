#!/usr/bin/env python3
"""FINITE-EXACT companion to THM-4581 (6c): for every residue r mod 2^K, (i) the pair chain from (0,1) driven by
v = r (absorbed by step K?) and (ii) direct integer orbits of n = 2^(K+40) + r and n + 1 (T^t(n) == T^t(n+1) for some
t <= K?).  Prints the exact unmerged fraction q(K) = #unmerged / 2^K for both, and checks they agree residue by residue."""
from coal_check import step
def T(x): return x >> 1 if x % 2 == 0 else (3 * x + 1) >> 1
for K in range(2, 21, 2):
    un_chain = un_int = disagree = 0
    big = 1 << (K + 40)
    for r in range(1 << K):
        # chain: parities of T^t(r) for t < K depend only on r mod 2^K
        k, N, v = 0, 1, r; ab = False
        for t in range(K):
            k, N = step(k, N, v & 1); v = T(v)
            if k == 0 and N == 0: ab = True; break
        n = big + r; a, b = n, n + 1; mi = False
        for t in range(K):
            a, b = T(a), T(b)
            if a == b: mi = True; break
        un_chain += (not ab); un_int += (not mi); disagree += (ab != mi)
    print(f"K={K:2d}: unmerged fraction chain = {un_chain}/{1<<K} = {un_chain/(1<<K):.5f}; integers = {un_int/(1<<K):.5f}; residue disagreements = {disagree}", flush=True)
