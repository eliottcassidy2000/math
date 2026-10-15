#!/usr/bin/env python3
"""Audit B: THM-4627 (2)(a) periodicity check.  M_h(K) = {n >= 1 : exists k <= K with T^k(n) = T^k(n+h) and equal
odd-step counts up to k} (T iterated K times, no stopping at 1).  Claim: M_h(K) is a union of residue classes mod 2^K.
Test: compute membership for every n <= 2^20 by direct iteration and compare with the membership of the class
representative r = n mod 2^K (r in [2^K, 2^(K+1)) used as the representative to avoid r = 0).  Also report the density
mu_h(K) and the fraction of ODD classes in M_h(K), against the bound 1 - 2(1 - mu_h(K))."""
import numpy as np

def member(n, h, K):
    x, y, a, b = n, n + h, 0, 0
    for k in range(1, K + 1):
        if x & 1: x = (3 * x + 1) >> 1; a += 1
        else: x >>= 1
        if y & 1: y = (3 * y + 1) >> 1; b += 1
        else: y >>= 1
        if x == y and a == b:
            return True
    return False

NMAX = 1 << 18
for h in (1, 2, 3):
    for K in (8, 12, 14):
        M = 1 << K
        cls = np.array([member(M + r, h, K) for r in range(M)])   # representative M + r of the class r mod 2^K
        bad = sum(1 for n in range(1, NMAX + 1) if member(n, h, K) != cls[n % M])
        mu = cls.mean()
        odd = cls[1::2].mean()
        print(f"h = {h}, K = {K:2d}: mu_h(K) = {mu:.5f}; odd-class fraction {odd:.5f} >= 1 - 2(1 - mu) = {1 - 2 * (1 - mu):.5f}: "
              f"{odd >= 1 - 2 * (1 - mu) - 1e-12}; n <= 2^18 disagreeing with their class: {bad}", flush=True)
