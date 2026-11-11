#!/usr/bin/env python3
# Auditor "hypotheses" (S92 audit) - own core for the F1 census. Written before reading the session scripts.
# psi(0)=0, psi(2h)=psi(h)+h, psi(2h+1)=psi(h)+floor(h/2)+1; x=2psi+1, y=2h+1.
# F1: N = x1 x2 + 2^s y1 y2 - [s=0] over ORDERED triples (h1,h2,s).
import numpy as np
from numba import njit, prange


@njit(cache=True)
def build_psi(n):
    psi = np.zeros(n, dtype=np.int64)
    for h in range(1, n):
        k = h >> 1
        if h & 1:
            psi[h] = psi[k] + (k >> 1) + 1
        else:
            psi[h] = psi[k] + k
    return psi


@njit(cache=True)
def gcd(a, b):
    while b:
        a, b = b, a % b
    return a


def states(n):
    psi = build_psi(n)
    h = np.arange(n, dtype=np.int64)
    return 2 * psi + 1, 2 * h + 1, psi


def hmax_for(X):
    # N >= x2 + y2 - 1 >= (h2+1) + (2h2+1) - 1 = 3 h2 + 1  (h1 = 0, s = 0), so h2 <= (X-1)/3
    return (X - 1) // 3 + 2
