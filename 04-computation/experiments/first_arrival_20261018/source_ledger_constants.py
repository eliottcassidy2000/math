#!/usr/bin/env python3
"""Exact source-level constants of the ROOT transport (mac-mini-2026-10-18-sixcolour; THM-4690 statement 4).

Sources: every odd m <= 4N carries K_N(m) = w(m) + thr(m)[m < 2N] and is aggregated on its seed x = seed(m) (strip
B(x) = 4x + 1). So c_0(x) = sum_j K_N(B^j(x)) over the B-chain x, Bx, B^2 x, ... <= 4N. Since B^j(x) = x + j (mod 3),
the chain residues cycle with period 3 (the Solfeggio three-lift: B^3(x) = x + 3 mod 9). Below 2N the weight is
kappa(r) = 0, -3, 3 for r = 0, 1, 2 (zero sum over a period); in [2N, 4N] it is lambda(r) = -1, -1, 2; at most one chain
element lies in [2N, 4N]. Hence c_0(x) depends only on (x mod 3, L mod 3, whether the chain ends in [2N, 4N]), where L is
the number of chain elements below 2N, and the x-ranges with a given L have length 2N/4^L (ending in [2N, 4N]) and
4N/4^L (not), with the seed density 3/8 of the integers.

Derived constants (exact rationals; the 7 is 63 = 4^3 - 1 divided by 9):
  E_0(N)  = sum_x c_0(x)^2                                  = (34/7) N + O(log N)
  SIB(N)  = sum_(x = 1 mod 8, 2x+1 <= 4N) c_0(x) c_0(2x+1)  = (5/7) N + O(log N)
Siblings x (a = 2) and 2x+1 (a = 1) have the same image P(x) = P(2x+1) and literal clocks differing by exactly 1, so the
pair enters every stopping set together (unless one of them lies in it) with clock lag 1 and opposite signs.
Part 1 derives the constants from the residue tables; part 2 checks them against brute-force sums.
"""
from fractions import Fraction as Fr
import numpy as np

kappa = [0, -3, 3]
lam = [-1, -1, 2]

def P(r, l):                                   # prefix sum of kappa over l = L mod 3 chain residues starting at r
    return sum(kappa[(r + j) % 3] for j in range(l))

def c0(r, L, ends):                            # source coefficient of a seed x = r (mod 3) with chain parameter L
    return P(r, L % 3) + (lam[(r + L) % 3] if ends else 0)

def geo(cond):                                 # sum over L >= 0 with cond(L) of 4^-L, exactly (period 3 in L)
    per = [Fr(1, 4 ** l) for l in range(3)]
    return sum(per[l] for l in range(3) if cond(l)) / (1 - Fr(1, 64))

def avg(f):
    return Fr(sum(f(r) for r in range(3)), 3)

# --- part 1: derivation -------------------------------------------------------------------------------------------
# E_0: seeds have density 3/8; L = 0 occurs only with 'ends' (x in [2N, 4N]); for L >= 1 both kinds occur.
E0 = Fr(3, 8) * (sum(2 * geo(lambda l, L=L: l == L) * avg(lambda r, L=L: c0(r, L, True) ** 2) for L in range(3))
                 + sum(4 * (geo(lambda l, L=L: l == L) - (1 if L == 0 else 0)) * avg(lambda r, L=L: c0(r, L, False) ** 2)
                       for L in range(3)))
# SIB: x = 1 (mod 8) has density 1/8, x mod 3 = r, and 2x + 1 = (2r + 1) mod 3. With y = 3x + 1 the sibling has
# 3(2x+1) + 1 = 2y + 2, so: x ending in [2N, 4N] at level L  <->  2x+1 not ending at level L        (case A, length 2N/4^L)
#                           x not ending at level L (L >= 1) <->  2x+1 ending at level L - 1        (case B, length 4N/4^L)
# L = 0 with 'ends' would need 2x + 1 > 4N, so case A starts at L = 1.
def r2(r): return (2 * r + 1) % 3
SIB = Fr(1, 8) * (sum(2 * (geo(lambda l, L=L: l == L) - (1 if L == 0 else 0)) * avg(lambda r, L=L: c0(r, L, True) * c0(r2(r), L, False))
                      for L in range(3))
                  + sum(4 * (geo(lambda l, L=L: l == L) - (1 if L == 0 else 0)) * avg(lambda r, L=L: c0(r, L, False) * c0(r2(r), L + 2, True))
                        for L in range(3)))
print('derived: E_0/N ->', E0, '=', float(E0), ';  SIB/N ->', SIB, '=', float(SIB))
assert E0 == Fr(34, 7) and SIB == Fr(5, 7)
# note: c0(r2(r), L + 2, True) uses L - 1 = L + 2 (mod 3) for the sibling's level; c0 depends on L only mod 3

# --- part 2: brute force --------------------------------------------------------------------------------------------
def c0_array(N):
    M = 4 * N
    m = np.arange(1, M + 1, 2, dtype=np.int64)
    r = m % 3
    K = np.where(r == 2, 2, -1) + np.where(m < 2 * N, np.where(r == 0, 1, np.where(r == 1, -2, 1)), 0)
    s = m.copy()
    while True:
        msk = (s & 7) == 5
        if not msk.any(): break
        s[msk] = (s[msk] - 1) >> 2
    c = np.zeros((M + 1) // 2, dtype=np.int64)
    np.add.at(c, (s - 1) // 2, K)
    return c

worst_e, worst_s = 0.0, 0.0
for N in [1, 2, 3, 10, 100, 1000, 4096, 10007, 65536, 100003, 262143, 1000003, 2 ** 21, 4000037]:
    c = c0_array(N)
    e0 = int((c * c).sum())
    x = np.arange(1, 4 * N + 1, 8, dtype=np.int64)          # x = 1 (mod 8)
    x = x[2 * x + 1 <= 4 * N]
    sib = int((c[(x - 1) // 2] * c[x]).sum())              # index of 2x+1 is x
    de, ds = e0 - Fr(34, 7) * N, sib - Fr(5, 7) * N
    lg = max(1.0, np.log2(N))
    worst_e, worst_s = max(worst_e, abs(float(de)) / lg), max(worst_s, abs(float(ds)) / lg)
    print('N = %8d  E_0 = %10d  E_0 - 34N/7 = %+9.2f   SIB = %9d  SIB - 5N/7 = %+8.2f   |c_0| <= %d' % (N, e0, float(de), sib, float(ds), int(abs(c).max())))
print('max |E_0 - 34N/7| / log2 N = %.3f ; max |SIB - 5N/7| / log2 N = %.3f' % (worst_e, worst_s))
