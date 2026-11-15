#!/usr/bin/env python3
"""Replay for the reverse-complement collision theorem (RC_COLLISIONS.md).

G(u) = c u + k p(u) + phi (u mod 2),  p(0)=0, p(2v+1)=p(v), p(2v)=p(v)+v+(v mod 2).
rc_a(x) = 2^a - 1 - rev_a(x) (reverse-complement in an a-bit window).

1. p(rc_a(x)) = p(x) and x - rc_a(x) = rev_a(x) - rev_a(rc_a(x)) for every a <= 16, x < 2^a.
2. For four rows and L = 18: every family member (u, u') below 2^L is an exact collision; the member count
   equals the closed-form count; the share of all unordered collisions is reported.
3. Sharpness: on every family prefix, N(l, w) * w / 2^l >= (c/k)/2 (ordered pairs).
Runs in about a minute with numpy.
"""
import numpy as np
from fractions import Fraction as Fr
from collections import Counter

def rev(x, a):
    r = 0
    for i in range(a): r = (r << 1) | ((x >> i) & 1)
    return r

def ptable(n):
    p = np.zeros(n, dtype=np.int64)
    for u in range(1, n):
        v = u >> 1
        p[u] = p[v] if u & 1 else p[v] + v + (v & 1)
    return p

P = ptable(1 << 20)
# 1. the two identities
for a in range(1, 17):
    for x in range(1 << a):
        y = (1 << a) - 1 - rev(x, a)
        assert P[y] == P[x], (a, x)
        assert x - y == rev(x, a) - rev(y, a), (a, x)
print('1. p(rc_a x) = p(x) and x - rc_a x = rev_a x - rev_a rc_a x hold for all a <= 16, x < 2^a')

def family(c, k, L):
    """members (u', u, l, w) of the RC family below 2^L, by direct construction"""
    m = 1 + Fr(c, k); out = []
    for j in range(1, L):
        if (Fr(c, k) * 2 ** (j - 1)).denominator != 1: continue      # needs k | c 2^(j-1)
        top = int(m * 2 ** j)
        for a in range(2, L - j + 1):
            for x in range((1 << (a - 1)) | 1, 1 << a, 2):
                y = (1 << a) - 1 - rev(x, a)
                for r in range(1 << j):
                    w = top - 1 - r; t = rev(r, j)
                    u = (((w << a) | x) << j) | t; v = (((w << a) | y) << j) | t
                    if u < (1 << L): out.append((v, u, j + a, w))
    return out

def closed_count(c, k, L):
    m = 1 + Fr(c, k); tot = 0
    for j in range(1, L):
        if (Fr(c, k) * 2 ** (j - 1)).denominator != 1: continue
        top = int(m * 2 ** j)
        for a in range(2, L - j + 1):
            Z = 2 ** (L - j - a)                         # need w <= Z - 1, w = top - 1 - r
            good = max(0, min(2 ** j, Z - (top - 2 ** j)))
            tot += 2 ** (a - 2) * good
    return tot

L = 18
for (c, k, phi) in [(18, 3, 2), (2, 1, 0), (12, 3, 0), (189, 11, 0)]:
    N = 1 << L
    U = np.arange(N, dtype=np.int64); Gv = c * U + k * P[:N] + phi * (U & 1)
    _, cnt = np.unique(Gv, return_counts=True); allc = int((cnt * (cnt - 1) // 2).sum())
    fam = family(c, k, L)
    bad = sum(1 for v, u, l, w in fam if Gv[u] != Gv[v])
    assert bad == 0 and len(set((v, u) for v, u, l, w in fam)) == len(fam) == closed_count(c, k, L)
    m = 1 + Fr(c, k)
    lim = Fr(1, 2 ** (1 + (int(m) - 1).bit_length())) if m.denominator == 1 else None
    print('2. %dU+%dp+%d[odd], L=%d: %d family members, all exact collisions, = closed count; '
          '%.1f%% of the %d unordered collisions; members/2^L = %.4f (limit %s)'
          % (c, k, phi, L, len(fam), 100 * len(fam) / max(allc, 1), allc, len(fam) / N, lim))
    if fam:
        per = Counter((l, w) for v, u, l, w in fam)
        worst = min(Fr(2 * n * w, 2 ** l) for (l, w), n in per.items())
        assert worst >= Fr(c, k) / 2
        print('3.   sharpness: min over family prefixes of N_fam(l,w) w / 2^l (ordered) = %.4f >= c/(2k) = %s'
              % (float(worst), Fr(c, 2 * k)))
