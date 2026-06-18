#!/usr/bin/env python3
"""Shared exact interval / lonely-set / gap helpers for LRC14 verification."""
from fractions import Fraction as F
from math import gcd
from functools import reduce

# ----------------------------------------------------------------------
# Exact interval arithmetic on [0,1): list of (a,b) with 0<=a<b<=1, disjoint, sorted.
# ----------------------------------------------------------------------

def normalize(intervals):
    """Merge a list of (a,b) Fraction intervals into disjoint sorted form, clipped to [0,1]."""
    cl = []
    for a, b in intervals:
        if a < 0: a = F(0)
        if b > 1: b = F(1)
        if a < b:
            cl.append((a, b))
    if not cl:
        return []
    cl.sort()
    out = [list(cl[0])]
    for a, b in cl[1:]:
        if a <= out[-1][1]:
            if b > out[-1][1]:
                out[-1][1] = b
        else:
            out.append([a, b])
    return [(a, b) for a, b in out]

def measure(intervals):
    return sum((b - a for a, b in intervals), F(0))

def complement(intervals):
    """Complement within [0,1) of a normalized interval list."""
    intervals = normalize(intervals)
    out = []
    prev = F(0)
    for a, b in intervals:
        if a > prev:
            out.append((prev, a))
        prev = max(prev, b)
    if prev < 1:
        out.append((prev, F(1)))
    return out

def intersect(A, B):
    """Intersection of two normalized interval lists."""
    A = normalize(A); B = normalize(B)
    out = []
    i = j = 0
    while i < len(A) and j < len(B):
        lo = max(A[i][0], B[j][0])
        hi = min(A[i][1], B[j][1])
        if lo < hi:
            out.append((lo, hi))
        if A[i][1] < B[j][1]:
            i += 1
        else:
            j += 1
    return out

# ----------------------------------------------------------------------
# Danger set D_v = { tau in [0,1) : ||v tau|| < 1/14 } as exact intervals.
#   ||v tau|| < 1/14  <=>  v tau is within 1/14 of an integer k
#   <=>  tau in ( (k - 1/14)/v , (k + 1/14)/v )  for k = 0,1,...,v-1 (and wrap).
# Each tooth has width 2/(14 v) = 1/(7 v). There are v teeth. Total = v/(7v) = 1/7.
# We use a general half-width h (default 1/14) so we can test variants.
# ----------------------------------------------------------------------

def danger(v, h=F(1, 14)):
    """Danger set of speed v with half-width h on the v*tau axis. Returns intervals in [0,1)."""
    out = []
    for k in range(0, v + 1):
        a = (F(k) - h) / v
        b = (F(k) + h) / v
        out.append((a, b))  # normalize clips to [0,1]; k=0 left part and k=v right part wrap
    return normalize(out)

def lonely_set(C, h=F(1, 14)):
    """G_C = complement of union of danger sets of speeds in C."""
    U = []
    for v in C:
        U += danger(v, h)
    U = normalize(U)
    return complement(U)

def lonely_measure(C, h=F(1, 14)):
    return measure(lonely_set(C, h))

# ----------------------------------------------------------------------
# Exact gap M(S) tool (verbatim from the prompt, validated).
# ----------------------------------------------------------------------

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r

def g(S, t):
    return min(nrm(v * t) for v in S)

def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2 * k + 1, 2 * v) <= F(1, 2):
            C.add(F(2 * k + 1, 2 * v)); k += 1
    for i in range(len(S)):
        for j in range(i + 1, len(S)):
            for d in (S[i] + S[j], S[j] - S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2):
                        C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C

def M(S):
    b = F(0); at = None
    for t in cand(S):
        v = g(S, t)
        if v > b:
            b = v; at = t
    return b, at
