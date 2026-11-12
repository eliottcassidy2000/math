#!/usr/bin/env python3
"""Verifier: extended invariant-plane scan (beyond the lane's single-block scope).
Letters = all products of length <= 2 of the 134 invertible blocks of order <= 5, and of length <= 4 of the 21
invertible blocks of order <= 4.  For every rational invariant plane (integer left eigenvalue lam, exact check),
record the words acting unimodularly on it (|det| = |lam|) and whether two such words fail to commute.
"""
import math
import numpy as np
import v_alg as A

def int_roots_fast(m):
    t, c, d = A.tr(m), A.c2(m), A.det(m)
    out = set()
    for r in np.roots([1.0, -float(t), float(c), -float(d)]):
        if abs(r.imag) < 1e-6 * max(1.0, abs(r)):
            x = int(round(r.real))
            for lam in (x - 1, x, x + 1):
                if lam ** 3 - t * lam ** 2 + c * lam - d == 0:
                    out.add(lam)
    return out

def planes_of(m):
    res = []
    for lam in int_roots_fast(m):
        B = [m[0] - lam, m[1], m[2], m[3], m[4] - lam, m[5], m[6], m[7], m[8] - lam]
        cols = [(B[0], B[3], B[6]), (B[1], B[4], B[7]), (B[2], B[5], B[8])]
        n = None
        for i in range(3):
            for j in range(i + 1, 3):
                c = A.cross(cols[i], cols[j])
                if c != (0, 0, 0):
                    n = A.prim(c)
        if n is not None:
            res.append((lam, n))
    return res

def scan(words, label, orb):
    pl = {}
    for m in words:
        for lam, n in planes_of(m):
            pl.setdefault(n, []).append((m, lam))
    shared = {n: v for n, v in pl.items() if len(v) >= 2}
    print(f"{label}: {len(words)} distinct invertible words; with a rational invariant plane: {sum(len(v) for v in pl.values())} (word, plane) pairs; "
          f"planes shared by >= 2 words: {len(shared)}")
    for n, v in sorted(shared.items(), key=lambda kv: -len(kv[1])):
        uni = [m for m, lam in v if abs(A.det(m)) == abs(lam)]
        nonc = any(A.mul(a, b) != A.mul(b, a) for a in uni for b in uni)
        on = sum(1 for s in orb if sum(x * y for x, y in zip(n, s)) == 0)
        print(f"   normal {n}: {len(v)} words; acting unimodularly on the plane: {len(uni)}; non-commuting unimodular pair: {nonc}; orbit states (h <= 2000, order <= 4 orbit) on the plane: {on}")

lib = {k: A.load(k) for k in (4, 5)}
orb = A.small_orbit(lib[4], 2000)
inv5 = [m for m in lib[5] if A.det(m) != 0]
w = set(inv5) | set(A.mul(a, b) for a in inv5 for b in inv5)
scan(w, "order <= 5, length <= 2", orb)
inv4 = [m for m in lib[4] if A.det(m) != 0]
w = set(inv4); cur = set(inv4)
for _ in range(3):
    cur = set(A.mul(a, b) for a in cur for b in inv4)
    w |= cur
scan(w, "order <= 4, length <= 4", orb)
