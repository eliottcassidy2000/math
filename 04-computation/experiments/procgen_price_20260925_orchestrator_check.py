#!/usr/bin/env python3
"""Orchestrator's independent re-implementation of the G_L construction (THM-4475), written from the note's prose
description only (not the lane's code).  For L = 8, 12, 16 and N = 10^6 it checks: no stuck n; every n <= N/2
descends within L steps under the final bits; flip density on pairs <= N/4; B-rescues only at n = 507 (mod 512)."""
import sys
from array import array
def run(L, N):
    MAXP = (N * 2**(L//2+3)) // 2 + 10   # generous pair-index cap for paths
    bit = bytearray(MAXP)   # 0 free, 1 frozen-0, 2 frozen-1
    def step(x):
        i = (x + 1) // 2
        flipped = (bit[i] == 2)
        if x % 2 == 1:
            return (x - 1) // 2 if flipped else (3 * x + 1) // 2
        else:
            return 3 * x // 2 if flipped else x // 2
    def path(n):
        # default path; returns list of points before descent and descent flag
        pts = [n]; x = n
        for j in range(1, L + 1):
            x = step(x)
            if x < n: return pts, True
            pts.append(x)
        return pts, False
    def freeze(pts):
        for x in pts:
            i = (x + 1) // 2
            if bit[i] == 0: bit[i] = 1
    T = lambda x: x // 2 if x % 2 == 0 else (3 * x + 1) // 2
    stats = {'A': 0, 'F': 0, 'B': 0, 'stuck': 0}
    Bres = []
    for n in range(3, N + 1):
        pts, ok = path(n)
        if ok: freeze(pts); continue
        def try_flip(v):
            i = (v + 1) // 2
            if bit[i] != 0: return False
            bit[i] = 2
            p2, ok2 = path(n)
            if ok2: freeze(p2); return True
            bit[i] = 0; return False
        def optA(): return n % 2 == 1 and try_flip(n)
        def optF():
            # first point w of D(n) (pos>=1) with w = T(v), v odd going up, w = 3 mod 4, w < 2n
            for k in range(1, len(pts)):
                v, w = pts[k-1], pts[k]
                if v % 2 == 1 and bit[(v+1)//2] != 2 and w == (3*v+1)//2 and w % 4 == 3 and w < 2*n:
                    return try_flip(w)
            return False
        def optB():
            t = T(n)
            return t % 2 == 1 and try_flip(t)
        order = (('A', optA), ('F', optF), ('B', optB)) if n % 8 == 3 else (('F', optF), ('B', optB), ('A', optA))
        for name, fn in order:
            if fn():
                stats[name] += 1
                if name == 'B': Bres.append(n)
                break
        else:
            stats['stuck'] += 1
    # verification: every n <= N/2 descends within L using final bits (all pairs on paths frozen)
    bad = 0
    for n in range(3, N // 2):
        _, ok = path(n)
        if not ok: bad += 1
    nflip = sum(1 for i in range(1, N // 4) if bit[i] == 2)
    return stats, bad, nflip / (N // 4), Bres
for L in (8, 12, 16):
    N = 1_000_000
    stats, bad, dens, Bres = run(L, N)
    Bmod = sorted({n % 512 for n in Bres})
    print(f"L={L}: rescues {stats}; n<=N/2 failing L-step descent: {bad}; flip density (pairs <= N/4): {dens:.5f}; B-rescues mod 512: {Bmod}")
    assert stats['stuck'] == 0 and bad == 0 and Bmod == [507]
print("ALL CHECKS PASSED")
