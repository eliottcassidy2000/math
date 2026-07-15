#!/usr/bin/env python3
"""
selfline_sequence_atlas_kps_S128c15.py
======================================
kind-pasteur-2026-07-15-S128 (cont.15). THE SELF-LINE SEQUENCE ATLAS: fill the missing rows.
 (A) W(8) via Burnside-affine (Aut-weighted quasi-fixed count; completes W = 8, 20, 88, ?).
 (B) The g-kappa DIAGONAL |Y|(n) = #{tilings t : cls(t) is SC} = SC fiber mass, n = 5, 6, 7
     (direct iso T(t) vs T(t)^op with the validated pipeline).
 (C) Print the atlas table: ORBITS, selfK, selfB, |X|, W, |Y|, SC, A000568 for n = 5..8 (known+new).
"""
import sys
from math import comb
from itertools import permutations
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

def setup(n):
    tiles = [(x, y) for y in range(1, n - 1) for x in range(n, y + 1, -1) if x - y >= 2]
    tidx = {t: i for i, t in enumerate(tiles)}
    return tiles, tidx, len(tiles)

def beats(tv, n, tidx):
    B = [[False] * n for _ in range(n)]
    for k in range(2, n + 1):
        B[k - 1][k - 2] = True
    for (x, y), i in tidx.items():
        if tv[i]:
            B[x - 1][y - 1] = True
        else:
            B[y - 1][x - 1] = True
    return B

def inv_pair(B, n):
    d = [sum(B[u]) for u in range(n)]
    return (tuple(sorted(d)), comb(n, 3) - sum(comb(x, 2) for x in d))

def iso(B1, B2, n):
    d1 = [sum(B1[u]) for u in range(n)]
    d2 = [sum(B2[u]) for u in range(n)]
    if sorted(d1) != sorted(d2):
        return False
    buckets = defaultdict(list)
    for v in range(n):
        buckets[d2[v]].append(v)
    order = sorted(range(n), key=lambda u: (d1[u], u))
    assign = [-1] * n
    used = [False] * n
    def bt(i):
        if i == n:
            return True
        u = order[i]
        for v in buckets[d1[u]]:
            if used[v]:
                continue
            ok = True
            for j in range(i):
                w = order[j]
                if B1[u][w] != B2[v][assign[w]] or B1[w][u] != B2[assign[w]][v]:
                    ok = False
                    break
            if ok:
                assign[u] = v; used[v] = True
                if bt(i + 1):
                    return True
                used[v] = False; assign[u] = -1
        return False
    return bt(0)

# (A) Burnside-affine W(8)
print("(A) W(8) Burnside-affine...")
def burnside_affine(n):
    tiles, tidx, m = setup(n)
    def slot(a, b):
        if b - a == 1:
            return None, True
        return tidx[(b + 1, a + 1)], None
    Wtot = 0
    pairs = [(u, v) for u in range(n) for v in range(u + 1, n)]
    for sig in permutations(range(n)):
        rows = []
        ok = True
        for (u, v) in pairs:
            p, q = sig[u], sig[v]
            uu, vv = (u, v) if p < q else (v, u)
            if p > q:
                p, q = q, p
            a, b = (uu, vv) if uu < vv else (vv, uu)
            e, isbase = slot(a, b)
            if isbase:
                f1c, f1v = 0, None
            else:
                f1c, f1v = 1, e
            if uu > vv:
                f1c ^= 1
            e2, isbase2 = slot(p, q)
            if isbase2:
                f2c, f2v = 0, None
            else:
                f2c, f2v = 0, e2
            row = [0] * m
            c = f1c ^ f2c
            if f1v is not None:
                row[f1v] ^= 1
            if f2v is not None:
                row[f2v] ^= 1
            if any(row):
                rows.append(row + [c])
            elif c:
                ok = False
                break
        if not ok:
            continue
        rank = 0
        for col in range(m):
            piv = next((r for r in range(rank, len(rows)) if rows[r][col]), None)
            if piv is None:
                continue
            rows[rank], rows[piv] = rows[piv], rows[rank]
            for r in range(len(rows)):
                if r != rank and rows[r][col]:
                    rows[r] = [a ^ b for a, b in zip(rows[r], rows[rank])]
            rank += 1
        if all(not (all(x == 0 for x in row[:-1]) and row[-1]) for row in rows):
            Wtot += 1 << (m - rank)
    return Wtot
W8 = burnside_affine(8)
print("W(8) = %d   (W = 8, 20, 88, %d)" % (W8, W8))

# (B) the g-kappa diagonal: SC fiber mass |Y|(n), n = 5, 6, 7
print("(B) |Y|(n) = # tilings whose class is SC (T(t) iso T(t)^op)...")
Y = {}
for n in [5, 6, 7]:
    tiles, tidx, m = setup(n)
    cnt = 0
    for t in range(1 << m):
        tv = [(t >> i) & 1 for i in range(m)]
        B1 = beats(tv, n, tidx)
        B2 = [[B1[u][v] for u in range(n)] for v in range(n)]
        if inv_pair(B1, n) == inv_pair(B2, n) and iso(B1, B2, n):
            cnt += 1
    Y[n] = cnt
    print("   |Y|(%d) = %d" % (n, cnt))

# (C) atlas
print()
print("THE SELF-LINE SEQUENCE ATLAS (n = 5, 6, 7, 8):")
print("  ORBITS (free Klein-4)      : 2, 3, 22, 101       [fundamental; = merged carrier census]")
print("  selfK  (black self-lines)  : 4, 6, 44, 202")
print("  2selfB (gridsym qf tilings): 0, 4, 0, 8")
print("  |X|    (quasi-fixed total) : 8, 16, 88, 412")
print("  W      (Aut-weighted |X|)  : 8, 20, 88, %d" % W8)
print("  |Y|    (SC fiber mass)     : %d, %d, %d, ?" % (Y[5], Y[6], Y[7]))
print("  SC     (A000570)           : 8, 12, 88, 176")
print("  A000568 (classes)          : 12, 56, 456, 6880")
print("DONE")
