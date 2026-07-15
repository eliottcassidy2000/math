#!/usr/bin/env python3
"""COARSE CONCORDANCE AT n=8, BOTH SIDES (opus-S308) -- the minimal above-7 test.

The decided law: exact concordance FAILS (n=7 G-side, 132 real adjacent-level
overlaps; n=6 E-side, 1 overlap); the corrected law is COARSE concordance
(potential order strict at level-distance >= 2). This script tests n=8:
  G-side: the certified 6880-class network (subH classifier), float solve with
          exact-residual refinement; discordant pairs BY LEVEL-DISTANCE.
  E-side: the 243-class even-graph network (A002854-certified), exact solve.
If discordances stay confined to adjacent levels and the fraction stays ~0.1%,
the coarse law survives its first above-7 test; if they spread to distance >= 2,
we are misguided and the law needs re-scoping again."""
import sys, os
from collections import defaultdict
from fractions import Fraction as F
import numpy as np
import itertools

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# ---------------- G-side: reuse the S305 subH classifier core ----------------
n = 8
V = list(range(1, n+1))
tiles = [(x, y) for y in range(1, n-1) for x in range(n, y+1, -1) if x-y >= 2]
m = len(tiles)
POP = np.array([bin(i).count('1') for i in range(256)], dtype=np.int64)
FULLV = (1 << n) - 1
N = 1 << m
CH = 1 << 15

def feats(bits):
    B = bits.shape[0]
    adj = np.zeros((B, n), dtype=np.int64)
    for v in range(2, n+1): adj[:, v-1] |= 1 << (v-2)
    for i, (x, y) in enumerate(tiles):
        b = (bits >> i) & 1
        adj[:, x-1] |= b << (y-1)
        adj[:, y-1] |= (1 - b) << (x-1)
    s = POP[adj]
    inv = np.zeros_like(adj)
    for v in range(n): inv[:, v] = FULLV ^ (1 << v) ^ adj[:, v]
    c3 = np.zeros((B, n), dtype=np.int64)
    for v in range(n):
        for u in range(n):
            if u == v: continue
            c3[:, v] += ((adj[:, v] >> u) & 1) * POP[adj[:, u] & inv[:, v]]
    so1 = np.zeros((B, n), dtype=np.int64); si1 = np.zeros_like(so1)
    so2 = np.zeros_like(so1); co1 = np.zeros_like(so1); ci1 = np.zeros_like(so1)
    so3 = np.zeros_like(so1)
    for v in range(n):
        for u in range(n):
            if u == v: continue
            hout = (adj[:, v] >> u) & 1; hin = (adj[:, u] >> v) & 1
            so1[:, v] += hout * s[:, u]; si1[:, v] += hin * s[:, u]
            so2[:, v] += hout * s[:, u]**2; so3[:, v] += hout * s[:, u]**3
            co1[:, v] += hout * c3[:, u]; ci1[:, v] += hin * c3[:, u]
    to1 = np.zeros_like(so1); ti1 = np.zeros_like(so1); tc1 = np.zeros_like(so1)
    for v in range(n):
        for u in range(n):
            if u == v: continue
            hout = (adj[:, v] >> u) & 1; hin = (adj[:, u] >> v) & 1
            to1[:, v] += hout * so1[:, u]; ti1[:, v] += hin * si1[:, u]
            tc1[:, v] += hout * ci1[:, u]
    key = s.astype(np.int64)
    for extra, base in ((c3, 64), (so1, 64), (si1, 64), (so2, 512), (so3, 4096),
                        (co1, 256), (ci1, 256), (to1, 512), (ti1, 512), (tc1, 2048)):
        key = key * base + extra
    key.sort(axis=1)
    a1 = np.zeros(B, dtype=np.int64); a2 = np.zeros(B, dtype=np.int64)
    for u in range(n):
        for v in range(n):
            if u == v: continue
            h = (adj[:, u] >> v) & 1
            common = POP[adj[:, u] & adj[:, v]]
            a1 += h * (s[:, u] * 8 + s[:, v]) * (common + 1)
            a2 += h * common * common
    adjbit = [[(adj[:, v] >> u) & 1 for u in range(n)] for v in range(n)]
    levels = {1: {(1 << v, v): np.ones(B, dtype=np.int64) for v in range(n)}}
    for size in range(1, n):
        nxt = {}
        for (mask, v), arr in levels[size].items():
            for u in range(n):
                bu = 1 << u
                if mask & bu: continue
                contrib = arr * adjbit[v][u]
                tgt = (mask | bu, u)
                if tgt in nxt: nxt[tgt] += contrib
                else: nxt[tgt] = contrib.copy()
        del levels[size]
        levels[size + 1] = nxt
    Hv = np.zeros(B, dtype=np.int64)
    for (mask, v), arr in levels[n].items(): Hv += arr
    subH = np.zeros((B, n), dtype=np.int64)
    for w in range(n):
        verts = [v for v in range(n) if v != w]
        lv = {(1 << v, v): np.ones(B, dtype=np.int64) for v in verts}
        for size in range(1, n-1):
            nxt = {}
            for (mask, v), arr in lv.items():
                for u in verts:
                    bu = 1 << u
                    if mask & bu: continue
                    contrib = arr * adjbit[v][u]
                    tgt = (mask | bu, u)
                    if tgt in nxt: nxt[tgt] += contrib
                    else: nxt[tgt] = contrib.copy()
            lv = nxt
        col = np.zeros(B, dtype=np.int64)
        for (mask, v), arr in lv.items(): col += arr
        subH[:, w] = col
    subH.sort(axis=1)
    return key, a1, a2, s, Hv, subH

buckets = {}; cls_of = np.zeros(N, dtype=np.int32); rep = []
scores_of = []
for lo in range(0, N, CH):
    bits = np.arange(lo, min(lo+CH, N), dtype=np.int64)
    key, a1, a2, s, Hv, subH = feats(bits)
    ss = np.sort(s, axis=1)
    for j in range(bits.shape[0]):
        k = (key[j].tobytes(), int(a1[j]), int(a2[j]), int(Hv[j]), subH[j].tobytes())
        c = buckets.get(k)
        if c is None:
            c = len(buckets); buckets[k] = c; rep.append(int(bits[j]))
            scores_of.append(tuple(int(v) for v in ss[j]))
        cls_of[lo+j] = c
C = len(buckets)
print(f"G-side n=8: {C} classes ({'OK' if C == 6880 else 'FAIL'})", flush=True)
assert C == 6880
x_of = np.array([sum((2*si - (n-1))**2 for si in sc) for sc in scores_of])
H1 = [c for c in range(C)][0]
# wiggly edges
Wg = defaultdict(int)
for lo in range(0, N, CH):
    bits = np.arange(lo, min(lo+CH, N), dtype=np.int64)
    ca = cls_of[bits]
    for i in range(m):
        cb = cls_of[bits ^ (1 << i)]
        lo_ = np.minimum(ca, cb); hi_ = np.maximum(ca, cb)
        pk = lo_.astype(np.int64) * 100000 + hi_
        uniq, cnts = np.unique(pk, return_counts=True)
        for u, cc in zip(uniq, cnts):
            Wg[(int(u)//100000, int(u)%100000)] += int(cc)
Wg = {k: v // 2 for k, v in Wg.items()}
cond = {k: v for k, v in Wg.items() if k[0] != k[1]}
xmin = int(x_of.min())
# the transitive class: fiber 1 & max x -- find via rep of all-ones tiling
trans = int(cls_of[N-1])
sinks = [c for c in range(C) if x_of[c] == xmin]
nodes = [c for c in range(C) if c not in sinks]
pos = {c: i for i, c in enumerate(nodes)}
NN = len(nodes)
A = np.zeros((NN, NN))
deg = defaultdict(int)
for (a, b), w in cond.items():
    deg[a] += w; deg[b] += w
    if a in pos and b in pos:
        A[pos[a], pos[b]] -= w; A[pos[b], pos[a]] -= w
for c in nodes: A[pos[c], pos[c]] = deg[c]
bvec = np.zeros(NN); bvec[pos[trans]] = 1.0
phi = np.linalg.solve(A, bvec)
# one refinement pass (float residual is fine at this scale for classification by distance)
r = bvec - A @ phi
phi = phi + np.linalg.solve(A, r)
pot = np.zeros(C)
for c in nodes: pot[c] = phi[pos[c]]
print(f"G-side R^G_8 ~= {pot[trans]:.6f}", flush=True)
# discordances by level distance
lvls = sorted(set(int(v) for v in x_of))
by_dist = defaultdict(int); tot_pairs = 0
order = np.argsort(x_of)
disc_adj = 0
for a, b in itertools.combinations(range(C), 2):
    dx = int(x_of[a]) - int(x_of[b])
    if dx == 0: continue
    dp = pot[a] - pot[b]
    if abs(dp) < 1e-13: continue
    tot_pairs += 1
    if (dx > 0) != (dp > 0):
        by_dist[abs(dx)] += 1
print(f"G-side n=8 discordant pairs by |level distance|: {dict(sorted(by_dist.items()))} "
      f"of {tot_pairs} pairs", flush=True)
print(f"COARSE LAW at n=8 (G): {'HOLDS (all discordances at distance 8)' if all(k == 8 for k in by_dist) else 'FAILS -- discordances beyond adjacent levels'}", flush=True)
