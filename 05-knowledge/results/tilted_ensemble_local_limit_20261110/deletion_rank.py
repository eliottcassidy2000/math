#!/usr/bin/env python3
"""Exploration: how does the pair form beta_R change when a vertex is deleted?
For every class R of order m (gentourng) and every vertex u: rank of beta_R, of beta_{R-u}, of the restriction of
beta_R to V - u, and of the difference (beta_R restricted) + beta_{R-u} over F_2."""
import subprocess, sys
from collections import Counter

def beta(n, out):
    full = (1 << n) - 1
    inn = [sum(1 << w for w in range(n) if out[w] >> v & 1) for v in range(n)]
    f = [[0] * n for _ in range(1 << n)]; g = [[0] * n for _ in range(1 << n)]
    for S in range(1, full + 1):
        for v in range(n):
            if not S >> v & 1: continue
            R = S & ~(1 << v)
            if R == 0: f[S][v] = g[S][v] = 1; continue
            a = 0; mm = R & inn[v]
            while mm: w = (mm & -mm).bit_length() - 1; mm &= mm - 1; a += f[R][w]
            b = 0; mm = R & out[v]
            while mm: w = (mm & -mm).bit_length() - 1; mm &= mm - 1; b += g[R][w]
            f[S][v] = a; g[S][v] = b
    rows = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            rest = full & ~(1 << i) & ~(1 << j); mu = 0; S0 = rest
            while True:
                S = S0 | (1 << j); mu += f[S][j] * g[full ^ S][i]
                S = S0 | (1 << i); mu += f[S][i] * g[full ^ S][j]
                if S0 == 0: break
                S0 = (S0 - 1) & rest
            assert mu % 2 == 0
            if (mu >> 1) & 1: rows[i] |= 1 << j; rows[j] |= 1 << i
    return rows

def rank(rows):
    rr = list(rows); rk = 0; n = len(rr)
    for c in range(n):
        p = next((i for i in range(rk, n) if rr[i] >> c & 1), None)
        if p is None: continue
        rr[p], rr[rk] = rr[rk], rr[p]
        for i in range(n):
            if i != rk and rr[i] >> c & 1: rr[i] ^= rr[rk]
        rk += 1
    return rk

def delete(out, u):
    n = len(out); keep = [x for x in range(n) if x != u]; pos = {x: k for k, x in enumerate(keep)}
    return [sum(1 << pos[w] for w in keep if out[x] >> w & 1) for x in keep]

def restrict(rows, u):
    n = len(rows); keep = [x for x in range(n) if x != u]; pos = {x: k for k, x in enumerate(keep)}
    return [sum(1 << pos[w] for w in keep if rows[x] >> w & 1) for x in keep]

m = int(sys.argv[1])
res = subprocess.run(['/opt/homebrew/bin/gentourng', '-q', str(m)], capture_output=True, text=True).stdout.split()
stat = Counter(); diffrank = Counter(); drop = Counter()
for line in res:
    out = [0] * m; k = 0
    for i in range(m):
        for j in range(i + 1, m):
            if line[k] == '1': out[i] |= 1 << j
            else: out[j] |= 1 << i
            k += 1
    b = beta(m, out); rb = rank(b)
    for u in range(m):
        bu = beta(m - 1, delete(out, u)); br = restrict(b, u)
        d = [x ^ y for x, y in zip(bu, br)]
        diffrank[rank(d)] += 1; drop[(rb, rank(bu))] += 1; stat[(rb, rank(br))] += 1
print("order", m, "classes", len(res))
print("  rank of the difference  beta_R|_(V-u) + beta_(R-u):", dict(sorted(diffrank.items())))
print("  (rank beta_R, rank beta_(R-u)):", dict(sorted(drop.items())))
