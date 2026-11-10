#!/usr/bin/env python3
"""AUDIT 2, late item L4 (auditor's own code).  Effect of deleting a vertex on the pair form beta_R.
beta_R is computed here from (4.4): beta_R(i,j) = parity of the number of directed paths between i and j (either
direction) on an even number of vertices -- a different route from the packet's deletion_rank.py (two-path-cover
matrix).  (4.4) itself is checked against mu/2 mod 2 on every class of orders 2-9 by x1_cube.c.)
For every class R of order m (gentourng representative) and every vertex u:
  rank over F_2 of  beta_R|_(V-u) + beta_(R-u),  and the pair (rank beta_R, rank beta_(R-u)).
Usage: python3 x7_deletion_rank.py m [m ...]"""
import subprocess, sys
from collections import Counter
G = '/opt/homebrew/bin/gentourng'

def beta_rows(n, out):
    inn = [0] * n
    for v in range(n):
        for w in range(n):
            if out[v] >> w & 1: inn[w] |= 1 << v
    P = [0] * n                                   # P[s] bit v = parity of the number of s -> v paths on an even number of vertices
    for s in range(n):
        par = [0] * (1 << n); par[1 << s] = 1 << s
        for S in range(1, 1 << n):
            if not S >> s & 1 or S == 1 << s: continue
            acc = 0; vv = S & ~(1 << s)
            while vv:
                v = (vv & -vv).bit_length() - 1; vv &= vv - 1
                if bin(par[S ^ (1 << v)] & inn[v]).count('1') & 1: acc |= 1 << v
            par[S] = acc
            if bin(S).count('1') % 2 == 0: P[s] ^= acc
    rows = [0] * n
    for i in range(n):
        for j in range(n):
            if i != j and ((P[i] >> j) ^ (P[j] >> i)) & 1: rows[i] |= 1 << j
    return rows

def rank(rows):
    rr = list(rows); rk = 0
    for c in range(len(rr)):
        p = next((i for i in range(rk, len(rr)) if rr[i] >> c & 1), None)
        if p is None: continue
        rr[p], rr[rk] = rr[rk], rr[p]
        for i in range(len(rr)):
            if i != rk and rr[i] >> c & 1: rr[i] ^= rr[rk]
        rk += 1
    return rk

def squeeze(masks, u):                             # drop vertex u from a list of bitmask rows, renumbering
    keep = [x for x in range(len(masks)) if x != u]
    return [sum(1 << k for k, w in enumerate(keep) if masks[x] >> w & 1) for x in keep]

for m in map(int, sys.argv[1:]):
    lines = subprocess.run([G, '-q', str(m)], capture_output=True, text=True).stdout.split()
    diffrank = Counter(); pairs = Counter(); sym = 0
    for line in lines:
        out = [0] * m; k = 0
        for i in range(m):
            for j in range(i + 1, m):
                if line[k] == '1': out[i] |= 1 << j
                else: out[j] |= 1 << i
                k += 1
        b = beta_rows(m, out); rb = rank(b)
        sym += any(((b[i] >> j) ^ (b[j] >> i)) & 1 for i in range(m) for j in range(m))
        for u in range(m):
            bu = beta_rows(m - 1, squeeze(out, u)); br = squeeze(b, u)
            diffrank[rank([x ^ y for x, y in zip(bu, br)])] += 1; pairs[(rb, rank(bu))] += 1
    print("order %d, %d classes, %d pairs (R, u); asymmetric beta: %d" % (m, len(lines), len(lines) * m, sym))
    print("  rank of beta_R|_(V-u) + beta_(R-u):", dict(sorted(diffrank.items())))
    print("  (rank beta_R, rank beta_(R-u)):", dict(sorted(pairs.items())))
