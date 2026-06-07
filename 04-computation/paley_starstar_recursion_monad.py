#!/usr/bin/env python3
"""
paley_starstar_recursion_monad.py
monad-explorer-2026-06-07 (deep-research, 5th session)

GOAL: find/verify the first-return recursion behind the number-theory-free identity
   (**)  S_k := sum_{sigma : EVEN-SERIES pattern of path [0..2k]} mu(0,sigma) = (-1)^k C_k
and corroborate the new structural findings:
  (1) even-series pattern COUNT = OEIS A215257(k+1)   (indecomposable deque-sortable perms)
  (2) the DRT-universal algebraic engine  M^2 = J - p I,  M 1 = 0  (M[x,y]=chi(y-x)).

We brute-enumerate set partitions of {0..2k}, keep EVEN-SERIES ones (connected, no loop,
no bridge, every cycle-space flow-line of even multiplicity), and break down S_k by:
  - cycle rank m,
  - the block-size multiset of block(0)  (how many times the start index's block is revisited),
  - a 'first-return' index r = smallest j>0 with block(j)==block(0).
This data is meant to expose the recursion (Catalan first-return) for S_k.
"""
import sys, math
from collections import defaultdict
import numpy as np


def set_partitions(c):
    c = list(c)
    if len(c) == 1:
        yield [c]; return
    f = c[0]
    for sm in set_partitions(c[1:]):
        for i, s in enumerate(sm):
            yield sm[:i] + [[f] + s] + sm[i + 1:]
        yield [[f]] + sm


def mu_partition(blocks):
    m = 1
    for B in blocks:
        b = len(B); m *= ((-1) ** (b - 1)) * math.factorial(b - 1)
    return m


def catalan(k):
    return math.comb(2 * k, k) // (k + 1)


def flow_lines(edges, nb):
    E = len(edges)
    Bm = np.zeros((nb, E))
    for ei, (u, v) in enumerate(edges):
        Bm[v, ei] += 1.0; Bm[u, ei] -= 1.0
    _, s, vh = np.linalg.svd(Bm)
    rank = int((s > 1e-9).sum()); m = E - rank
    if m == 0:
        return [("Z",)] * E, 0
    ns = vh[rank:]
    lines = []
    for e in range(E):
        v = ns[:, e]
        if np.max(np.abs(v)) < 1e-7:
            lines.append(("Z",)); continue
        v = v / np.max(np.abs(v))
        for x in v:
            if abs(x) > 1e-7:
                if x < 0: v = -v
                break
        lines.append(tuple(round(float(x), 6) for x in v))
    return lines, m


def even_series_info(blocks, L):
    pos2blk = {}
    for bi, B in enumerate(blocks):
        for pos in B: pos2blk[pos] = bi
    edges = [(pos2blk[i], pos2blk[i + 1]) for i in range(L)]
    nb = len(blocks)
    if any(u == v for (u, v) in edges):
        return None
    adj = defaultdict(list)
    for (u, v) in edges:
        adj[u].append(v); adj[v].append(u)
    seen = {0}; stk = [0]
    while stk:
        x = stk.pop()
        for w in adj[x]:
            if w not in seen: seen.add(w); stk.append(w)
    if len(seen) != nb:
        return None
    lines, m = flow_lines(edges, nb)
    if m == 0 or any(l == ("Z",) for l in lines):
        return None
    grp = defaultdict(int)
    for l in lines: grp[l] += 1
    if not all(c % 2 == 0 for c in grp.values()):
        return None
    # first return of index 0
    b0 = pos2blk[0]
    r = next((j for j in range(1, L + 1) if pos2blk[j] == b0), None)
    blk0_size = sum(1 for pos in range(L + 1) if pos2blk[pos] == b0)
    return m, nb, r, blk0_size


KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 5
print("=" * 70)
print("(**) S_k = sum_{even-series} mu(0,sigma);  target (-1)^k C_k")
print(f"{'k':>2} {'S_k':>8} {'(-1)^kC_k':>10} {'match':>6} {'#es':>6}  A215257?")
A215257 = {1:1, 2:3, 3:13, 4:67, 5:383, 6:2345}
alldata = {}
for k in range(1, KMAX + 1):
    L = 2 * k
    S = 0; cnt = 0
    by_m = defaultdict(int)             # S_k contribution by cycle rank
    by_r = defaultdict(int)             # by first-return index of 0
    by_b0 = defaultdict(int)            # by block(0) size
    cnt_by_r = defaultdict(int)
    for blocks in set_partitions(range(L + 1)):
        info = even_series_info(blocks, L)
        if info is None: continue
        m, nb, r, b0 = info
        w = mu_partition(blocks)
        S += w; cnt += 1
        by_m[m] += w; by_r[r] += w; by_b0[b0] += w
        cnt_by_r[r] += 1
    tgt = (-1) ** k * catalan(k)
    ok = (S == tgt); okcnt = (cnt == A215257.get(k))
    print(f"{k:>2} {S:>8} {tgt:>10} {str(ok):>6} {cnt:>6}  count==A215257(k+1): {okcnt}")
    alldata[k] = (by_m, by_r, by_b0, cnt_by_r)

print("=" * 70)
print("Breakdown of S_k by CYCLE RANK m  (m=k bigon-trees ... m=1 single 2k-cycle):")
for k in range(1, KMAX + 1):
    by_m = alldata[k][0]
    print(f"  k={k}: " + ", ".join(f"m={m}:{by_m[m]:+d}" for m in sorted(by_m)))
print("Breakdown of S_k by FIRST-RETURN index r of block(0):")
for k in range(1, KMAX + 1):
    by_r = alldata[k][1]
    print(f"  k={k}: " + ", ".join(f"r={r}:{by_r[r]:+d}" for r in sorted(by_r, key=lambda x:(x is None,x))))
print("Breakdown of S_k by block(0) SIZE (revisits of start):")
for k in range(1, KMAX + 1):
    by_b0 = alldata[k][2]
    print(f"  k={k}: " + ", ".join(f"|b0|={b}:{by_b0[b]:+d}" for b in sorted(by_b0)))
print("COUNTS by first-return r (unsigned, structure check):")
for k in range(1, KMAX + 1):
    cr = alldata[k][3]
    print(f"  k={k}: " + ", ".join(f"r={r}:{cr[r]}" for r in sorted(cr, key=lambda x:(x is None,x))))
