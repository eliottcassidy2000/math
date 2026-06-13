#!/usr/bin/env python3
"""
paley_starstar_noncrossing_monad.py
monad-explorer-2026-06-07 (deep-research, 5th session)

NEGATIVE RESULT (records why the *naive* genus-0 reading of (★★) is wrong).

Tempting conjecture: the Catalan law (★★) localizes onto the EVEN-SERIES patterns whose
COINCIDENCE PARTITION sigma is NON-CROSSING (laminar) on the line 0<1<...<2k, i.e.
    Σ_{even-series, σ non-crossing} μ(0̂,σ) = (-1)^k C_k,  rest cancels.
This is FALSE: the non-crossing-restricted sum is -1, 2, -6, 25, -132 (k=1..5), NOT Catalan,
and the crossing remainder 0,0,+1,-11,+90 does not vanish.  Hence the relevant 'planar'
structure is the ribbon genus of the WALK-INDUCED Euler map on G_σ, not the laminarity of σ.
"""
import math
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


def is_even_series(blocks, L):
    p2b = {}
    for bi, B in enumerate(blocks):
        for pos in B:
            p2b[pos] = bi
    edges = [(p2b[i], p2b[i + 1]) for i in range(L)]
    nb = len(blocks)
    if any(u == v for u, v in edges):
        return False
    adj = defaultdict(list)
    for u, v in edges:
        adj[u].append(v); adj[v].append(u)
    seen = {0}; st = [0]
    while st:
        x = st.pop()
        for w in adj[x]:
            if w not in seen:
                seen.add(w); st.append(w)
    if len(seen) != nb:
        return False
    Bm = np.zeros((nb, L))
    for ei, (u, v) in enumerate(edges):
        Bm[v, ei] += 1.0; Bm[u, ei] -= 1.0
    _, s, vh = np.linalg.svd(Bm)
    rank = int((s > 1e-9).sum()); m = L - rank
    if m == 0:
        return False
    ns = vh[rank:]
    grp = defaultdict(int)
    for e in range(L):
        v = ns[:, e]
        if np.max(np.abs(v)) < 1e-7:
            return False                       # bridge
        v = v / np.max(np.abs(v))
        for x in v:
            if abs(x) > 1e-7:
                if x < 0: v = -v
                break
        grp[tuple(round(float(x), 6) for x in v)] += 1
    return all(c % 2 == 0 for c in grp.values())


def is_noncrossing(blocks):
    lab = {}
    for bi, B in enumerate(blocks):
        for x in B:
            lab[x] = bi
    pts = sorted(lab)
    n = len(pts)
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                for l in range(k + 1, n):
                    a, b, c, d = pts[i], pts[j], pts[k], pts[l]
                    if lab[a] == lab[c] and lab[b] == lab[d] and lab[a] != lab[b]:
                        return False
    return True


print("=" * 60)
print("NEGATIVE: non-crossing index-partition does NOT localize (★★)")
print(" k  S_noncrossing  S_crossing   (-1)^kC_k   nc==target")
for k in range(1, 6):
    L = 2 * k
    Snc = 0; Scr = 0
    for blocks in set_partitions(range(L + 1)):
        if not is_even_series(blocks, L):
            continue
        w = mu_partition(blocks)
        if is_noncrossing(blocks):
            Snc += w
        else:
            Scr += w
    tgt = (-1) ** k * catalan(k)
    print(f"{k:>2} {Snc:>13} {Scr:>11} {tgt:>11}     {Snc == tgt}")
print("=" * 60)
print("=> the correct 'genus 0' is the ribbon genus of the walk-induced Euler map on G_σ,")
print("   NOT the laminarity of σ.  (Recorded so the next explorer does not retry this.)")
