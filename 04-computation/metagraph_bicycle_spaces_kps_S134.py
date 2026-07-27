#!/usr/bin/env python3
r"""
metagraph_bicycle_spaces_kps_S134.py
(kind-pasteur-2026-07-26-S134; companion to THM-2467)

The bicycle spaces of the star-flip GF(2) split -- klein-S399's
"named, never computed" top pick (TOURNAMENT-INVARIANT-ZOO-ATLAS
SS II.e #1).  Two ambient graphs, per THM-1405's conventions:

  K_n with the base path as spanning tree (the CLAUDE.md-era
  Cut (+) Cycle split), and the TILE GRAPH K_n minus the base path
  (THM-1405's incidence structure, dim Cut = n-1,
  dim Cycle = m-n+1, m = C(n-1,2)).

LAWS (this run, exact GF(2) ranks):
  1. dim bicycle(K_n) = (n-2) * [n even], verified n = 4..18.
     Proof sketch: delta(S) has even degree everywhere iff n even
     and |S| even; even-size cuts span a codimension-1 subspace of
     the cut space.  For odd n the split is direct -- confirming and
     quantifying the zoo table's "direct iff n odd" (which is about
     K_n; see law 2 for why it must not be read on the tile graph).
  2. dim bicycle(tile graph) in {0, 1} for all 4 <= n <= 30, and
     = 1  iff  n = 2, 3, 6, 9, 10 (mod 12)   (palindromic residue
     set, closed under n -> 12-n).  A mod-9 fit matching n <= 18
     BREAKS at n = 21 (the MISTAKE-055 motif: five-point fits lie);
     the mod-12 law survives 11 nonzero points and 16 zeros.

Checks are hard failures.
Reproduction: python 04-computation/metagraph_bicycle_spaces_kps_S134.py
"""
import numpy as np
import collections


def fail(msg):
    raise SystemExit("CHECK FAILED: " + msg)


def gf2_rank(M):
    if M.size == 0:
        return 0
    M = M.copy() % 2
    r = 0
    rows, cols = M.shape
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if M[i, c]:
                piv = i
                break
        if piv is None:
            continue
        M[[r, piv]] = M[[piv, r]]
        for i in range(rows):
            if i != r and M[i, c]:
                M[i] ^= M[r]
        r += 1
    return r


def bicycle_dim(n, full):
    if full:
        edges = [(i, j) for i in range(1, n + 1) for j in range(i + 1, n + 1)]
    else:
        edges = [(i, j) for i in range(1, n + 1) for j in range(i + 1, n + 1)
                 if j - i >= 2]
    m = len(edges)
    eidx = {e: k for k, e in enumerate(edges)}
    C = np.zeros((n, m), dtype=np.int8)
    for v in range(1, n + 1):
        for e, k in eidx.items():
            if v in e:
                C[v - 1, k] = 1
    adj = collections.defaultdict(list)
    for (a, b) in edges:
        adj[a].append(b)
        adj[b].append(a)
    parent = {1: None}
    dq = collections.deque([1])
    while dq:
        u = dq.popleft()
        for w in adj[u]:
            if w not in parent:
                parent[w] = u
                dq.append(w)
    tree = set((min(u, w), max(u, w)) for w, u in parent.items()
               if u is not None)

    def ptr(x):
        out = []
        while parent[x] is not None:
            out.append((min(x, parent[x]), max(x, parent[x])))
            x = parent[x]
        return out

    Zr = []
    for e in edges:
        if e in tree:
            continue
        cyc = set(ptr(e[0])) ^ set(ptr(e[1]))
        cyc.add(e)
        row = np.zeros(m, dtype=np.int8)
        for ee in cyc:
            row[eidx[ee]] = 1
        Zr.append(row)
    if not Zr:
        return 0
    Z = np.array(Zr, dtype=np.int8)
    return gf2_rank(C) + gf2_rank(Z) - gf2_rank(np.vstack([C, Z]))


print("n : bicycle(K_n)  bicycle(tile)")
kn_ok = tile_ok = True
tile_dims = {}
for n in range(4, 31):
    bk = bicycle_dim(n, True) if n <= 18 else None
    bt = bicycle_dim(n, False)
    tile_dims[n] = bt
    if bk is not None:
        expect = (n - 2) if n % 2 == 0 else 0
        if bk != expect:
            fail(f"K_n law at n={n}: {bk} != {expect}")
    if bt not in (0, 1):
        fail(f"tile dim > 1 at n={n}")
    expect_t = 1 if n % 12 in (2, 3, 6, 9, 10) else 0
    if bt != expect_t:
        fail(f"tile mod-12 law at n={n}: {bt} != {expect_t}")
    print(f"{n:2d}:      {bk if bk is not None else '-'}            {bt}")

# extended verification for the PROOF upgrade: mod-12 law to n = 60
for n in range(31, 61):
    bt = bicycle_dim(n, False)
    exp = 1 if n % 12 in (2, 3, 6, 9, 10) else 0
    if bt != exp:
        fail(f"mod-12 law extended check at n={n}")
print("mod-12 law verified n = 5..60 (56 points): PASS")

# hostile control: the mod-9 fit that matches through n=18 must break
mod9 = [n for n in range(6, 31) if n % 9 in (0, 1, 5, 6)]
actual = [n for n in range(6, 31) if tile_dims[n] > 0]
if mod9 == actual:
    fail("mod-9 fit unexpectedly survives")
print("hostile control: mod-9 fit breaks (first at n=19/21):",
      sorted(set(mod9) ^ set(actual))[:4])
print("ALL CHECKS PASSED")
