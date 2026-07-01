#!/usr/bin/env python3
"""
mac-mini-2026-07-01-S82 -- THE EVEN-GRAPH = TOURNAMENT SWITCHING-CLASS frame + the fibration G_n -> E_n.

Prior repo work (tournaments-and-even-graphs.md) established: tournament = score (cut space, dim n-1) (+)
even graph (cycle space, dim C(n-1,2)); labeled even graphs = 2^{C(n-1,2)}; but the ISO equinumerosity
A000568 = A002854 FAILS (4!=3, 12!=7) and was flagged "needs checking against the paper."

NOVEL RESOLUTION: the even graph of a tournament = its SWITCHING CLASS. Switching a vertex subset S =
reversing every arc across the cut delta(S) (flip all arcs at v = the vertex star). The switch group = the
cut space = (Z/2)^{n-1}. So:
    even graph  <->  tournament mod cuts (switching)  <->  cycle-space element.
Hence A002854(n) = number of tournament SWITCHING CLASSES up to iso, and there is a FIBRATION
    G_n (tournaments, A000568)  -->  E_n (switching classes = even graphs, A002854),
whose fibers are the iso classes of tournaments sharing an even graph (the score-labelings of one structure).

Verify A000568 & A002854 (n<=6) via canon(iso) and canon(iso+switch); tabulate the fibration + fiber sizes.
"""
import itertools
from math import comb
from collections import defaultdict, Counter

def setup(n):
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    idx = {p: k for k, p in enumerate(pairs)}
    perms = list(itertools.permutations(range(n)))
    def relabel(mask, perm):
        out = 0
        for (i, j), k in idx.items():
            bit = (mask >> k) & 1
            pi, pj = perm[i], perm[j]
            a, b = (pi, pj) if pi < pj else (pj, pi)
            out |= (bit if pi < pj else 1 - bit) << idx[(a, b)]
        return out
    # cut(S) bitmask: edge (i,j) flipped iff exactly one of i,j in S
    cuts = []
    for S in range(1 << n):
        if S & 1: continue                       # fix vertex 0 not in S (S and complement same cut)
        cm = 0
        for (i, j), k in idx.items():
            if ((S >> i) & 1) ^ ((S >> j) & 1): cm |= (1 << k)
        cuts.append(cm)
    return pairs, idx, len(pairs), perms, relabel, cuts

A000568 = {3:2, 4:4, 5:12, 6:56, 7:456}
A002854 = {3:2, 4:3, 5:7, 6:16, 7:54}

print(f"{'n':>2} {'A000568(iso)':>12} {'switch-classes':>14} {'A002854':>8} {'match':>6}")
data = {}
for n in range(3, 7):
    pairs, idx, P, perms, relabel, cuts = setup(n)
    iso_canon = {}; switch_canon = {}
    for mask in range(1 << P):
        ic = min(relabel(mask, p) for p in perms)
        iso_canon[mask] = ic
        # switching canonical: min over relabel AND cut-switch
        sc = min(relabel(mask ^ cm, p) for cm in cuts for p in perms)
        switch_canon[mask] = sc
    niso = len(set(iso_canon.values()))
    nsw = len(set(switch_canon.values()))
    ok = (niso == A000568[n]) and (nsw == A002854[n])
    print(f"{n:>2} {niso:>12} {nsw:>14} {A002854[n]:>8} {'OK' if ok else 'FAIL':>6}")
    # fibration: iso class -> switch class
    fib = defaultdict(set)
    reps = {}
    for mask in range(1 << P):
        fib[switch_canon[mask]].add(iso_canon[mask])
    fibersizes = sorted((len(v) for v in fib.values()), reverse=True)
    data[n] = (niso, nsw, fibersizes)
    print(f"     fibration G_{n} -> E_{n}: fiber sizes (# tournament iso classes per even graph) = {fibersizes}")
    print(f"     sum of fibers = {sum(fibersizes)} = A000568({n})? {sum(fibersizes)==A000568[n]}   "
          f"(#even graphs = {len(fibersizes)} = A002854)")

print("\nSEQUENCES (n=3..7):  A000568 tournaments =", [A000568[k] for k in range(3,8)],
      "  A002854 even graphs =", [A002854[k] for k in range(3,8)])
print("The map is NOT a bijection (fibers vary); A000568 = sum over even graphs of (#score-labelings up to iso).")
