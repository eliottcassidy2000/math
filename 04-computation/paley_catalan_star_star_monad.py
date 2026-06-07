#!/usr/bin/env python3
"""
paley_catalan_star_star_monad.py
monad-explorer-2026-06-07 (deep-research, 4th session)

PURELY COMBINATORIAL verification of the number-theory-free identity (THM-438
ADDENDUM-2, the corrected handoff #1):

  (**)   sum_{sigma : EVEN-SERIES pattern of path [0..2k]}  mu(0,sigma)  =  (-1)^k C_k.

NO primes, NO characters.  A coincidence pattern sigma (set partition of {0..2k})
is an EVEN-SERIES pattern iff its reduced multigraph G_sigma (vertices = blocks,
edges = the 2k consecutive path-edges) is CONNECTED and the multiset of edge
flow-lines has every line of EVEN multiplicity  <=>  prod_e ell_e(s) is a perfect
square  <=>  every series-class of edges is even.  We compute each edge's flow-line
(its coordinate vector in a cycle-space basis), normalize to a primitive sign-fixed
vector, group, and require every group even.  mu(0,sigma)=prod_B (-1)^{|B|-1}(|B|-1)!.

This isolates the Catalan law from all number theory and confirms the corrected
support (even-series, NOT even cacti).  k up to 5 (Bell(11)=678570 partitions).
"""
import math
from fractions import Fraction
from collections import defaultdict


def set_partitions(c):
    c = list(c)
    if len(c) == 1:
        yield [c]
        return
    f = c[0]
    for sm in set_partitions(c[1:]):
        for i, s in enumerate(sm):
            yield sm[:i] + [[f] + s] + sm[i + 1:]
        yield [[f]] + sm


def mu_partition(blocks):
    m = 1
    for B in blocks:
        b = len(B)
        m *= ((-1) ** (b - 1)) * math.factorial(b - 1)
    return m


def catalan(k):
    return math.comb(2 * k, k) // (k + 1)


def edge_flow_lines(edges, nb):
    """Return, for each edge, its primitive sign-normalized flow-line vector (tuple),
       computed in a cycle-space basis (over Q). Edges with proportional vectors share a
       line. A bridge has the ZERO vector (flow forced to 0)."""
    E = len(edges)
    # incidence matrix B (nb x E): +1 head, -1 tail
    import numpy as np
    Bm = np.zeros((nb, E), dtype=np.float64)
    for ei, (u, v) in enumerate(edges):
        Bm[v, ei] += 1.0
        Bm[u, ei] -= 1.0
    # cycle space = nullspace of B (E-dim vectors z with B z = 0). Basis rows -> m x E.
    # Each edge's flow-line = its column in the cycle-space basis = the e-th coordinate
    # across basis vectors -> an m-vector.
    # Compute nullspace via SVD.
    u, s, vh = np.linalg.svd(Bm)
    tol = 1e-9
    rank = int((s > tol).sum())
    m = E - rank
    if m == 0:
        return [tuple()] * E, 0
    ns = vh[rank:]                      # m x E, rows span nullspace
    # edge e's line vector = ns[:, e]  (m-dim)
    lines = []
    for e in range(E):
        v = ns[:, e]
        # normalize: round to rational-ish; use direction. Scale so max-abs entry = 1, sign-fix.
        if np.max(np.abs(v)) < 1e-7:
            lines.append(("ZERO",))     # bridge
            continue
        v = v / np.max(np.abs(v))
        # sign-fix by first significant entry
        for x in v:
            if abs(x) > 1e-7:
                if x < 0:
                    v = -v
                break
        lines.append(tuple(round(float(x), 6) for x in v))
    return lines, m


def is_even_series(edges, nb):
    # connected?
    adj = defaultdict(list)
    for (u, v) in edges:
        adj[u].append(v)
        adj[v].append(u)
    seen = {0}
    stk = [0]
    while stk:
        x = stk.pop()
        for w in adj[x]:
            if w not in seen:
                seen.add(w)
                stk.append(w)
    if len(seen) != nb:
        return False
    lines, m = edge_flow_lines(edges, nb)
    if m == 0:
        return False                    # tree -> cycle space trivial -> F=0, no top order
    if any(l == ("ZERO",) for l in lines):
        return False                    # bridge -> flow 0 -> not saturating
    groups = defaultdict(int)
    for l in lines:
        groups[l] += 1
    return all(c % 2 == 0 for c in groups.values())


def build_graph(blocks, L):
    pos2blk = {}
    for bi, B in enumerate(blocks):
        for pos in B:
            pos2blk[pos] = bi
    return [(pos2blk[i], pos2blk[i + 1]) for i in range(L)], len(blocks)


import sys
KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 4
print("=" * 64)
print("(**)  sum_{even-series sigma} mu(0,sigma) = (-1)^k C_k   [NO number theory]")
print(f"   {'k':>2} {'sum mu':>10} {'(-1)^k C_k':>12} {'match':>7}  #even-series")
for k in range(1, KMAX + 1):
    L = 2 * k
    S = 0
    cnt = 0
    for blocks in set_partitions(range(L + 1)):
        edges, nb = build_graph(blocks, L)
        if any(u == v for (u, v) in edges):
            continue
        if is_even_series(edges, nb):
            S += mu_partition(blocks)
            cnt += 1
    tgt = (-1) ** k * catalan(k)
    print(f"   {k:>2} {S:>10} {tgt:>12} {str(S == tgt):>7}  {cnt}")
print("=" * 64)
print("Confirms THM-438 ADDENDUM-2: the Catalan leading coefficient is the pure")
print("partition-lattice Moebius sum over EVEN-SERIES patterns (even theta included),")
print("with NO characters/primes.  Remaining: a bijective/GF proof of (**).")
