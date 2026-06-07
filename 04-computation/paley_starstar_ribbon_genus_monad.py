#!/usr/bin/env python3
"""
paley_starstar_ribbon_genus_monad.py
monad-explorer-2026-06-07 (deep-research, 6th session)

THE FLAGGED NEXT STEP (THM-438 ADDENDUM-3 point 4): the naive laminar/non-crossing
localization of (★★) was REFUTED.  The session flagged the CORRECT object: the
RIBBON GENUS of the WALK-INDUCED Euler map on G_σ.  This script builds that ribbon
map and tests whether (★★) localizes by genus:

   (★★)   S_k = Σ_{even-series σ} μ(0̂,σ) = (−1)^k C_k.

Conjecture to test:  Σ_{even-series σ, walk-genus(σ)=0} μ(0̂,σ) = (−1)^k C_k,
and the genus≥1 part sums to 0.

The walk on positions 0,1,…,2k traverses edges e_i=(block(i−1),block(i)), each edge
ONCE (an Eulerian trail of G_σ).  The trail induces a rotation system (cyclic order
of darts at each vertex = visit order).  Faces = orbits of α∘σ_rot (α = edge flip);
genus g = (2 − V + E − F)/2.  Open trails (block(0)≠block(2k)) are closed with one
phantom wrap edge.

NO number theory.
"""
import math, sys
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
        b = len(B)
        m *= ((-1) ** (b - 1)) * math.factorial(b - 1)
    return m


def catalan(k):
    return math.comb(2 * k, k) // (k + 1)


def edge_flow_lines(edges, nb):
    E = len(edges)
    Bm = np.zeros((nb, E))
    for ei, (u, v) in enumerate(edges):
        Bm[v, ei] += 1; Bm[u, ei] -= 1
    u, s, vh = np.linalg.svd(Bm)
    tol = 1e-9; rank = int((s > tol).sum()); m = E - rank
    if m == 0:
        return [tuple()] * E, 0
    ns = vh[rank:]; lines = []
    for e in range(E):
        v = ns[:, e]
        if np.max(np.abs(v)) < 1e-7:
            lines.append(("ZERO",)); continue
        v = v / np.max(np.abs(v))
        for x in v:
            if abs(x) > 1e-7:
                if x < 0:
                    v = -v
                break
        lines.append(tuple(round(float(x), 6) for x in v))
    return lines, m


def is_even_series(edges, nb):
    adj = defaultdict(list)
    for (u, v) in edges:
        adj[u].append(v); adj[v].append(u)
    seen = {0}; stk = [0]
    while stk:
        x = stk.pop()
        for w in adj[x]:
            if w not in seen:
                seen.add(w); stk.append(w)
    if len(seen) != nb:
        return False
    lines, m = edge_flow_lines(edges, nb)
    if m == 0:
        return False
    if any(l == ("ZERO",) for l in lines):
        return False
    groups = defaultdict(int)
    for l in lines:
        groups[l] += 1
    return all(c % 2 == 0 for c in groups.values())


def build_graph(blocks, L):
    pos2blk = {}
    for bi, B in enumerate(blocks):
        for pos in B:
            pos2blk[pos] = bi
    return [(pos2blk[i], pos2blk[i + 1]) for i in range(L)], len(blocks), pos2blk


def walk_genus(pos2blk, L):
    """Genus of the walk-induced ribbon map.
    Walk: positions 0..L, c_j=pos2blk[j].  Steps i=1..L: edge e_i from c_{i-1} to c_i.
    If open (c_0!=c_L) add a phantom wrap step L+1 from c_L to c_0.
    Darts: each step s has tail dart T_s (at its start vertex) and head dart H_s (at its
    end vertex). alpha: T_s <-> H_s.
    Rotation at vertex v = visit order: the walk, in cyclic time, passes corners; at the
    moment the walk LEAVES v on step s (tail T_s) it just ARRIVED via H_{s-1} (head of
    previous step).  We order darts at v by their step-time; rotation links them cyclically.
    Faces = cycles of perm  next(d) = sigma_rot(alpha(d)).  g=(2-V+E-F)/2.
    Returns (g, V, E, F).
    """
    c = [pos2blk[j] for j in range(L + 1)]
    steps = []  # (start_vertex, end_vertex)
    for i in range(1, L + 1):
        steps.append((c[i - 1], c[i]))
    closed = (c[0] == c[L])
    if not closed:
        steps.append((c[L], c[0]))  # phantom wrap
    N = len(steps)  # number of edges/steps in the closed walk
    # darts: for step index s (0..N-1): tail dart ('T',s) at steps[s][0], head ('H',s) at steps[s][1]
    # alpha swaps ('T',s)<->('H',s)
    def alpha(d):
        typ, s = d
        return ('H', s) if typ == 'T' else ('T', s)
    # rotation: gather darts per vertex with a time key, then cyclic-link in time order.
    # time of leaving on step s (tail) ~ time s (start of step s).
    # time of arriving on step s (head) ~ time s (end of step s) = just before step s+1 leaves.
    # We order corners by the cyclic walk: a "corner" at vertex v is (arrive head H_{s}, leave tail T_{s+1}).
    # Equivalent rotation: at v, sort darts; use the natural interleaving from the closed walk.
    # Build rotation as: for the closed walk, the sequence of darts in traversal is
    #   T_0, H_0, T_1, H_1, ..., T_{N-1}, H_{N-1}, (back to T_0)
    # At each vertex, take its darts in the order they appear in this global cyclic sequence,
    # and make that the rotation (cyclic next).  This is the canonical "thickened walk" ribbon.
    global_seq = []
    for s in range(N):
        global_seq.append(('T', s))
        global_seq.append(('H', s))
    # per-vertex ordered dart list
    vert_darts = defaultdict(list)
    pos_in_seq = {}
    for idx, d in enumerate(global_seq):
        # vertex of dart d:
        typ, s = d
        v = steps[s][0] if typ == 'T' else steps[s][1]
        vert_darts[v].append(d)
        pos_in_seq[d] = idx
    # rotation next: within each vertex's ordered list, cyclic successor
    sigma_next = {}
    for v, dl in vert_darts.items():
        for i, d in enumerate(dl):
            sigma_next[d] = dl[(i + 1) % len(dl)]
    # face permutation: phi(d) = sigma_next(alpha(d))
    def phi(d):
        return sigma_next[alpha(d)]
    # count cycles of phi
    alld = list(pos_in_seq.keys())
    seen = set(); F = 0
    for d in alld:
        if d in seen:
            continue
        F += 1
        x = d
        while x not in seen:
            seen.add(x); x = phi(x)
    V = len(set(c))
    E = N
    chi = V - E + F
    g = (2 - chi) // 2
    return g, V, E, F, closed


def main():
    KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    print("=" * 70)
    print("(★★) genus localization test: walk-induced ribbon genus on G_σ")
    print("  Conjecture: Σ_{even-series, genus 0} μ = (−1)^k C_k ; genus≥1 sums to 0")
    print("=" * 70)
    for k in range(1, KMAX + 1):
        L = 2 * k
        by_genus = defaultdict(int)          # genus -> sum of mu
        by_genus_cnt = defaultdict(int)
        S = 0
        for blocks in set_partitions(range(L + 1)):
            edges, nb, p2b = build_graph(blocks, L)
            if any(u == v for (u, v) in edges):
                continue
            if not is_even_series(edges, nb):
                continue
            mu = mu_partition(blocks)
            g, V, E, F, closed = walk_genus(p2b, L)
            by_genus[g] += mu
            by_genus_cnt[g] += 1
            S += mu
        tgt = (-1) ** k * catalan(k)
        gmax = max(by_genus.keys())
        g0 = by_genus.get(0, 0)
        print(f"\nk={k}: S_k={S}  target (−1)^kC_k={tgt}  {'OK' if S==tgt else 'MISMATCH'}")
        for g in range(0, gmax + 1):
            print(f"    genus {g}: Σμ={by_genus.get(g,0):>6}   (#patterns={by_genus_cnt.get(g,0)})")
        print(f"    >>> genus-0 Σμ = {g0}   target = {tgt}   "
              f"{'LOCALIZES' if g0==tgt else 'does NOT localize'}")
        higher = S - g0
        print(f"    >>> genus≥1 Σμ = {higher}   {'(=0, cancels)' if higher==0 else '(nonzero)'}")
    print("=" * 70)


if __name__ == "__main__":
    main()
