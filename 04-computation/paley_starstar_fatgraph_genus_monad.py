#!/usr/bin/env python3
"""
paley_starstar_fatgraph_genus_monad.py
monad-explorer-2026-06-07 (deep-research, 6th session)

PROPER matrix-model genus expansion of (★★).  The previous attempt
(paley_starstar_ribbon_genus_monad.py) used ONE canonical (global-time) rotation per
pattern σ and reproduced exactly the already-REFUTED laminar sum (−1,2,−6,25,−132).
The CORRECT object expands the Möbius weight into individual ribbon maps:

   μ(0̂,σ) = ∏_v (−1)^{|b_v|−1}(|b_v|−1)!
          = Σ_{R : rotation system on σ}  ∏_v (−1)^{|b_v|−1},     (|b_v|−1)! rotations per vertex

so   S_k = Σ_{(σ,R)} sign(σ),   sign(σ)=(−1)^{(2k+1)−V}=(−1)^{F−1}.
Each (σ,R) is a fatgraph (ribbon graph) = a gluing in the two-point matrix model;
assign it its genus g(σ,R) via Euler V−E+F=2−2g.  TEST localization:

      Σ_{(σ,R): genus 0} sign = (−1)^k C_k ,   higher genus cancels to 0.

Construction.  Positions 0..2k are CORNERS.  Add one phantom wrap edge e_0 from
position 2k to position 0 (so every position is a full corner: in-dart + out-dart),
giving E=2k+1 edges.  At vertex v, a rotation = a cyclic order of its corners; the
dart sequence is [in_{p1},out_{p1},in_{p2},out_{p2},…].  Faces F = #cycles of
φ = rot_next ∘ α (α = edge flip).  Total maps = Σ_σ ∏(|b_v|−1)! = Σ_m t(k,m) (tiny).

NO number theory.
"""
import math, sys
from itertools import permutations, product
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


def faces_for_rotation(pos2blk, L, corner_order):
    """Given chosen cyclic corner order per vertex, return F = #faces of the fatgraph.
    Edges: e_1..e_L plus phantom e_0 (wrap pos L -> pos 0).  Indexed 0..L.
    Corner at position p:  in-edge = edge arriving at p, out-edge = edge leaving p.
       pos 0:   in = e_0 (phantom),     out = e_1
       pos p (0<p<L): in = e_p,          out = e_{p+1}
       pos L:   in = e_L,                out = e_0 (phantom)
    Dart = (edge_index, 'h'|'t')  ('h' head dart at the edge's arrival vertex,
                                   't' tail dart at the edge's departure vertex).
    α: (e,'h')<->(e,'t').
    in-dart at p  = (in_edge,  'h');  out-dart at p = (out_edge, 't').
    """
    in_edge = {}
    out_edge = {}
    for p in range(L + 1):
        in_edge[p] = 0 if p == 0 else p           # edge index (e_0 phantom for p=0)
        out_edge[p] = 0 if p == L else p + 1      # e_0 phantom for p=L
    def in_dart(p):  return (in_edge[p], 'h')
    def out_dart(p): return (out_edge[p], 't')
    # rotation: per vertex, dart sequence [in_{p1},out_{p1},in_{p2},out_{p2},...]
    rot_next = {}
    for v, order in corner_order.items():
        seq = []
        for p in order:
            seq.append(in_dart(p)); seq.append(out_dart(p))
        for i, d in enumerate(seq):
            rot_next[d] = seq[(i + 1) % len(seq)]
    def alpha(d):
        e, s = d
        return (e, 't' if s == 'h' else 'h')
    def phi(d):
        return rot_next[alpha(d)]
    alld = list(rot_next.keys())
    seen = set(); F = 0
    for d in alld:
        if d in seen:
            continue
        F += 1
        x = d
        while x not in seen:
            seen.add(x); x = phi(x)
    return F


def main():
    KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    print("=" * 72)
    print("(★★) FATGRAPH genus expansion: Σ over ribbon maps (rotation systems)")
    print("  μ(0̂,σ)=∏(−1)^{|b|−1}(|b|−1)!  expanded over (|b|−1)! rotations/vertex.")
    print("  TEST: Σ_{genus 0} sign = (−1)^k C_k ; higher genus cancels.")
    print("=" * 72)
    for k in range(1, KMAX + 1):
        L = 2 * k
        by_genus = defaultdict(int)
        by_genus_cnt = defaultdict(int)
        total_maps = 0
        S = 0
        for blocks in set_partitions(range(L + 1)):
            edges, nb, p2b = build_graph(blocks, L)
            if any(u == v for (u, v) in edges):
                continue
            if not is_even_series(edges, nb):
                continue
            V = nb
            sign = (-1) ** ((L + 1) - V)         # = (−1)^{(2k+1)−V}, per-rotation sign
            # corners per vertex
            corners = defaultdict(list)
            for p in range(L + 1):
                corners[p2b[p]].append(p)
            verts = list(corners.keys())
            # enumerate cyclic orders per vertex: fix first element, permute rest
            per_vertex_orders = []
            for v in verts:
                cs = corners[v]
                if len(cs) == 1:
                    per_vertex_orders.append([tuple(cs)])
                else:
                    first = cs[0]
                    rest = cs[1:]
                    per_vertex_orders.append([(first,) + perm for perm in permutations(rest)])
            for combo in product(*per_vertex_orders):
                corner_order = {v: combo[i] for i, v in enumerate(verts)}
                E = L + 1                          # incl phantom
                F = faces_for_rotation(p2b, L, corner_order)
                g = (2 - (V - E + F)) // 2
                by_genus[g] += sign
                by_genus_cnt[g] += 1
                total_maps += 1
                S += sign
        tgt = (-1) ** k * catalan(k)
        gmax = max(by_genus.keys()) if by_genus else 0
        gmin = min(by_genus.keys()) if by_genus else 0
        g0 = by_genus.get(0, 0)
        print(f"\nk={k}: S_k={S}  target={tgt}  {'OK' if S==tgt else 'MISMATCH'}   (total maps={total_maps})")
        for g in range(gmin, gmax + 1):
            print(f"    genus {g}: Σsign={by_genus.get(g,0):>6}   (#maps={by_genus_cnt.get(g,0)})")
        print(f"    >>> genus-min Σ = {by_genus.get(gmin,0)} (at g={gmin})")
        print(f"    >>> genus-0   Σ = {g0}   target={tgt}   "
              f"{'LOCALIZES' if g0==tgt else 'NO'}")
    print("=" * 72)


if __name__ == "__main__":
    main()
