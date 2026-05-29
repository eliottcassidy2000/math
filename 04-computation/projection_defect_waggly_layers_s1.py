#!/usr/bin/env python3
"""
projection_defect_waggly_layers_s1.py

kind-pasteur-2026-05-29-S1

Measure how two quotient lenses on the same tiling hypercube disagree:

  Q_m -> merged tournament class G_n/Z_2
  Q_m -> even graph class E_n

Earlier bridge work measured this mostly for d=1 tile flips.  This script
extends the measurement to every waggly Hamming layer d=1,...,m for n<=6.

Definitions follow AGENTS.md:
- Base path is n-1 -> n-2 -> ... -> 0.
- A tile is a nonconsecutive pair hi > lo+1.
- bit 0 means hi -> lo, bit 1 means lo -> hi.
- bit 1 also means XOR the corresponding fundamental cycle into the even graph.
"""

from collections import Counter, defaultdict
from itertools import permutations
from math import comb, factorial, sqrt


def tile_pairs(n):
    """Explorer order: lo increases, hi decreases."""
    return [(hi, lo) for lo in range(n - 2) for hi in range(n - 1, lo + 1, -1)]


def all_edges(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def tiling_to_tournament(bits, n, tiles):
    A = [[0] * n for _ in range(n)]

    # Base path: high vertex beats the next lower vertex.
    for lo in range(n - 1):
        A[lo + 1][lo] = 1

    for k, (hi, lo) in enumerate(tiles):
        if bits & (1 << k):
            A[lo][hi] = 1
        else:
            A[hi][lo] = 1
    return A


def complement_tournament(A):
    n = len(A)
    return [[0 if i == j else 1 - A[i][j] for j in range(n)] for i in range(n)]


def tournament_canon(A, perms):
    n = len(A)
    best = None
    for p in perms:
        key = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n) if i != j)
        if best is None or key < best:
            best = key
    return best


def merged_tournament_canon(A, perms):
    c = tournament_canon(A, perms)
    cc = tournament_canon(complement_tournament(A), perms)
    return min(c, cc)


def fundamental_cycles(n, tiles):
    edges = all_edges(n)
    edge_idx = {e: i for i, e in enumerate(edges)}
    cycles = []
    for hi, lo in tiles:
        lo2, hi2 = min(lo, hi), max(lo, hi)
        mask = 0
        for v in range(lo2, hi2):
            mask ^= 1 << edge_idx[(v, v + 1)]
        mask ^= 1 << edge_idx[(lo2, hi2)]
        cycles.append(mask)
    return edges, edge_idx, cycles


def tiling_to_even(bits, cycles):
    g = 0
    for k, cyc in enumerate(cycles):
        if bits & (1 << k):
            g ^= cyc
    return g


def even_graph_canon(g_bits, n, edges, edge_idx, perms):
    best = None
    for p in perms:
        nb = 0
        for k, (i, j) in enumerate(edges):
            if not (g_bits & (1 << k)):
                continue
            a, b = sorted((p[i], p[j]))
            nb |= 1 << edge_idx[(a, b)]
        if best is None or nb < best:
            best = nb
    return best


def graph_degree_sequence(g_bits, n, edges):
    deg = [0] * n
    for k, (i, j) in enumerate(edges):
        if g_bits & (1 << k):
            deg[i] += 1
            deg[j] += 1
    return tuple(sorted(deg))


def hamiltonian_paths(A):
    n = len(A)
    dp = {(1 << v, v): 1 for v in range(n)}
    full = (1 << n) - 1
    for S in range(1 << n):
        for v in range(n):
            val = dp.get((S, v), 0)
            if not val:
                continue
            for u in range(n):
                if S & (1 << u):
                    continue
                if A[v][u]:
                    key = (S | (1 << u), u)
                    dp[key] = dp.get(key, 0) + val
    return sum(dp.get((full, v), 0) for v in range(n))


def phi(counter):
    p00 = counter["silent_both"]
    p10 = counter["tournament_only"]
    p01 = counter["even_only"]
    p11 = counter["joint"]
    denom = (p10 + p11) * (p00 + p01) * (p01 + p11) * (p00 + p10)
    if denom == 0:
        return 0.0
    return (p11 * p00 - p10 * p01) / sqrt(denom)


def summarize_classes(n, t_class, e_class, H_by_t, even_deg_by_class):
    t_fibers = Counter(t_class.values())
    e_fibers = Counter(e_class.values())
    return {
        "merged_tournament_classes": len(t_fibers),
        "even_graph_classes": len(e_fibers),
        "tournament_fibers": sorted(t_fibers.values()),
        "even_fibers": sorted(e_fibers.values()),
        "H_values": sorted(set(H_by_t.values())),
        "even_degree_sequences": sorted(set(even_deg_by_class.values())),
    }


def analyze_n(n):
    tiles = tile_pairs(n)
    m = len(tiles)
    N = 1 << m
    perms = list(permutations(range(n)))
    edges, edge_idx, cycles = fundamental_cycles(n, tiles)

    t_canon_to_id = {}
    e_canon_to_id = {}
    t_class = {}
    e_class = {}
    H_by_t = {}
    even_deg_by_class = {}

    for bits in range(N):
        A = tiling_to_tournament(bits, n, tiles)
        tc = merged_tournament_canon(A, perms)
        if tc not in t_canon_to_id:
            t_canon_to_id[tc] = len(t_canon_to_id)
            H_by_t[t_canon_to_id[tc]] = hamiltonian_paths(A)
        t_class[bits] = t_canon_to_id[tc]

        eg = tiling_to_even(bits, cycles)
        ec = even_graph_canon(eg, n, edges, edge_idx, perms)
        if ec not in e_canon_to_id:
            e_canon_to_id[ec] = len(e_canon_to_id)
            even_deg_by_class[e_canon_to_id[ec]] = graph_degree_sequence(eg, n, edges)
        e_class[bits] = e_canon_to_id[ec]

    by_d = defaultdict(Counter)
    masks_by_d = defaultdict(list)
    for mask in range(1, N):
        masks_by_d[mask.bit_count()].append(mask)

    for bits in range(N):
        for d, masks in masks_by_d.items():
            c = by_d[d]
            for mask in masks:
                other = bits ^ mask
                if other <= bits:
                    continue
                t_changed = t_class[bits] != t_class[other]
                e_changed = e_class[bits] != e_class[other]
                if t_changed and e_changed:
                    c["joint"] += 1
                elif t_changed:
                    c["tournament_only"] += 1
                elif e_changed:
                    c["even_only"] += 1
                else:
                    c["silent_both"] += 1

    return m, summarize_classes(n, t_class, e_class, H_by_t, even_deg_by_class), by_d


def print_table(n, m, summary, by_d):
    print("\n" + "=" * 88)
    print(f"n={n}, m={m}, tilings=2^{m}={1 << m}")
    print(
        "classes: merged tournament={merged_tournament_classes}, even graph={even_graph_classes}".format(
            **summary
        )
    )
    print(f"H-values on merged tournament classes: {summary['H_values']}")
    print(f"even graph degree sequences: {summary['even_degree_sequences']}")
    print("\nd  lines    silent%   T-only%   E-only%   joint%    union-joint%  phi     T/E defect")
    print("-- -------  --------  --------  --------  --------  ------------  ------  ----------")
    for d in range(1, m + 1):
        c = by_d[d]
        total = sum(c.values())
        silent = 100 * c["silent_both"] / total
        t_only = 100 * c["tournament_only"] / total
        e_only = 100 * c["even_only"] / total
        joint = 100 * c["joint"] / total
        union = c["joint"] + c["tournament_only"] + c["even_only"]
        jacc = 100 * c["joint"] / union if union else 0.0
        defect = (c["tournament_only"] - c["even_only"]) / total
        print(
            f"{d:2d} {total:7d}  {silent:8.2f}  {t_only:8.2f}  {e_only:8.2f}  "
            f"{joint:8.2f}  {jacc:12.2f}  {phi(c):6.3f}  {defect:+10.4f}"
        )

    cumulative = Counter()
    for d in range(1, m + 1):
        cumulative.update(by_d[d])
    total = sum(cumulative.values())
    union = cumulative["joint"] + cumulative["tournament_only"] + cumulative["even_only"]
    print("\nall layers:")
    print(f"  total waggly lines: {total} = C(2^{m}, 2)")
    print(f"  silent in both quotients: {cumulative['silent_both']} ({100*cumulative['silent_both']/total:.2f}%)")
    print(f"  tournament-only defects: {cumulative['tournament_only']} ({100*cumulative['tournament_only']/total:.2f}%)")
    print(f"  even-only defects:       {cumulative['even_only']} ({100*cumulative['even_only']/total:.2f}%)")
    print(f"  joint changes:           {cumulative['joint']} ({100*cumulative['joint']/total:.2f}%)")
    print(f"  joint share of non-silent union: {100*cumulative['joint']/union if union else 0:.2f}%")


def main():
    print("PROJECTION DEFECTS ACROSS WAGGLY LAYERS")
    print("kind-pasteur-2026-05-29-S1")
    print("T = merged tournament class, E = even graph class")
    print("T-only/E-only are the two quotient-commutator defects.")
    for n in range(3, 7):
        m, summary, by_d = analyze_n(n)
        print_table(n, m, summary, by_d)


if __name__ == "__main__":
    main()
