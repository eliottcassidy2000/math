#!/usr/bin/env python3
"""
goodcut_projection_defect_s14.py

opus-2026-05-29-S14

Cross S13's good-cut bucket coordinate with kind-pasteur S1/S2's
projection-defect lens:

  Q_m -> G_n/Z_2       merged tournament quotient
  Q_m -> E_n           even-graph quotient
  Q_m -> {0,2,...}     good-cut bucket g

Exact for n=3..6.  The goal is to see whether the tournament-only vs
even-only defect polarity is visible in the interval-cover coordinate.
"""

from collections import Counter, defaultdict
from itertools import permutations

from projection_defect_structured_moves_s2 import structured_masks
from projection_defect_waggly_layers_s1 import (
    even_graph_canon,
    fundamental_cycles,
    graph_degree_sequence,
    hamiltonian_paths,
    merged_tournament_canon,
    tile_pairs,
    tiling_to_even,
    tiling_to_tournament,
)


def good_cut_count(bits, tiles):
    cuts = set()
    for k, (hi, lo) in enumerate(tiles):
        if bits & (1 << k):
            cuts.update(range(lo + 1, hi + 1))
    return len(cuts)


def defect_label(t_changed, e_changed):
    if t_changed and e_changed:
        return "joint"
    if t_changed:
        return "tournament_only"
    if e_changed:
        return "even_only"
    return "silent_both"


def pct(x, total):
    return 100 * x / total if total else 0.0


def defect(counter):
    total = sum(v for k, v in counter.items() if k in LABELS)
    return (counter["tournament_only"] - counter["even_only"]) / total if total else 0.0


LABELS = ("silent_both", "tournament_only", "even_only", "joint")


def build_context(n):
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

    g = {bits: good_cut_count(bits, tiles) for bits in range(N)}
    return {
        "n": n,
        "m": m,
        "N": N,
        "tiles": tiles,
        "t_class": t_class,
        "e_class": e_class,
        "g": g,
        "t_classes": len(t_canon_to_id),
        "e_classes": len(e_canon_to_id),
        "H_values": sorted(set(H_by_t.values())),
        "even_degrees": sorted(set(even_deg_by_class.values())),
    }


def add_line(counter, ctx, bits, other):
    t_changed = ctx["t_class"][bits] != ctx["t_class"][other]
    e_changed = ctx["e_class"][bits] != ctx["e_class"][other]
    label = defect_label(t_changed, e_changed)
    dg = abs(ctx["g"][other] - ctx["g"][bits])
    counter[label] += 1
    counter["abs_dg_sum"] += dg
    counter["max_abs_dg"] = max(counter["max_abs_dg"], dg)


def line_total(counter):
    return sum(counter[label] for label in LABELS)


def print_counter_row(name, counter):
    total = line_total(counter)
    mean_dg = counter["abs_dg_sum"] / total if total else 0.0
    print(
        f"{name:42s} {total:7d} {pct(counter['silent_both'], total):8.2f} "
        f"{pct(counter['tournament_only'], total):8.2f} "
        f"{pct(counter['even_only'], total):8.2f} "
        f"{pct(counter['joint'], total):8.2f} "
        f"{defect(counter):+9.4f} {mean_dg:8.3f} {counter['max_abs_dg']:4d}"
    )


def analyze_single_tiles(ctx):
    by_range = defaultdict(Counter)
    by_abs_dg = defaultdict(Counter)
    tiles = ctx["tiles"]

    for bits in range(ctx["N"]):
        for k, (hi, lo) in enumerate(tiles):
            other = bits ^ (1 << k)
            if other <= bits:
                continue
            range_key = hi - lo
            dg = abs(ctx["g"][other] - ctx["g"][bits])
            add_line(by_range[range_key], ctx, bits, other)
            add_line(by_abs_dg[dg], ctx, bits, other)

    print("\nsingle-tile lines by tile range")
    print(HEADER)
    for r in sorted(by_range):
        print_counter_row(f"range={r}", by_range[r])

    print("\nsingle-tile lines by |Delta g|")
    print(HEADER)
    for dg in sorted(by_abs_dg):
        print_counter_row(f"|Delta g|={dg}", by_abs_dg[dg])


def analyze_structured_moves(ctx):
    individual = []
    for label, mask in structured_masks(ctx["n"], ctx["tiles"]):
        counter = Counter()
        for bits in range(ctx["N"]):
            other = bits ^ mask
            if other <= bits:
                continue
            add_line(counter, ctx, bits, other)
        individual.append((label, mask.bit_count(), counter))

    print("\nstructured moves with largest tournament-only polarity")
    print(HEADER)
    for label, distance, counter in sorted(individual, key=lambda x: defect(x[2]), reverse=True)[:8]:
        print_counter_row(f"d={distance}:{label}", counter)

    print("\nstructured moves with largest even-only polarity")
    print(HEADER)
    for label, distance, counter in sorted(individual, key=lambda x: defect(x[2]))[:8]:
        print_counter_row(f"d={distance}:{label}", counter)


HEADER = (
    "bucket/profile                            lines  silent%  T-only%  E-only%   joint%    T-Edef  mean|dg| max"
)


def analyze_n(n):
    ctx = build_context(n)
    print("\n" + "=" * 100)
    print(
        f"n={n}, m={ctx['m']}, tilings={ctx['N']}, "
        f"T-classes={ctx['t_classes']}, E-classes={ctx['e_classes']}"
    )
    print(f"H-values: {ctx['H_values']}")
    print(f"even degree sequences: {ctx['even_degrees']}")
    analyze_single_tiles(ctx)
    analyze_structured_moves(ctx)


def main():
    print("GOOD-CUT BUCKETS VS PROJECTION DEFECTS")
    print("opus-2026-05-29-S14")
    print("T-Edef = (tournament_only - even_only) / lines; dg = change in good-cut bucket.")
    for n in range(3, 7):
        analyze_n(n)


if __name__ == "__main__":
    main()
