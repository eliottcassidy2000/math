#!/usr/bin/env python3
"""
projection_defect_structured_moves_s2.py

kind-pasteur-2026-05-29-S2

Follow-up to projection_defect_waggly_layers_s1.py.

The S1 scan averaged over whole Hamming layers. This script asks whether
geometrically meaningful move families have distinct projection-defect
signatures under the two quotient maps:

  Q_m -> merged tournament class G_n/Z_2
  Q_m -> even graph class E_n

Structured move families:
- single tile, grouped by range hi-lo
- range flip: all tiles with fixed hi-lo
- upper strip: all tiles with fixed upper vertex hi
- lower strip: all tiles with fixed lower vertex lo
- vertex star: all tiles incident to a fixed vertex
- full complement tiling: all tiles
"""

from collections import Counter, defaultdict
from itertools import permutations

from projection_defect_waggly_layers_s1 import (
    even_graph_canon,
    fundamental_cycles,
    graph_degree_sequence,
    hamiltonian_paths,
    merged_tournament_canon,
    phi,
    tile_pairs,
    tiling_to_even,
    tiling_to_tournament,
)


def pct(x, total):
    return 100 * x / total if total else 0.0


def defect_label(t_changed, e_changed):
    if t_changed and e_changed:
        return "joint"
    if t_changed:
        return "tournament_only"
    if e_changed:
        return "even_only"
    return "silent_both"


def build_classes(n):
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

    return {
        "n": n,
        "m": m,
        "N": N,
        "tiles": tiles,
        "t_class": t_class,
        "e_class": e_class,
        "H_values": sorted(set(H_by_t.values())),
        "even_degree_sequences": sorted(set(even_deg_by_class.values())),
        "merged_tournament_classes": len(t_canon_to_id),
        "even_graph_classes": len(e_canon_to_id),
    }


def mask_for(tiles, predicate):
    mask = 0
    for k, tile in enumerate(tiles):
        if predicate(tile):
            mask |= 1 << k
    return mask


def structured_masks(n, tiles):
    m = len(tiles)
    masks = []

    for k, (hi, lo) in enumerate(tiles):
        masks.append((f"single_tile:range={hi-lo}:tile=({hi},{lo})", 1 << k))

    for r in range(2, n):
        mask = mask_for(tiles, lambda tile, r=r: tile[0] - tile[1] == r)
        if mask:
            masks.append((f"range_flip:r={r}", mask))

    for hi in range(n):
        mask = mask_for(tiles, lambda tile, hi=hi: tile[0] == hi)
        if mask:
            masks.append((f"upper_strip:hi={hi}", mask))

    for lo in range(n):
        mask = mask_for(tiles, lambda tile, lo=lo: tile[1] == lo)
        if mask:
            masks.append((f"lower_strip:lo={lo}", mask))

    for v in range(n):
        mask = mask_for(tiles, lambda tile, v=v: tile[0] == v or tile[1] == v)
        if mask:
            masks.append((f"vertex_star:v={v}", mask))

    masks.append(("full_complement_tiling", (1 << m) - 1))

    # Different labels can intentionally produce the same mask at small n; keep
    # each geometric interpretation, but remove exact duplicate labels.
    seen = set()
    unique = []
    for label, mask in masks:
        key = (label, mask)
        if mask and key not in seen:
            seen.add(key)
            unique.append((label, mask))
    return unique


def profile_for_mask(ctx, mask):
    counter = Counter()
    N = ctx["N"]
    t_class = ctx["t_class"]
    e_class = ctx["e_class"]

    for bits in range(N):
        other = bits ^ mask
        if other <= bits:
            continue
        t_changed = t_class[bits] != t_class[other]
        e_changed = e_class[bits] != e_class[other]
        counter[defect_label(t_changed, e_changed)] += 1
    return counter


def family_key(label):
    return label.split(":", 1)[0]


def summarize_counter(counter):
    total = sum(counter.values())
    union = counter["joint"] + counter["tournament_only"] + counter["even_only"]
    return {
        "total": total,
        "silent": pct(counter["silent_both"], total),
        "t_only": pct(counter["tournament_only"], total),
        "e_only": pct(counter["even_only"], total),
        "joint": pct(counter["joint"], total),
        "union_joint": pct(counter["joint"], union),
        "phi": phi(counter),
        "defect": (counter["tournament_only"] - counter["even_only"]) / total if total else 0.0,
    }


def print_row(name, counter):
    summary = summarize_counter(counter)
    print(
        f"{name:28s} {summary['total']:7d} {summary['silent']:8.2f} "
        f"{summary['t_only']:8.2f} {summary['e_only']:8.2f} "
        f"{summary['joint']:8.2f} {summary['union_joint']:8.2f} "
        f"{summary['phi']:7.3f} {summary['defect']:+9.4f}"
    )


def analyze_n(n):
    ctx = build_classes(n)
    print("\n" + "=" * 100)
    print(
        f"n={n}, m={ctx['m']}, tilings={ctx['N']}, "
        f"T-classes={ctx['merged_tournament_classes']}, "
        f"E-classes={ctx['even_graph_classes']}"
    )
    print(f"H-values: {ctx['H_values']}")
    print(f"even degree sequences: {ctx['even_degree_sequences']}")
    print(
        "\nfamily/move                    lines  silent%  T-only%  E-only%   joint%  J/union     phi    T-Edef"
    )
    print("-" * 100)

    family_totals = defaultdict(Counter)
    individual = []

    for label, mask in structured_masks(n, ctx["tiles"]):
        counter = profile_for_mask(ctx, mask)
        family_totals[family_key(label)].update(counter)
        individual.append((label, mask.bit_count(), counter))

    for family in sorted(family_totals):
        print_row(family, family_totals[family])

    print("\nmost tournament-only biased individual moves:")
    ranked = sorted(
        individual,
        key=lambda item: summarize_counter(item[2])["defect"],
        reverse=True,
    )
    for label, distance, counter in ranked[:8]:
        summary = summarize_counter(counter)
        print(
            f"  d={distance:2d} {label:34s} "
            f"defect={summary['defect']:+.4f} joint={summary['joint']:.2f}%"
        )

    print("\nmost even-only biased individual moves:")
    for label, distance, counter in ranked[-8:]:
        summary = summarize_counter(counter)
        print(
            f"  d={distance:2d} {label:34s} "
            f"defect={summary['defect']:+.4f} joint={summary['joint']:.2f}%"
        )


def main():
    print("STRUCTURED PROJECTION-DEFECT PROFILES")
    print("kind-pasteur-2026-05-29-S2")
    print("Defect = (tournament_only - even_only) / lines.")
    for n in range(3, 7):
        analyze_n(n)


if __name__ == "__main__":
    main()
