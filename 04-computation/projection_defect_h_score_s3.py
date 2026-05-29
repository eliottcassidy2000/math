#!/usr/bin/env python3
"""
projection_defect_h_score_s3.py

kind-pasteur-2026-05-29-S3

Follow-up to projection_defect_structured_moves_s2.py.

Question: do endpoint/boundary structured moves look tournament-only because
they travel farther in tournament invariants such as H and score profile?

For selected structured moves, this script records the same four projection
defect outcomes and augments each line with:
- |Delta H| between the two merged tournament classes
- labeled score-vector L1 distance
- sorted score-sequence L1 distance

This is deliberately small-case exact work for n=5,6.
"""

from collections import Counter, defaultdict
from itertools import permutations

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
from projection_defect_structured_moves_s2 import (
    defect_label,
    mask_for,
    pct,
    summarize_counter,
)


def score_vector(A):
    return tuple(sum(row) for row in A)


def l1(xs, ys):
    return sum(abs(x - y) for x, y in zip(xs, ys))


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
    score_by_bits = {}
    sorted_score_by_bits = {}
    even_deg_by_class = {}

    for bits in range(N):
        A = tiling_to_tournament(bits, n, tiles)
        score = score_vector(A)
        score_by_bits[bits] = score
        sorted_score_by_bits[bits] = tuple(sorted(score))

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

    H_by_bits = {bits: H_by_t[t_class[bits]] for bits in range(N)}
    return {
        "n": n,
        "m": m,
        "N": N,
        "tiles": tiles,
        "t_class": t_class,
        "e_class": e_class,
        "H_by_bits": H_by_bits,
        "score_by_bits": score_by_bits,
        "sorted_score_by_bits": sorted_score_by_bits,
        "merged_tournament_classes": len(t_canon_to_id),
        "even_graph_classes": len(e_canon_to_id),
    }


def selected_move_groups(n, tiles):
    groups = []

    def add(label, predicate):
        mask = mask_for(tiles, predicate)
        if mask:
            groups.append((label, [mask]))

    def add_single_tile_family(r):
        masks = [1 << k for k, tile in enumerate(tiles) if tile[0] - tile[1] == r]
        if masks:
            groups.append((f"single_tiles:range={r}", masks))

    add("endpoint_star:v=0", lambda tile: tile[0] == 0 or tile[1] == 0)
    add(f"endpoint_star:v={n-1}", lambda tile: tile[0] == n - 1 or tile[1] == n - 1)
    for v in sorted({n // 2, (n - 1) // 2}):
        add(f"middle_star:v={v}", lambda tile, v=v: tile[0] == v or tile[1] == v)

    for r in range(2, n):
        add(f"range_flip:r={r}", lambda tile, r=r: tile[0] - tile[1] == r)

    for r in (2, 3, n - 1):
        if 2 <= r < n:
            add_single_tile_family(r)

    add("full_complement_tiling", lambda tile: True)
    return groups


def profile_masks(ctx, masks):
    by_outcome = defaultdict(lambda: {"count": 0, "dH": 0, "score": 0, "sorted_score": 0})
    counter = Counter()
    dH_hist = Counter()
    N = ctx["N"]

    for mask in masks:
        for bits in range(N):
            other = bits ^ mask
            if other <= bits:
                continue
            t_changed = ctx["t_class"][bits] != ctx["t_class"][other]
            e_changed = ctx["e_class"][bits] != ctx["e_class"][other]
            outcome = defect_label(t_changed, e_changed)
            dH = abs(ctx["H_by_bits"][bits] - ctx["H_by_bits"][other])
            score_dist = l1(ctx["score_by_bits"][bits], ctx["score_by_bits"][other])
            sorted_score_dist = l1(
                ctx["sorted_score_by_bits"][bits],
                ctx["sorted_score_by_bits"][other],
            )

            counter[outcome] += 1
            dH_hist[(outcome, dH)] += 1
            row = by_outcome[outcome]
            row["count"] += 1
            row["dH"] += dH
            row["score"] += score_dist
            row["sorted_score"] += sorted_score_dist

    return counter, by_outcome, dH_hist


def mean(row, key):
    return row[key] / row["count"] if row["count"] else 0.0


def distance_label(masks):
    distances = sorted({mask.bit_count() for mask in masks})
    if len(distances) == 1:
        return str(distances[0])
    return ",".join(str(d) for d in distances)


def print_profile(label, masks, counter, by_outcome, dH_hist):
    total = sum(counter.values())
    summary = summarize_counter(counter)
    print(
        f"\n{label}  d={distance_label(masks)}  masks={len(masks)}  lines={total}  "
        f"defect={summary['defect']:+.4f}  joint={summary['joint']:.2f}%"
    )
    print("  outcome           pct     mean|dH|  mean scoreL1  mean sortedScoreL1")
    for outcome in ("silent_both", "tournament_only", "even_only", "joint"):
        row = by_outcome[outcome]
        print(
            f"  {outcome:15s} {pct(counter[outcome], total):6.2f} "
            f"{mean(row, 'dH'):10.2f} {mean(row, 'score'):13.2f} "
            f"{mean(row, 'sorted_score'):19.2f}"
        )

    top_dh = dH_hist.most_common(5)
    rendered = ", ".join(f"{outcome}:dH={dH} x{count}" for (outcome, dH), count in top_dh)
    print(f"  common outcome/H bins: {rendered}")


def analyze_n(n):
    ctx = build_context(n)
    print("\n" + "=" * 100)
    print(
        f"n={n}, m={ctx['m']}, tilings={ctx['N']}, "
        f"T-classes={ctx['merged_tournament_classes']}, E-classes={ctx['even_graph_classes']}"
    )
    for label, masks in selected_move_groups(n, ctx["tiles"]):
        counter, by_outcome, dH_hist = profile_masks(ctx, masks)
        print_profile(label, masks, counter, by_outcome, dH_hist)


def main():
    print("PROJECTION DEFECTS VS H/SCORE GRADIENTS")
    print("kind-pasteur-2026-05-29-S3")
    print("Exact structured-move profiles for n=5,6.")
    for n in (5, 6):
        analyze_n(n)


if __name__ == "__main__":
    main()
