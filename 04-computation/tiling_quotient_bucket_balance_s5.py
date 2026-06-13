#!/usr/bin/env python3
"""
tiling_quotient_bucket_balance_s5.py

kind-pasteur-2026-05-29-S5

Formalize and test the tiling-bucket balance constraints.

Setup:
  Q_m is the tiling hypercube with m=C(n-1,2) non-base-path tiles.
  A quotient map q: Q_m -> B partitions tilings into buckets.
  For any set M of nonzero masks, every tiling has |M| incident lines
  t -- t xor mask.

Bucket balance theorem:
  2*self_b(M) + sum_{c != b} w_{bc}(M) = |q^{-1}(b)| * |M|.

This script checks the theorem for two quotient maps:
  1. merged tournament class Q_m -> G_n/Z_2
  2. even graph class Q_m -> E_n

It also decomposes merged-tournament cross-lines into the class-level
SC-SC spine, SC-NS ribs, and NS-NS sea.
"""

from collections import Counter, defaultdict
from itertools import permutations

from projection_defect_waggly_layers_s1 import (
    complement_tournament,
    even_graph_canon,
    fundamental_cycles,
    graph_degree_sequence,
    hamiltonian_paths,
    merged_tournament_canon,
    tile_pairs,
    tiling_to_even,
    tiling_to_tournament,
    tournament_canon,
)


def is_self_complementary(A, perms):
    return tournament_canon(A, perms) == tournament_canon(complement_tournament(A), perms)


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
    t_info = {}
    e_info = {}

    for bits in range(N):
        A = tiling_to_tournament(bits, n, tiles)

        tc = merged_tournament_canon(A, perms)
        if tc not in t_canon_to_id:
            tid = len(t_canon_to_id)
            t_canon_to_id[tc] = tid
            t_info[tid] = {
                "H": hamiltonian_paths(A),
                "SC": is_self_complementary(A, perms),
            }
        t_class[bits] = t_canon_to_id[tc]

        eg = tiling_to_even(bits, cycles)
        ec = even_graph_canon(eg, n, edges, edge_idx, perms)
        if ec not in e_canon_to_id:
            eid = len(e_canon_to_id)
            e_canon_to_id[ec] = eid
            e_info[eid] = {"degseq": graph_degree_sequence(eg, n, edges)}
        e_class[bits] = e_canon_to_id[ec]

    return {
        "n": n,
        "m": m,
        "N": N,
        "t_class": t_class,
        "e_class": e_class,
        "t_info": t_info,
        "e_info": e_info,
    }


def masks_at_distance(m, d):
    return [mask for mask in range(1, 1 << m) if mask.bit_count() == d]


def selected_mask_families(m):
    families = []

    def add(label, masks):
        key = tuple(masks)
        if masks and all(key != old_key for _, _, old_key in families):
            families.append((label, masks, key))

    add("d=1 wiggly", masks_at_distance(m, 1))
    if 1 < m // 2 < m:
        add(f"d={m // 2} middle", masks_at_distance(m, m // 2))
    if 1 < (m + 1) // 2 < m and (m + 1) // 2 != m // 2:
        add(f"d={(m + 1) // 2} middle", masks_at_distance(m, (m + 1) // 2))
    add("d=m complement-tiling", masks_at_distance(m, m))
    add("all waggly", list(range(1, 1 << m)))
    return [(label, masks) for label, masks, _ in families]


def line_stats(N, q_class, info, masks, categorize=False):
    sizes = Counter(q_class.values())
    self_lines = Counter()
    cross_lines = Counter()
    cross_kind = Counter()
    incident = Counter()
    same_neighbors = [0] * N

    for mask in masks:
        for bits in range(N):
            other = bits ^ mask
            if other <= bits:
                continue
            b1 = q_class[bits]
            b2 = q_class[other]
            if b1 == b2:
                self_lines[b1] += 1
                same_neighbors[bits] += 1
                same_neighbors[other] += 1
                continue

            lo, hi = sorted((b1, b2))
            cross_lines[(lo, hi)] += 1
            incident[b1] += 1
            incident[b2] += 1

            if categorize:
                sc1 = info[b1]["SC"]
                sc2 = info[b2]["SC"]
                if sc1 and sc2:
                    cross_kind["spine_SC_SC"] += 1
                elif sc1 != sc2:
                    cross_kind["ribs_SC_NS"] += 1
                else:
                    cross_kind["sea_NS_NS"] += 1

    degree = len(masks)
    violations = []
    parity_violations = []
    spreads = {}
    escape_rows = []
    for bucket, size in sorted(sizes.items()):
        lhs = 2 * self_lines[bucket] + incident[bucket]
        rhs = size * degree
        if lhs != rhs:
            violations.append((bucket, lhs, rhs))
        if degree == 1 and incident[bucket] % 2 != size % 2:
            parity_violations.append(bucket)

        vals = [same_neighbors[bits] for bits, b in q_class.items() if b == bucket]
        spreads[bucket] = max(vals) - min(vals) if vals else 0
        escape = incident[bucket] / rhs if rhs else 0.0
        neutral = (2 * self_lines[bucket]) / rhs if rhs else 0.0
        escape_rows.append((escape, neutral, bucket, size, self_lines[bucket], incident[bucket]))

    total_lines = sum(self_lines.values()) + sum(cross_lines.values())
    expected_lines = N * degree // 2
    return {
        "sizes": sizes,
        "self_lines": self_lines,
        "cross_lines": cross_lines,
        "cross_kind": cross_kind,
        "incident": incident,
        "degree": degree,
        "total_lines": total_lines,
        "expected_lines": expected_lines,
        "violations": violations,
        "parity_violations": parity_violations,
        "spreads": spreads,
        "escape_rows": escape_rows,
    }


def bucket_label(bucket, info, quotient):
    if quotient == "T":
        typ = "SC" if info[bucket]["SC"] else "NS"
        return f"b={bucket} {typ} H={info[bucket]['H']}"
    return f"b={bucket} degseq={info[bucket]['degseq']}"


def print_stats(ctx, quotient, label, stats, info):
    sizes = stats["sizes"]
    total_buckets = len(sizes)
    odd_buckets = sum(1 for size in sizes.values() if size % 2)
    self_total = sum(stats["self_lines"].values())
    cross_total = sum(stats["cross_lines"].values())
    non_equitable = sum(1 for spread in stats["spreads"].values() if spread)
    max_spread = max(stats["spreads"].values()) if stats["spreads"] else 0

    print(f"\n  [{quotient}] {label}")
    print(
        f"    buckets={total_buckets}, size-range={min(sizes.values())}..{max(sizes.values())}, "
        f"odd-buckets={odd_buckets}"
    )
    print(
        f"    lines={stats['total_lines']} expected={stats['expected_lines']}, "
        f"self={self_total}, cross={cross_total}, violations={len(stats['violations'])}"
    )
    print(
        f"    internal-neighbor non-equitable buckets={non_equitable}/{total_buckets}, "
        f"max spread={max_spread}"
    )

    if stats["degree"] == 1:
        print(
            f"    complement parity violations={len(stats['parity_violations'])} "
            "(incident_cross(b) == bucket_size(b) mod 2)"
        )

    if quotient == "T":
        kinds = stats["cross_kind"]
        print(
            "    cross-line geometry: "
            f"spine={kinds['spine_SC_SC']}, ribs={kinds['ribs_SC_NS']}, sea={kinds['sea_NS_NS']}"
        )

    rows = sorted(stats["escape_rows"])
    print("    most trapped buckets:")
    for escape, neutral, bucket, size, self_l, inc in rows[:3]:
        print(
            f"      {bucket_label(bucket, info, quotient)} size={size} "
            f"escape={escape:.3f} neutral={neutral:.3f} self_lines={self_l} cross_half={inc}"
        )
    print("    most escaping buckets:")
    for escape, neutral, bucket, size, self_l, inc in rows[-3:][::-1]:
        print(
            f"      {bucket_label(bucket, info, quotient)} size={size} "
            f"escape={escape:.3f} neutral={neutral:.3f} self_lines={self_l} cross_half={inc}"
        )


def analyze_n(n):
    ctx = build_context(n)
    print("\n" + "=" * 96)
    print(f"n={n}, m={ctx['m']}, tilings={ctx['N']}")
    print(
        f"merged tournament buckets={len(set(ctx['t_class'].values()))}, "
        f"even graph buckets={len(set(ctx['e_class'].values()))}"
    )

    for label, masks in selected_mask_families(ctx["m"]):
        t_stats = line_stats(ctx["N"], ctx["t_class"], ctx["t_info"], masks, categorize=True)
        e_stats = line_stats(ctx["N"], ctx["e_class"], ctx["e_info"], masks, categorize=False)
        print_stats(ctx, "T", label, t_stats, ctx["t_info"])
        print_stats(ctx, "E", label, e_stats, ctx["e_info"])


def main():
    print("MERGED TILING BUCKET CONSTRAINTS")
    print("kind-pasteur-2026-05-29-S5")
    print("T = merged tournament quotient G_n/Z_2; E = even graph quotient E_n.")
    print("Balance equation: 2*self_b + incident_cross_b = |bucket_b| * |mask_family|.")
    for n in range(3, 7):
        analyze_n(n)


if __name__ == "__main__":
    main()
