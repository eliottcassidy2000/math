#!/usr/bin/env python3
"""
endpoint_transfer_witnesses_s95.py

Mechanism probe for endpoint-transfer GF(2) rank.

The previous scripts established:
  - tournament and complement-merged endpoint transfers have full row rank
  - even-graph endpoint transfer has rank defects

This script asks what kind of rank phenomenon this is:
  - Do rows have private child columns?
  - Does the odd-support bipartite graph have a full matching?
  - Are even-graph defects already visible as support/Hall failures, or only
    as algebraic GF(2) cancellations?
  - What do the first even-graph row dependencies look like?
"""

import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(__file__))

import endpoint_transfer_bucket_recursion_s95 as tourn
import even_graph_endpoint_transfer_s95 as even


def row_bits(rows):
    masks = []
    for row in rows:
        mask = 0
        for j, weight in row.items():
            if weight % 2:
                mask |= 1 << j
        masks.append(mask)
    return masks


def popcount(x):
    return x.bit_count()


def gf2_rank_and_dependencies(masks):
    basis = {}
    dependencies = []
    pivot_cols = []
    pivot_combos = []

    for i, row in enumerate(masks):
        x = row
        combo = 1 << i
        while x:
            pivot = x.bit_length() - 1
            if pivot in basis:
                bx, bc = basis[pivot]
                x ^= bx
                combo ^= bc
            else:
                basis[pivot] = (x, combo)
                pivot_cols.append(pivot)
                pivot_combos.append(combo)
                break
        if x == 0:
            dependencies.append(combo)

    return {
        "rank": len(basis),
        "pivot_cols": pivot_cols,
        "pivot_combos": pivot_combos,
        "dependencies": dependencies,
    }


def max_matching(masks, ncols):
    col_match = [-1] * ncols

    def dfs(i, seen):
        x = masks[i]
        while x:
            lsb = x & -x
            j = lsb.bit_length() - 1
            x ^= lsb
            if seen[j]:
                continue
            seen[j] = True
            if col_match[j] == -1 or dfs(col_match[j], seen):
                col_match[j] = i
                return True
        return False

    matched = 0
    for i in sorted(range(len(masks)), key=lambda r: popcount(masks[r])):
        if dfs(i, [False] * ncols):
            matched += 1
    return matched


def private_column_stats(masks, ncols):
    supports = [0] * ncols
    owner = [-1] * ncols
    for i, mask in enumerate(masks):
        x = mask
        while x:
            lsb = x & -x
            j = lsb.bit_length() - 1
            x ^= lsb
            supports[j] += 1
            owner[j] = i

    private_by_row = Counter()
    for j, support in enumerate(supports):
        if support == 1:
            private_by_row[owner[j]] += 1

    return {
        "private_cols": sum(1 for s in supports if s == 1),
        "rows_with_private": len(private_by_row),
        "max_private_for_row": max(private_by_row.values()) if private_by_row else 0,
        "column_support_spectrum": Counter(supports),
        "private_by_row": private_by_row,
    }


def combo_indices(combo):
    rows = []
    i = 0
    while combo:
        if combo & 1:
            rows.append(i)
        combo >>= 1
        i += 1
    return rows


def tournament_row_label(level, idx, merged=False):
    if merged:
        node = level["merged"][idx]
        return f"{idx}:{node['type']}:M={node['M']}:members={node['members']}"
    can = level["cans"][idx]
    info = level["classes"][can]
    return f"{idx}:H={info['H']}:F={info['F']}:aut={info['aut']}:SC={info['sc']}"


def even_row_label(level, idx):
    info = level["classes"][idx]
    return f"{idx}:fiber={info['fiber']}:edges={info['edges']}:deg={info['degree_seq']}"


def summarize_matrix(name, parent_level, child_count, rows, label_fn):
    masks = row_bits(rows)
    r = gf2_rank_and_dependencies(masks)
    matching = max_matching(masks, child_count)
    priv = private_column_stats(masks, child_count)
    row_weights = [popcount(m) for m in masks]

    print(f"  {name}: rows={len(rows)}, cols={child_count}")
    print(
        "    rank/matching/private: "
        f"rank={r['rank']}, matching={matching}, "
        f"private_cols={priv['private_cols']}, rows_with_private={priv['rows_with_private']}"
    )
    print(
        "    odd row weights: "
        f"min={min(row_weights) if row_weights else 0}, "
        f"max={max(row_weights) if row_weights else 0}, "
        f"avg={sum(row_weights)/len(row_weights) if row_weights else 0:.2f}, "
        f"spectrum={dict(sorted(Counter(row_weights).items()))}"
    )
    print(
        "    child odd-column support spectrum: "
        f"{dict(sorted(priv['column_support_spectrum'].items()))}"
    )

    if r["dependencies"]:
        print(f"    dependencies={len(r['dependencies'])}; first dependency rows:")
        for combo in r["dependencies"][:3]:
            inds = combo_indices(combo)
            print(f"      size={len(inds)} -> {[label_fn(parent_level, i) for i in inds[:8]]}")
    else:
        pivot_sample = r["pivot_cols"][:10]
        print(f"    no dependencies; first pivot child columns={pivot_sample}")

    return {
        "rank": r["rank"],
        "matching": matching,
        "private_cols": priv["private_cols"],
        "rows_with_private": priv["rows_with_private"],
        "dependencies": len(r["dependencies"]),
    }


def main():
    print("=" * 78)
    print("ENDPOINT TRANSFER WITNESSES S95")
    print("=" * 78)
    print("Odd-support, matching, private-column, and dependency probes.")

    t_levels = {n: tourn.build_level(n) for n in range(2, 8)}
    e_levels = {n: even.build_even_level(n) for n in range(2, 8)}

    sequences = {
        "t_rank": [],
        "t_matching": [],
        "t_private_rows": [],
        "m_rank": [],
        "m_matching": [],
        "m_private_rows": [],
        "e_rank": [],
        "e_matching": [],
        "e_private_rows": [],
        "e_dependencies": [],
    }

    for n in range(2, 7):
        print("\n" + "-" * 78)
        print(f"n={n}->{n+1}")
        td = tourn.analyze_transfer(t_levels[n], t_levels[n + 1])
        ed = even.analyze_transfer(e_levels[n], e_levels[n + 1])

        ts = summarize_matrix(
            "tournament classes",
            t_levels[n],
            len(t_levels[n + 1]["cans"]),
            td["class_rows"],
            lambda level, i: tournament_row_label(level, i, merged=False),
        )
        ms = summarize_matrix(
            "merged tournament nodes",
            t_levels[n],
            len(t_levels[n + 1]["merged"]),
            td["merged_rows"],
            lambda level, i: tournament_row_label(level, i, merged=True),
        )
        es = summarize_matrix(
            "even-graph classes",
            e_levels[n],
            len(e_levels[n + 1]["classes"]),
            ed["rows"],
            even_row_label,
        )

        sequences["t_rank"].append(ts["rank"])
        sequences["t_matching"].append(ts["matching"])
        sequences["t_private_rows"].append(ts["rows_with_private"])
        sequences["m_rank"].append(ms["rank"])
        sequences["m_matching"].append(ms["matching"])
        sequences["m_private_rows"].append(ms["rows_with_private"])
        sequences["e_rank"].append(es["rank"])
        sequences["e_matching"].append(es["matching"])
        sequences["e_private_rows"].append(es["rows_with_private"])
        sequences["e_dependencies"].append(es["dependencies"])

    print("\n" + "=" * 78)
    print("SEQUENCES")
    print("=" * 78)
    for key, seq in sequences.items():
        print(f"  {key}: {seq}")

    print("\n" + "=" * 78)
    print("INTERPRETATION")
    print("=" * 78)
    print("""
1. A full odd-support matching is weaker than GF(2) independence.
   If the even-graph quotient has full matching but deficient rank at a level,
   then Hall support is not the obstruction; the obstruction is algebraic
   cancellation among rows.

2. Private child columns are a candidate proof mechanism.
   If every tournament parent row had a private odd child column after a
   suitable ordering, full rank would follow from triangularity. The raw
   private-column counts show how close that direct mechanism is.

3. Even-graph dependencies are explicit collapsed memories.
   The dependency rows printed above are parent even-graph classes whose
   endpoint-extension parity shadows cancel. These are the first concrete
   witnesses for what the cycle-space quotient forgets.
""")


if __name__ == "__main__":
    main()
