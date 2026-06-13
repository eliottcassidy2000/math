#!/usr/bin/env python3
"""
endpoint_collision_geometry_s95.py

Geometric probe for the merged endpoint-transfer SC collision columns.

The earlier collision script found that every non-private self-complementary
child column through 6->7 has odd support exactly 3.  This script asks whether
those three parent nodes form a recognizable object in the parent merged
arc-flip metagraph: a clique, path, independent triple, or mixed residual.
"""

from collections import Counter, defaultdict

import endpoint_private_goodcut_s95 as epg
import endpoint_transfer_bucket_recursion_s95 as et


def odd_column_owners(rows, ncols):
    owners = [[] for _ in range(ncols)]
    for i, row in enumerate(rows):
        for j, weight in row.items():
            if weight % 2:
                owners[j].append(i)
    return owners


def profile_bucket(profile):
    if len(profile) != 1:
        return "mixed"
    return next(iter(profile))


def flip_arc(A, i, j):
    B = [row[:] for row in A]
    B[i][j], B[j][i] = B[j][i], B[i][j]
    return B


def build_full_merged_metagraph(level):
    """Merged arc-reversal graph on tournament isomorphism classes."""
    n = level["n"]
    can_to_cid = {can: cid for cid, can in enumerate(level["cans"])}
    class_edges = set()
    merged_edges = set()
    edge_color = {}

    for cid, can in enumerate(level["cans"]):
        A = level["classes"][can]["rep"]
        src_mid = level["cid_to_mid"][cid]
        for i in range(n):
            for j in range(i + 1, n):
                dst_can = et.canonical(flip_arc(A, i, j))
                dst_cid = can_to_cid[dst_can]
                if cid != dst_cid:
                    class_edges.add(tuple(sorted((cid, dst_cid))))
                dst_mid = level["cid_to_mid"][dst_cid]
                if src_mid != dst_mid:
                    edge = tuple(sorted((src_mid, dst_mid)))
                    merged_edges.add(edge)
                    kinds = tuple(
                        sorted(
                            (
                                level["merged"][src_mid]["type"],
                                level["merged"][dst_mid]["type"],
                            )
                        )
                    )
                    color = "-".join(kinds)
                    if edge in edge_color and edge_color[edge] != color:
                        edge_color[edge] = "mixed"
                    else:
                        edge_color[edge] = color

    adjacency = defaultdict(set)
    for a, b in merged_edges:
        adjacency[a].add(b)
        adjacency[b].add(a)

    return {
        "class_edges": class_edges,
        "merged_edges": merged_edges,
        "adjacency": adjacency,
        "edge_color": edge_color,
    }


def h_value(level, mid):
    node = level["merged"][mid]
    return tuple(
        level["classes"][level["cans"][cid]]["H"] for cid in node["members"]
    )


def induced_edge_data(graph, owners):
    edges = []
    colors = []
    for a_i in range(len(owners)):
        for b_i in range(a_i + 1, len(owners)):
            a = owners[a_i]
            b = owners[b_i]
            edge = tuple(sorted((a, b)))
            if b in graph["adjacency"][a]:
                edges.append(edge)
                colors.append(graph["edge_color"].get(edge, "?"))
    return edges, colors


def peel_hyperedges(triples):
    """Greedy leaf peeling for a 3-uniform incidence hypergraph."""
    remaining = set(range(len(triples)))
    row_to_cols = defaultdict(set)
    for idx, triple in enumerate(triples):
        for row in triple:
            row_to_cols[row].add(idx)

    order = []
    while remaining:
        leaf = None
        for col in sorted(remaining):
            unique_rows = [
                row for row in triples[col]
                if len(row_to_cols[row] & remaining) == 1
            ]
            if unique_rows:
                leaf = (col, unique_rows[0])
                break
        if leaf is None:
            break
        col, row = leaf
        order.append((col, row))
        remaining.remove(col)

    return order, remaining


def summarize_transition(n, parent, child):
    transfer = et.analyze_transfer(parent, child)
    rows = transfer["merged_rows"]
    owners = odd_column_owners(rows, len(child["merged"]))
    parent_g = epg.class_goodcut_profiles(parent, merged=True)
    child_g = epg.class_goodcut_profiles(child, merged=True)
    parent_graph = build_full_merged_metagraph(parent)
    child_top_bucket = child["n"] - 1

    edge_count_spectrum = Counter()
    edge_color_spectrum = Counter()
    owner_g_tuples = Counter()
    owner_kind_tuples = Counter()
    owner_h_tuple_spectrum = Counter()
    child_h_values = Counter()
    child_h_mod8 = Counter()
    support_spectrum = Counter()
    triples = []
    top_bucket_failures = []
    examples = []

    for mid, node in enumerate(child["merged"]):
        if node["type"] != "SC":
            continue
        support = len(owners[mid])
        if support == 1:
            continue
        owner_list = owners[mid]
        support_spectrum[support] += 1
        if support == 3:
            triples.append(tuple(owner_list))
        bucket = profile_bucket(child_g[mid])
        if bucket != child_top_bucket:
            top_bucket_failures.append((mid, bucket, owners[mid]))

        edges, colors = induced_edge_data(parent_graph, owner_list)
        edge_count_spectrum[len(edges)] += 1
        edge_color_spectrum[tuple(sorted(colors))] += 1
        owner_g_tuples[tuple(sorted(profile_bucket(parent_g[o]) for o in owner_list))] += 1
        owner_kind_tuples[
            tuple(sorted(parent["merged"][o]["type"] for o in owner_list))
        ] += 1
        owner_h_tuple_spectrum[
            tuple(sorted(h_value(parent, o) for o in owner_list))
        ] += 1

        child_h = h_value(child, mid)[0]
        child_h_values[child_h] += 1
        child_h_mod8[child_h % 8] += 1

        if len(examples) < 16:
            examples.append(
                {
                    "child": mid,
                    "child_h": child_h,
                    "child_g": bucket,
                    "owners": owner_list,
                    "owner_g": [profile_bucket(parent_g[o]) for o in owner_list],
                    "owner_kind": [parent["merged"][o]["type"] for o in owner_list],
                    "owner_h": [h_value(parent, o) for o in owner_list],
                    "edges": edges,
                    "edge_colors": colors,
                }
            )

    owner_degree = Counter()
    for triple in triples:
        for owner in triple:
            owner_degree[owner] += 1
    pair_intersections = Counter()
    for i, a in enumerate(triples):
        aset = set(a)
        for b in triples[i + 1:]:
            pair_intersections[len(aset & set(b))] += 1
    peel_order, peel_core = peel_hyperedges(triples)

    print("\n" + "-" * 78)
    print(f"n={n}->{n+1}")
    print(
        f"  parent merged graph: V={len(parent['merged'])}, "
        f"E={len(parent_graph['merged_edges'])}"
    )
    print(f"  non-private SC support spectrum={dict(sorted(support_spectrum.items()))}")
    print(f"  owner induced edge-count spectrum={dict(sorted(edge_count_spectrum.items()))}")
    print(f"  owner edge-color tuples={dict(sorted(edge_color_spectrum.items()))}")
    print(f"  owner good-cut tuples={dict(sorted(owner_g_tuples.items()))}")
    print(f"  owner kind tuples={dict(sorted(owner_kind_tuples.items()))}")
    print(f"  collision owner-degree spectrum={dict(sorted(Counter(owner_degree.values()).items()))}")
    print(f"  collision pair-intersection spectrum={dict(sorted(pair_intersections.items()))}")
    print(f"  leaf peeling: removed={len(peel_order)}, core={sorted(peel_core)}")
    print(f"  child H mod 8={dict(sorted(child_h_mod8.items()))}")
    print(f"  child H values={dict(sorted(child_h_values.items()))}")
    print(f"  top-bucket failures={top_bucket_failures}")

    if examples:
        print("  examples:")
        for ex in examples:
            print(
                f"    child={ex['child']}:H={ex['child_h']}:g={ex['child_g']} "
                f"owners={ex['owners']}"
            )
            print(
                f"      owner_g={ex['owner_g']} owner_kind={ex['owner_kind']} "
                f"owner_H={ex['owner_h']}"
            )
            print(f"      parent_edges={ex['edges']} colors={ex['edge_colors']}")

    return {
        "support_spectrum": support_spectrum,
        "edge_count_spectrum": edge_count_spectrum,
        "edge_color_spectrum": edge_color_spectrum,
        "owner_g_tuples": owner_g_tuples,
        "owner_kind_tuples": owner_kind_tuples,
        "owner_degree_spectrum": Counter(owner_degree.values()),
        "pair_intersections": pair_intersections,
        "peel_removed": len(peel_order),
        "peel_core": peel_core,
        "top_bucket_failures": top_bucket_failures,
    }


def main():
    print("=" * 78)
    print("ENDPOINT COLLISION GEOMETRY S95")
    print("=" * 78)
    print("Do support-3 SC collisions form visible triples in the parent metagraph?")

    levels = {n: et.build_level(n) for n in range(2, 8)}
    summaries = []
    for n in range(2, 7):
        summaries.append(summarize_transition(n, levels[n], levels[n + 1]))

    print("\n" + "=" * 78)
    print("SEQUENCES")
    print("=" * 78)
    print(
        "  edge_count_spectrum: "
        f"{[dict(sorted(s['edge_count_spectrum'].items())) for s in summaries]}"
    )
    print(
        "  owner_degree_spectrum: "
        f"{[dict(sorted(s['owner_degree_spectrum'].items())) for s in summaries]}"
    )
    print(
        "  pair_intersections: "
        f"{[dict(sorted(s['pair_intersections'].items())) for s in summaries]}"
    )
    print(f"  peel_removed: {[s['peel_removed'] for s in summaries]}")
    print(f"  peel_core: {[sorted(s['peel_core']) for s in summaries]}")
    print(
        "  top_bucket_failures: "
        f"{[s['top_bucket_failures'] for s in summaries]}"
    )

    print("\n" + "=" * 78)
    print("INTERPRETATION")
    print("=" * 78)
    print(
        "A support-3 endpoint collision is not automatically a literal triangle "
        "in the parent merged arc-flip graph.  The induced edge count records "
        "whether the ternary residual is visible as adjacency geometry or only "
        "as endpoint-transfer incidence."
    )


if __name__ == "__main__":
    main()
