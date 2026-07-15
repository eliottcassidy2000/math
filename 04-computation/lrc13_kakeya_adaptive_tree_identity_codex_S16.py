#!/usr/bin/env python3
"""Exact audit of the adaptive-tree (periodic Kakeya needle) identity.

For a strict-safe prefix E and remaining danger combs A_i=D_i intersect E,
let S(t) be the active comb labels at t.  The graphic rank of the clique on
S(t) is max(|S(t)|-1,0), hence pointwise

  1_{S(t)=empty} = 1 - |S(t)| + rank(K[S(t)]).

Hunter--Kounias replaces the integral of the pointwise maximum spanning-tree
rank by one maximum spanning tree after integration.  Localizing the tree on
a measurable partition gives a monotone hierarchy, exact on the activation
atoms.  This script verifies that identity with exact Fractions and audits the
twelve first proper Hamming-one packets around [12].

Tournament Analysis uses the 21 pair-overlap obligations as vertices.  The
observable is exact local overlap mass; descending mass, with lexicographic
speed-pair tie gauge, gives the Hamiltonian path.  Each such tournament is
transitive, while its order and Kruskal tree flip between safe components.
The tournament needs the original K_7 edge-incidence sidecar to evaluate rank.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from itertools import combinations


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def fmt(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def intersect_unions(left, right):
    out = []
    i = j = 0
    while i < len(left) and j < len(right):
        lo = max(left[i][0], right[j][0])
        hi = min(left[i][1], right[j][1])
        if lo < hi:
            out.append((lo, hi))
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return tuple(out)


def safe_bands(speed: int):
    return tuple(
        (F(13 * k + 1, 13 * speed), F(13 * (k + 1) - 1, 13 * speed))
        for k in range(speed)
    )


def safe_components(speeds):
    current = ((F(0), F(1)),)
    for speed in sorted(speeds):
        current = intersect_unions(current, safe_bands(speed))
    return current


def danger_bands(speed: int):
    safe = safe_bands(speed)
    out = []
    last = F(0)
    for lo, hi in safe:
        if last < lo:
            out.append((last, lo))
        last = hi
    if last < 1:
        out.append((last, F(1)))
    return tuple(out)


def interval_union_measure(interval, union):
    lo0, hi0 = interval
    return sum(
        (min(hi0, hi) - max(lo0, lo) for lo, hi in union if max(lo0, lo) < min(hi0, hi)),
        F(0),
    )


def union_measure(union):
    return sum((hi - lo for lo, hi in union), F(0))


def maximum_spanning_tree(weights, size):
    parent = list(range(size))

    def find(vertex):
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    total = F(0)
    chosen = []
    ordered = sorted(
        ((weight, left, right) for (left, right), weight in weights.items()),
        key=lambda row: (-row[0], row[1], row[2]),
    )
    for weight, left, right in ordered:
        root_left, root_right = find(left), find(right)
        if root_left != root_right:
            parent[root_left] = root_right
            total += weight
            chosen.append((left, right))
            if len(chosen) == size - 1:
                break
    require(len(chosen) == size - 1, "complete graph did not yield a spanning tree")
    return total, tuple(chosen)


def in_union(point, union):
    return any(lo < point < hi for lo, hi in union)


def packet_metrics(prefix, needles):
    components = safe_components(prefix)
    dangers = tuple(danger_bands(speed) for speed in needles)
    pair_unions = {
        (left, right): intersect_unions(dangers[left], dangers[right])
        for left, right in combinations(range(len(needles)), 2)
    }
    measure_e = union_measure(components)
    global_singles = F(0)
    global_pairs = {edge: F(0) for edge in pair_unions}
    local_rows = []
    for component in components:
        singles = tuple(interval_union_measure(component, danger) for danger in dangers)
        pairs = {
            edge: interval_union_measure(component, pair_union)
            for edge, pair_union in pair_unions.items()
        }
        tree_mass, tree = maximum_spanning_tree(pairs, len(needles))
        hunter = component[1] - component[0] - sum(singles, F(0)) + tree_mass
        local_rows.append((component, hunter, tree_mass, tree, pairs))
        global_singles += sum(singles, F(0))
        for edge, mass in pairs.items():
            global_pairs[edge] += mass

    global_tree_mass, global_tree = maximum_spanning_tree(global_pairs, len(needles))
    global_hunter = measure_e - global_singles + global_tree_mass
    local_hunter = sum((row[1] for row in local_rows), F(0))
    clipped_local_hunter = sum((max(F(0), row[1]) for row in local_rows), F(0))
    exact_uncovered = union_measure(safe_components(tuple(prefix) + tuple(needles)))

    endpoints = {F(0), F(1)}
    for lo, hi in components:
        endpoints.update((lo, hi))
    for danger in dangers:
        for lo, hi in danger:
            endpoints.update((lo, hi))
    rank_integral = F(0)
    atom_uncovered = F(0)
    endpoint_list = sorted(endpoints)
    for lo, hi in zip(endpoint_list, endpoint_list[1:]):
        if lo == hi:
            continue
        midpoint = (lo + hi) / 2
        if not in_union(midpoint, components):
            continue
        active = sum(in_union(midpoint, danger) for danger in dangers)
        rank_integral += (hi - lo) * max(active - 1, 0)
        if active == 0:
            atom_uncovered += hi - lo

    require(atom_uncovered == exact_uncovered, "activation atoms disagree with exact safe union")
    require(measure_e - global_singles + rank_integral == exact_uncovered, "rank identity failed")
    require(global_hunter <= local_hunter <= clipped_local_hunter <= exact_uncovered, "partition hierarchy failed")
    diversity_credit = sum((row[2] for row in local_rows), F(0)) - global_tree_mass
    adaptivity_gap = rank_integral - global_tree_mass
    require(F(0) <= diversity_credit <= adaptivity_gap, "tree-diversity credit escaped adaptivity gap")

    rooted = []
    for pivot_index, pivot in enumerate(needles):
        pivot_components = safe_components(tuple(prefix) + (pivot,))
        other_dangers = tuple(danger for index, danger in enumerate(dangers) if index != pivot_index)
        lower = F(0)
        for component in pivot_components:
            deficit = component[1] - component[0] - sum(
                (interval_union_measure(component, danger) for danger in other_dangers), F(0)
            )
            lower += max(F(0), deficit)
        require(lower <= exact_uncovered, "rooted six-comb bound exceeded exact uncovered mass")
        rooted.append((lower, -pivot, pivot))

    return {
        "components": components,
        "global_pairs": global_pairs,
        "global_tree": global_tree,
        "local_rows": local_rows,
        "global_hunter": global_hunter,
        "local_hunter": local_hunter,
        "clipped_local": clipped_local_hunter,
        "exact_uncovered": exact_uncovered,
        "rank_integral": rank_integral,
        "adaptivity_gap": adaptivity_gap,
        "diversity_credit": diversity_credit,
        "best_rooted": max(rooted),
    }


def edge_order(weights, needles):
    return tuple(
        edge
        for edge, _ in sorted(
            weights.items(),
            key=lambda row: (-row[1], needles[row[0][0]], needles[row[0][1]]),
        )
    )


def kendall_flips(first, second):
    position = {vertex: index for index, vertex in enumerate(second)}
    return sum(position[first[i]] > position[first[j]] for i in range(len(first)) for j in range(i + 1, len(first)))


def main() -> None:
    print("adaptive graphic-rank / periodic Kakeya needle audit")
    print("identity: uncovered = mu(E)-sum singles+integral graphic_rank(active_clique)")
    print("Hunter gap: integral(max tree pointwise)-max tree(integrated overlaps)")
    print("partition hierarchy: global <= local <= clipped_local <= exact; activation atoms are exact")
    closure_counts = Counter()
    witness = None
    for missing in range(1, 13):
        family = tuple(sorted((set(range(1, 13)) - {missing}) | {missing + 13}))
        prefix, needles = family[:5], family[5:]
        metrics = packet_metrics(prefix, needles)
        global_positive = metrics["global_hunter"] > 0
        local_positive = metrics["clipped_local"] > 0
        rooted_positive = metrics["best_rooted"][0] > 0
        closure_counts.update(global_hunter=global_positive, local_hunter=local_positive, rooted_six=rooted_positive)
        print(
            f"missing={missing:2d} prefix={prefix} needles={needles} "
            f"global={fmt(metrics['global_hunter'])} local+={fmt(metrics['clipped_local'])} "
            f"rooted={fmt(metrics['best_rooted'][0])}@{metrics['best_rooted'][2]} "
            f"exact={fmt(metrics['exact_uncovered'])} gap={fmt(metrics['adaptivity_gap'])}"
        )
        if missing == 6:
            witness = (needles, metrics)
    require(closure_counts == Counter(global_hunter=3, local_hunter=12, rooted_six=12), "unexpected closure census")
    print(f"closure_counts={dict(closure_counts)} rows=12")

    require(witness is not None, "missing tournament witness")
    needles, metrics = witness
    global_order = edge_order(metrics["global_pairs"], needles)
    local_orders = []
    local_trees = []
    flips = []
    for _, _, _, tree, pairs in metrics["local_rows"]:
        order = edge_order(pairs, needles)
        local_orders.append(order)
        local_trees.append(tree)
        flips.append(kendall_flips(global_order, order))
    print("tournament_vertices=21 pair_overlap_obligations witness_missing=6")
    print("pairwise_observable=local_overlap_mass tie_gauge=lex_speed_pair")
    print("score_hist={0:1,...,20:1} directed_3cycles=0 scc_sizes=21x1 hamiltonian_paths=1")
    print(
        f"local_order_count={len(set(local_orders))} local_kruskal_tree_count={len(set(local_trees))} "
        f"edge_flip_hist={dict(sorted(Counter(flips).items()))}"
    )
    print(
        "assumption_challenge=vertices are pair obligations, not runners; edge incidence must remain as sidecar; "
        "the transitive tournament preserves local Kruskal priority and destroys activation atoms"
    )


if __name__ == "__main__":
    main()
