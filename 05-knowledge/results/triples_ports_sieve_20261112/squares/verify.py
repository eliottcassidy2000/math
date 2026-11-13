"""Exact degree-forcing and theta-port certificates for the first 25 Q_n.

Python standard library only. No removable assertions. Normal and -O runs
must write identical receipts. The independent path census does not call the
forcing or theta routines. No graph-isomorphism quotient is used.
"""

from collections import Counter, deque
from itertools import combinations, permutations, product
from math import comb, isqrt
from pathlib import Path
import argparse
import json


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def square_edges(n):
    return tuple((a, b) for a in range(1, n + 1)
                 for b in range(a + 1, n + 1)
                 if isqrt(a + b) ** 2 == a + b)


def components(n, edges):
    adjacency = {v: [] for v in range(1, n + 1)}
    for a, b in edges:
        adjacency[a].append(b)
        adjacency[b].append(a)
    unseen = set(adjacency)
    out = []
    while unseen:
        root = min(unseen)
        unseen.remove(root)
        pending = [root]
        found = []
        while pending:
            v = pending.pop()
            found.append(v)
            for w in adjacency[v]:
                if w in unseen:
                    unseen.remove(w)
                    pending.append(w)
        out.append(sorted(found))
    return out


def cycle_exists(n, edges):
    parent = list(range(n + 1))

    def find(v):
        while parent[v] != v:
            v = parent[v]
        return v

    for a, b in edges:
        a, b = find(a), find(b)
        if a == b:
            return True
        parent[a] = b
    return False


def degree_forcing(n, endpoints):
    """Fixed deterministic order; both endpoints have demand 1, others 2."""
    state = {e: None for e in square_edges(n)}
    target = {v: 1 if v in endpoints else 2 for v in range(1, n + 1)}
    trace = []
    while True:
        changed = False
        for v in target:
            forced = [e for e, bit in state.items() if v in e and bit == 1]
            undecided = [e for e, bit in state.items() if v in e and bit is None]
            if len(forced) > target[v] or len(forced) + len(undecided) < target[v]:
                trace.append(["degree_contradiction", v, len(forced),
                              len(undecided), target[v]])
                return "degree", state, trace
            bit = (0 if len(forced) == target[v]
                   else 1 if len(forced) + len(undecided) == target[v]
                   else None)
            if bit is not None and undecided:
                for edge in undecided:
                    state[edge] = bit
                trace.append(["force" if bit else "delete", v, undecided])
                changed = True
        forced_edges = [e for e, bit in state.items() if bit == 1]
        if cycle_exists(n, forced_edges):
            trace.append(["forced_cycle"])
            return "cycle", state, trace
        if not changed:
            return ("unique" if all(bit is not None for bit in state.values())
                    else "open"), state, trace


def ordered_path(n, edges, endpoints):
    adjacency = {v: [] for v in range(1, n + 1)}
    for a, b in edges:
        adjacency[a].append(b)
        adjacency[b].append(a)
    check(all(len(adjacency[v]) == (1 if v in endpoints else 2)
              for v in adjacency), "candidate degree mismatch")
    if len(components(n, edges)) != 1:
        return None
    walk = [min(endpoints)]
    previous = None
    while True:
        options = [w for w in adjacency[walk[-1]] if w != previous]
        if not options:
            break
        check(len(options) == 1, "path unexpectedly branches")
        previous = walk[-1]
        walk.append(options[0])
    check(len(walk) == n and len(set(walk)) == n, "invalid path traversal")
    return tuple(walk)


def independent_paths(n):
    """Vertex DFS with square-shell adjacency, independent of edge forcing."""
    if n == 1:
        return [(1,)]
    adjacency = {v: set() for v in range(1, n + 1)}
    for root in range(2, isqrt(2 * n - 1) + 1):
        for a in range(1, n + 1):
            b = root * root - a
            if a < b <= n:
                adjacency[a].add(b)
                adjacency[b].add(a)
    leaves = [v for v in adjacency if len(adjacency[v]) == 1]
    if len(leaves) > 2 or any(not adjacency[v] for v in adjacency):
        return []
    starts = [min(leaves)] if leaves else list(adjacency)
    found = set()

    def walk(path, visited):
        if len(path) == n:
            found.add(min(tuple(path), tuple(reversed(path))))
            return
        for w in sorted(adjacency[path[-1]]):
            if not (visited >> (w - 1)) & 1:
                walk(path + [w], visited | (1 << (w - 1)))

    for start in starts:
        walk([start], 1 << (start - 1))
    return sorted(found)


def theta_arms(edges):
    adjacency = {}
    for a, b in edges:
        adjacency.setdefault(a, []).append(b)
        adjacency.setdefault(b, []).append(a)
    hubs = sorted(v for v in adjacency if len(adjacency[v]) == 3)
    check(len(hubs) == 2, "residual graph does not have two trivalent hubs")
    check(all(v in hubs or len(adjacency[v]) == 2 for v in adjacency),
          "residual graph has another junction")
    left, right = hubs
    arms = []
    for first in sorted(adjacency[left]):
        arm = [left, first]
        while arm[-1] != right:
            options = [v for v in adjacency[arm[-1]] if v != arm[-2]]
            check(len(options) == 1, "bad theta interior")
            check(options[0] not in arm, "theta arm loops")
            arm.append(options[0])
        arms.append(arm)
    arms.sort(key=lambda p: (len(p), p))
    covered = [tuple(sorted(e)) for p in arms for e in zip(p, p[1:])]
    check(len(covered) == len(set(covered)) and set(covered) == set(edges),
          "theta arms do not partition residual edges")
    return hubs, arms


def binomial(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def theta_coefficient(lengths, a, b):
    odd = sum(length % 2 for length in lengths)
    even = len(lengths) - odd
    twice_i, twice_j = a + b - even, a - b + even
    if twice_i % 2 or twice_j % 2:
        return 0
    return binomial(odd, twice_i // 2) * binomial(even, twice_j // 2)


def quotient_tree(n, forced, chosen):
    blocks = components(n, forced)
    owner = {v: i for i, block in enumerate(blocks) for v in block}
    parent = list(range(len(blocks)))

    def find(v):
        while parent[v] != v:
            v = parent[v]
        return v

    for a, b in chosen:
        a, b = find(owner[a]), find(owner[b])
        if a == b:
            return False
        parent[a] = b
    return len(chosen) == len(blocks) - 1


def theta_certificate(n, endpoints):
    status, state, trace = degree_forcing(n, endpoints)
    check(status == "open", "requested theta instance is already decided")
    forced = [e for e, bit in state.items() if bit == 1]
    unresolved = [e for e, bit in state.items() if bit is None]
    hubs, arms = theta_arms(unresolved)
    target = {v: 1 if v in endpoints else 2 for v in range(1, n + 1)}
    residual = {v: target[v] - sum(v in e for e in forced)
                for v in range(1, n + 1)}
    internal = {v for arm in arms for v in arm[1:-1]}
    check(all(residual[v] == 1 for v in internal), "internal demand not 1")
    check(all(residual[v] == 0 for v in target
              if v not in internal and v not in hubs), "hidden demand")
    completions = []
    for bits in product((0, 1), repeat=3):
        chosen = sorted(tuple(sorted(edge)) for bit, arm in zip(bits, arms)
                        for i, edge in enumerate(zip(arm, arm[1:]))
                        if bit ^ (i % 2))
        if not all(sum(v in e for e in chosen) == residual[v] for v in target):
            continue
        support = sorted(forced + chosen)
        path = ordered_path(n, support, endpoints)
        check(quotient_tree(n, forced, chosen) == (path is not None),
              "quotient tree disagrees with full graph")
        completions.append({"first_edge_bits": bits, "selected_residual": chosen,
                            "components": components(n, support), "path": path})
    lengths = [len(arm) - 1 for arm in arms]
    expected = theta_coefficient(lengths, residual[hubs[0]], residual[hubs[1]])
    check(len(completions) == expected, "theta coefficient mismatch")
    return {"n": n, "endpoints": endpoints, "hubs": hubs, "arms": arms,
            "lengths": lengths, "hub_demands": [residual[v] for v in hubs],
            "residual_demands": [[v, residual[v]] for v in sorted(set(hubs) | internal)],
            "forced_edges": forced, "forced_components": components(n, forced),
            "deleted_edges": [e for e, bit in state.items() if bit == 0],
            "unresolved_edges": unresolved, "forcing_trace": trace,
            "degree_completions": len(completions),
            "hamilton_completions": sum(c["path"] is not None for c in completions),
            "completions": completions}


def brute_theta_controls():
    """Enumerate every edge subset, not merely the alternating-arm choices."""
    graph_count = subset_count = 0
    for lengths in product(range(1, 5), repeat=3):
        if lengths.count(1) > 1:  # Preserve the declared simple-graph universe.
            continue
        edges = []
        internal = []
        next_vertex = 2
        for length in lengths:
            middle = list(range(next_vertex, next_vertex + length - 1))
            next_vertex += length - 1
            internal.extend(middle)
            arm = [0] + middle + [1]
            edges.extend(zip(arm, arm[1:]))
        counts = Counter()
        for mask in range(1 << len(edges)):
            degree = [0] * next_vertex
            for i, (a, b) in enumerate(edges):
                if mask >> i & 1:
                    degree[a] += 1
                    degree[b] += 1
            if all(degree[v] == 1 for v in internal):
                counts[degree[0], degree[1]] += 1
            subset_count += 1
        for a, b in product(range(4), repeat=2):
            check(counts[a, b] == theta_coefficient(lengths, a, b),
                  f"independent edge-subset theta mismatch {lengths,a,b}")
        graph_count += 1
    return {"simple_ordered_length_triples_1_to_4": graph_count,
            "edge_subsets_checked": subset_count,
            "hub_demand_pairs_per_graph": 16}


def johnson_controls():
    """All parity words with 1..6 arms, all demands, all degree-completions."""
    parity_words = sectors = states = 0
    for arm_count in range(1, 7):
        for parity in product((0, 1), repeat=arm_count):
            # Length 2 for even, 3 for odd avoids coincident direct edges.
            lengths = [3 if p else 2 for p in parity]
            groups = {}
            for bits in product((0, 1), repeat=arm_count):
                a = sum(bits)
                b = sum(bit if p else 1 - bit for bit, p in zip(bits, parity))
                groups.setdefault((a, b), []).append(bits)
            for (a, b), configurations in groups.items():
                odd, even = sum(parity), arm_count - sum(parity)
                i, j = (a + b - even) // 2, (a - b + even) // 2
                check(len(configurations) == theta_coefficient(lengths, a, b),
                      "multi-arm coefficient mismatch")
                expected_degree = i * (odd - i) + j * (even - j)
                expected_diameter = min(i, odd - i) + min(j, even - j)
                adjacency = {x: [] for x in configurations}
                for x, y in combinations(configurations, 2):
                    changed = [k for k in range(arm_count) if x[k] != y[k]]
                    if len(changed) == 2 and parity[changed[0]] == parity[changed[1]]:
                        adjacency[x].append(y)
                        adjacency[y].append(x)
                check(all(len(v) == expected_degree for v in adjacency.values()),
                      "Johnson product degree mismatch")
                root = configurations[0]
                distance, queue = {root: 0}, deque([root])
                while queue:
                    x = queue.popleft()
                    for y in adjacency[x]:
                        if y not in distance:
                            distance[y] = distance[x] + 1
                            queue.append(y)
                check(len(distance) == len(configurations), "Johnson graph disconnected")
                check(max(distance.values()) == expected_diameter,
                      "Johnson product diameter mismatch")
                sectors += 1
                states += len(configurations)
            parity_words += 1
    return {"arm_range": [1, 6], "parity_words": parity_words,
            "feasible_demand_sectors": sectors, "configurations": states}


def relocation_controls():
    """Every labelled degree-(1,1,2,...) support on 4..7 vertices."""
    totals = Counter()
    by_size = []

    def test_support(n, support):
        degree = Counter(v for edge in support for v in edge)
        endpoints = sorted(v for v in range(1, n + 1) if degree[v] == 1)
        if len(endpoints) != 2 or any(degree[v] not in (1, 2) for v in range(1, n + 1)):
            return 0
        blocks = components(n, support)
        owner = {v: i for i, block in enumerate(blocks) for v in block}
        full_path = ordered_path(n, support, endpoints)
        for b, t in (endpoints, endpoints[::-1]):
            path = full_path
            if path and path[0] != b:
                path = path[::-1]
            for edge in support:
                if any(v in endpoints for v in edge):
                    continue
                for a, v in (edge, edge[::-1]):
                    if path:
                        predicted = path.index(a) < path.index(v)
                        case = "path_retained" if predicted else "path_closes_cycle"
                    elif len(blocks) == 2 and owner[a] != owner[b]:
                        predicted, case = True, "one_cycle_opened"
                    else:
                        predicted, case = False, "other_cycle_remains"
                    moved = [e for e in support if e != edge] + [(v, n + 1), (b, n + 1)]
                    actual = ordered_path(n + 1, moved, (a, t)) is not None
                    check(actual == predicted, "endpoint relocation criterion mismatch")
                    totals[case] += 1
        return 1

    for n in range(4, 8):
        complete = list(combinations(range(1, n + 1), 2))
        supports = sum(test_support(n, edges) for edges in combinations(complete, n - 1))
        by_size.append({"n": n, "degree_supports": supports,
                        "edge_subsets_tested": comb(len(complete), n - 1)})
    # The least order supporting two cycles plus a nonempty endpoint path.
    multi_cycle = [(1, 2), (2, 3), (1, 3), (4, 5), (5, 6), (4, 6), (7, 8)]
    check(test_support(8, multi_cycle) == 1, "multiple-cycle hostile rejected")
    # The first25 residual candidates realize both positive cases at the SAME birth.
    source = theta_certificate(25, (10, 18))
    target = theta_certificate(26, (13, 18))
    relocated = []
    for candidate in source["completions"]:
        support = sorted(source["forced_edges"] + candidate["selected_residual"])
        check((13, 23) in support, "birth relocation source edge absent")
        moved = [edge for edge in support if edge != (13, 23)] + [(10, 26), (23, 26)]
        relocated.append(ordered_path(26, moved, (13, 18)))
    check(set(relocated) == {tuple(c["path"]) for c in target["completions"]},
          "actual Q25 to Q26 relocation mismatch")
    arithmetic = 0
    for m in range(4, 81):
        for r in range((m * m - 2 * m - 4) // 2 + 1):
            z, b, v, a = m * m - 1 - r, r + 1, 2 * m + 2 + r, m * m - 2 * m - 2 - r
            if a == v:
                continue
            check(len({a, b, v}) == 3 and min(a, b, v) >= 1 and max(a, b, v) < z,
                  "arithmetic relocation vertices invalid")
            check(a + v == b + z == m * m and v + z == (m + 1) ** 2,
                  "arithmetic relocation square identities fail")
            arithmetic += 1
    return {"complete_support_universe": by_size, "classification_checks": dict(sorted(totals.items())),
            "two_triangle_plus_edge_hostile": multi_cycle,
            "actual_birth": {"source": 25, "target": 26, "a": 13, "b": 10, "v": 23,
                             "t": 18, "z": 26, "relocated_paths": relocated},
            "arithmetic_gadgets_m4_to_80": arithmetic}


def directed_cycle_opening_controls():
    """All loopless directed degree supports, then their tournament-realizable subuniverse."""
    totals = Counter()
    tournament_totals = Counter()
    rows = []

    def inspect(n, successor, b, t):
        path = [b]
        while path[-1] != t:
            path.append(successor[path[-1]])
            check(len(path) <= n, "degree support path loops")
        remaining = set(range(1, n + 1)) - set(path)
        cycles = []
        while remaining:
            root = min(remaining)
            cycle = [root]
            remaining.remove(root)
            w = successor[root]
            while w != root:
                check(w in remaining, "directed cycle decomposition failed")
                remaining.remove(w)
                cycle.append(w)
                w = successor[w]
            cycles.append(cycle)
        tournament_realizable = all(len(cycle) >= 3 for cycle in cycles)
        for a, v in successor.items():
            if len({a, b, v, t}) != 4:
                continue
            moved = dict(successor)
            moved[a], moved[n + 1] = n + 1, b
            walk = [v]
            while walk[-1] in moved and moved[walk[-1]] not in walk:
                walk.append(moved[walk[-1]])
            actual = len(walk) == n + 1 and walk[-1] == t
            expected = len(cycles) == 1 and a not in path
            check(actual == expected, "directed cycle opening criterion mismatch")
            reason = ("one_cycle_opened" if expected else
                      "arc_on_main_path" if a in path else "other_cycles_remain")
            totals[reason] += 1
            if tournament_realizable:
                tournament_totals[reason] += 1
            if actual:
                # Explicitly verify the inverse map on the target support.
                inverse = dict(moved)
                del inverse[n + 1]
                inverse[a] = v
                check(inverse == successor, "directed inverse fails")
        return tournament_realizable

    for n in range(4, 8):
        supports = realizable = 0
        for b, t in permutations(range(1, n + 1), 2):
            domain = [v for v in range(1, n + 1) if v != t]
            codomain = [v for v in range(1, n + 1) if v != b]
            for image in permutations(codomain):
                if any(a == v for a, v in zip(domain, image)):
                    continue
                supports += 1
                realizable += inspect(n, dict(zip(domain, image)), b, t)
        rows.append({"n": n, "loopless_degree_supports": supports,
                     "tournament_realizable_supports": realizable})
    inspect(8, {1: 2, 2: 3, 3: 1, 4: 5, 5: 6, 6: 4, 7: 8}, 7, 8)
    # Every tournament completion of this fixed source and prescribed two new arcs.
    source_arcs = {(1, 2), (2, 3), (3, 1), (4, 5)} | {
        (x, y) for x in (1, 2, 3) for y in (4, 5)}
    base_paths = [p for p in permutations(range(1, 6))
                  if all(edge in source_arcs for edge in zip(p, p[1:]))]
    check(len(base_paths) == 3, "minimal tournament host source count")
    completions = []
    for bits in product((0, 1), repeat=3):
        arcs = source_arcs | {(3, 6), (6, 4)}
        arcs |= {(6, v) if bit else (v, 6) for v, bit in zip((1, 2, 5), bits)}
        paths = [p for p in permutations(range(1, 7))
                 if all(edge in arcs for edge in zip(p, p[1:]))]
        selected = [p for p in paths if p[0] == 1 and p[-1] == 5
                    and p.index(6) > 0 and p.index(6) < 5
                    and p[p.index(6) - 1:p.index(6) + 2] == (3, 6, 4)]
        check(selected == [(1, 2, 3, 6, 4, 5)], "conditioned tournament host count")
        completions.append({"new_vertex_beats_1_2_5": bits,
                            "total_hamilton_paths": len(paths),
                            "conditioned_hamilton_paths": len(selected)})
    return {"complete_support_universe": rows,
            "all_loopless_classification_checks": dict(sorted(totals.items())),
            "tournament_realizable_classification_checks": dict(sorted(tournament_totals.items())),
            "multiple_cycle_hostile": "two directed triangles plus7->8 at order8",
            "minimal_tournament_host": {"source_arcs": sorted(source_arcs),
                                        "source_hamilton_paths": base_paths,
                                        "fixed_new_arcs": [[3, 6], [6, 4]],
                                        "all_eight_extensions": completions}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="receipt.json")
    args = parser.parse_args()
    rows = []
    census = {1: [(1,)]}
    theta_cases = []
    for n in range(2, 26):
        census[n] = independent_paths(n)
        by_endpoints = {}
        for path in census[n]:
            by_endpoints.setdefault(tuple(sorted((path[0], path[-1]))), set()).add(path)
        counts = Counter()
        pair_rows = []
        for endpoints in combinations(range(1, n + 1), 2):
            status, state, _ = degree_forcing(n, endpoints)
            counts[status] += 1
            expected = by_endpoints.get(endpoints, set())
            if status in ("degree", "cycle"):
                found = set()
            elif status == "unique":
                path = ordered_path(n, [e for e, bit in state.items() if bit == 1], endpoints)
                check(path is not None, "fully forced forest is not Hamiltonian")
                found = {path}
            else:
                case = theta_certificate(n, endpoints)
                theta_cases.append(case)
                found = {tuple(c["path"]) for c in case["completions"] if c["path"]}
            check(found == expected, f"independent DFS disagrees at {n,endpoints}")
            pair_rows.append({"endpoints": endpoints, "status": status,
                              "hamilton_paths": len(found)})
        rows.append({"n": n, "edge_count": len(square_edges(n)),
                     "forcing_statuses": dict(sorted(counts.items())),
                     "hamilton_paths": len(census[n]), "endpoint_cases": pair_rows})
    check([(c["n"], c["endpoints"]) for c in theta_cases]
          == [(25, (10, 18)), (25, (18, 22))], "first25 open kernel list changed")
    check({n: len(paths) for n, paths in census.items() if paths}
          == {1: 1, 15: 1, 16: 1, 17: 1, 23: 3, 25: 10}, "inherited census mismatch")
    hostile = theta_certificate(26, (13, 18))
    first = theta_cases[0]
    for coordinate in ("hubs", "arms", "lengths", "hub_demands",
                       "residual_demands", "unresolved_edges"):
        check(first[coordinate] == hostile[coordinate], "theta kernel not identical")
    check(first["hamilton_completions"] == 1 and hostile["hamilton_completions"] == 2,
          "missing outside-partition hostile failure")
    independent26 = independent_paths(26)
    check(len(independent26) == 12, "Q26 inherited control census mismatch")
    check({tuple(c["path"]) for c in hostile["completions"]}
          == {p for p in independent26 if tuple(sorted((p[0], p[-1]))) == (13, 18)},
          "Q26 hostile mismatch independent DFS")
    receipt = {"status": "FINITE-EXACT", "universe": "Q_n, vertices 1..n; distinct a,b adjacent iff a+b is square",
               "census_range": [1, 25], "endpoint_instances": sum(comb(n, 2) for n in range(2, 26)),
               "ordering": "ascending vertices, ascending lexicographic edges; cycle check after each forcing sweep",
               "singleton_Q1_paths": [[1]], "rows": rows, "theta_cases": theta_cases,
               "outside_partition_hostile": hostile,
               "independent_paths": [{"n": n, "paths": paths} for n, paths in census.items()],
               "independent_Q26_paths": independent26,
               "theta_edge_subset_controls": brute_theta_controls(),
               "johnson_product_controls": johnson_controls(),
               "endpoint_relocation_controls": relocation_controls(),
               "directed_cycle_opening_controls": directed_cycle_opening_controls()}
    destination = Path(__file__).with_name(args.output)
    destination.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"endpoint_instances": receipt["endpoint_instances"],
                      "first25_theta_cases": len(theta_cases),
                      "theta_edge_subset_controls": receipt["theta_edge_subset_controls"],
                      "johnson_product_controls": receipt["johnson_product_controls"],
                      "undirected_relocation_checks": sum(receipt["endpoint_relocation_controls"]["classification_checks"].values()),
                      "directed_relocation_checks": sum(receipt["directed_cycle_opening_controls"]["all_loopless_classification_checks"].values()),
                      "tournament_relocation_checks": sum(receipt["directed_cycle_opening_controls"]["tournament_realizable_classification_checks"].values()),
                      "eight_tournament_extension_totals": [row["total_hamilton_paths"] for row in receipt["directed_cycle_opening_controls"]["minimal_tournament_host"]["all_eight_extensions"]],
                      "output": str(destination)}, indent=2))


if __name__ == "__main__":
    main()
