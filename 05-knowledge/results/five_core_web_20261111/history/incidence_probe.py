#!/usr/bin/env python3
"""Small exact controls for triangle--arc incidence; no inherited census rerun."""
import itertools as it
import json
from pathlib import Path
import networkx as nx


def rows_from_code(n, code):
    rows = [0] * n
    for k, (u, v) in enumerate(it.combinations(range(n), 2)):
        a, b = (u, v) if code >> k & 1 else (v, u)
        rows[a] |= 1 << b
    return rows


def path_count(rows):
    n = len(rows)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            choices = rows[v] & ~mask
            while choices:
                bit = choices & -choices
                w = bit.bit_length() - 1
                dp[mask | bit][w] += dp[mask][v]
                choices ^= bit
    return sum(dp[-1])


def triangles(rows):
    return [frozenset(q) for q in it.combinations(range(len(rows)), 3)
            if all(sum(rows[v] >> w & 1 for w in q) == 1 for v in q)]


def overlap(rows):
    ts = triangles(rows)
    graph = nx.Graph()
    graph.add_nodes_from(range(len(ts)))
    graph.add_edges_from((i, j) for i, j in it.combinations(range(len(ts)), 2)
                         if len(ts[i] & ts[j]) == 2)
    return ts, graph


def check_incidence(rows):
    ts, graph = overlap(rows)
    for clique in nx.find_cliques(graph):
        if len(clique) >= 2:
            common = set.intersection(*(set(ts[v]) for v in clique))
            if len(common) != 2:
                raise RuntimeError((rows, clique, common))
    for v in graph:
        cliques = [c for c in nx.find_cliques(graph) if v in c and len(c) >= 2]
        if len(cliques) > 3:
            raise RuntimeError((rows, v, cliques))


def odd_cycles(rows):
    n = len(rows)
    masks = []
    for size in range(3, n + 1, 2):
        for subset in it.combinations(range(n), size):
            first = subset[0]
            for tail in it.permutations(subset[1:]):
                cycle = (first,) + tail
                if all(rows[u] >> v & 1 for u, v in zip(cycle, cycle[1:] + cycle[:1])):
                    masks.append(sum(1 << v for v in subset))
    return masks


def packing_polynomial(rows):
    states = {(0, 0): 1}
    for mask in odd_cycles(rows):
        for (used, degree), count in list(states.items()):
            if not used & mask:
                key = (used | mask, degree + 1)
                states[key] = states.get(key, 0) + count
    out = [0] * (1 + len(rows) // 3)
    for (_, degree), count in states.items():
        out[degree] += count
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def chain(n):
    rows = [0] * n
    for u, v in it.combinations(range(n), 2):
        a, b = (u, v) if v == u + 1 else (v, u)
        rows[a] |= 1 << b
    return rows


def blowup(sizes):
    blocks = []
    start = 0
    for size in sizes:
        blocks.append(list(range(start, start + size)))
        start += size
    rows = [0] * start
    for block in blocks:
        for u, v in it.combinations(block, 2):
            rows[u] |= 1 << v
    for i in range(3):
        for u in blocks[i]:
            for v in blocks[(i + 1) % 3]:
                rows[u] |= 1 << v
    return rows


def main():
    checked = 0
    for n in range(1, 6):
        for code in range(1 << (n * (n - 1) // 2)):
            check_incidence(rows_from_code(n, code))
            checked += 1
    families = []
    for n in range(3, 9):
        item = {"n": n}
        for name, rows in [("chain", chain(n)), ("book", blowup((1, 1, n - 2)))]:
            ts, graph = overlap(rows)
            expected = nx.path_graph(n - 2) if name == "chain" else nx.complete_graph(n - 2)
            if not nx.is_isomorphic(graph, expected):
                raise RuntimeError((name, n, "overlap mismatch"))
            poly = packing_polynomial(rows)
            h = path_count(rows)
            if sum(c * 2**k for k, c in enumerate(poly)) != h:
                raise RuntimeError((name, n, "OCF mismatch"))
            item[name] = {"rows": rows, "H": h, "packing_polynomial": poly}
        families.append(item)
    hamming = []
    for sizes in [(1, 2, 2), (2, 2, 2), (1, 2, 3)]:
        rows = blowup(sizes)
        ts, graph = overlap(rows)
        expected = nx.cartesian_product(nx.cartesian_product(nx.complete_graph(sizes[0]),
                    nx.complete_graph(sizes[1])), nx.complete_graph(sizes[2]))
        if not nx.is_isomorphic(graph, expected):
            raise RuntimeError((sizes, "Hamming mismatch"))
        hamming.append({"sizes": sizes, "rows": rows, "H": path_count(rows),
                        "triangles": len(ts), "overlap_edges": graph.number_of_edges()})
    q13 = [18, 8, 3, 5, 14]
    q13_triangles, q13_graph = overlap(q13)
    if not nx.is_isomorphic(q13_graph, nx.path_graph(4)):
        raise RuntimeError("H13 is not golden P4")
    switches = []
    for cut in [{0, 1}, {1, 2}]:
        base = [r | 32 for r in q13] + [0]
        rows = [r | (64 if v not in cut else 0) for v, r in enumerate(base)]
        rows.append(sum(1 << v for v in cut))
        switches.append({"cut": sorted(cut), "rows": rows, "H": path_count(rows)})
    golden_loss = {"historical_source": "THM-1780-H-leaves-the-spectral-ladder-at-n6.md",
                   "H13_rows": q13, "H13_triangles": [sorted(t) for t in q13_triangles],
                   "H13_packing_polynomial": packing_polynomial(q13),
                   "H17_rows": chain(6), "H17_packing_polynomial": packing_polynomial(chain(6)),
                   "unit_switch_source": "THM-5117-physical-cut-switches-reversal-cancellation-and-primitive-sparse-cubes.md",
                   "unit_switches": switches}
    # Search for the first arc-overlap collision with unequal H, by labelled code.
    # Orders 1..5 exhaustive. Order 6 stops immediately on its first witness.
    seen = {}
    collision = None
    checked_collision = 0
    for n in range(1, 7):
        for code in range(1 << (n * (n - 1) // 2)):
            rows = rows_from_code(n, code)
            ts, graph = overlap(rows)
            h = path_count(rows)
            key = nx.weisfeiler_lehman_graph_hash(graph)
            for old_rows, old_h, old_graph in seen.get(key, []):
                if old_h != h and nx.is_isomorphic(graph, old_graph):
                    collision = {"rows1": old_rows, "H1": old_h, "rows2": rows, "H2": h,
                                 "overlap_edges1": list(old_graph.edges()),
                                 "overlap_edges2": list(graph.edges()),
                                 "triangle_count": len(ts), "last_order": n, "last_code": code}
                    break
            checked_collision += 1
            if collision:
                break
            if not any(old_h == h and nx.is_isomorphic(graph, old_graph)
                       for _, old_h, old_graph in seen.get(key, [])):
                seen.setdefault(key, []).append((rows, h, graph))
        if collision:
            break
    out = {"scope": "All labelled tournaments through five for clique controls; three explicit blowups; chains/books n=3..8; collision search stops at first witness.",
           "incidence_controls": checked, "families": families, "hamming_controls": hamming,
           "golden_loss": golden_loss,
           "collision_search_checked": checked_collision, "collision": collision}
    destination = Path(__file__).with_name("incidence_receipt.json")
    destination.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
