"""Independent small exact audit of ordered path/graph/tournament coordinates.

Universe: every path-normalized tournament for n=1..4; every labelled
three-edge simple graph on [4]. No inherited filters and no external packages.
Run from any directory: python <this-file> [--output <json-file>].
"""
import argparse
from collections import Counter, defaultdict
from itertools import combinations, permutations
import json
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def tournament(n, backward):
    return {(j, i) if (i, j) in backward else (i, j)
            for i, j in combinations(range(1, n + 1), 2)}


def canon(n, arcs):
    return min(tuple(int((p[i], p[j]) in arcs)
                     for i, j in combinations(range(n), 2))
               for p in permutations(range(1, n + 1)))


def profiles(n, arcs):
    pp = list(permutations(range(1, n + 1)))
    h = sum(all((p[i], p[i + 1]) in arcs for i in range(n - 1)) for p in pp)
    aut = sum({(p[i - 1], p[j - 1]) for i, j in arcs} == arcs for p in pp)
    scores = sorted(sum(a == i for a, b in arcs) for i in range(1, n + 1))
    cyclic = sum(all(sum((u, v) in arcs for v in q if v != u) == 1 for u in q)
                 for q in combinations(range(1, n + 1), 3))
    require(cyclic == n * (n - 1) * (n - 2) // 6 - sum(s * (s - 1) // 2 for s in scores),
            f'triangle count and score formula disagree: n={n}, arcs={sorted(arcs)}')
    return h, aut, scores, cyclic


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    rows = []
    checks = Counter()
    all_counts = {}
    for n in range(1, 5):
        free = [(i, j) for i, j in combinations(range(1, n + 1), 2) if j > i + 1]
        fibres = defaultdict(list)
        seen_g = set()
        seen_u = set()
        for mask in range(1 << len(free)):
            backward = {edge for k, edge in enumerate(free) if mask & (1 << k)}
            arcs = tournament(n, backward)
            graph = {(i, j - 1) for i, j in backward}
            u = tournament(n - 1, graph)
            require({(i, j + 1) for i, j in graph} == backward,
                    f'graph inverse failed: n={n}, mask={mask}')
            require({(i, j) for i, j in combinations(range(1, n), 2) if (j, i) in u} == graph,
                    f'tournament inverse failed: n={n}, mask={mask}')
            seen_g.add(tuple(sorted(graph)))
            seen_u.add(tuple(sorted(u)))
            for i, j, k in combinations(range(1, n + 1), 3):
                bits = (int((i, j) in backward), int((j, k) in backward), int((i, k) in backward))
                directed_cycle = all(sum((a, b) in arcs for b in (i, j, k) if b != a) == 1
                                     for a in (i, j, k))
                require(directed_cycle == (bits in ((0, 0, 1), (1, 1, 0))),
                        f'triangle bit predicate failed: n={n}, mask={mask}, triple={(i,j,k)}')
                checks['triangle_predicate_cases'] += 1
            reflected_graph = {tuple(sorted((n - i, n - j))) for i, j in graph}
            reversed_converse = {(n + 1 - j, n + 1 - i) for i, j in arcs}
            recovered_reflection = {(i, j - 1) for i, j in free if (j, i) in reversed_converse}
            require(reflected_graph == recovered_reflection,
                    f'mirror map failed: n={n}, mask={mask}')
            checks['mirror_cases'] += 1
            h, aut, scores, cyclic = profiles(n, arcs)
            fibres[canon(n, arcs)].append((h, aut))
            require((cyclic == 0) == (not graph),
                    f'transitivity predicate failed: n={n}, mask={mask}')
            checks['transitivity_cases'] += 1
            reached = {n}
            while True:
                grown = reached | {b for a, b in arcs if a in reached}
                if grown == reached:
                    break
                reached = grown
            interval_covers = all(any(a <= k <= b for a, b in graph) for k in range(1, n))
            require((len(reached) == n) == interval_covers,
                    f'strong interval-cover predicate failed: n={n}, mask={mask}')
            checks['strong_interval_cover_cases'] += 1
            if n == 4:
                class_names = {(0, 1, 2, 3): 'transitive', (0, 2, 2, 2): 'C3 -> sink',
                               (1, 1, 1, 3): 'source -> C3', (1, 1, 2, 2): 'strong'}
                rows.append({'mask': mask, 'bits_x13_x14_x24': ''.join(str(int(edge in backward)) for edge in free),
                             'graph_edges': sorted(graph), 'tournament_type': class_names[tuple(scores)],
                             'outdegrees_sorted': scores, 'directed_triangles': cyclic,
                             'H': h, 'Aut': aut, 'U_type': 'C3' if profiles(3, u)[3] else 'transitive'})
        require(len(seen_g) == len(seen_u) == 1 << len(free),
                f'target bijection cardinality failed: n={n}')
        for fibre in fibres.values():
            require(all(h == len(fibre) * aut for h, aut in fibre),
                    f'H/Aut fibre formula failed: n={n}, fibre={fibre}')
            checks['path_orbit_fibres'] += 1
        all_counts[n] = {'path_marked_classes': 1 << len(free), 'unmarked_classes': len(fibres)}

    graph_classes = Counter()
    edge_universe = set(combinations(range(1, 5), 2))
    for edges in combinations(sorted(edge_universe), 3):
        degree = sorted(sum(i in e for e in edges) for i in range(1, 5))
        kind = {(1, 1, 2, 2): 'P4', (1, 1, 1, 3): 'K1,3', (0, 2, 2, 2): 'K3 + K1'}[tuple(degree)]
        graph_classes[kind] += 1
        complement_degree = sorted(sum(i in e for e in edge_universe - set(edges)) for i in range(1, 5))
        require(complement_degree == sorted(3 - d for d in degree),
                f'complement degree formula failed: edges={edges}')
    require(graph_classes == {'P4': 12, 'K1,3': 4, 'K3 + K1': 4},
            f'three-edge graph census failed: counts={graph_classes}')
    single_edges = [r for r in rows if len(r['graph_edges']) == 1]
    require(len({r['tournament_type'] for r in single_edges}) == 3,
            'single-edge order-loss hostile failed')
    for n, expected in ((1, 1), (2, 1), (3, 2), (4, 4)):
        require(all_counts[n]['unmarked_classes'] == expected,
                f'unmarked tournament census failed: n={n}')
    result = {'status': 'FINITE-EXACT; all explicit checks passed (also active under python -O)',
              'universe': 'All path-normalized tournaments n=1..4; all 20 labelled 3-edge graphs on [4]',
              'filters': 'None beyond the stated universes',
              'all_counts': all_counts, 'checks': dict(checks), 'n4_table': rows,
              'three_edge_graph_counts': dict(graph_classes),
              'minimal_order_loss_hostile': single_edges}
    rendered = json.dumps(result, indent=2) + '\n'
    if args.output:
        args.output.write_text(rendered, encoding='utf-8')
    print(rendered)


if __name__ == '__main__':
    main()
