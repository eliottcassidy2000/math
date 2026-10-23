"""Independent exhaustive small insertion/cut-energy controls; stdlib only."""
from collections import Counter
from itertools import combinations, permutations
import json

checks = 0


def require(test, why):
    global checks
    checks += 1
    if not test:
        raise RuntimeError(why)


def literal_arcs(graph, order):
    arcs = set()
    for i, u in enumerate(order):
        for v in order[i+1:]:
            arcs.add((v, u) if frozenset((u, v)) in graph else (u, v))
    return arcs


def count_cycles(arcs, order):
    return sum(((u, v) in arcs and (v, w) in arcs and (w, u) in arcs)
               or ((v, u) in arcs and (w, v) in arcs and (u, w) in arcs)
               for u, v, w in combinations(order, 3))


def cyclic_triangles(graph, order):
    # Independent literal arcs rather than the main checker's bit formula.
    return count_cycles(literal_arcs(graph, order), order)


def scores(arcs, vertices):
    return {u: sum((u, v) in arcs for v in vertices if v != u) for u in vertices}


graphs = orders = insertions = 0
neighborhood_insertions = 0
for n in range(1, 5):
    pairs = list(combinations(range(n), 2))
    for mask in range(1 << len(pairs)):
        graph = {frozenset(pairs[k]) for k in range(len(pairs)) if mask >> k & 1}
        universal = graph | {frozenset((v, n)) for v in range(n)}
        graphs += 1
        for order in permutations(range(n)):
            old_arcs = literal_arcs(graph, order)
            old = count_cycles(old_arcs, order)
            old_scores = scores(old_arcs, order)
            orders += 1
            for k in range(n + 1):
                new = order[:k] + (n,) + order[k:]
                crossing = sum(frozenset((u, v)) in graph for u in order[:k]
                               for v in order[k:])
                require(cyclic_triangles(graph, new) - old == crossing,
                        ('isolated insertion cut', n, mask, order, k))
                require(cyclic_triangles(universal, new) - old == k*(n-k)-crossing,
                        ('universal insertion cocut', n, mask, order, k))
                insertions += 1
                for neighborhood_mask in range(1 << n):
                    neighborhood = {v for v in range(n) if neighborhood_mask >> v & 1}
                    new_graph = graph | {frozenset((v, n)) for v in neighborhood}
                    new_arcs = literal_arcs(new_graph, new)
                    O = (set(order[:k]) & neighborhood) | (set(order[k:]) - neighborhood)
                    require(O == {v for v in order if (n, v) in new_arcs},
                            ('actual addressed outneighbors', n, mask, order, k, neighborhood_mask))
                    crossing_directed = sum((u, v) in old_arcs for u in O
                                            for v in set(order) - O)
                    by_scores = sum(old_scores[u] for u in O) - len(O)*(len(O)-1)//2
                    require(count_cycles(new_arcs, new) - old == crossing_directed == by_scores,
                            ('general insertion increment', n, mask, order, k, neighborhood_mask))
                    predicted_scores = {u: old_scores[u] + (u not in O) for u in order}
                    predicted_scores[n] = len(O)
                    require(scores(new_arcs, new) == predicted_scores,
                            ('owner-indexed score transition', n, mask, order, k, neighborhood_mask))
                    neighborhood_insertions += 1

empty = Counter(cyclic_triangles(set(), p) for p in permutations(range(3)))
edge = Counter(cyclic_triangles({frozenset((0, 1))}, p)
               for p in permutations(range(3)))
require(empty == {0: 6}, 'empty plus isolate polynomial')
require(edge == {0: 4, 1: 2}, 'edge plus isolate polynomial')
require(cyclic_triangles({frozenset((0, 2))}, (0, 1, 2)) == 1,
        'source selected as single outneighbor')
require(cyclic_triangles({frozenset((1, 2))}, (0, 1, 2)) == 0,
        'sink selected as single outneighbor')
print(json.dumps({'status': 'PASS', 'checks': checks, 'all_graphs_orders_through': 4,
                  'graphs': graphs, 'orders': orders, 'insertion_positions': insertions,
                  'all_neighborhood_insertions': neighborhood_insertions,
                  'minimal_indexed_score_hostile': {'old_scores': [1, 0], 'outneighbor_set_size': 1,
                                                   'increments': [1, 0]},
                  'complement_sheet_hostile': {'empty3': dict(empty), 'edge_plus_isolate': dict(edge)}}))
