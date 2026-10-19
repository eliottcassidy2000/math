#!/usr/bin/env python3
"""Independent exact replay of the four-vertex, three-color classification.

Standard library only. No imports from classify_rank2.py or exact_utils.py.
Ranks use vertex-set reachability. Rainbow trees use precomputed six-bit edge
masks. Width-two components use DSU over coordinate-value incidence buckets.
Run from any directory: python verify_rank2_independent.py
"""
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations, permutations, product
from pathlib import Path
import json

EDGES = tuple(combinations(range(4), 2))
ALL_VERTICES = (1 << 4) - 1


def component_count(edge_mask):
    """Connected-component count by Boolean vertex-set closure, not edge DSU."""
    adjacency = [0] * 4
    for i, (a, b) in enumerate(EDGES):
        if edge_mask & (1 << i):
            adjacency[a] |= 1 << b
            adjacency[b] |= 1 << a
    unvisited = ALL_VERTICES
    count = 0
    while unvisited:
        reached = unvisited & -unvisited
        while True:
            expanded = reached
            for vertex in range(4):
                if reached & (1 << vertex):
                    expanded |= adjacency[vertex]
            if expanded == reached:
                break
            reached = expanded
        unvisited &= ~reached
        count += 1
    return count


def edge_ids(mask):
    return tuple(i for i in range(6) if mask & (1 << i))


def components_of_states(states):
    """Join precisely the states sharing at least one colored contact type."""
    parent = list(range(len(states)))
    size = [1] * len(states)

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i, j):
        i, j = find(i), find(j)
        if i != j:
            if size[i] < size[j]:
                i, j = j, i
            parent[j] = i
            size[i] += size[j]

    first = {}
    for index, state in enumerate(states):
        for color, edge in enumerate(state):
            key = color, edge
            if key in first:
                union(index, first[key])
            else:
                first[key] = index
    buckets = defaultdict(list)
    for i, state in enumerate(states):
        buckets[find(i)].append(state)
    return tuple(sorted(tuple(sorted(bucket)) for bucket in buckets.values()))


def verify():
    nonempty_masks = tuple(range(1, 1 << 6))
    ranks = {mask: 4 - component_count(mask) for mask in nonempty_masks}
    supports = tuple(mask for mask in nonempty_masks if ranks[mask] <= 2)
    # Cross-check the rank computation against an explicit combinatorial list.
    explicit = {1 << i for i in range(6)}
    explicit.update((1 << i) | (1 << j) for i, j in combinations(range(6), 2))
    for vertices in combinations(range(4), 3):
        explicit.add(sum(1 << EDGES.index(edge)
                         for edge in combinations(vertices, 2)))
    assert set(supports) == explicit
    assert len(supports) == 25
    assert Counter(mask.bit_count() for mask in supports) == {1: 6, 2: 15, 3: 4}

    tree_masks = {mask for mask in nonempty_masks
                  if mask.bit_count() == 3 and component_count(mask) == 1}
    assert len(tree_masks) == 16
    expected = set()
    expected_tree = {}
    for tree in tree_masks:
        pair_supports = tuple(tree ^ (1 << edge) for edge in edge_ids(tree))
        for triple in permutations(pair_supports):
            expected.add(triple)
            expected_tree[triple] = tree
    assert len(expected) == 96

    outcomes = Counter()
    states_histogram = Counter()
    component_histogram = Counter()
    disconnected = set()
    witnesses = []
    tree_shape_counts = Counter()
    support_edges = {mask: edge_ids(mask) for mask in supports}
    for triple in product(supports, repeat=3):
        states = tuple(state for state in product(*(support_edges[m] for m in triple))
                       if ((1 << state[0]) | (1 << state[1]) | (1 << state[2]))
                       in tree_masks)
        states_histogram[len(states)] += 1
        components = components_of_states(states)
        component_histogram[len(components)] += 1
        if not states:
            outcomes['empty'] += 1
            assert triple not in expected
        elif len(components) == 1:
            outcomes['connected_nonempty'] += 1
            assert triple not in expected
        else:
            outcomes['disconnected_nonempty'] += 1
            disconnected.add(triple)
            assert triple in expected, ('unexpected obstruction', triple, states)
            assert len(states) == 2
            assert tuple(map(len, components)) == (1, 1)
            assert all(a != b for a, b in zip(*states))
            union = triple[0] | triple[1] | triple[2]
            assert union == expected_tree[triple]
            assert union in tree_masks
            assert len(set(triple)) == 3
            assert all(mask.bit_count() == 2 for mask in triple)
            degrees = Counter(vertex for edge in edge_ids(union) for vertex in EDGES[edge])
            shape = 'star' if max(degrees.values()) == 3 else 'path'
            tree_shape_counts[shape] += 1
            witnesses.append({'support_masks': list(triple),
                              'tree_mask': union,
                              'states': [list(state) for state in states]})
    assert disconnected == expected, ('missing predicted patterns', expected - disconnected)
    assert sum(outcomes.values()) == 25 ** 3
    assert outcomes['disconnected_nonempty'] == 96
    assert tree_shape_counts == {'star': 24, 'path': 72}

    # Relate this independent result to the previously staged, compact receipt.
    directory = Path(__file__).resolve().parent
    previous_receipt = json.loads((directory / 'classification.json').read_text())
    assert previous_receipt == {
        'supports': 25,
        'ordered_support_triples': 15625,
        'disconnected_by_type': {'(3, (2, 2, 2), 2)': 96},
        'unexpected': 0,
    }
    inspected_files = ('THREE_FREE_COLOR_CLASSIFICATION.md', 'classify_rank2.py',
                       'exact_utils.py', 'classification.json')
    source_hashes = {name: sha256((directory / name).read_bytes()).hexdigest()
                     for name in inspected_files}
    return {
        'verified': True,
        'checker': 'verify_rank2_independent.py',
        'algorithm': {
            'support_rank': 'Boolean vertex-set reachability',
            'tree_test': 'precomputed connected three-edge bit masks',
            'state_components': 'DSU of shared coordinate-value incidence buckets',
            'obstruction_test': 'set equality with all tree-pair permutations',
            'imports_original_checker_or_helpers': False,
        },
        'edge_index': [list(edge) for edge in EDGES],
        'nonempty_graph_masks_checked': len(nonempty_masks),
        'rank_at_most_two_supports': len(supports),
        'support_size_counts': dict(sorted(Counter(m.bit_count() for m in supports).items())),
        'ordered_support_triples': 25 ** 3,
        'outcome_counts': dict(sorted(outcomes.items())),
        'state_count_histogram': dict(sorted(states_histogram.items())),
        'component_count_histogram': dict(sorted(component_histogram.items())),
        'labeled_tree_count': len(tree_masks),
        'expected_disconnected_patterns': len(expected),
        'disconnected_tree_shape_counts': dict(sorted(tree_shape_counts.items())),
        'all_disconnected_cases_have_two_singleton_components': True,
        'source_sha256': source_hashes,
        'disconnected_witnesses': witnesses,
    }


if __name__ == '__main__':
    print(json.dumps(verify(), sort_keys=True, indent=2))
