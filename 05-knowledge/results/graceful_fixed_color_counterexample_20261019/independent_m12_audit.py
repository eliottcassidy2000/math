"""Independent exhaustive audit of the explicit m=12 fixed-color witness.

Uses a hard-coded transcription of the mathematical statement and a disjoint-set
cycle test, not the search mapping, stored support lists, or supplied verifier.
Prints a JSON receipt; does not modify the supplied witness files.
"""

import itertools
import json
from pathlib import Path


M = 12
FIXED = {
    1: (5, 6), 2: (1, 3), 3: (8, 11), 4: (3, 7), 5: (0, 5),
    7: (4, 11), 8: (1, 9), 10: (0, 10), 12: (0, 12),
}
FREE = (6, 9, 11)


def forest_components(edges, vertex_count):
    parent = list(range(vertex_count))

    def root(v):
        while parent[v] != v:
            v = parent[v]
        return v

    acyclic = True
    for a, b in edges:
        ra, rb = root(a), root(b)
        if ra == rb:
            acyclic = False
        else:
            parent[ra] = rb
    components = {}
    for v in range(vertex_count):
        components.setdefault(root(v), []).append(v)
    return acyclic, sorted(components.values())


assert set(FIXED).isdisjoint(FREE)
assert set(FIXED) | set(FREE) == set(range(1, M + 1))
assert all(0 <= a < b <= M and b - a == d for d, (a, b) in FIXED.items())
is_forest, components = forest_components(FIXED.values(), M + 1)
assert is_forest
assert components == [[0, 5, 6, 10, 12], [1, 3, 7, 9], [2], [4, 8, 11]]
component_of = {v: chr(65 + i) for i, comp in enumerate(components) for v in comp}
choices = {d: [(a, a + d) for a in range(M + 1 - d)] for d in FREE}
contacts = {}
loops = {}
for d in FREE:
    contacts[d] = {}
    loops[d] = []
    for edge in choices[d]:
        a, b = map(component_of.__getitem__, edge)
        if a == b:
            loops[d].append(edge)
        else:
            contacts[d].setdefault(''.join(sorted((a, b))), []).append(edge)
assert {d: sorted(types) for d, types in contacts.items()} == {
    6: ['AD', 'CD'], 9: ['AB', 'CD'], 11: ['AB', 'AD']
}

states = []
tested = 0
rejected = 0
for candidate in itertools.product(*(choices[d] for d in FREE)):
    tested += 1
    all_edges = list(FIXED.values()) + list(candidate)
    assert len(all_edges) == M
    assert sorted(abs(a - b) for a, b in all_edges) == list(range(1, M + 1))
    acyclic, parts = forest_components(all_edges, M + 1)
    if acyclic and len(parts) == 1:
        states.append(candidate)
    else:
        rejected += 1
assert tested == 7 * 4 * 2 == 56
expected = [
    ((2, 8), (0, 9), (0, 11)),
    ((2, 8), (1, 10), (0, 11)),
    ((2, 8), (3, 12), (0, 11)),
    ((4, 10), (2, 11), (1, 12)),
    ((5, 11), (2, 11), (1, 12)),
]
assert states == expected
distances = [[sum(e != f for e, f in zip(s, t)) for t in states] for s in states]
width_two_edges = [(i, j) for i in range(len(states)) for j in range(i + 1, len(states))
                   if distances[i][j] <= 2]
assert width_two_edges == [(0, 1), (0, 2), (1, 2), (3, 4)]
_, state_components = forest_components(width_two_edges, len(states))
assert state_components == [[0, 1, 2], [3, 4]]
assert all(distances[i][j] == 3 for i in range(3) for j in (3, 4))

# Compare the independently transcribed witness to the stored data only after
# the mathematical assertions above have passed.
stored = json.loads(Path(__file__).with_name('witness.json').read_text())
assert stored['m'] == M and tuple(stored['free']) == FREE
assert {int(d): tuple(e) for d, e in stored['fixed_forest'].items()} == FIXED

print(json.dumps({
    'status': 'PASS',
    'method': 'Hard-coded witness; exhaustive choices; disjoint-set cycle/connectivity test',
    'm': M,
    'free_colors': FREE,
    'fixed_components': dict(zip('ABCD', components)),
    'nonloop_realizations_by_color_and_type': contacts,
    'contracted_loops_by_color': loops,
    'candidate_count': tested,
    'rejected_candidate_count': rejected,
    'graceful_tree_count': len(states),
    'states_in_free_color_order': states,
    'color_hamming_distance_matrix': distances,
    'width_two_simple_graph_edges': width_two_edges,
    'width_two_components': state_components,
    'width_two_component_sizes': list(map(len, state_components)),
    'width_two_graph': 'K3 disjoint union K2',
    'stored_witness_matches_transcription': True,
}, indent=2, sort_keys=True))
