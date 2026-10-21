"""Independent exact controls; never import the candidate checker."""
from hashlib import sha256
from itertools import combinations, permutations
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / 'session23b-duality'
EXPECTED = {
    'PROOF.md': 'af297f77570dc82643e1453eaccc3f7b6976001e4dc016996a07a01d96d43729',
    'check.py': 'bc531f4099b72dce0671c692158b06d46870a8175d9a7ee68cba8e25af12f695',
}
checks = 0
def require(condition, message):
    global checks
    checks += 1
    if not condition:
        raise RuntimeError(message)

actual = {name: sha256((SOURCE / name).read_bytes()).hexdigest() for name in EXPECTED}
require(actual == EXPECTED, 'Reviewed source drift: review must be rebound explicitly')
vertices = range(4)
edges = tuple(combinations(vertices, 2))
edge_index = {edge: i for i, edge in enumerate(edges)}
def edge_mask(es):
    return sum(1 << edge_index[tuple(sorted(e))] for e in es)
def support(mask):
    return [e for i, e in enumerate(edges) if (mask >> i) & 1]
def even_degrees(mask):
    return all(sum(v in e for e in support(mask)) % 2 == 0 for v in vertices)

# Derive the cycle space directly from all 64 edge subsets, not from faces/quads.
space = tuple(mask for mask in range(64) if even_degrees(mask))
cycles = tuple(mask for mask in space if mask)
require(len(space) == 8, 'Dimension-three binary cycle space')
require(sorted(mask.bit_count() for mask in cycles) == [3,3,3,3,4,4,4], 'Cycle types')
for a in space:
    for b in space:
        require(a ^ b in space, 'Cycle-space closure')

# Construct edge points from affine-tetrahedron differences, and all kernels.
x = (1, 2, 4, 0)
omega = 7
coordinates = {e: x[e[0]] ^ x[e[1]] for e in edges}
points = frozenset(range(1,8))
dot = lambda a, b: (a & b).bit_count() & 1
kernels = {a: frozenset(v for v in points if dot(a, v) == 0) for a in points}
require(len(set(kernels.values())) == 7, 'Seven distinct edge-plane lines')
for u, v in combinations(points, 2):
    require(sum({u,v} <= L for L in kernels.values()) == 1, 'Unique edge-plane join')
def dual(mask):
    if mask.bit_count() == 3:
        return frozenset(coordinates[e] for e in support(mask))
    return frozenset({omega} | {coordinates[e] for e in edges if e not in support(mask)})
require(set(map(dual, cycles)) == set(kernels.values()), 'Combinatorial dual bijection')

cycle_lines = set()
for a, b in combinations(cycles, 2):
    cycle_lines.add(frozenset((a,b,a ^ b)))
require(len(cycle_lines) == 7, 'Seven cycle-plane lines')
pencil_point = {}
for triple in combinations(cycles, 3):
    meet = set.intersection(*(set(dual(c)) for c in triple))
    require((frozenset(triple) in cycle_lines) == (len(meet) == 1), 'All 35 triples: dual incidence')
    if frozenset(triple) in cycle_lines:
        pencil_point[frozenset(triple)] = next(iter(meet))
require(set(pencil_point.values()) == points, 'Pencils are bijective with edge points')
for C in cycles:
    for line, point in pencil_point.items():
        require((C in line) == (point in dual(C)), 'All 49 flags: incidence reversal')

words = {a: edge_mask(e for e in edges if dot(a, coordinates[e] ^ omega)) for a in range(8)}
require(set(words.values()) == set(space), 'Coordinate cycle-space isomorphism')
for a in range(8):
    for b in range(8):
        require(words[a] ^ words[b] == words[a ^ b], 'All 64 coordinate linearity equations')
for a in points:
    require(dual(words[a]) == kernels[a], 'Coordinate and combinatorial dualities agree')
    require((words[a].bit_count() == 3) == bool(dot(a,omega)), 'Face/quad coordinate partition')

# Compute every induced linear edge action and its CONTRAGREDIENT dual action.
actions = []
wrong_standard_action_failures = 0
for p in permutations(vertices):
    def act_edge(e):
        return tuple(sorted((p[e[0]],p[e[1]])))
    def act_cycle(mask):
        return edge_mask(act_edge(e) for e in support(mask))
    cols = [coordinates[act_edge((i,3))] for i in range(3)]
    def T(v):
        result = 0
        for i in range(3):
            if (v >> i) & 1:
                result ^= cols[i]
        return result
    require({T(v) for v in range(8)} == set(range(8)), 'Edge action invertible')
    require(T(omega) == omega, 'Natural action fixes omega')
    for e in edges:
        require(T(coordinates[e]) == coordinates[act_edge(e)], 'Edge coordinate action')
    point_action = {v: T(v) for v in points}
    cycle_action = {C: act_cycle(C) for C in cycles}
    actions.append((point_action, cycle_action))
    for a in points:
        dual_candidates = [b for b in points if all(dot(b,T(v)) == dot(a,v) for v in range(8))]
        require(len(dual_candidates) == 1, 'Unique contragredient covector')
        b = dual_candidates[0]
        require(words[b] == act_cycle(words[a]), 'Coordinate cycle action is contragredient')
        require(kernels[b] == frozenset(T(v) for v in kernels[a]), 'Coordinate kernel action')
        require(dual(act_cycle(words[a])) == frozenset(T(v) for v in dual(words[a])), 'Duality equivariance')
        wrong_standard_action_failures += words[T(a)] != act_cycle(words[a])
    for line, point in pencil_point.items():
        require(pencil_point[frozenset(act_cycle(C) for C in line)] == T(point), 'Pencil map equivariance')

def orbit_sizes(domain, which):
    unseen = set(domain)
    sizes = []
    while unseen:
        seed = min(unseen)
        orbit = {action[which][seed] for action in actions}
        unseen -= orbit
        sizes.append(len(orbit))
    return sorted(sizes)
edge_orbits = orbit_sizes(points, 0)
cycle_orbits = orbit_sizes(cycles, 1)
require(edge_orbits == [1,6] and cycle_orbits == [3,4], 'Independently computed point-orbit obstruction')
require(wrong_standard_action_failures > 0, 'Hostile: do not use ordinary rather than dual covector action')

# Essential translation/complement negative controls.
untranslated = {a: edge_mask(e for e in edges if dot(a, coordinates[e])) for a in points}
require(sum(not even_degrees(mask) for mask in untranslated.values()) == 4, 'Hostile: omitting omega translation produces four noncycles')
quad_cycles = [C for C in cycles if C.bit_count() == 4]
require(all(frozenset(coordinates[e] for e in support(C)) not in kernels.values() for C in quad_cycles), 'Hostile: quad supports themselves are never edge-plane lines')

# Scope control for the one tournament sentence in the bridge, not flow theory.
strong_count = 0
for orientation in range(64):
    arcs = {e if (orientation >> i) & 1 else e[::-1] for i,e in enumerate(edges)}
    def reaches_all(start):
        seen = {start}
        while True:
            more = seen | {v for u,v in arcs if u in seen}
            if more == seen:
                return len(seen) == 4
            seen = more
    if not all(reaches_all(v) for v in vertices):
        continue
    strong_count += 1
    directed_supports = set()
    for length in (3,4):
        for ordered in permutations(vertices, length):
            cycle_arcs = [(ordered[i],ordered[(i+1)%length]) for i in range(length)]
            if all(arc in arcs for arc in cycle_arcs):
                directed_supports.add(edge_mask(cycle_arcs))
    require(len(directed_supports) == 3, 'Strong tournament has three directed-cycle supports')
    require(frozenset(directed_supports) not in cycle_lines, 'Strong directed supports noncollinear')
    require(not set.intersection(*(set(dual(C)) for C in directed_supports)), 'Corresponding three dual lines nonconcurrent')
require(strong_count == 24, 'All strong order-four tournaments controlled')

print(json.dumps({
    'status': 'PASS',
    'checks': checks,
    'source_hashes': actual,
    'independent_cycle_space_size': len(space),
    'cycle_line_count': len(cycle_lines),
    'triple_incidence_cases': 35,
    'flag_incidence_cases': 49,
    'tetrahedral_permutations': len(actions),
    'computed_edge_point_orbits': edge_orbits,
    'computed_cycle_point_orbits': cycle_orbits,
    'wrong_standard_covector_action_failures': wrong_standard_action_failures,
    'strong_tournaments_checked': strong_count,
}, indent=2, sort_keys=True))
