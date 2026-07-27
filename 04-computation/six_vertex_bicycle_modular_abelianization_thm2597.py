#!/usr/bin/env python3
"""Exact, dependency-free companion for THM-2597."""

from collections import deque
from itertools import combinations, permutations


def require(condition, message):
    if not condition:
        raise RuntimeError("CHECK FAILED: " + message)


V = tuple(range(1, 7))
ALL_EDGES = tuple(combinations(V, 2))
H_EDGES = frozenset(e for e in ALL_EDGES if e[1] - e[0] >= 2)
TILES = ((6, 1), (5, 1), (4, 1), (3, 1), (6, 2),
         (5, 2), (4, 2), (6, 3), (5, 3), (6, 4))
TILE_EDGES = tuple((y, x) for x, y in TILES)
require(frozenset(TILE_EDGES) == H_EDGES, "repo tile order")


def cut(s):
    s = set(s)
    return frozenset(e for e in H_EDGES if ((e[0] in s) != (e[1] in s)))


def even_degrees(es):
    return all(sum(v in e for e in es) % 2 == 0 for v in V)


nonzero_bicycles = set()
cut_witnesses = {}
for bits in range(1 << 6):
    s = frozenset(v for v in V if bits & (1 << (v - 1)))
    b = cut(s)
    if b and even_degrees(b):
        nonzero_bicycles.add(b)
        cut_witnesses.setdefault(b, []).append(s)

require(len(nonzero_bicycles) == 1, "unique nonzero n=6 bicycle")
bicycle = next(iter(nonzero_bicycles))
EXPECTED = frozenset({(1, 3), (3, 5), (2, 5), (2, 4), (4, 6), (1, 6)})
require(bicycle == EXPECTED, "bicycle support")
require({frozenset({1, 4, 5}), frozenset({2, 3, 6})}.issubset(
    set(cut_witnesses[bicycle])), "complementary cut witnesses")
require(all(sum(v in e for e in bicycle) == 2 for v in V), "support is 2-regular")
print("unique nonzero bicycle: delta({1,4,5})=delta({2,3,6})=C6")


def star(v):
    return frozenset(e for e in H_EDGES if v in e)


def xor_edges(*sets):
    out = set()
    for es in sets:
        out.symmetric_difference_update(es)
    return frozenset(out)


require(xor_edges(star(1), star(4), star(5)) == bicycle,
        "first star-flip word")
require(xor_edges(star(2), star(3), star(6)) == bicycle,
        "complement star-flip word")
print("star words: star1+star4+star5 = star2+star3+star6")


CYCLE = (1, 3, 5, 2, 4, 6)
require(frozenset(tuple(sorted((CYCLE[i], CYCLE[(i + 1) % 6])))
                  for i in range(6)) == bicycle, "displayed cycle order")


# Partial-cube realization C6 = Q3 minus an antipodal pair.
CUBE = {
    1: (0, 0, 0), 3: (1, 0, 0), 5: (1, 1, 0),
    2: (1, 1, 1), 4: (0, 1, 1), 6: (0, 0, 1),
}
require(set(CUBE.values()) == {
    (0, 0, 0), (1, 0, 0), (1, 1, 0),
    (1, 1, 1), (0, 1, 1), (0, 0, 1),
}, "cube image")


def graph_distance(source, target):
    q = deque([(source, 0)])
    seen = {source}
    while q:
        x, d = q.popleft()
        if x == target:
            return d
        for e in bicycle:
            if x not in e:
                continue
            y = e[0] if e[1] == x else e[1]
            if y not in seen:
                seen.add(y)
                q.append((y, d + 1))
    raise RuntimeError("disconnected")


for x, y in combinations(V, 2):
    hamming = sum(a != b for a, b in zip(CUBE[x], CUBE[y]))
    require(graph_distance(x, y) == hamming, "isometric cube embedding")

theta = []
for coordinate in range(3):
    cls = frozenset(e for e in bicycle
                    if CUBE[e[0]][coordinate] != CUBE[e[1]][coordinate])
    theta.append(cls)
require(theta == [
    frozenset({(1, 3), (2, 4)}),
    frozenset({(3, 5), (4, 6)}),
    frozenset({(2, 5), (1, 6)}),
], "three Theta classes")
print("partial cube: Q3 minus {010,101}; Theta classes have sizes 2,2,2")


# Perfect-matching and hafnian/Pfaffian boundary.
LEFT = (1, 4, 5)
RIGHT = (2, 3, 6)
FROZEN = frozenset({(1, 2), (3, 4), (5, 6)})
K33 = frozenset(tuple(sorted((a, b))) for a in LEFT for b in RIGHT)
require(bicycle == K33 - FROZEN, "K3,3 minus frozen matching")


def perfect_matchings(vertices, edges):
    vertices = tuple(sorted(vertices))
    edges = set(edges)
    if not vertices:
        return (frozenset(),)
    a = vertices[0]
    out = []
    for b in vertices[1:]:
        e = tuple(sorted((a, b)))
        if e not in edges:
            continue
        rest = tuple(v for v in vertices if v not in (a, b))
        for m in perfect_matchings(rest, edges):
            out.append(frozenset(set(m) | {e}))
    return tuple(out)


cycle_matchings = set(perfect_matchings(V, bicycle))
M_PLUS = frozenset({(1, 3), (2, 5), (4, 6)})
M_MINUS = frozenset({(3, 5), (2, 4), (1, 6)})
require(cycle_matchings == {M_PLUS, M_MINUS}, "two C6 matchings")
require(len(perfect_matchings(V, K33)) == 6, "six K3,3 matchings")
require(M_PLUS.isdisjoint(M_MINUS) and M_PLUS | M_MINUS == bicycle,
        "alternating matching decomposition")
require(all(len(m & t) == 1 for m in (M_PLUS, M_MINUS) for t in theta),
        "each matching crosses every Theta class once")
print("matching boundary: C6 has 2 derangement matchings; + frozen matching gives K3,3 with 6")
print("two monochromatic edge colors realize the d=2,n=6 inherited-color construction")


# The rotations realize C2 x C3 = C6 on the support cycle.
def permutation_from_shift(k):
    return tuple(CYCLE[(CYCLE.index(v) + k) % 6] for v in V)


def compose(p, q):
    return tuple(p[q[v - 1] - 1] for v in V)


def ppower(p, n):
    out = V
    for _ in range(n):
        out = compose(p, out)
    return out


rho2 = permutation_from_shift(3)
rho3 = permutation_from_shift(2)
require(rho2 == (2, 1, 4, 3, 6, 5), "order-two rotation")
require(rho3 == (5, 6, 2, 1, 4, 3), "order-three rotation")
require(ppower(rho2, 2) == V and rho2 != V, "rho2 order 2")
require(ppower(rho3, 3) == V and rho3 != V, "rho3 order 3")
require(compose(rho2, rho3) == compose(rho3, rho2), "commuting rotations")
generated = {compose(ppower(rho2, a), ppower(rho3, b))
             for a in range(2) for b in range(3)}
require(len(generated) == 6, "C2 x C3 support action")


def permute_edge(e, p):
    return tuple(sorted((p[e[0] - 1], p[e[1] - 1])))


for p in generated:
    require(frozenset(permute_edge(e, p) for e in bicycle) == bicycle,
            "support rotation")
for cls in theta:
    require(frozenset(permute_edge(e, rho2) for e in cls) == cls,
            "rho2 preserves each Theta class")
    require(all(permute_edge(e, rho2) != e for e in cls),
            "rho2 swaps the two edges in each Theta class")
rho3_theta = [frozenset(permute_edge(e, rho3) for e in cls) for cls in theta]
require(set(rho3_theta) == set(theta) and all(rho3_theta[i] != theta[i] for i in range(3)),
        "rho3 cycles the Theta classes")
print("support rotations: C2 swaps each Theta pair; C3 cycles the three classes; group C6")


# Tournament meaning in the fixed-path tile convention.
def tournament(mask):
    arcs = set()
    tile_index = {e: i for i, e in enumerate(TILE_EDGES)}
    for i, j in ALL_EDGES:
        if (i, j) in tile_index and mask & (1 << tile_index[(i, j)]):
            arcs.add((i, j))
        else:
            arcs.add((j, i))
    return frozenset(arcs)


def relabel(t, p):
    return frozenset((p[a - 1], p[b - 1]) for a, b in t)


def encode(t):
    return sum((1 << k) for k, (i, j) in enumerate(ALL_EDGES) if (i, j) in t)


PERMS = tuple(permutations(V))


def canonical(t):
    return min(encode(relabel(t, p)) for p in PERMS)


mask = sum(1 << TILE_EDGES.index(e) for e in bicycle)
require(mask == 873 and (mask ^ 1023) == 150, "tile/complement masks")
tb = tournament(mask)
tc = tournament(mask ^ 1023)
require(canonical(tb) == canonical(tc), "bicycle mask is self-loop endpoint")

gs_map = tuple(TILES.index((7 - y, 7 - x)) for x, y in TILES)


def grid_symmetric(x):
    return all(((x >> i) & 1) == ((x >> gs_map[i]) & 1)
               for i in range(10))


small_endpoints = []
for x in range(1 << 10):
    y = x ^ 1023
    if x < y and grid_symmetric(x) and canonical(tournament(x)) == canonical(tournament(y)):
        small_endpoints.append(x)
require(small_endpoints == [19, 150], "two n=6 blue self-loop lines")

scores = tuple(sum(a == v for a, _ in tb) for v in V)
require(scores == (2, 3, 2, 3, 2, 3), "bicycle tournament scores")
hamilton = sum(all((p[i], p[i + 1]) in tb for i in range(5)) for p in PERMS)
require(hamilton == 45, "Hamiltonian path count")
auts = tuple(p for p in PERMS if relabel(tb, p) == tb)
aut_generator = (3, 4, 5, 6, 1, 2)  # (1 3 5)(2 4 6)
require(len(auts) == 3 and aut_generator in auts, "tournament automorphism C3")
require(rho2 not in auts and rho3 not in auts, "support rotations are not tournament automorphisms")
iso_to_complement = (1, 3, 4, 6, 2, 5)  # (2 3 4 6 5), inverse of (2 5 6 4 3)
require(relabel(tb, iso_to_complement) == tc, "explicit complement isomorphism")
print("tournament: tile mask 873/complement 150, scores (2,3)^3, H=45, Aut=C3")
print("blue self-loop smaller endpoints at n=6: [19, 150]")


# Graceful-labeling hostile: partial cube does not imply graceful.
graceful_cycle = []
for labels in permutations(range(7), 6):
    lab = dict(zip(V, labels))
    if {abs(lab[a] - lab[b]) for a, b in bicycle} == set(range(1, 7)):
        graceful_cycle.append(labels)
require(not graceful_cycle, "C6 is not graceful")
path_labels = dict(zip(CYCLE, (0, 5, 1, 4, 2, 3)))
path_edges = frozenset(tuple(sorted((CYCLE[i], CYCLE[i + 1]))) for i in range(5))
require({abs(path_labels[a] - path_labels[b]) for a, b in path_edges} == set(range(1, 6)),
        "displayed graceful P6")
print("graceful hostile: exhaustive 7P6 gives no C6 labeling; delete one edge and P6 is graceful")
print("ALL CHECKS PASSED")
