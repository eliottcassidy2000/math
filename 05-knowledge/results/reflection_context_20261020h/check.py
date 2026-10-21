"""Exact marked-reflection and palette-kernel certificates; no third-party deps."""
from itertools import combinations, permutations, product
import json
from pathlib import Path

checks = 0


def require(value, message):
    global checks
    checks += 1
    if not value:
        raise RuntimeError(message)


def relabel(edges, p):
    return frozenset(tuple(sorted((p[i], p[j]))) for i, j in edges)


def opposite(arcs):
    return frozenset((j, i) for i, j in arcs)


def lift(graph, n=5):
    return frozenset((j, i) if (i, j-1) in graph and j > i+1 else (i, j)
                     for i, j in combinations(range(n), 2))


def directed_relabel(arcs, p):
    return frozenset((p[i], p[j]) for i, j in arcs)


def strong(arcs, n):
    for start in range(n):
        seen = {start}
        while True:
            new = seen | {j for i, j in arcs if i in seen}
            if new == seen:
                break
            seen = new
        if len(seen) != n:
            return False
    return True


def encode(d):
    a, b, c = d
    return frozenset(e for e, bit in [((0, 3), a == 0), ((1, 2), a == 2),
                                      ((0, 1), b), ((2, 3), c)] if bit)


def decode(g):
    return (1 + int((1, 2) in g) - int((0, 3) in g),
            int((0, 1) in g), int((2, 3) in g))


domain = [d for d in product(range(3), range(2), range(2))
          if d not in [(0, 0, 0), (2, 1, 1)]]
C = lambda d: (2-d[0], 1-d[1], 1-d[2])
Q = lambda d: (d[0], d[2], d[1])
F = lambda d: C(Q(d))
rho = (3, 2, 1, 0)
gamma = (1, 0, 3, 2)
mask = frozenset([(0, 1), (2, 3)])
switch = lambda g: relabel(g, gamma) ^ mask
fixed = {name: sum(op(d) == d for d in domain)
         for name, op in [("C", C), ("Q", Q), ("F", F)]}
require(fixed == {"C": 0, "Q": 4, "F": 2}, "domain fixed points")
require(len({encode(d) for d in domain}) == 10, "injective encoding")
rows = []
for d in domain:
    g, a = encode(d), lift(encode(d))
    require(decode(g) == d, "source recovery")
    require(relabel(g, rho) == encode(Q(d)), "horizontal reflection")
    require(switch(g) == encode(C(d)), "complement of divisor")
    require(switch(relabel(g, rho)) == encode(F(d)), "vertical reflection")
    require(switch(switch(g)) == g, "involution")
    require(switch(relabel(g, rho)) == relabel(switch(g), rho), "commutation")
    require(lift(relabel(g, rho)) == directed_relabel(opposite(a), tuple(reversed(range(5)))),
            "literal converse and path reversal")
    actual_strong = strong(a, 5)
    require(actual_strong == (d[0] == 0 or d == (1, 1, 1)), "strong interval test")
    if actual_strong:
        require(not strong(lift(switch(g)), 5), "strongness hostile")
    paths = [p for p in permutations(range(5)) if all((p[i], p[i+1]) in a for i in range(4))]
    rows.append({"divisor_exponents": d, "graph_edges": sorted(g), "arcs": sorted(a),
                 "strong": actual_strong, "hamilton_paths": len(paths),
                 "f": d[0]-d[1]-d[2], "g": d[1]-d[2]})

# Independent complete ordered-graph universes, including the size-two degeneracy.
target_rows = []
for m in range(2, 7):
    edges = list(combinations(range(m), 2))
    universe = frozenset(edges)
    rev = tuple(reversed(range(m)))
    counts = [0, 0, 0]
    for bits in range(1 << len(edges)):
        g = frozenset(e for i, e in enumerate(edges) if (bits >> i) & 1)
        r = relabel(g, rev)
        counts[0] += g == universe-g
        counts[1] += g == r
        counts[2] += g == universe-r
        require((0, m-1) in r if (0, m-1) in g else (0, m-1) not in r,
                "fixed outside edge")
        require(g != universe-g and g != universe-r, "complement fixed-point obstruction")
        # All three operations under the marked-path trijection, independent arc route.
        if m <= 4:
            n = m+1
            a = lift(g, n)
            path_edges = frozenset((i, i+1) for i in range(n-1))
            free_flip = path_edges | opposite(a-path_edges)
            require(free_flip == lift(universe-g, n), "nonpath flip")
            all_reverse = tuple(reversed(range(n)))
            require(directed_relabel(opposite(a), all_reverse) == lift(r, n), "marked converse")
            path_only = (a-path_edges) | opposite(path_edges)
            require(directed_relabel(path_only, all_reverse) == lift(universe-r, n), "path-only flip")
    target_rows.append({"vertices": m, "graphs": 1 << len(edges), "fixed_K_R_KR": counts})

palettes = list(combinations(range(6), 3))
K = [[len(set(a) & set(b)) for b in palettes] for a in palettes]
mm = lambda a, b: [[sum(x*y for x, y in zip(row, col)) for col in zip(*b)] for row in a]
power = K
palette_rows = []
for r in range(1, 6):
    uniform = (30**r - 30*6**(r-1)) // 20
    require(20*uniform == 30**r - 30*6**(r-1), "integer uniform term")
    for i in range(20):
        require(sum(power[i]) == 30**r, "palette row sum")
        for j in range(20):
            require(power[i][j] == 6**(r-1)*K[i][j] + uniform, "palette chain formula")
            if r >= 2:
                require(power[i][j] > 0, "all separating support restored")
    palette_rows.append({"length": r, "min_unscaled_entry": min(map(min, power)),
                         "max_unscaled_entry": max(map(max, power)), "row_sum": 30**r})
    power = mm(power, K)
require(K[0][-1] == 0 and mm(K, K)[0][-1] == 36, "complementary palette hostile")

receipt = {"checks": checks, "domain_fixed_points": fixed, "marked_rows": rows,
           "target_graph_universes": target_rows, "palette_chains": palette_rows,
           "scope": "Exact finite controls for the all-order proofs in README.md; no Collatz coverage or forcing claim."}
Path(__file__).with_name("receipt.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
print(json.dumps({"checks": checks, "fixed": fixed, "palette_chains": palette_rows}))
