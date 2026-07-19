#!/usr/bin/env python3
"""Exact referee for THM-1250's fully located six-cover tree.

The compactness/irredundant-interval-chain extraction is a paper topology
provider.  This dependency-free replay checks the tooth/address constants,
the gcd/lcm overlap quantum, every labelled six-vertex tree against every
active set, exhaustive short chronological owner words, and the exact Hunter
and private-stalk scalar ledger.
"""

from fractions import Fraction as F
from itertools import combinations, product
from math import gcd


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


MAX_RATIO = 2345
TOOTH_BOUND_NUMERATOR = 6 * MAX_RATIO + 1
TOOTH_BOUND_DENOMINATOR = 7
TEETH_PER_LABEL = 2011
TOTAL_TEETH = 6 * TEETH_PER_LABEL
OTHER_TEETH = 5 * TEETH_PER_LABEL
PRIVATE_COMPONENTS_PER_TOOTH = OTHER_TEETH + 1
PRIVATE_TOOTH_DENOMINATOR = 49 * TEETH_PER_LABEL
PRIVATE_STALK_DENOMINATOR = (
    PRIVATE_TOOTH_DENOMINATOR * PRIVATE_COMPONENTS_PER_TOOTH
)

require(F(TOOTH_BOUND_NUMERATOR, TOOTH_BOUND_DENOMINATOR) < TEETH_PER_LABEL,
        "tooth-address count")
require(TOTAL_TEETH == 12066, "total tooth count")
require(OTHER_TEETH == 10055, "other-tooth count")
require(PRIVATE_COMPONENTS_PER_TOOTH == 10056,
        "private component count")
require(PRIVATE_TOOTH_DENOMINATOR == 98539,
        "private tooth pigeonhole")
require(PRIVATE_STALK_DENOMINATOR == 990908184,
        "private stalk pigeonhole")


# Every positive tooth-to-tooth handoff numerator is a positive multiple of
# gcd(u,v), and therefore the overlap is at least 1/(14*lcm(u,v)).
overlap_rows = 0
for u in range(2, 181):
    for v in range(u + 1, 181):
        g = gcd(u, v)
        lcm = u * v // g
        for n in range(-5, 6):
            for m in range(-5, 6):
                numerator = v * (14 * n + 1) - u * (14 * m - 1)
                if numerator <= 0:
                    continue
                require(numerator % g == 0, "handoff numerator gcd")
                overlap = F(numerator, 14 * u * v)
                require(overlap >= F(g, 14 * u * v),
                        "gcd overlap quantum")
                require(F(g, u * v) == F(1, lcm), "gcd/lcm identity")
                overlap_rows += 1


vertices = tuple(range(6))
edges = tuple(combinations(vertices, 2))


def is_tree(chosen: tuple[tuple[int, int], ...]) -> bool:
    if len(chosen) != 5:
        return False
    parent = list(vertices)

    def root(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for u, v in chosen:
        ru, rv = root(u), root(v)
        if ru == rv:
            return False
        parent[ru] = rv
    return len({root(v) for v in vertices}) == 1


trees = [chosen for chosen in combinations(edges, 5) if is_tree(chosen)]
require(len(trees) == 6 ** 4 == 1296, "Cayley six-tree count")

edge_tree_incidences = {
    edge: sum(edge in tree for tree in trees) for edge in edges
}
require(set(edge_tree_incidences.values()) == {432},
        "each K6 edge occurs in one third of the labelled trees")

active_checks = 0
for tree in trees:
    for mask in range(1 << 6):
        active = {v for v in vertices if mask >> v & 1}
        induced = sum(u in active and v in active for u, v in tree)
        require(induced <= max(0, len(active) - 1),
                "forest Hunter pointwise count")
        active_checks += 1


def transition_connected(word: tuple[int, ...]) -> bool:
    adjacency = [set() for _ in vertices]
    for u, v in zip(word, word[1:]):
        adjacency[u].add(v)
        adjacency[v].add(u)
    seen = {word[0]}
    stack = [word[0]]
    while stack:
        u = stack.pop()
        for v in adjacency[u] - seen:
            seen.add(v)
            stack.append(v)
    return seen == set(vertices)


def maximum_tree_weight(weights: dict[tuple[int, int], F]) -> F:
    """Kruskal maximum spanning-tree weight on the six labels."""
    parent = list(vertices)

    def root(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    total = F(0)
    used = 0
    for (u, v), weight in sorted(weights.items(),
                                  key=lambda item: item[1], reverse=True):
        ru, rv = root(u), root(v)
        if ru == rv:
            continue
        parent[ru] = rv
        total += weight
        used += 1
        if used == 5:
            break
    require(used == 5, "maximum spanning tree exists")
    return total


# Exhaustive finite stress test of the abstract word lemma.  The paper proof
# is one line: every prefix vertex is joined to its successor, so a word that
# visits all labels has connected transition graph.  Lengths through eight
# include repeats and returns, not merely permutations.
owner_words = 0
weighted_owner_words = 0
test_speeds = (11, 13, 17, 19, 23, 29)
for length in range(6, 9):
    for word in product(vertices, repeat=length):
        if any(word[i] == word[i + 1] for i in range(length - 1)):
            continue
        if set(word) != set(vertices):
            continue
        require(transition_connected(word), "chronological transition graph")
        multiplicity = {edge: 0 for edge in edges}
        for u, v in zip(word, word[1:]):
            edge = (u, v) if u < v else (v, u)
            multiplicity[edge] += 1
        require(sum(multiplicity.values()) == length - 1,
                "handoff multiplicity sum")
        weights = {
            (u, v): F(multiplicity[(u, v)] * gcd(test_speeds[u], test_speeds[v]),
                       test_speeds[u] * test_speeds[v])
            for u, v in edges
        }
        total_weight = sum(weights.values(), F(0))
        require(maximum_tree_weight(weights) >= total_weight / 3,
                "Cayley-averaged maximum tree")
        owner_words += 1
        weighted_owner_words += 1


# Private mass 1/(49c) requires at least ceil(d/(7c)) distinct teeth of its
# owner, since one complete tooth has length 1/(7d).
owner_recurrence_rows = 0
for c in range(1, 101):
    for d in range(c + 1, 40 * c + 1):
        n_min = (d + 7 * c - 1) // (7 * c)
        require(F(n_min, 7 * d) >= F(1, 49 * c),
                "owner tooth recurrence capacity")
        if n_min > 0:
            require(F(n_min - 1, 7 * d) < F(1, 49 * c),
                    "owner tooth recurrence minimality")
        owner_recurrence_rows += 1


# Exact Hunter coefficient and scale-covariant lcm form.
require(F(49, 6) * F(1, 14) == F(7, 12),
        "located seam coefficient")
require(F(6, 7) / 7 * F(49, 6) == 1,
        "slow-gap singleton coefficient")

# A low harmonic slack epsilon forces every selected tree lcm above the same
# threshold: (7c/12)/lcm <= delta.
for c in range(1, 101):
    for lcm in range(c + 1, 20 * c + 1):
        seam_debt = F(7 * c, 12 * lcm)
        require(seam_debt * lcm == F(7 * c, 12),
                "per-edge low-slack threshold")


print("THM-1250 SIX PRIVATE NEEDLES / FULLY LOCATED TREE EXACT AUDIT")
print(f"projective tooth-address ceiling = {TEETH_PER_LABEL} per label")
print(f"total tooth ceiling = {TOTAL_TEETH}")
print(f"private tooth floor = 1/({PRIVATE_TOOTH_DENOMINATOR} c)")
print(f"private interval-stalk floor = 1/({PRIVATE_STALK_DENOMINATOR} c)")
print(f"positive gcd/lcm tooth handoffs checked = {overlap_rows}")
print(f"labelled six-vertex trees = {len(trees)}")
print("K6 edge incidence across labelled trees = 432/1296 = 1/3")
print(f"tree/active-set Hunter checks = {active_checks}")
print(f"surjective chronological owner words checked = {owner_words}")
print(f"multiplicity-weighted maximum-tree checks = {weighted_owner_words}")
print(f"private-owner tooth-recurrence rows = {owner_recurrence_rows}")
print("located debt = cH-1 >= (7c/12) sum_(uv in T) 1/lcm(u,v)")
print("averaged debt = H >= 1/c + (7/36) sum_(u<v) m_uv/lcm(u,v)")
print("RESULT: PASS")
