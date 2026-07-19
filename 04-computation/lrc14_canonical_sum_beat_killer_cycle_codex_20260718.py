#!/usr/bin/env python3
"""Exact replay for THM-1193's canonical sum-beat killer cycle."""

from fractions import Fraction as F
from itertools import combinations, permutations


def frac_dist(x: F) -> F:
    r = x % 1
    return min(r, 1 - r)


def nearest_tie_up(num: int, den: int) -> int:
    """Nearest integer to num/den, with a half-integer tie rounded upward."""
    return (2 * num + den) // (2 * den)


def canonical_point(a: int, b: int, k: int) -> tuple[int, int, F]:
    q = a + b
    p = nearest_tie_up(q * (2 * k + 1), 2 * a)
    return p, q, F(p, q)


def slow_gap(a: int, k: int) -> tuple[F, F]:
    return F(14 * k + 1, 14 * a), F(14 * k + 13, 14 * a)


def closed_teeth(a: int, lo: F, hi: F) -> list[tuple[F, F]]:
    radius = F(1, 14 * a)
    first = (a * lo).numerator // (a * lo).denominator - 2
    last = (a * hi).numerator // (a * hi).denominator + 2
    ans = []
    for m in range(first, last + 1):
        left = max(lo, F(m, a) - radius)
        right = min(hi, F(m, a) + radius)
        if left <= right:
            ans.append((left, right))
    return ans


def survivor(a: int, bs: tuple[int, ...], k: int) -> F:
    lo, hi = slow_gap(a, k)
    pieces = sorted(piece for b in bs for piece in closed_teeth(b, lo, hi))
    covered = F(0)
    left, right = pieces[0]
    for next_left, next_right in pieces[1:]:
        if next_left <= right:
            right = max(right, next_right)
        else:
            covered += right - left
            left, right = next_left, next_right
    covered += right - left
    return hi - lo - covered


print("THM-1193 exact replay: canonical sum-beat killer cycle")

# Lemma 1 over a substantial deterministic exact bank.
rows = 0
minimum_margin = None
for a in range(1, 31):
    for b in range(a + 1, 12 * a + 1):
        for k in range(a):
            p, q, t = canonical_point(a, b, k)
            lo, hi = slow_gap(a, k)
            assert lo < t < hi
            da = frac_dist(a * t)
            db = frac_dist(b * t)
            assert da == db
            assert da >= F(b, 2 * q) > F(1, 4) > F(1, 14)
            margin = da - F(1, 4)
            minimum_margin = margin if minimum_margin is None else min(minimum_margin, margin)
            rows += 1
print(f"canonical interior-beat checks={rows}")
print(f"minimum observed margin above 1/4={minimum_margin}")

# Exact guardrail: all six canonical obligations fire but the interval survives.
a = 3
k = 1
bs = (4, 5, 6, 7, 9, 11)
points = tuple(canonical_point(a, b, k) for b in bs)
assert tuple((p, q) for p, q, _ in points) == (
    (4, 7), (4, 8), (5, 9), (5, 10), (6, 12), (7, 14)
)
killers = []
for i, (_, _, t) in enumerate(points):
    row = tuple(j for j, b in enumerate(bs) if j != i and frac_dist(b * t) < F(1, 14))
    assert row
    killers.append(row)
killers = tuple(killers)
assert killers == ((3,), (0, 2), (4,), (0, 2), (0, 2), (0, 2))
least_map = tuple(row[0] for row in killers)
assert least_map[0] == 3 and least_map[3] == 0
remaining = survivor(a, bs, k)
assert remaining == F(879, 10780) > 0
print("guardrail a=3 k=1 bs=4,5,6,7,9,11")
print("canonical (p,q)=" + str(tuple((p, q) for p, q, _ in points)))
print("killer rows=" + str(killers) + "; least-map cycle=0->3->0")
print(f"closed-tooth survivor={remaining}")

# Tournament quotient from reciprocal canonical-beat clearances.
clearance = [
    [F(1, 14) - frac_dist(bs[j] * points[i][2]) for j in range(6)]
    for i in range(6)
]
arcs = set()
ties = []
for i, j in combinations(range(6), 2):
    observable = clearance[i][j] - clearance[j][i]
    if observable > 0:
        arcs.add((i, j))
    elif observable < 0:
        arcs.add((j, i))
    else:
        arcs.add((i, j))
        ties.append((i, j))
scores = tuple(sum((i, j) in arcs for j in range(6) if i != j) for i in range(6))
triangles = sum(
    ((i, j) in arcs and (j, ell) in arcs and (ell, i) in arcs)
    or ((i, ell) in arcs and (ell, j) in arcs and (j, i) in arcs)
    for i, j, ell in combinations(range(6), 3)
)
paths = tuple(
    order
    for order in permutations(range(6))
    if all((order[z], order[z + 1]) in arcs for z in range(5))
)
reach = [[i == j or (i, j) in arcs for j in range(6)] for i in range(6)]
for mid in range(6):
    for i in range(6):
        for j in range(6):
            reach[i][j] = reach[i][j] or (reach[i][mid] and reach[mid][j])
seen = set()
sccs = []
for i in range(6):
    if i not in seen:
        component = {j for j in range(6) if reach[i][j] and reach[j][i]}
        seen.update(component)
        sccs.append(tuple(sorted(component)))
assert tuple(sorted(scores)) == (1, 2, 2, 2, 3, 5)
assert triangles == 4
assert tuple(sorted(map(len, sccs), reverse=True)) == (5, 1)
assert len(paths) == 15 and paths[0] == (1, 0, 3, 2, 4, 5)
print("tournament observable=kappa_ij-kappa_ji; zero gauge=increasing label")
print("score multiset=1,2,2,2,3,5; directed triangles=4; SCC sizes=5,1")
print("Hamiltonian paths=15; tie Hamiltonian path=1,0,3,2,4,5")
print(f"observable ties={len(ties)}")
print("ALL EXACT CHECKS PASSED")
