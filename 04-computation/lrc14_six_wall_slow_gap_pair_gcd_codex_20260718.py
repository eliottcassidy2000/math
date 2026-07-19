#!/usr/bin/env python3
"""Exact replay for THM-1179.

The proof is analytic.  This dependency-free replay checks its finite banks,
pointwise majorants, Fano incidence identities, rational constants, sharp
periodic-remainder envelope, and an endpoint-exact guardrail packet.
"""

from fractions import Fraction as F
from itertools import combinations
from math import comb, gcd


def fold14(r: int) -> int:
    r %= 14
    return r * (14 - r)


def pair_density(a: int, b: int) -> F:
    """Global measure of D_a intersect D_b at radius 1/14."""
    g = gcd(a, b)
    a //= g
    b //= g
    if a > b:
        a, b = b, a
    return F(
        4 * a * b + fold14(a + b) - fold14(b - a),
        196 * a * b,
    )


def chi7(v: int) -> int:
    while v % 7 == 0:
        v //= 7
    return 1 if v % 7 in (1, 2, 4) else -1


def q6(c: int) -> F:
    return F(c) - F(1, 3) * comb(c, 2)


def induced_density(n: int, edges: tuple[tuple[int, int], ...]) -> F:
    """max_A |E(A)|/(|A|-1), over vertex sets of size at least two."""
    edge_set = {tuple(sorted(e)) for e in edges}
    ans = F(0)
    for mask in range(1 << n):
        vertices = [i for i in range(n) if mask >> i & 1]
        if len(vertices) < 2:
            continue
        inside = sum(tuple(sorted(e)) in edge_set for e in combinations(vertices, 2))
        ans = max(ans, F(inside, len(vertices) - 1))
    return ans


def danger_pieces(s: int) -> list[tuple[F, F]]:
    """Measure-equivalent pieces of D_s in [0,1]."""
    radius = F(1, 14 * s)
    pieces: list[tuple[F, F]] = [(F(0), radius), (F(1) - radius, F(1))]
    for m in range(1, s):
        centre = F(m, s)
        pieces.append((centre - radius, centre + radius))
    pieces.sort()
    return pieces


def intersect_pieces(
    left: list[tuple[F, F]], right: list[tuple[F, F]]
) -> list[tuple[F, F]]:
    out: list[tuple[F, F]] = []
    i = j = 0
    while i < len(left) and j < len(right):
        lo = max(left[i][0], right[j][0])
        hi = min(left[i][1], right[j][1])
        if lo < hi:
            out.append((lo, hi))
        if left[i][1] < right[j][1]:
            i += 1
        else:
            j += 1
    return out


def interval_measure(pieces: list[tuple[F, F]], lo: F, hi: F) -> F:
    return sum(
        (min(hi, y) - max(lo, x) for x, y in pieces if max(lo, x) < min(hi, y)),
        F(0),
    )


def ceil_fraction(x: F) -> int:
    return -((-x.numerator) // x.denominator)


def lifted_teeth(s: int, lo: F, hi: F) -> list[tuple[int, F, F]]:
    first = (s * lo.numerator) // lo.denominator - 2
    last = ceil_fraction(s * hi) + 2
    out: list[tuple[int, F, F]] = []
    for m in range(first, last + 1):
        left = max(lo, F(14 * m - 1, 14 * s))
        right = min(hi, F(14 * m + 1, 14 * s))
        if left < right:
            out.append((m, left, right))
    return out


def lifted_pair_measure(
    left: list[tuple[int, F, F]], right: list[tuple[int, F, F]]
) -> F:
    return sum(
        (
            min(y, v) - max(x, u)
            for _, x, y in left
            for _, u, v in right
            if max(x, u) < min(y, v)
        ),
        F(0),
    )


print("THM-1179 exact replay: six-wall slow-gap pair/gcd obstruction")

# 1. The optimal six-comb quadratic majorant and graph-density variants.
q_values = tuple(q6(c) for c in range(7))
assert q_values == (F(0), F(1), F(5, 3), F(2), F(2), F(5, 3), F(1))
assert all(q6(c) >= 1 for c in range(1, 7))
k6 = tuple(combinations(range(6), 2))
path6 = tuple((i, i + 1) for i in range(5))
assert induced_density(6, k6) == 3
assert induced_density(6, path6) == 1
chi_density_rows = []
for p in range(7):
    q = 6 - p
    if p > q:
        continue
    edges = tuple(combinations(range(p), 2)) + tuple(combinations(range(p, 6), 2))
    density = induced_density(6, edges)
    assert density == F(q, 2)
    chi_density_rows.append((p, q, len(edges), density))
print("quadratic Q6 values:", ", ".join(map(str, q_values)))
print("graph densities: K6=3, path=1, chi rows=" + str(chi_density_rows))

# 2. The collective pair floor and all displayed complete-graph constants.
r6_floor = F(comb(6, 3), 24 * comb(4, 1))
assert r6_floor == F(5, 24)
union_ceiling = F(6, 7) - F(1, 3) * r6_floor
assert union_ceiling == F(397, 504)
assert union_ceiling * F(7, 6) == F(397, 432)
assert F(196, 13) * (F(5, 28) - F(18, 49)) == F(-37, 13)
# Symbolically: (196/13)(5/28-(18/49)delta)=(35-72delta)/13.
assert F(196, 13) * F(5, 28) == F(35, 13)
assert F(196, 13) * F(18, 49) == F(72, 13)
print("six-pair floor R6=5/24; union ceiling=397/504; period ratio=397/432")
print("mixed branch: a*sum(1/gij) >= (35-72*delta)/13")

# 3. Pair-density range and the new chi7 same-colour floor.
small_upper = []
for a, b in ((1, 2), (1, 3), (1, 4)):
    small_upper.append((a, b, pair_density(a, b)))
assert small_upper == [(1, 2, F(1, 14)), (1, 3, F(1, 21)), (1, 4, F(1, 28))]
# For reduced ab>=5, the folded correction is at most 49<10ab, proving rho<1/14.
assert 49 < 10 * 5

same_colour_bank = []
same_colour_minima = []
for b in range(2, 20):
    row = []
    for a in range(1, b):
        if gcd(a, b) == 1 and chi7(a) == chi7(b):
            value = pair_density(a, b)
            row.append((value, a))
            same_colour_bank.append((value, a, b))
    if row:
        value, a = min(row)
        same_colour_minima.append((b, value, a, len(row)))
assert len(same_colour_bank) == 58
assert min(value for value, _, _ in same_colour_bank) == F(1, 77)
assert [(a, b) for value, a, b in same_colour_bank if value == F(1, 77)] == [
    (1, 11),
    (2, 11),
]
tail_at_20 = F(1, 49) - F(1, 7 * 20)
assert tail_at_20 == F(13, 980) > F(1, 77)
print("same-chi7 reduced bank: 58 rows through b'=19")
print("same-chi7 floor=1/77 at reduced ratios 1:11 and 2:11")
print("tail b'>=20 starts at 13/980 > 1/77")

cross_nonseam_bank = []
for b in range(2, 17):
    for a in range(1, b):
        if gcd(a, b) == 1 and chi7(a) != chi7(b) and (a + b) % 14 != 0:
            cross_nonseam_bank.append((pair_density(a, b), a, b))
assert len(cross_nonseam_bank) == 33
assert min(value for value, _, _ in cross_nonseam_bank) == F(1, 84)
assert [(a, b) for value, a, b in cross_nonseam_bank if value == F(1, 84)] == [
    (1, 12)
]
cross_nonseam_tail = F(1, 49) - F(1, 7 * 17)
assert cross_nonseam_tail == F(10, 833) > F(1, 84)
print("nonseam cross-chi7 floor=1/84 at reduced ratio 1:12; 33-row bank, tail 10/833")

# Same-colour graph corollary constants.
for p, q, edge_count, _ in chi_density_rows:
    assert F(196, 13) * F(6 * edge_count, 539) == F(24 * edge_count, 143)
    assert F(196, 13) * F(3 * q, 49) == F(12 * q, 13)
assert F(196, 13) * F(1, 98) == F(2, 13)
assert F(196, 13) * F(6, 49) == F(24, 13)
print("chi7 branch: a*sum_same(1/gij) >= 24E/143-(12q/13)*delta")

# 4. Sharp periodic remainder formula, plus direct interval geometry.
envelope_checks = 0
for density in (F(1, 91), F(1, 77), F(1, 49), F(1, 14)):
    for j in range(1001):
        r = F(j, 1000)
        universal_lower = max(F(0), r - (1 - density))
        affine_lower = density * r - density * (1 - density)
        assert universal_lower >= affine_lower
        envelope_checks += 1

pieces = {s: danger_pieces(s) for s in range(1, 29)}
intervals = set()
for i in range(0, 25):
    for width in (1, 2, 5, 13, 29):
        j = i + width
        if j <= 84:
            intervals.add((F(i, 84), F(j, 84)))
intervals.update(
    {
        (F(1, 196), F(13, 196)),
        (F(29, 196), F(41, 196)),
        (F(5, 28), F(3, 7)),
        (F(0), F(1)),
    }
)
single_interval_checks = 0
pair_mean_checks = 0
pair_interval_checks = 0
for s in range(1, 29):
    for lo, hi in intervals:
        value = interval_measure(pieces[s], lo, hi)
        assert value <= (hi - lo) / 7 + F(6, 49 * s)
        single_interval_checks += 1
for a, b in combinations(range(1, 29), 2):
    overlap = intersect_pieces(pieces[a], pieces[b])
    rho = pair_density(a, b)
    assert interval_measure(overlap, F(0), F(1)) == rho
    pair_mean_checks += 1
    eta = rho * (1 - rho)
    g = gcd(a, b)
    for lo, hi in intervals:
        value = interval_measure(overlap, lo, hi)
        assert value >= (hi - lo) * rho - eta / g
        pair_interval_checks += 1
print(
    "periodic discrepancy checks:",
    envelope_checks,
    "envelopes,",
    single_interval_checks,
    "single intervals,",
    pair_mean_checks,
    "pair means,",
    pair_interval_checks,
    "pair intervals",
)

# 5. Carrier-rooted Fano identity and its exact line-period constants.
# Nonzero F_2^3 vectors are 1,...,7; a line is {x,y,x xor y}.
fano_lines = sorted(
    {
        tuple(sorted((x, y, x ^ y)))
        for x in range(1, 8)
        for y in range(x + 1, 8)
        if x ^ y
    }
)
assert len(fano_lines) == 7
carrier = 1
through = [line for line in fano_lines if carrier in line]
off = [line for line in fano_lines if carrier not in line]
assert (len(through), len(off)) == (3, 4)
killers = tuple(v for v in range(1, 8) if v != carrier)
fano_mask_checks = 0
for mask in range(1 << 6):
    active = {killers[i] for i in range(6) if mask >> i & 1}
    c = len(active)
    line_sum = F(0)
    for line in fano_lines:
        c_line = len(active.intersection(line))
        line_sum += F(c_line, 3) - F(comb(c_line, 2), 3)
    assert line_sum == q6(c)
    fano_mask_checks += 1

nu_three = F(1, 7) - F(1, 72)
e_three = nu_three * (1 - 3 * nu_three)
nu_pair = F(2, 21) - F(1, 273)
e_pair = nu_pair * (1 - 3 * nu_pair)
nu_same_pair = F(2, 21) - F(1, 231)
e_same_pair = nu_same_pair * (1 - 3 * nu_same_pair)
assert nu_three == F(65, 504) and e_three == F(6695, 84672)
assert nu_pair == F(25, 273) and e_pair == F(550, 8281)
assert nu_same_pair == F(1, 11) and e_same_pair == F(8, 121)
fano_all_line_constant = 4 * e_three + 3 * e_pair
fano_all_line_ratio = F(588, 107) * fano_all_line_constant
assert fano_all_line_constant == F(263465, 511056)
assert fano_all_line_ratio == F(1844255, 650988)
fano_chi_constant = 4 * e_three + 2 * e_same_pair + e_pair
fano_chi_ratio = F(588, 107) * fano_chi_constant
assert fano_chi_constant == F(222893927, 432864432)
assert fano_chi_ratio == F(222893927, 78769548)
print("Fano incidence masks:", fano_mask_checks, "through/off lines=3/4")
print("Fano line errors: triple=6695/84672, pair=550/8281, same-chi pair=8/121")
print("all-line period ratio:", fano_all_line_ratio, "=", float(fano_all_line_ratio))
print("chi-matched period ratio:", fano_chi_ratio, "=", float(fano_chi_ratio))

# 6. Exact zero-local-pair guardrail: endpoint errors cannot be discarded.
a = 14
k = 2
gap_lo = F(14 * k + 1, 14 * a)
gap_hi = F(14 * k + 13, 14 * a)
speeds = tuple(range(18, 24))
teeth = {s: lifted_teeth(s, gap_lo, gap_hi) for s in speeds}
single_mass = sum(
    (sum((right - left for _, left, right in teeth[s]), F(0)) for s in speeds),
    F(0),
)
local_pair_mass = sum(
    (lifted_pair_measure(teeth[x], teeth[y]) for x, y in combinations(speeds, 2)),
    F(0),
)
global_pair_mass = sum((pair_density(x, y) for x, y in combinations(speeds, 2)), F(0))
endpoint_tax = sum(
    (
        pair_density(x, y) * (1 - pair_density(x, y)) / gcd(x, y)
        for x, y in combinations(speeds, 2)
    ),
    F(0),
)
assert (gap_lo, gap_hi, gap_hi - gap_lo) == (F(29, 196), F(41, 196), F(3, 49))
assert single_mass == F(850237, 16959096)
assert local_pair_mass == 0
assert global_pair_mass == F(286151, 921690)
assert (gap_hi - gap_lo) * global_pair_mass - endpoint_tax < 0
occupied = sorted((left, right, s) for s in speeds for _, left, right in teeth[s])
assert tuple(s for _, _, s in occupied) == (20, 19, 18, 23, 22, 21, 20, 19)
assert all(right < next_left for (_, right, _), (next_left, _, _) in zip(occupied, occupied[1:]))
print(
    "zero-pair guardrail: a=14, k=2, b=18..23, G=[29/196,41/196],",
    "single mass=850237/16959096, local pair mass=0",
)

# 7. Tournament audit required by the repository methodology.
# O(i,j)=a/b_i-a/b_j, positive sign oriented from i to j.
scores = [0] * 6
ties = 0
for i, j in combinations(range(6), 2):
    observable = F(a, speeds[i]) - F(a, speeds[j])
    if observable > 0:
        scores[i] += 1
    elif observable < 0:
        scores[j] += 1
    else:
        ties += 1
assert sorted(scores) == list(range(6)) and ties == 0
print(
    "runner tournament: score histogram=0,1,2,3,4,5; cycles=0;",
    "SCCs=1^6; Hamilton paths=1; ties=0; endpoint-hypergraph verdict lost",
)

print("ALL EXACT CHECKS PASSED")
