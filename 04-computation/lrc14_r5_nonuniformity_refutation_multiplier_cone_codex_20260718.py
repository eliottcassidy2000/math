#!/usr/bin/env python3
"""Exact audit of the THM-1141 four-comb nonuniformity target.

This file has two logically independent jobs.

1.  It gives a rational counterexample to

        longest survivor gap >= (4/3) * mean survivor gap

    even when the interval is a maximal component of an eight-runner core and
    the four integer moduli are distinct and legal.

2.  It checks the finite Z/13 lemma behind a replacement analytic gate:
    every set of at most four residues has, after a nonzero multiplier, a
    cyclic gap of at least seven.  The accompanying note turns that finite
    fact into the cone B >= 15 max(A,84).

All arrangement endpoints and comparisons below use Fraction arithmetic.
The script prints data only and has no repository side effects.
"""

from fractions import Fraction as F
from itertools import combinations
from math import gcd


def norm_dist(x: F) -> F:
    r = x % 1
    return min(r, 1 - r)


def remove_comb(intervals: list[tuple[F, F]], k: int) -> list[tuple[F, F]]:
    """Remove D_k={t: ||kt||<1/14} from a rational interval list."""
    out: list[tuple[F, F]] = []
    for a, b in intervals:
        ak, bk = a * k, b * k
        j_lo = ak.numerator // ak.denominator - 1
        j_hi = bk.numerator // bk.denominator + 2
        cur = a
        for j in range(j_lo, j_hi + 1):
            left = F(14 * j - 1, 14 * k)
            right = F(14 * j + 1, 14 * k)
            if right <= a or left >= b:
                continue
            left, right = max(left, a), min(right, b)
            if left > cur:
                out.append((cur, left))
            cur = max(cur, right)
        if cur < b:
            out.append((cur, b))
    return out


def core_safe_components(P: tuple[int, ...]) -> list[tuple[F, F]]:
    """Reconstruct all maximal core-safe components from exact wall endpoints."""
    breakpoints = {F(0), F(1)}
    for p in P:
        for j in range(p + 1):
            for sign in (-1, 1):
                wall = F(j, p) + sign * F(1, 14 * p)
                if 0 <= wall <= 1:
                    breakpoints.add(wall)
    atomic = []
    ordered = sorted(breakpoints)
    for a, b in zip(ordered, ordered[1:]):
        if b <= a:
            continue
        midpoint = (a + b) / 2
        if all(norm_dist(p * midpoint) >= F(1, 14) for p in P):
            atomic.append((a, b))
    merged: list[tuple[F, F]] = []
    for a, b in atomic:
        if merged and merged[-1][1] == a:
            merged[-1] = (merged[-1][0], b)
        else:
            merged.append((a, b))
    return merged


def cyclic_gaps(points: tuple[int, ...], p: int = 13) -> tuple[int, ...]:
    s = sorted(points)
    return tuple(s[i + 1] - s[i] for i in range(len(s) - 1)) + (
        s[0] + p - s[-1],
    )


def affine_image(points: tuple[int, ...], u: int, v: int) -> tuple[int, ...]:
    return tuple(sorted((u * x + v) % 13 for x in points))


def multiplier_best(points: tuple[int, ...]) -> tuple[int, int, tuple[int, ...]]:
    rows = []
    for u in range(1, 13):
        image = affine_image(points, u, 0)
        gaps = cyclic_gaps(image)
        rows.append((max(gaps), -u, gaps))
    largest, neg_u, gaps = max(rows)
    return largest, -neg_u, gaps


print("FOUR-COMB NONUNIFORMITY AUDIT (exact rational arithmetic)")
print()
print("[1] Exact counterexample on a maximal core-safe component")

P = (1, 2, 3, 4, 5, 6, 7, 8)
J = (F(1, 14), F(13, 112))
killers = (108, 109, 110, 111)

midpoint = sum(J, F(0)) / 2
assert J in core_safe_components(P)
assert all(norm_dist(p * midpoint) >= F(1, 14) for p in P)
assert norm_dist(P[0] * J[0]) == F(1, 14)
assert norm_dist(P[-1] * J[1]) == F(1, 14)
assert killers[0] > 13 * max(P)

survivors = [J]
for k in killers:
    survivors = remove_comb(survivors, k)

lengths = [b - a for a, b in survivors]
mass = sum(lengths, F(0))
mean = mass / len(lengths)
longest = max(lengths)
ratio = longest / mean
metric = 7 * killers[-1] * longest

assert len(survivors) == 5
assert mass == F(955, 37296)
assert mean == F(191, 37296)
assert longest == F(319, 55944)
assert ratio == F(638, 573)
assert ratio < F(4, 3)
assert metric == F(319, 72) > 1

print(f"  core P             = {P}")
print(f"  maximal component  = [{J[0]}, {J[1]}]")
print(f"  distinct killers   = {killers} (legal base: {killers[0]} > 13*{max(P)})")
print(f"  survivor lengths   = {' '.join(str(x) for x in lengths)}")
print(f"  survivor mass      = {mass}")
print(f"  component count    = {len(lengths)}")
print(f"  actual mean gap    = {mean}")
print(f"  longest gap        = {longest}")
print(f"  longest / mean     = {ratio} = {float(ratio):.12f} < 4/3")
print(f"  sharp metric       = 7*111*L = {metric} = {float(metric):.12f} > 1")
print()

print("[2] Four-residue multiplier lemma on Z/13")
remaining = set(combinations(range(13), 4))
orbit_rows = []
covered = set()
while remaining:
    representative = min(remaining)
    orbit = {
        affine_image(representative, u, v)
        for u in range(1, 13)
        for v in range(13)
    }
    assert not (covered & orbit)
    covered.update(orbit)
    remaining.difference_update(orbit)
    largest, witness, gap_word = multiplier_best(representative)
    assert all(multiplier_best(points)[0] == largest for points in orbit)
    orbit_rows.append(
        (representative, len(orbit), witness, gap_word, largest)
    )

assert len(orbit_rows) == 7
assert sum(row[1] for row in orbit_rows) == 715
assert len(covered) == 715
assert min(row[4] for row in orbit_rows) == 7

print("  representative  orbit  witness  witness gap word  best largest gap")
for representative, size, witness, gap_word, largest in orbit_rows:
    print(
        "  %-14s %5d %8d  %-16s %d"
        % (
            "".join(format(x, "x") for x in representative),
            size,
            witness,
            " ".join(map(str, gap_word)),
            largest,
        )
    )

cardinality_minima = []
sharp_four_sets = 0
for cardinality in range(1, 5):
    minimum = 13
    for points in combinations(range(13), cardinality):
        largest, _, _ = multiplier_best(points)
        minimum = min(minimum, largest)
        if cardinality == 4 and largest == 7:
            sharp_four_sets += 1
    cardinality_minima.append(minimum)
assert cardinality_minima == [13, 12, 9, 7]
assert sharp_four_sets == 52
print(f"  minima for cardinalities 1..4 = {cardinality_minima}")
print(f"  labelled sharp four-sets       = {sharp_four_sets}")
print()

print("[3] Exact constants for the replacement multiplier cone")
epsilon_at_M84 = F(14, 365 * 84)
y = F(28, 365)
initial_vertical_width = F(7, 13) - F(1, 7)
common_vertical_width = initial_vertical_width - y
chosen_arc = F(11, 73)
crossing_budget = 15 * y
core_floor = F(1, 13) - 12 * epsilon_at_M84

assert epsilon_at_M84 == F(1, 2190) < F(1, 2184)
assert core_floor > F(1, 14)
assert y < F(1, 13)
assert common_vertical_width == F(10592, 33215) > chosen_arc
assert crossing_budget == 1 + chosen_arc
assert 11 * 7 > 73

print(f"  epsilon at M=84             = {epsilon_at_M84}")
print(f"  core floor                  = {core_floor} > 1/14")
print(f"  two-wall drift y            = {y} < 1/13")
print(f"  seven-grid-gap safe width   = {initial_vertical_width}")
print(f"  common width after drift    = {common_vertical_width} > {chosen_arc}")
print(f"  B|I| at B/M=15             = {crossing_budget} = 1 + {chosen_arc}")
print(f"  needle component length     = {chosen_arc}/B > 1/(7*k4) when B<=k4")
print("  midpoint corollary          = k1 >= max(1260, 7*(Delta+1))")
print()


def best_fixed_chart(k1: int, k4: int) -> tuple[int, int]:
    """Return (B,M) maximizing B/max(max(B-k1,k4-B),84)."""
    return max(
        (
            (B, max(B - k1, k4 - B, 84))
            for B in range(k1, k4 + 1)
        ),
        key=lambda row: F(row[0], row[1]),
    )


def best_fixed_chart_closed_form(k1: int, k4: int) -> tuple[int, int]:
    delta = k4 - k1
    if delta <= 84:
        return k4, 84
    if delta <= 168:
        return k1 + 84, 84
    half = (delta + 1) // 2
    return k1 + half, half


# Exhaustively guard the piecewise chart optimizer on a finite box containing
# every transition and both midpoint parities.
for k1_test in range(1, 201):
    for k4_test in range(k1_test + 3, 401):
        brute = best_fixed_chart(k1_test, k4_test)
        closed = best_fixed_chart_closed_form(k1_test, k4_test)
        assert F(*brute) == F(*closed)

print("[4] Exact four-comb separated-ratio gate Q4")
print("  survivor mass lower bound A0 = 3*ell/7 - (6/49)*sum_i(1/k_i)")
print("  component upper bound     C0 = ell*sum_i(k_i) + 39/7")
print("  Q4 = ell*(21*K-7*sum_i(k_i)) - 6*K*sum_i(1/k_i) - 39")
print("     = ell*(14*K-7*(k1+k2+k3))")
print("       - 6*K*(1/k1+1/k2+1/k3) - 45")
print("  identity: 7*K*A0-C0 = Q4/7; hence Q4>0 gives L>1/(7*K)")
print("  exact best cone chart by Delta=k4-k1:")
print("    Delta<=84:       (B,M)=(k4,84), gate k4>=1260")
print("    85<=Delta<=168:  (B,M)=(k1+84,84), gate k1>=1176")
print("    Delta>=169:      B=k1+ceil(Delta/2), M=ceil(Delta/2)")
print("                     gate k1>=14*ceil(Delta/2)")
print()

# The longest-component atlas puts the first transfer in its saturated regime
# uniformly: ell*k1>13/7 at every legal base.
atlas_min = None
for core in combinations(range(1, 13), 8):
    ell = max(b - a for a, b in core_safe_components(core))
    legal_product = ell * (13 * max(core) + 1)
    row = (legal_product, core, ell)
    if atlas_min is None or row < atlas_min:
        atlas_min = row
assert atlas_min == (
    F(72, 35),
    (1, 2, 6, 7, 8, 9, 10, 11),
    F(1, 70),
)

print("[5] Gate-union test and exact residual")
print("  core-atlas minimum ell*(13*max(P)+1) = 72/35 > 13/7")
print("  attained by P=(1,2,6,7,8,9,10,11), ell=1/70")


def transfer_phi_pair(x_num: int, x_den: int) -> tuple[int, int]:
    """Exact Phi(x), for x_num/x_den >= 1, as a reduced pair."""
    assert x_num >= x_den
    if 7 * x_num >= 13 * x_den:
        return (6, 7)
    num, den = 7 * x_num - x_den, 14 * x_den
    common = gcd(num, den)
    return (num // common, den // common)


def exact_transfer_gate(a: int, b: int, c: int, d: int) -> bool:
    """THM-1137 transfer from saturated c1=6/7 through three ratios."""
    num, den = 6, 7
    for old, new in ((a, b), (b, c), (c, d)):
        x_num, x_den = new * num, old * den
        if x_num < x_den:
            return False
        num, den = transfer_phi_pair(x_num, x_den)
    return 7 * num > den


# A transparent common-ratio corollary improves the corrected 7/3 coarse
# recursion.  The exact optimal common-ratio threshold in this transfer scheme
# is the positive root of 6r^3-r^2-2r-28=0 (about 1.79711); 9/5 is rational.
common_ratio = F(9, 5)
transfer_ledger = [F(6, 7)]
transfer_inputs = []
for _ in range(3):
    x = common_ratio * transfer_ledger[-1]
    transfer_inputs.append(x)
    assert x >= 1
    transfer_ledger.append(min(F(6, 7), (x - F(1, 7)) / 2))
assert transfer_inputs == [F(54, 35), F(63, 50), F(3519, 3500)]
assert transfer_ledger == [F(6, 7), F(7, 10), F(391, 700), F(3019, 7000)]

print("  corrected exact transfer Phi(x)=min(6/7,(x-1/7)/2)")
print("  saturated start c1=6/7; adjacent ratio 9/5 gives")
print("    inputs 54/35, 63/50, 3519/3500")
print("    outputs 7/10, 391/700, 3019/7000 > 1/7")
print("  optimal common-ratio threshold solves 6r^3-r^2-2r-28=0")
print("  (r=1.797111878...; 9/5 is a transparent rational corollary)")

# Census primitive ratio rays.  A ray is counted once after division by the
# common gcd.  These are asymptotic gates: strict cone, strict Q4 slope, and
# the corrected exact scale-free transfer.  The finite span<=30 bank disappears
# along every nonconstant ray and is handled separately in the exact residual
# predicate in the note.
ray_cap = 100
ray_count = cone_count = q4_count = transfer_count = union_count = 0
cone_boundary_count = q4_boundary_count = residual_count = 0
first_residual = None
for a, b, c, d in combinations(range(1, ray_cap + 1), 4):
    if gcd(gcd(a, b), gcd(c, d)) != 1:
        continue
    ray_count += 1
    cone = 8 * a > 7 * d
    q4_tail = 2 * d > a + b + c
    transfer = exact_transfer_gate(a, b, c, d)
    cone_count += int(cone)
    q4_count += int(q4_tail)
    transfer_count += int(transfer)
    cone_boundary_count += int(8 * a == 7 * d)
    q4_boundary_count += int(2 * d == a + b + c)
    covered_ray = cone or q4_tail or transfer
    union_count += int(covered_ray)
    if not covered_ray:
        residual_count += 1
        if first_residual is None:
            first_residual = (a, b, c, d)

assert first_residual == (3, 4, 5, 6)
assert (ray_count, union_count, residual_count) == (3646069, 3054412, 591657)

print(f"  primitive rays with 1<=a1<...<a4<={ray_cap}: {ray_count}")
print(f"    strict multiplier-cone rays = {cone_count}")
print(f"    strict Q4-tail rays          = {q4_count}")
print(f"    exact-transfer rays          = {transfer_count}")
print(f"    union                        = {union_count}")
print(f"    residual                     = {residual_count} ({residual_count/ray_count:.6f})")
print(f"    cone equality rays           = {cone_boundary_count}")
print(f"    Q4-slope equality rays       = {q4_boundary_count}")
print(f"    first residual ray           = {first_residual}")

# The first residual is an infinite legal family and sits simultaneously on
# the Q4-slope wall and outside the exact-transfer domain.
for m in (53, 54, 101):
    r = (3 * m, 4 * m, 5 * m, 6 * m)
    assert r[0] > 13 * 12
    assert r[-1] - r[0] > 30
    B, M = best_fixed_chart(r[0], r[-1])
    assert B < 15 * M
    assert 14 * r[-1] - 7 * sum(r[:-1]) == 0
    q4_value = -6 * r[-1] * sum((F(1, k) for k in r[:-1]), F(0)) - 45
    assert q4_value == -F(366, 5)
    assert not exact_transfer_gate(*r)

print("  infinite exact survivor       = m*(3,4,5,6), every m>=53")
print("    span=3m>30; cone fails; Q4=-366/5 (zero slope)")
print("    transfer: c1=6/7 -> c2=1/2, then (5/4)c2=5/8<1")
print("  asymptotic residual ratios (x,y,z)=(k1,k2,k3)/k4 satisfy")
print("    x<=7/8, x+y+z>=2, and the three-step exact Phi path fails")
print()
print("DONE")
