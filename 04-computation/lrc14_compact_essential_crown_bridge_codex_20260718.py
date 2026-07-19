#!/usr/bin/env python3
"""Exact audit for the compact strict-cover / essential-crown bridge.

This script supports the accompanying research note.  It does not claim an
LRC(14) proof.  Exact decisions use integers and fractions.Fraction.

The pair-sum-ruler evaluator is complete by the repository's THM-668/THM-1002:
at a global maximum below 1/2, a reduced denominator divides the sum of two
oppositely sloped active speeds (the equal-speed case is included).
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations
from math import gcd, isqrt


def divisors(number: int) -> set[int]:
    result: set[int] = set()
    for candidate in range(1, isqrt(number) + 1):
        if number % candidate == 0:
            result.add(candidate)
            result.add(number // candidate)
    return result


def clearance(speed: int, numerator: int, denominator: int) -> int:
    residue = speed * numerator % denominator
    return min(residue, denominator - residue)


def exact_M(values: tuple[int, ...]) -> tuple[F, tuple[int, int, int]]:
    rulers: set[int] = set()
    for index, left in enumerate(values):
        for right in values[index:]:
            rulers.update(divisors(left + right))

    best = F(0)
    owner = (0, 1, 0)
    for denominator in sorted(rulers):
        if denominator < 2:
            continue
        for numerator in range(1, denominator):
            if gcd(numerator, denominator) != 1:
                continue
            raw_clearance = min(
                clearance(speed, numerator, denominator) for speed in values
            )
            value = F(raw_clearance, denominator)
            if value > best:
                best = value
                owner = (numerator, denominator, raw_clearance)
    return best, owner


def covers(values: tuple[int, ...], upper: int) -> bool:
    return all(any(speed % modulus == 0 for speed in values) for modulus in range(2, upper + 1))


def bounded_witness(
    values: tuple[int, ...], upper_denominator: int, stop_above: F | None = None
) -> tuple[F, tuple[int, int, int]]:
    """A sound lower bound on M; optionally stop after a strict witness."""
    best = F(0)
    owner = (0, 1, 0)
    for denominator in range(2, upper_denominator + 1):
        for numerator in range(1, denominator):
            if gcd(numerator, denominator) != 1:
                continue
            raw_clearance = min(
                clearance(speed, numerator, denominator) for speed in values
            )
            value = F(raw_clearance, denominator)
            if value > best:
                best = value
                owner = (numerator, denominator, raw_clearance)
            if stop_above is not None and value > stop_above:
                return best, owner
    return best, owner


def danger_teeth(speed: int, denominator: int) -> tuple[tuple[F, F], ...]:
    half_width = F(1, denominator * speed)
    teeth = []
    for integer in range(speed + 1):
        left = max(F(0), F(integer, speed) - half_width)
        right = min(F(1), F(integer, speed) + half_width)
        if left < right:
            teeth.append((left, right))
    return tuple(teeth)


def multiplicity_profile(
    values: tuple[int, ...], denominator: int
) -> tuple[dict[int, F], dict[int, tuple[tuple[F, F], ...]]]:
    """Open-chamber multiplicities; endpoints have measure zero."""
    teeth = {speed: danger_teeth(speed, denominator) for speed in values}
    endpoints = {F(0), F(1)}
    for speed_teeth in teeth.values():
        for left, right in speed_teeth:
            endpoints.add(left)
            endpoints.add(right)

    histogram: defaultdict[int, F] = defaultdict(F)
    private: defaultdict[int, list[tuple[F, F]]] = defaultdict(list)
    ordered = sorted(endpoints)
    for left, right in zip(ordered, ordered[1:]):
        midpoint = (left + right) / 2
        active = tuple(
            speed
            for speed in values
            if any(a < midpoint < b for a, b in teeth[speed])
        )
        histogram[len(active)] += right - left
        if len(active) == 1:
            private[active[0]].append((left, right))
    return dict(histogram), {speed: tuple(parts) for speed, parts in private.items()}


def deleted_core_rows(values: tuple[int, ...], threshold_denominator: int):
    threshold = F(1, threshold_denominator)
    histogram, private = multiplicity_profile(values, threshold_denominator)
    rows = []
    for speed in values:
        core = tuple(value for value in values if value != speed)
        core_M, owner = exact_M(core)
        delta = core_M - threshold
        private_parts = private.get(speed, ())
        longest = max((right - left for left, right in private_parts), default=F(0))
        lipschitz_floor = 2 * delta / max(core) if delta > 0 else F(0)
        tooth_ceiling = F(2, threshold_denominator * speed)
        if delta > 0:
            assert lipschitz_floor <= longest < tooth_ceiling
        rows.append((speed, core_M, owner, lipschitz_floor, longest, tooth_ceiling))
    return histogram, private, tuple(rows)


print("COMPACT STRICT-COVER / ESSENTIAL-CROWN EXACT AUDIT")

# The exact all-loose crown which misses Cover14.
V0 = (1, 2, 3, 5, 7, 8, 9, 10, 11, 12, 17, 19, 104)
M0, owner0 = exact_M(V0)
assert M0 == F(8, 105)
assert gcd(*V0) == 1 and F(V0[-1], V0[-2]) < 13
assert covers(V0, 13) and not covers(V0, 14)
histogram0, private0, rows0 = deleted_core_rows(V0, 13)
assert len(private0) == 13
assert all(row[1] > F(1, 13) for row in rows0)

private_mass = histogram0[1]
triple_excess = sum(
    F(multiplicity - 2) * measure
    for multiplicity, measure in histogram0.items()
    if multiplicity >= 3
)
assert histogram0.get(0, F(0)) == 0
assert sum(F(multiplicity) * measure for multiplicity, measure in histogram0.items()) == 2
assert private_mass == triple_excess

print("V0 (the exact all-loose crown):")
print(f"  values={V0}")
print(f"  primitive=True compact_rho={F(V0[-1], V0[-2])} cover13=True cover14=False")
print(f"  M={M0} owner={owner0}; threshold=1/13")
print(f"  deleted_cores={len(rows0)} tight_deleted_cores=0 positive_private_owners={len(private0)}")
print(f"  minimum_deleted_core_M={min(row[1] for row in rows0)}")
print(f"  private_mass={private_mass}")
print(f"  weighted_triple_excess={triple_excess}")
sum_lipschitz_floors = sum(row[3] for row in rows0)
print(f"  sum_Lipschitz_component_floors={sum_lipschitz_floors}")
print(f"  private_mass_over_floor_sum={float(private_mass / sum_lipschitz_floors):.12f}")
print(f"  multiplicity_histogram={tuple(sorted(histogram0.items()))}")
print("  owner rows: v, M(V\\v), owner(a,q,val), Lipschitz floor, longest private chamber, tooth width")
for row in rows0:
    print(f"    {row}")

# The exact lift changes only an integer representative in the 8/105 chart,
# supplies Cover14, and creates a new global witness at 3/20.
V1 = (1, 2, 3, 5, 8, 9, 10, 11, 12, 17, 19, 104, 112)
M1, owner1 = exact_M(V1)
assert M1 == F(3, 20)
assert gcd(*V1) == 1 and F(V1[-1], V1[-2]) < 13 and covers(V1, 14)
residues0 = tuple(sorted(8 * value % 105 for value in V0))
residues1 = tuple(sorted(8 * value % 105 for value in V1))
assert residues0 == residues1
assert 112 == 7 + 105

print("V0 -> V1 cross-modulus lift:")
print(f"  replace 7 by 112=7+105; residues_at_8/105_equal={residues0 == residues1}")
print(f"  V1 primitive=True compact_rho={F(V1[-1], V1[-2])} cover14=True")
print(f"  V1 M={M1} owner={owner1}; the strict 1/13 crown disappears")

lift_rows = []
for lift_index in range(6):
    carrier = 112 + 210 * lift_index
    lifted = tuple(sorted((1, 2, 3, 5, 8, 9, 10, 11, 12, 17, 19, 104, carrier)))
    assert F(lifted[-1], lifted[-2]) < 13 and covers(lifted, 14)
    lifted_M, lifted_owner = exact_M(lifted)
    assert lifted_M == F(3, 20)
    lift_rows.append((carrier, F(lifted[-1], lifted[-2]), lifted_M, lifted_owner))
print(f"  all compact 14-carrier lifts 112+210k, k=0..5: {tuple(lift_rows)}")

# Farey-toothpick blocker: the d copies of the twelve base maximizers are
# n/(13d), 13 not dividing n.  A nonzero residue v mod 13d cannot be dangerous
# at all 12d points.  Brute-force the exact assertion on a deterministic box.
blocker_rows = 0
for scale in range(1, 31):
    modulus = 13 * scale
    base_maximizers = tuple(n for n in range(modulus) if n % 13)
    assert len(base_maximizers) == 12 * scale
    for speed_residue in range(1, modulus):
        assert any(
            clearance(speed_residue, numerator, modulus) >= scale
            for numerator in base_maximizers
        )
        blocker_rows += 1
print("Farey-toothpick AP-core blocker:")
print(f"  brute_force_nonzero_residue_rows={blocker_rows} scales=1..30")
print("  symbolic count: Q>=2 has b(Q)=1+2*floor((Q-1)/13) < 12Q/13")
print("  conclusion: M(d*[12] union {v})<1/13 implies 13d divides v")
print("  primitive + Cover14 then d=1 and 182 divides v; rho>=182/12>13")

# Exhaust the full single-substitution neighbourhood of the sharp compact
# boundary through candidate speed 500.  A bounded witness >1/13 soundly
# rejects a row; the few unrejected rows are then evaluated exactly.
boundary = (2, 4, 6, 8, 10, 12, 13, 14, 16, 18, 20, 22, 24)
single_substitution_legal = 0
single_substitution_survivors = []
for removed in boundary:
    remainder = set(boundary) - {removed}
    for added in range(1, 501):
        if added == removed or added in remainder:
            continue
        candidate = tuple(sorted(remainder | {added}))
        if len(candidate) != 13:
            continue
        if gcd(*candidate) != 1 or not covers(candidate, 14):
            continue
        if F(candidate[-1], candidate[-2]) >= 13:
            continue
        single_substitution_legal += 1
        lower_bound, _ = bounded_witness(candidate, 200, F(1, 13))
        if lower_bound <= F(1, 13):
            candidate_M, candidate_owner = exact_M(candidate)
            assert candidate_M == F(1, 13)
            assert removed == 13
            assert added % 26 == 13
            assert set(candidate) - {added} == set(range(2, 25, 2))
            single_substitution_survivors.append((removed, added, candidate_M, candidate_owner))

assert tuple(row[1] for row in single_substitution_survivors) == tuple(range(39, 300, 26))
print("Sharp-boundary single-substitution hard test:")
print(f"  legal primitive Cover14 compact rows (added<=500)={single_substitution_legal}")
print(f"  rows with no q<=200 witness above 1/13={len(single_substitution_survivors)}")
print(f"  exact survivors={tuple(single_substitution_survivors)}")
print("  every survivor is equality with the same regenerated core 2*[12]; below-threshold rows=0")

# Equality extraction is not a topological consequence of a minimal strict
# double cover.  These lower-n rows have every analogous side condition but
# every deletion is strictly above 1/n.
small_rows = (
    (3, (1, 3, 4)),
    (4, (1, 3, 4, 5)),
    (5, (1, 4, 5, 6, 7)),
    (6, (1, 2, 5, 6, 7, 8)),
    (7, (1, 4, 5, 6, 7, 11, 16)),
)
print("Lower-n counterexamples to a purely topological tight-deletion lemma:")
for denominator, values in small_rows:
    value_M, owner = exact_M(values)
    deleted = tuple(
        exact_M(values[:index] + values[index + 1 :])[0]
        for index in range(len(values))
    )
    assert gcd(*values) == 1
    assert covers(values, denominator + 1)
    assert F(values[-1], values[-2]) < denominator
    assert value_M < F(1, denominator) < min(deleted)
    print(
        f"  n={denominator} V={values} M={value_M} owner={owner} "
        f"min_deleted_M={min(deleted)} cover_2_to_{denominator + 1}=True compact=True"
    )

# Alternate-carrier / tournament audit.  A total order on private stalks is
# transitive by construction and therefore loses the multiplicity-3 gluing
# which the exact balance identity detects.
stalks = []
for owner, parts in private0.items():
    for left, right in parts:
        stalks.append(((left + right) / 2, owner, left, right))
stalks.sort()
stalk_count = len(stalks)
scores = tuple(range(stalk_count))
score_histogram = {score: 1 for score in scores}
print("Alternate carrier / tournament audit:")
print(f"  vertices=private_stalks count={stalk_count}")
print("  observable=cyclic midpoint order after cut at 0; tie path=(midpoint,owner)")
print(f"  tournament=transitive score_histogram_size={len(score_histogram)} directed_3_cycles=0 SCCs={stalk_count} Hamiltonian_paths=1")
print("  destroys=tooth labels, simultaneous >=3 overlap, Cover14 lift compatibility")
print("  faithful carrier=sheet/event word + private-stalk/triple-overlap chamber complex + owner-lift sidecar")

print("ALL EXACT CHECKS PASS")
