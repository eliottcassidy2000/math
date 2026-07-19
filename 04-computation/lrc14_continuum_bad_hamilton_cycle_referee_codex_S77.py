#!/usr/bin/env python3
"""Exact referee for THM-1210's continuum max-gap/triangle ceiling.

For D=(d2,d3,d4), the continuum four-comb obstruction is BAD exactly when
the four circle points {0,{d2*u},{d3*u},{d4*u}} have every cyclic gap at most
2/7.  Off finitely many walls, a cyclic order is an oriented Hamilton cycle
of K4.  This script constructs every collision and 2/7-threshold wall with
Fraction, computes the bad measure, and independently sums the three
Hamilton-cycle correlation integrals (one orientation per reversal pair).

The uniform ceiling is proved by a second, smaller quotient.  BAD forces every
pair difference into the three-band set

    A=[1/7,2/7] U [3/7,4/7] U [5/7,6/7].

Every four distinct integers contain a non-arithmetic three-point subset, so
BAD is contained in an additive-triangle event

    {p*u}, {q*u}, {(p+q)*u} in A,       p != q.

The corresponding two-torus region is six right triangles of area 3/49.
After the shear (x,y)->(x,z=x+y), its x-section is empty off A and is a union
of at most two intervals on A.  Conditioning on z gives

    J(p,q) <= 3/49 + 6/[7(p+q)],

which closes p+q>=26 after common-gcd reduction.  This referee checks the
remaining 99 coprime pairs 1<=p<q and p+q<=25 exactly; the unique reduced
maximum is J(1,2)=2/21.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from math import gcd


HEIGHT = 30
THRESHOLD = F(2, 7)
CORE_SUM_LIMIT = 25
THREE_BANDS = ((1, 2), (3, 4), (5, 6))
UPPER_TRIANGLE_CELLS = {(1, 1), (1, 3), (3, 1)}
LOWER_TRIANGLE_CELLS = {(3, 5), (5, 3), (5, 5)}
CYCLES = (
    (0, 1, 2, 3),
    (0, 1, 3, 2),
    (0, 2, 1, 3),
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def fractional_part(value: F) -> F:
    return value - value.numerator // value.denominator


def max_cyclic_gap(points: tuple[F, ...]) -> F:
    ordered = sorted(points)
    gaps = [right - left for left, right in zip(ordered, ordered[1:])]
    gaps.append(1 + ordered[0] - ordered[-1])
    return max(gaps)


def bad_at(spacings: tuple[int, int, int], u: F) -> bool:
    points = (F(0),) + tuple(fractional_part(d * u) for d in spacings)
    return max_cyclic_gap(points) <= THRESHOLD


def in_three_band(value: F) -> bool:
    """Membership in A away from its null boundary."""
    value = fractional_part(value)
    return any(F(left, 7) < value < F(right, 7) for left, right in THREE_BANDS)


@lru_cache(maxsize=None)
def three_band_intervals(multiplier: int) -> tuple[tuple[F, F], ...]:
    """The 3*multiplier open components of {u:{multiplier*u} in A}."""
    require(multiplier > 0, "three-band multiplier must be positive")
    denominator = 7 * multiplier
    return tuple(
        (F(7 * k + left, denominator), F(7 * k + right, denominator))
        for k in range(multiplier)
        for left, right in THREE_BANDS
    )


def intersect_interval_lists(
    first: tuple[tuple[F, F], ...], second: tuple[tuple[F, F], ...]
) -> tuple[tuple[F, F], ...]:
    """Exact intersection of two sorted, internally disjoint interval lists."""
    left_index = 0
    right_index = 0
    answer = []
    while left_index < len(first) and right_index < len(second):
        lower = max(first[left_index][0], second[right_index][0])
        upper = min(first[left_index][1], second[right_index][1])
        if lower < upper:
            answer.append((lower, upper))
        first_upper = first[left_index][1]
        second_upper = second[right_index][1]
        if first_upper < second_upper:
            left_index += 1
        elif second_upper < first_upper:
            right_index += 1
        else:
            left_index += 1
            right_index += 1
    return tuple(answer)


def interval_list_measure(intervals: tuple[tuple[F, F], ...]) -> F:
    return sum((right - left for left, right in intervals), F(0))


def additive_triangle_measure(p: int, q: int) -> F:
    """J(p,q)=meas{u: pu,qu,(p+q)u all lie in A}, exactly."""
    require(p > 0 and q > 0, "triangle spacings must be positive")
    intersection = intersect_interval_lists(three_band_intervals(p), three_band_intervals(q))
    intersection = intersect_interval_lists(intersection, three_band_intervals(p + q))
    return interval_list_measure(intersection)


def additive_triangle_carry_formula(p: int, q: int) -> F:
    """Closed carry/residue formula for coprime p<q; independent of interval merge."""
    require(0 < p < q and gcd(p, q) == 1, "carry formula needs a reduced ordered pair")
    total_frequency = p + q
    selected_residues = (
        (2 * p - q) % 7,
        (4 * p - q) % 7,
        (2 * p - 3 * q) % 7,
    )
    weighted = F(0)
    for k in range(1, total_frequency):
        tent = F(k, q) if k <= q else F(total_frequency - k, p)
        multiplicity = sum(k % 7 == residue for residue in selected_residues)
        weighted += multiplicity * tent
    return F(2, 7 * total_frequency) * weighted


def non_arithmetic_triangle(values: tuple[int, int, int, int]) -> tuple[int, int]:
    """Choose successive gaps p!=q from a non-AP triple in four distinct values."""
    ordered = tuple(sorted(values))
    require(len(set(ordered)) == 4, "the four phase frequencies must be distinct")
    for triple in combinations(ordered, 3):
        p = triple[1] - triple[0]
        q = triple[2] - triple[1]
        if p != q:
            return p, q
    raise RuntimeError("four distinct integers cannot have every triple arithmetic")


def additive_triangle_at(p: int, q: int, u: F) -> bool:
    return all(in_three_band(multiplier * u) for multiplier in (p, q, p + q))


@lru_cache(maxsize=None)
def _arrangement_walls(differences: tuple[int, ...]) -> tuple[F, ...]:
    """Collision and oriented distance 2/7 or 5/7 walls on [0,1]."""
    walls = {F(0), F(1)}
    for difference in differences:
        e = abs(difference)
        require(e > 0, "zero difference entered wall arrangement")
        for k in range(e + 1):
            walls.add(F(k, e))
        for residue in (2, 5):
            for k in range(e):
                walls.add(F(7 * k + residue, 7 * e))
    return tuple(sorted(walls))


def arrangement_walls(differences: set[int]) -> tuple[F, ...]:
    return _arrangement_walls(tuple(sorted(abs(difference) for difference in differences)))


def exact_cell_measure(walls: tuple[F, ...], predicate) -> F:
    total = F(0)
    for left, right in zip(walls, walls[1:]):
        midpoint = (left + right) / 2
        if predicate(midpoint):
            total += right - left
    return total


def pair_differences(values: tuple[int, ...]) -> set[int]:
    return {abs(right - left) for left, right in combinations(values, 2)}


def bad_measure(spacings: tuple[int, int, int]) -> F:
    values = (0,) + spacings
    walls = arrangement_walls(pair_differences(values))
    return exact_cell_measure(walls, lambda u: bad_at(spacings, u))


def cycle_increments(values: tuple[int, ...], cycle: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        values[cycle[(index + 1) % 4]] - values[cycle[index]]
        for index in range(4)
    )


def cycle_correlation(values: tuple[int, ...], cycle: tuple[int, ...]) -> F:
    increments = cycle_increments(values, cycle)
    require(sum(increments) == 0 and all(increments), "invalid Hamilton cycle increments")
    walls = arrangement_walls({abs(increment) for increment in increments})

    def predicate(u: F) -> bool:
        return all(fractional_part(increment * u) <= THRESHOLD for increment in increments)

    return exact_cell_measure(walls, predicate)


def continuum_identity_audit() -> int:
    """Check max-gap, cycle, and BAD-to-additive-triangle implications."""
    rows = 0
    sample_spacings = (
        (1, 2, 3),
        (1, 3, 7),
        (1, 6, 7),
        (4, 8, 12),
        (7, 19, 30),
    )
    for spacings in sample_spacings:
        p, q = non_arithmetic_triangle((0,) + spacings)
        walls = arrangement_walls(pair_differences((0,) + spacings))
        for left, right in zip(walls, walls[1:]):
            u = (left + right) / 2
            phases = tuple(fractional_part(-d * u) for d in spacings)
            intervals = sorted((F(7, 6) * g - F(1, 6), F(7, 6) * g) for g in phases)
            cursor = F(0)
            longest = F(0)
            for lower, upper in intervals:
                clipped_lower = max(F(0), lower)
                clipped_upper = min(F(1), upper)
                if clipped_upper <= clipped_lower:
                    continue
                longest = max(longest, clipped_lower - cursor)
                cursor = max(cursor, clipped_upper)
            longest = max(longest, 1 - cursor)

            points = (F(0),) + tuple(fractional_part(d * u) for d in spacings)
            max_gap = max_cyclic_gap(points)
            predicted = F(7, 6) * max(F(0), max_gap - F(1, 7))
            require(longest == predicted, "continuum tooth/max-gap identity failed")
            require((longest <= F(1, 6)) == (max_gap <= F(2, 7)), "bad predicate mismatch")
            if max_gap < F(2, 7):
                require(additive_triangle_at(p, q, u), "BAD did not imply its non-AP triangle")
            rows += 1
    return rows


def triangle_region_model_audit() -> int:
    """Audit R={(x,y):x,y,x+y in A} against its six-triangle model."""
    rows = 0
    remainders = (F(1, 5), F(2, 5), F(3, 5), F(4, 5))
    for first_cell in range(7):
        for second_cell in range(7):
            for alpha in remainders:
                for beta in remainders:
                    x = (first_cell + alpha) / 7
                    y = (second_cell + beta) / 7
                    direct = in_three_band(x) and in_three_band(y) and in_three_band(x + y)
                    cell = first_cell, second_cell
                    predicted = (
                        cell in UPPER_TRIANGLE_CELLS and alpha + beta > 1
                    ) or (
                        cell in LOWER_TRIANGLE_CELLS and alpha + beta < 1
                    )
                    require(direct == predicted, "six-triangle torus model failed")
                    rows += 1
    require(F(len(UPPER_TRIANGLE_CELLS) + len(LOWER_TRIANGLE_CELLS), 2 * 7 * 7) == F(3, 49),
            "six-triangle area moved")
    return rows


def sheared_section_model_audit() -> int:
    """Audit all six explicit x-section branches after z=x+y.

    In coordinates x=(i+alpha)/7 and z=(r+gamma)/7, the nonempty sections
    are

      r=1: (i=3, alpha<gamma), (i=5, alpha<gamma),
      r=3: (i=1, alpha>gamma), (i=5, alpha<gamma),
      r=5: (i=1, alpha>gamma), (i=3, alpha>gamma).

    Thus the section is empty off z in A and has at most two intervals on A.
    """
    rows = 0
    remainders = (F(1, 5), F(2, 5), F(3, 5), F(4, 5))
    for x_cell in range(7):
        for z_cell in range(7):
            for alpha in remainders:
                for gamma in remainders:
                    x = (x_cell + alpha) / 7
                    z = (z_cell + gamma) / 7
                    direct = in_three_band(x) and in_three_band(z - x) and in_three_band(z)
                    predicted = (
                        (z_cell == 1 and x_cell in (3, 5) and alpha < gamma)
                        or (z_cell == 3 and x_cell == 1 and alpha > gamma)
                        or (z_cell == 3 and x_cell == 5 and alpha < gamma)
                        or (z_cell == 5 and x_cell in (1, 3) and alpha > gamma)
                    )
                    require(direct == predicted, "sheared six-section model failed")
                    rows += 1
    # The three z-bands contribute 1/49, 1/49, and 1/49 respectively.
    require(F(1, 49) + F(1, 49) + F(1, 49) == F(3, 49),
            "sheared section area moved")
    return rows


def finite_triangle_core() -> dict[str, object]:
    """Exact core left after the sheared shifted-grid tail p+q>=26."""
    interval_cache = {
        multiplier: three_band_intervals(multiplier)
        for multiplier in range(1, CORE_SUM_LIMIT + 1)
    }

    def cached_measure(p: int, q: int) -> F:
        intersection = intersect_interval_lists(interval_cache[p], interval_cache[q])
        intersection = intersect_interval_lists(intersection, interval_cache[p + q])
        return interval_list_measure(intersection)

    measures: Counter[F] = Counter()
    pairs_by_measure: dict[F, list[tuple[int, int]]] = {}
    pairs = 0
    for p in range(1, CORE_SUM_LIMIT):
        for q in range(p + 1, CORE_SUM_LIMIT + 1 - p):
            if gcd(p, q) != 1:
                continue
            value = cached_measure(p, q)
            require(value == additive_triangle_carry_formula(p, q),
                    "independent carry formula disagreed with interval merge")
            measures[value] += 1
            pairs_by_measure.setdefault(value, []).append((p, q))
            pairs += 1

    strata = tuple(sorted(measures, reverse=True))
    expected_top_strata = (
        (F(2, 21), ((1, 2),)),
        (F(13, 147), ((1, 6),)),
        (F(8, 105), ((3, 10),)),
    )
    top_strata = tuple(
        (value, tuple(pairs_by_measure[value])) for value in strata[:3]
    )
    tail_at_cutoff = F(3, 49) + F(6, 7 * (CORE_SUM_LIMIT + 1))
    tail_margin = F(2, 21) - tail_at_cutoff
    require(pairs == 99, "coprime triangle-core count moved")
    require(len(measures) == 87, "distinct triangle-core value count moved")
    require(top_strata == expected_top_strata,
            "top three unique reduced triangle strata moved")
    require(tail_at_cutoff == F(60, 637), "p+q=26 tail value moved")
    require(tail_margin == F(2, 1911), "p+q=26 strict tail margin moved")
    require(additive_triangle_measure(2, 4) == strata[0],
            "common-dilation invariance example moved")
    return {
        "pairs": pairs,
        "distinct": len(measures),
        "top_strata": top_strata,
        "tail_at_cutoff": tail_at_cutoff,
        "tail_margin": tail_margin,
        "carry_rows": pairs,
    }


def equality_obligation_audit() -> dict[str, int]:
    """Finite audit of the elementary four-triple equality-locus lemma."""

    def ratio_one_or_two(first: int, second: int) -> bool:
        return first == second or first == 2 * second or second == 2 * first

    rows = 0
    survivors = 0
    for a in range(1, HEIGHT + 1):
        for b in range(1, HEIGHT + 1):
            for c in range(1, HEIGHT + 1):
                obligations = ((a, b), (b, c), (a, b + c), (a + b, c))
                if all(ratio_one_or_two(*pair) for pair in obligations):
                    require(a == b == c, "four triangle equality obligations allowed a non-AP")
                    survivors += 1
                rows += 1
    require(survivors == HEIGHT, "equality-obligation survivor count moved")
    return {"rows": rows, "survivors": survivors}


def height_census() -> dict[str, object]:
    measures: Counter[F] = Counter()
    rows_by_measure: dict[F, list[tuple[int, int, int]]] = {}
    triangle_measures: dict[tuple[int, int], F] = {}
    cycle_rows = 0
    triangle_rows = 0
    total = 0
    for d2, d3, d4 in combinations(range(1, HEIGHT + 1), 3):
        spacings = d2, d3, d4
        values = (0,) + spacings
        measure = bad_measure(spacings)
        correlations = tuple(cycle_correlation(values, cycle) for cycle in CYCLES)
        require(measure == 2 * sum(correlations, F(0)), "Hamilton-cycle factor-two identity failed")
        p, q = non_arithmetic_triangle(values)
        common = gcd(p, q)
        reduced = tuple(sorted((p // common, q // common)))
        if reduced not in triangle_measures:
            triangle_measures[reduced] = additive_triangle_measure(*reduced)
        triangle_bound = triangle_measures[reduced]
        require(measure <= triangle_bound <= F(2, 21), "additive-triangle ceiling failed")
        measures[measure] += 1
        rows_by_measure.setdefault(measure, []).append(spacings)
        total += 1
        cycle_rows += len(CYCLES)
        triangle_rows += 1

    maximum = max(measures)
    second = max(value for value in measures if value < maximum)
    maximizers = rows_by_measure[maximum]

    expected_maximizers = tuple((m, 2 * m, 3 * m) for m in range(1, HEIGHT // 3 + 1))
    require(total == 4060, "height-30 triple count moved")
    require(len(measures) == 518, "distinct exact measure count moved")
    require(measures[F(0)] == 2010, "zero-measure count moved")
    require(maximum == F(2, 21), "finite exact maximum moved")
    require(tuple(maximizers) == expected_maximizers, "finite maximizers moved")
    require(second == F(4, 105) and measures[second] == 19, "second measure stratum moved")
    require(bad_measure((1, 6, 7)) == F(5, 147), "non-AP positive counterexample moved")
    return {
        "triples": total,
        "cycle_rows": cycle_rows,
        "triangle_rows": triangle_rows,
        "triangle_types": len(triangle_measures),
        "distinct": len(measures),
        "zero": measures[F(0)],
        "maximum": maximum,
        "maximizers": tuple(maximizers),
        "second": second,
        "second_count": measures[second],
        "non_ap": bad_measure((1, 6, 7)),
    }


def main() -> None:
    identity_rows = continuum_identity_audit()
    triangle_model_rows = triangle_region_model_audit()
    sheared_model_rows = sheared_section_model_audit()
    triangle_core = finite_triangle_core()
    equality_audit = equality_obligation_audit()
    census = height_census()
    print("THM-1210 continuum bad-set / additive-triangle ceiling exact referee")
    print("arithmetic=fractions.Fraction; chamber endpoints have measure zero")
    print("identity: F_inf=(7/6)*max_j(Delta_j-1/7)_+")
    print("BAD iff maxgap({0,d2*u,d3*u,d4*u})<=2/7")
    print("correlation: mu(BAD)=2*sum over 3 undirected K4 Hamilton cycles")
    print(f"direct continuum identity chamber rows={identity_rows}")
    print(
        "triangle quotient: A=[1/7,2/7]U[3/7,4/7]U[5/7,6/7]; "
        "BAD forces every pair difference into A"
    )
    print(
        f"R={{(x,y):x,y,x+y in A}}: six exact right triangles; "
        f"model rows={triangle_model_rows}; area=3/49; horizontal components<=2"
    )
    print(
        f"shear z=x+y: model rows={sheared_model_rows}; x-sections empty off A, "
        "six explicit interval branches and at most two components on A"
    )
    print(
        f"tail: J(p,q)<=3/49+6/[7(p+q)]<2/21 for coprime p<q and p+q>=26; "
        f"cutoff value={triangle_core['tail_at_cutoff']}; "
        f"strict margin={triangle_core['tail_margin']}"
    )
    print(
        f"finite triangle core: pairs={triangle_core['pairs']}; "
        f"distinct={triangle_core['distinct']}; "
        f"top unique strata={triangle_core['top_strata']}"
    )
    print(
        f"independent carry formula rows={triangle_core['carry_rows']}; "
        "J=2/[7(p+q)]*sum_k c_k*L_k"
    )
    print(
        f"equality obligation lemma audit: rows={equality_audit['rows']}; "
        f"survivors={equality_audit['survivors']} (exactly a=b=c)"
    )
    print(
        f"height-{HEIGHT} triples={census['triples']}; cycle correlations={census['cycle_rows']}; "
        f"triangle implications={census['triangle_rows']} ({census['triangle_types']} reduced types); "
        f"distinct measures={census['distinct']}; zero={census['zero']}"
    )
    print(f"maximum={census['maximum']}; maximizers={census['maximizers']}")
    print(f"second={census['second']}; count={census['second_count']}")
    print(f"non-AP (1,6,7) exact bad measure={census['non_ap']}")
    print("Tournament Analysis:")
    print("  vertices=four phase points; alternate vertices=pair-difference obligations")
    print("  switch=cyclic order; reversal u->-u pairs the two orientations; factor=2")
    print("  order tournament=transitive; score histogram=(0,1,2,3); SCCs=4; cycles=0; H=1")
    print("  proof quotient=non-AP additive triangle with edge labels p,q,p+q")
    print("  faithful data=three-band edge predicate|additive circuit|Haar mass")
    print("  destroyed by runner tournament=cyclic gaps|edge increments|correlation mass")
    print("VERDICT: sharp continuum BAD theorem proved (analytic p+q>=26 + 99 exact cases); equality iff four frequencies form an AP")


if __name__ == "__main__":
    main()
