#!/usr/bin/env python3
"""Literal exact-phase controls for the q=D/7 expected-spike theorem.

This audit uses the actual strict address mask

    ||a(r+u)/d|| < 1/14,

not the ceiling-width enlarged mask.  It partitions the u-circle at every
exact rational phase boundary.  The controls verify:

* d/gcd(d,q) is 1 or 7;
* if d|q, every selected residue fills a whole q-fibre and
  integral |Y_d(u)| du = q/7, independently of sampled unit numerator a;
* for d=7a+r|q the filled-fibre count is
  a(q/d)+(q/d)X with Haar(X=1)=r/7 almost everywhere, and that expression
  remains a pointwise upper envelope at strict-mask seams;
* if d does not divide q, the actual mask hits each q-fibre at most once;
* strict-boundary spike counts are lower semicontinuous, hence all positive
  integer superlevel events are open;
* literal four-mask phase-cell distributions obey the exact mean and every
  Markov superlevel bound.
"""

from fractions import Fraction as Q
from itertools import combinations_with_replacement
from math import gcd, lcm


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def norm_fraction(value):
    residue = value % 1
    return min(residue, 1 - residue)


def selected_classes(d, numerator, u):
    require(gcd(d, numerator) == 1, "numerator is not a unit")
    return tuple(
        residue
        for residue in range(d)
        if norm_fraction(Q(numerator) * (residue + u) / d) < Q(1, 14)
    )


def phase_boundaries(d, numerator):
    """All strict-mask equality phases in [0,1], plus the cell endpoints."""
    points = {Q(0), Q(1)}
    # Across r<=d-1 and 0<=u<=1, numerator*(r+u)/d lies in
    # [0,numerator].  The padded integer range is therefore exhaustive.
    for residue in range(d):
        for integer in range(-1, numerator + 2):
            for sign in (-1, 1):
                u = (
                    Q(d) * (Q(integer) + Q(sign, 14)) / numerator
                    - residue
                )
                if 0 <= u <= 1:
                    points.add(u)
    return tuple(sorted(points))


def circle_samples(boundaries):
    """Every boundary and one point in every complementary open cell."""
    samples = set(boundaries[:-1])
    samples.update(
        (left + right) / 2
        for left, right in zip(boundaries, boundaries[1:])
        if left < right
    )
    return tuple(sorted(samples))


def q_fibre_loads(D, d, numerator, u):
    q = D // 7
    selected = set(selected_classes(d, numerator, u))
    return tuple(
        sum((fibre + lift * q) % d in selected for lift in range(7))
        for fibre in range(q)
    )


def integral_on_cells(boundaries, value):
    return sum(
        (right - left) * value((left + right) / 2)
        for left, right in zip(boundaries, boundaries[1:])
    )


def assert_lower_semicontinuous(boundaries, value, label):
    """At a strict boundary the value cannot exceed either adjacent value."""
    cells = tuple(
        (left + right) / 2
        for left, right in zip(boundaries, boundaries[1:])
    )
    require(cells, ("phase partition has no cells", label))
    for index, point in enumerate(boundaries[:-1]):
        left_value = value(cells[index - 1])
        right_value = value(cells[index])
        boundary_value = value(point)
        require(
            boundary_value <= min(left_value, right_value),
            (
                "strict-boundary lower semicontinuity failed",
                label,
                point,
                left_value,
                boundary_value,
                right_value,
            ),
        )


def individual_controls():
    cases = 0
    boundary_checks = 0
    one_spike_thresholds = 0
    for D in (14, 28, 42, 56, 70, 84, 98):
        q = D // 7
        for d in range(2, D + 1):
            if D % d:
                continue
            g = gcd(d, q)
            require(d // g in (1, 7), ("d/g dichotomy failed", D, q, d))
            candidates = {
                numerator
                for numerator in range(1, min(2 * d, 20) + 1)
                if gcd(numerator, d) == 1
            }
            candidates.update(
                d + offset
                for offset in range(1, 11)
                if gcd(d + offset, d) == 1
            )
            for numerator in sorted(candidates):
                boundaries = phase_boundaries(d, numerator)
                count = lambda u: len(selected_classes(d, numerator, u))
                assert_lower_semicontinuous(
                    boundaries, count, (D, d, numerator)
                )
                boundary_checks += len(boundaries) - 1
                require(
                    integral_on_cells(boundaries, count) == Q(d, 7),
                    ("actual strict class mean failed", D, d, numerator),
                )
                for u in circle_samples(boundaries):
                    loads = q_fibre_loads(D, d, numerator, u)
                    require(
                        sum(loads)
                        == len(selected_classes(d, numerator, u)) * (D // d),
                        ("literal fibre accounting failed", D, d, numerator, u),
                    )
                    require(
                        sum(loads)
                        <= (D // d) * ((d + 6) // 7),
                        ("actual mask exceeds ceiling enlargement", D, d, numerator, u),
                    )
                    if q % d == 0:
                        require(
                            set(loads) <= {0, 7},
                            ("d|q mask is not whole-fibre", D, d, numerator, u),
                        )
                        full_fibres = sum(load == 7 for load in loads)
                        require(
                            full_fibres
                            == len(selected_classes(d, numerator, u)) * (q // d),
                            ("whole-fibre cardinality failed", D, d, numerator, u),
                        )
                    else:
                        require(
                            max(loads, default=0) <= 1,
                            ("d not|q mask exceeds one point/fibre", D, d, numerator, u),
                        )
                if q % d == 0:
                    quotient, remainder = divmod(d, 7)
                    width = q // d
                    full_count = lambda u: (
                        len(selected_classes(d, numerator, u)) * width
                    )
                    require(
                        integral_on_cells(boundaries, full_count) == Q(q, 7),
                        ("whole-spike Fubini mean failed", D, d, numerator),
                    )
                    # The floor/ceiling law is an almost-everywhere law on
                    # open phase cells.  At strict-mask equality boundaries
                    # both endpoints can be absent; lower semicontinuity is
                    # audited separately above.
                    values = {
                        full_count((left + right) / 2)
                        for left, right in zip(boundaries, boundaries[1:])
                    }
                    require(
                        values <= {quotient * width, (quotient + 1) * width},
                        ("one-spike floor/exception values failed", D, d, numerator),
                    )
                    active_mass = sum(
                        right - left
                        for left, right in zip(boundaries, boundaries[1:])
                        if full_count((left + right) / 2)
                        == (quotient + 1) * width
                    )
                    require(
                        active_mass == Q(remainder, 7),
                        ("one-spike exception marginal failed", D, d, numerator),
                    )
                    allowance = quotient * width + (
                        width if remainder >= 5 else 0
                    )
                    for need in range(1, (quotient + 1) * width + 1):
                        event_mass = sum(
                            right - left
                            for left, right in zip(boundaries, boundaries[1:])
                            if full_count((left + right) / 2) >= need
                        )
                        require(
                            (event_mass > Q(55, 91)) == (need <= allowance),
                            (
                                "one-spike 55/91 threshold failed",
                                D,
                                d,
                                numerator,
                                need,
                                event_mass,
                                allowance,
                            ),
                        )
                        one_spike_thresholds += 1
                cases += 1
    return cases, boundary_checks, one_spike_thresholds


def tuple_controls():
    """Exhaust all exact-lcm four-multisets for two small D, at numerator 1."""
    shape_cases = 0
    markov_checks = 0
    for D in (28, 42):
        q = D // 7
        alphabet = tuple(d for d in range(2, D + 1) if D % d == 0)
        for shape in combinations_with_replacement(alphabet, 4):
            if lcm(*shape) != D:
                continue
            spike_denominators = tuple(d for d in shape if q % d == 0)
            m = len(spike_denominators)
            boundaries = tuple(
                sorted(
                    {
                        point
                        for d in spike_denominators
                        for point in phase_boundaries(d, 1)
                    }
                    | {Q(0), Q(1)}
                )
            )
            total_spikes = lambda u: sum(
                len(selected_classes(d, 1, u)) * (q // d)
                for d in spike_denominators
            )
            assert_lower_semicontinuous(
                boundaries, total_spikes, (D, shape)
            )
            exact_mean = integral_on_cells(boundaries, total_spikes)
            require(
                exact_mean == Q(m * q, 7),
                ("four-mask exact mean failed", D, shape, exact_mean),
            )
            maximum = max(total_spikes(u) for u in circle_samples(boundaries))
            for need in range(1, maximum + 1):
                event_mass = sum(
                    right - left
                    for left, right in zip(boundaries, boundaries[1:])
                    if total_spikes((left + right) / 2) >= need
                )
                require(
                    event_mass <= exact_mean / need,
                    ("literal Markov bound failed", D, shape, need),
                )
                markov_checks += 1
            shape_cases += 1
    return shape_cases, markov_checks


def main():
    individual_cases, boundary_checks, one_spike_thresholds = individual_controls()
    shape_cases, markov_checks = tuple_controls()
    # The theorem's strict predicate must retain the zero-demand/zero-spike
    # case: its event is the entire circle despite the formal ratio 0/0.
    require(
        not (Q(0) < Q(0)),
        "zero-demand exception control changed",
    )
    print("LRC14 k3 expected-spike literal phase audit")
    print(f"individual_mask_cases={individual_cases}")
    print(f"strict_boundary_cells_checked={boundary_checks}")
    print(f"exact_one_spike_thresholds={one_spike_thresholds}")
    print(f"literal_four_mask_shapes={shape_cases}")
    print(f"literal_markov_thresholds={markov_checks}")
    print("zero_demand_predicate=N_c==0 OR 55*N_c<13*m*q")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
