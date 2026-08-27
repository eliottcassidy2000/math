#!/usr/bin/env python3
"""Exact companion for THM-4258's three-sample W=0 observer.

The lattice census is a corollary computation: it imports the dependency-free
THM-4249 audit and applies THM-4253's proved norm-three profile deletion.  The
F4 recurrence audit below is independent standard-library arithmetic.

Nothing here evaluates the mixed maps at the attachment coordinates.  In
particular, it does not prove either finite marked-ratio set empty, close W=0,
or prove JC(2).
"""

from __future__ import annotations

from collections import Counter
from itertools import product

import jc23_w0_cyclic_projector_squeeze_thm4249_independent_audit as base


Vector = base.Vector


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def vector_sub(left: Vector, right: Vector) -> Vector:
    return tuple(
        base.e_sub(left[index], right[index]) for index in range(4)
    )  # type: ignore[return-value]


def observer_map(vector: Vector) -> Vector:
    """F_m=(T-1)m, whose T^j translates evaluate to m(D_j)."""

    return vector_sub(base.tau(vector), vector)


SOURCE_GRAM = (
    ((4, 0), base.ZERO, base.ZERO, base.ZERO),
    (base.ZERO, (6, 0), (-4, -2), (-2, 2)),
    (base.ZERO, (-2, 2), (6, 0), (2, -2)),
    (base.ZERO, (-4, -2), (4, 2), (4, 0)),
)


def observer_gram() -> tuple[tuple[base.E, ...], ...]:
    """Pull q back by T-1 on the M0 basis [f,g,h]."""

    basis: tuple[Vector, ...] = (
        (base.ZERO, base.ONE, base.ZERO, base.ZERO),
        (base.ZERO, base.ZERO, base.ONE, base.ZERO),
        (base.ZERO, base.ZERO, base.ZERO, base.ONE),
    )
    columns = tuple(observer_map(vector) for vector in basis)
    answer = []
    for left in columns:
        row = []
        for right in columns:
            value = base.ZERO
            for i in range(4):
                for j in range(4):
                    value = base.e_add(
                        value,
                        base.e_mul(
                            left[i],
                            base.e_mul(SOURCE_GRAM[i][j], base.e_conjugate(right[j])),
                        ),
                    )
            row.append(value)
        answer.append(tuple(row))
    return tuple(answer)


def verify_observer_gram() -> tuple[tuple[base.E, ...], ...]:
    pulled = observer_gram()
    expected = (
        ((18, 0), (-20, -10), (-10, 4)),
        ((-10, 10), (18, 0), (4, -10)),
        ((-14, -4), (14, 10), (12, 0)),
    )
    require(pulled == expected, "pulled observer Gram changed")

    # Tr(a*O) is generated over Z by Tr(a) and Tr(a*omega).
    for index in range(3):
        require(pulled[index][index][1] == 0, "observer Gram diagonal left Z")
        require(pulled[index][index][0] % 6 == 0,
                "observer Gram diagonal lost six-divisibility")
        for other in range(index + 1, 3):
            entry = pulled[index][other]
            require(base.e_trace(entry) % 6 == 0,
                    "observer trace ideal lost six-divisibility")
            require(base.e_trace(base.e_mul(entry, base.OMEGA)) % 6 == 0,
                    "observer omega-trace ideal lost six-divisibility")
    return pulled


def residual_sets() -> tuple[dict[int, set[Vector]], set[Vector], set[Vector]]:
    """Reconstruct THM-4249 residuals and THM-4253's live degree-42 shell."""

    shells = base.enumerate_shells()
    pre_4253: dict[int, set[Vector]] = {}
    for degree, shell in shells.items():
        a_zero = {vector for vector in shell if vector[0] == base.ZERO}
        visible_bound = {
            vector for vector in a_zero
            if vector[3] == base.ZERO or base.e_norm(vector[3]) >= 3
        }
        hidden_bound = {
            vector for vector in visible_bound
            if base.split_coordinates(vector)[2] not in (1, 2)
        }
        pre_4253[degree] = {
            vector for vector in hidden_bound
            if not (degree == 42 and vector[3] == base.ZERO)
        }

    require(
        {degree: len(rows) for degree, rows in pre_4253.items()}
        == {34: 4224, 42: 3168},
        "THM-4249 residual changed",
    )

    deleted_42 = {
        vector for vector in pre_4253[42]
        if base.split_coordinates(vector)[1:3] == (3, 13)
    }
    post_4253_42 = pre_4253[42] - deleted_42
    require(len(deleted_42) == 672, "THM-4253 deleted profile changed")
    require(len(post_4253_42) == 2496, "THM-4253 live shell changed")

    live = {34: pre_4253[34], 42: post_4253_42}
    return live, pre_4253[42], deleted_42


def degree_histogram(vectors: set[Vector]) -> tuple[list[Vector], Counter[int]]:
    representatives, sizes = base.orbit_partition(vectors)
    require(sizes == Counter({24: len(representatives)}), "orbit size changed")
    histogram = Counter(base.full_degree(observer_map(row)) for row in representatives)

    # q((T-1)m) is constant under target units and source T-translates.
    for representative in representatives:
        expected = base.full_degree(observer_map(representative))
        current = representative
        for _ in range(12):
            for unit in base.UNITS:
                scaled = base.vector_scale(unit, current)
                require(
                    base.full_degree(observer_map(scaled)) == expected,
                    "observer degree escaped its source-target orbit",
                )
            current = base.tau(current)
    return representatives, histogram


def verify_lattice_observer() -> tuple[
    Counter[int], Counter[int], Counter[int], Counter[int], Counter[int]
]:
    live, pre_4253_42, deleted_42 = residual_sets()

    # P0(T)=(T-omega^2)(T^2+omega) annihilates exactly the a_u=0 lane.
    p_zero = base.polynomial_multiply(
        [base.e_neg(base.OMEGA2), base.ONE],
        [base.OMEGA, base.ZERO, base.ONE],
    )
    require(
        p_zero == [(-1, 0), base.OMEGA, base.e_neg(base.OMEGA2), base.ONE],
        "P0 coefficients changed",
    )
    for vectors in live.values():
        for vector in vectors:
            require(
                base.apply_tau_polynomial(p_zero, vector) == base.ZERO_VECTOR,
                "P0 failed on the a_u=0 residual",
            )

    representatives: dict[int, list[Vector]] = {}
    histograms: dict[int, Counter[int]] = {}
    for degree in (34, 42):
        representatives[degree], histograms[degree] = degree_histogram(live[degree])

    expected_34 = Counter({
        30: 8, 36: 4, 42: 10, 48: 6, 54: 7, 60: 8, 66: 21, 72: 6,
        78: 19, 84: 13, 90: 16, 96: 10, 102: 20, 108: 9, 114: 17,
        120: 2,
    })
    expected_42 = Counter({
        36: 1, 42: 2, 54: 8, 66: 6, 72: 6, 78: 3, 84: 4, 90: 14,
        96: 4, 102: 10, 108: 9, 114: 4, 120: 4, 126: 12, 132: 4,
        138: 10, 144: 2, 150: 1,
    })
    require(histograms == {34: expected_34, 42: expected_42},
            "post-sieve observer histogram changed")

    _, pre_histogram_42 = degree_histogram(pre_4253_42)
    _, deleted_histogram_42 = degree_histogram(deleted_42)
    expected_pre_42 = Counter({
        30: 2, 36: 3, 42: 2, 48: 4, 54: 10, 66: 8, 72: 8, 78: 3,
        84: 4, 90: 14, 96: 4, 102: 12, 108: 11, 114: 4, 120: 6,
        126: 16, 132: 4, 138: 12, 144: 4, 150: 1,
    })
    expected_deleted_42 = Counter({
        30: 2, 36: 2, 48: 4, 54: 2, 66: 2, 72: 2, 102: 2, 108: 2,
        120: 2, 126: 4, 138: 2, 144: 2,
    })
    require(pre_histogram_42 == expected_pre_42,
            "pre-THM-4253 histogram changed")
    require(deleted_histogram_42 == expected_deleted_42,
            "deleted-profile histogram changed")
    require(pre_histogram_42 == expected_42 + expected_deleted_42,
            "pre/post histogram subtraction failed")

    # Weight each map orbit by its already-proved live marked-ratio count.
    pi = base.e_sub(base.OMEGA2, base.ONE)
    excluded_small = base.torsion_kernel(pi) | base.torsion_kernel((2, 0))
    incidence_histograms: dict[int, Counter[int]] = {34: Counter(), 42: Counter()}
    ratio_count_profiles: dict[int, Counter[int]] = {34: Counter(), 42: Counter()}
    for degree in (34, 42):
        for vector in representatives[degree]:
            annihilator = base.e_mul(pi, vector[3])
            raw_count = sum(
                orbit.isdisjoint(excluded_small)
                for orbit in base.torsion_orbits(base.torsion_kernel(annihilator))
            )
            # THM-4249 removes the common ratio 1/3 from every degree-42 map.
            live_count = raw_count if degree == 34 else raw_count - 1
            require(live_count > 0, "live ratio count became nonpositive")
            degree_f = base.full_degree(observer_map(vector))
            incidence_histograms[degree][degree_f] += live_count
            ratio_count_profiles[degree][live_count] += 1

    expected_incidence_34 = Counter({
        30: 12, 36: 12, 42: 18, 48: 30, 54: 33, 60: 36, 66: 99,
        72: 42, 78: 90, 84: 84, 90: 114, 96: 60, 102: 105, 108: 60,
        114: 63, 120: 6,
    })
    expected_incidence_42 = Counter({
        36: 3, 42: 6, 54: 30, 66: 24, 72: 42, 78: 9, 84: 36, 90: 90,
        96: 30, 102: 60, 108: 69, 114: 30, 120: 36, 126: 84, 132: 30,
        138: 60, 144: 6, 150: 3,
    })
    require(
        incidence_histograms == {34: expected_incidence_34, 42: expected_incidence_42},
        "incidence-weighted observer histogram changed",
    )
    require(sum(expected_incidence_34.values()) == 864,
            "degree-34 incidence total changed")
    require(sum(expected_incidence_42.values()) == 648,
            "degree-42 incidence total changed")
    require(ratio_count_profiles == {
        34: Counter({1: 36, 3: 52, 6: 32, 7: 24, 9: 24, 12: 8}),
        42: Counter({3: 24, 4: 36, 9: 32, 12: 12}),
    }, "map-ratio profile changed")

    return (
        expected_34, expected_42, expected_deleted_42,
        expected_incidence_34, expected_incidence_42,
    )


# F4 is represented by a+b*w as a+2b, with w^2=w+1 in characteristic two.
F4_ZERO = 0
F4_ONE = 1
F4_OMEGA = 2
F4_OMEGA2 = 3


def f4_add(left: int, right: int) -> int:
    return left ^ right


def f4_mul(left: int, right: int) -> int:
    a, b = left & 1, (left >> 1) & 1
    c, d = right & 1, (right >> 1) & 1
    return ((a * c + b * d) & 1) | (((a * d + b * c + b * d) & 1) << 1)


def f4_pow(value: int, exponent: int) -> int:
    answer = F4_ONE
    for _ in range(exponent):
        answer = f4_mul(answer, value)
    return answer


def f4_scale(bit: int, value: int) -> int:
    return value if bit & 1 else F4_ZERO


def recurrence_sequence(c0: int, c1: int, c2: int) -> tuple[int, ...]:
    """delta_j=(omega^2)^j(c0+j*c1+binom(j,2)c2)."""

    answer = []
    for index in range(12):
        epsilon = f4_add(
            c0,
            f4_add(
                f4_scale(index, c1),
                f4_scale(index * (index - 1) // 2, c2),
            ),
        )
        answer.append(f4_mul(f4_pow(F4_OMEGA2, index), epsilon))
    return tuple(answer)


def least_period(sequence: tuple[int, ...]) -> int:
    for period in (1, 2, 3, 4, 6, 12):
        if all(sequence[index] == sequence[(index + period) % 12]
               for index in range(12)):
            return period
    raise RuntimeError("sequence lost period twelve")


def verify_f4_recurrence() -> tuple[Counter[int], tuple[int, ...]]:
    require(f4_mul(F4_OMEGA, F4_OMEGA) == F4_OMEGA2,
            "F4 omega square changed")
    require(f4_mul(F4_OMEGA2, F4_OMEGA2) == F4_OMEGA,
            "F4 omega-fourth-power changed")

    sequences = {
        recurrence_sequence(c0, c1, c2)
        for c0, c1, c2 in product(range(4), repeat=3)
    }
    require(len(sequences) == 64, "three coefficients lost injectivity")

    prefixes_three = Counter(sequence[:3] for sequence in sequences)
    prefixes_two = Counter(sequence[:2] for sequence in sequences)
    require(len(prefixes_three) == 64 and set(prefixes_three.values()) == {1},
            "three samples are not complete")
    require(len(prefixes_two) == 16 and set(prefixes_two.values()) == {4},
            "two-sample fibre profile changed")

    for sequence in sequences:
        total = F4_ZERO
        for value in sequence:
            total = f4_add(total, value)
        require(total == F4_ZERO, "difference cycle lost telescoping sum")

        # P0(S)=S^3-omega^2*S^2+omega*S-1; signs vanish in char two.
        for index in range(12):
            recurrence = f4_add(
                sequence[(index + 3) % 12],
                f4_add(
                    f4_mul(F4_OMEGA2, sequence[(index + 2) % 12]),
                    f4_add(
                        f4_mul(F4_OMEGA, sequence[(index + 1) % 12]),
                        sequence[index],
                    ),
                ),
            )
            require(recurrence == F4_ZERO, "cubic recurrence failed")

        # Twisting epsilon_j=omega^j*delta_j gives Delta^3 epsilon=0.
        epsilon = tuple(
            f4_mul(f4_pow(F4_OMEGA, index), sequence[index])
            for index in range(12)
        )
        difference = epsilon
        for _ in range(3):
            difference = tuple(
                f4_add(difference[(index + 1) % 12], difference[index])
                for index in range(12)
            )
        require(set(difference) == {F4_ZERO}, "twisted third difference failed")

    hostile = recurrence_sequence(F4_ZERO, F4_ZERO, F4_ONE)
    require(hostile[0] == hostile[1] == F4_ZERO and hostile[2] != F4_ZERO,
            "two-sample hostile disappeared")
    periods = Counter(least_period(sequence) for sequence in sequences)
    return periods, hostile


def format_counter(counter: Counter[int]) -> str:
    return "{" + ",".join(f"{key}:{counter[key]}" for key in sorted(counter)) + "}"


def main() -> None:
    pulled_gram = verify_observer_gram()
    (histogram_34, histogram_42, deleted_42,
     incidence_34, incidence_42) = verify_lattice_observer()
    periods, hostile = verify_f4_recurrence()

    print("THM-4258 W=0 three-sample attachment observer exact certificate")
    print("dependencies=THM-4249_residual+THM-4253_norm3_profile_deletion")
    print("operator=P0(T)=(T-omega^2)(T^2+omega) annihilates a_u=0 PASS")
    print("attachment_recurrence=delta[j+3]-omega^2*delta[j+2]+omega*delta[j+1]-delta[j]=O")
    print("three_consecutive_zero_deltas=all_twelve_zero_deltas IFF")
    print("observer=F_m=(T-1)m samples=(F_m,T*F_m,T^2*F_m)at_Q0")
    print(f"observer_gram_M0={pulled_gram} q(F_m)_in_6Z=PASS spectral_bound=(2+sqrt(3))*q(m)")
    print(f"degree34_map_orbits={sum(histogram_34.values())} observer_degree_hist={format_counter(histogram_34)} max=120")
    print(f"degree42_live_map_orbits={sum(histogram_42.values())} observer_degree_hist={format_counter(histogram_42)} max=150")
    print(f"degree42_deleted_orbits={sum(deleted_42.values())} observer_degree_hist={format_counter(deleted_42)}")
    print(f"degree34_live_incidences={sum(incidence_34.values())} incidence_degree_hist={format_counter(incidence_34)}")
    print(f"degree42_live_incidences={sum(incidence_42.values())} incidence_degree_hist={format_counter(incidence_42)}")
    print("uniform_observer_table_rows naive=1512*12=18144 three_sample=1512*3=4536 reduction=13608 minimality=NOT_CLAIMED")
    print("two_torsion_condition=after_visible_and_hidden_projector_collapse_only")
    print(f"F4_abstract_recurrence_patterns=64 telescoping_sum=PASS period_profile={format_counter(periods)}")
    print(f"two_sample_hostile_delta={hostile} first_nonzero_index=2")
    print("abstract_64_pattern_surjectivity_from_geometric_maps=NOT_CLAIMED")
    print("mixed_attachment_values=UNEVALUATED normalization_dictionary=NOT_IN_THIS_THEOREM")
    print("verdict=EXACT_OBSERVER_COMPRESSION JC2_OPEN W0_OPEN")


if __name__ == "__main__":
    main()
