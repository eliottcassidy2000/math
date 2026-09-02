#!/usr/bin/env python3
"""Exact half-turn occupation-depth dual for the LRC(14) minority branch.

For S={2h} union W with twelve distinct positive odd tails, pair the first
half of every anchor-safe component with its translate by 1/2.  If a(t) and
b(t) are the two tail-danger multiplicities, then no tail contributes to
both, hence q(t)=min(a(t),b(t)) lies in {0,...,6}.

Any polynomial p(q)=sum c_j binom(q,j) with p(0)<=1 and p(r)<=0 for
r=1,...,6 is therefore a pointwise lower dual for the number of uncovered
members of the pair.  This program checks the universal quadratic and cubic
duals and evaluates them by an exact rational wall sweep.

The computation is a certificate-format probe, not a proof of LRC(14).
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import ceil, comb, floor, gcd


@dataclass(frozen=True)
class Interval:
    left: Fraction
    right: Fraction


def anchor_component(h: int, k: int) -> Interval:
    return Interval(
        Fraction(14 * k + 1, 28 * h),
        Fraction(14 * k + 13, 28 * h),
    )


def danger_intervals_meeting(w: int, target: Interval) -> tuple[Interval, ...]:
    """Open components of ||w t||<1/14 meeting target's interior."""
    lo = floor(w * target.left - Fraction(1, 14))
    hi = ceil(w * target.right + Fraction(1, 14))
    answer = []
    for n in range(lo - 1, hi + 2):
        tooth = Interval(
            Fraction(14 * n - 1, 14 * w),
            Fraction(14 * n + 1, 14 * w),
        )
        if tooth.left < target.right and target.left < tooth.right:
            answer.append(tooth)
    return tuple(answer)


def clipped(interval: Interval, target: Interval, shift: Fraction) -> Interval | None:
    left = max(target.left, interval.left - shift)
    right = min(target.right, interval.right - shift)
    if left >= right:
        return None
    return Interval(left, right)


def dual_value(q: int, coefficients: tuple[Fraction, ...]) -> Fraction:
    return sum(
        (coefficient * comb(q, degree) for degree, coefficient in enumerate(coefficients)),
        start=Fraction(),
    )


QUADRATIC = (Fraction(1), Fraction(-1), Fraction(1, 3))
CUBIC = (Fraction(1), Fraction(-1), Fraction(1), Fraction(-1))
EXACT_ZERO = tuple(Fraction((-1) ** degree) for degree in range(7))


def audit_pointwise_duals() -> None:
    quadratic_values = tuple(dual_value(q, QUADRATIC) for q in range(7))
    cubic_values = tuple(dual_value(q, CUBIC) for q in range(7))
    exact_values = tuple(dual_value(q, EXACT_ZERO) for q in range(7))
    target = (1, 0, 0, 0, 0, 0, 0)

    assert quadratic_values == (
        Fraction(1),
        Fraction(0),
        Fraction(-2, 3),
        Fraction(-1),
        Fraction(-1),
        Fraction(-2, 3),
        Fraction(0),
    )
    assert cubic_values == (
        Fraction(1),
        Fraction(0),
        Fraction(0),
        Fraction(0),
        Fraction(-1),
        Fraction(-4),
        Fraction(-10),
    )
    assert exact_values == target
    assert all(
        cubic_values[q] == -Fraction((q - 1) * (q - 2) * (q - 3), 6)
        for q in range(7)
    )
    assert all(
        cubic_values[q] >= Fraction(1) - Fraction(11, 6) * q
        for q in range(7)
    )

    # If the normalization p(0)=1 and p(1)=0 is retained, a quadratic of
    # the form 1-q+c*C(q,2) must have c<=1/3: q=6 is the sharp constraint.
    assert dual_value(6, (Fraction(1), Fraction(-1), Fraction(1, 3))) == 0
    assert dual_value(6, (Fraction(1), Fraction(-1), Fraction(1, 3) + Fraction(1, 1000))) > 0

    print("pointwise_q=0..6")
    print(" quadratic=" + ",".join(map(str, quadratic_values)))
    print(" cubic=" + ",".join(map(str, cubic_values)))
    print(" exact_degree6=" + ",".join(map(str, exact_values)))
    print(" quadratic_factor=(q-1)(q-6)/6")
    print(" cubic_factor=-(q-1)(q-2)(q-3)/6")
    print(" cubic_qmass_identity=H0-H4-4H5-10H6")
    print(" cubic_sharp_chord=1-(11/6)q equality_at_q=0,6")
    print(" quadratic_coefficient_1_over_3=SHARP_AT_q6")


def abstract_marginal_extremizer() -> None:
    """Exact abstract obstruction and current thresholds.

    Normalize the full half-turn base [0,1/2) to probability one.  With
    probability 2/7 choose a uniform six-label subset for the low sheet and
    put its complement on the high sheet; otherwise all labels are absent.
    This has the exact one-sheet marginals of twelve odd danger combs, but it
    is deliberately only an abstract coupling, not a realizable clock row.
    """
    absent = Fraction(5, 7)
    balanced = Fraction(2, 7)
    per_label_low = balanced * Fraction(6, 12)
    per_label_high = balanced * Fraction(6, 12)
    expected_q = balanced * 6
    cubic = absent * dual_value(0, CUBIC) + balanced * dual_value(6, CUBIC)

    assert per_label_low == per_label_high == Fraction(1, 7)
    assert expected_q == Fraction(12, 7)
    assert cubic == Fraction(-15, 7)

    # The same balanced state can be placed wholly inside an abstract core of
    # nu-mass 6/7, with q=0 on the other 4/7 of that core and the deleted
    # anchor strip (nu-mass 1/7) empty.  This retains all full-base marginals.
    core_absent = Fraction(4, 7)
    core_cubic = core_absent + balanced * dual_value(6, CUBIC)
    assert core_absent + balanced == Fraction(6, 7)
    assert core_cubic == Fraction(-16, 7)

    # On the full probability base, E(a+b)=24/7.  On the anchor-safe core,
    # use nu=2 dt, a subprobability measure of mass 6/7; deleting anchor-bad
    # points can only reduce the total tail load from 24/7.
    full_current_threshold = Fraction(180, 77)
    core_subprobability_current_threshold = Fraction(192, 77)
    core_physical_current_threshold = Fraction(96, 77)
    core_conditional_current_threshold = Fraction(32, 11)
    assert (
        -Fraction(15, 7) + Fraction(11, 12) * full_current_threshold == 0
    )
    assert (
        -Fraction(16, 7)
        + Fraction(11, 12) * core_subprobability_current_threshold
        == 0
    )
    assert core_physical_current_threshold * 2 == core_subprobability_current_threshold
    assert (
        core_subprobability_current_threshold / Fraction(6, 7)
        == core_conditional_current_threshold
    )

    print("abstract_exact_marginal_extremizer")
    print(" full_base_normalization=probability_one_on_[0,1/2)")
    print(" state_masses=q0:5/7,q6:2/7")
    print(" each_label_low_marginal=1/7 each_label_high_marginal=1/7")
    print(" cubic_expectation=-15/7 SHARP_FROM_MARGINALS")
    print(" abstract_anchor_core_state_masses=q0:4/7,q6:2/7 removed:1/7")
    print(" abstract_anchor_core_cubic_nu=-16/7 SHARP_FROM_MARGINALS")
    print(" cubic_positive_iff=H0>H4+4H5+10H6")
    print(" full_base_mean_current_sufficient_strictly_gt=180/77")
    print(" anchor_core_nu_equals_2dt_mass=6/7 tail_load_at_most=24/7")
    print(" anchor_core_nu_current_sufficient_strictly_gt=192/77")
    print(" anchor_core_dt_current_sufficient_strictly_gt=96/77")
    print(" anchor_core_conditional_mean_current_sufficient_strictly_gt=32/11")


def profile(h: int, speeds: tuple[int, ...], label: str) -> dict[str, Fraction | int]:
    assert len(speeds) == len(set(speeds)) == 12
    assert all(w > 0 and w % 2 == 1 for w in speeds)

    q_masses = [Fraction() for _ in range(7)]
    free_mass = Fraction()
    both_zero_mass = Fraction()
    total_mass = Fraction()
    tail_load = Fraction()
    absolute_current = Fraction()

    for k in range(h):
        base = anchor_component(h, k)
        partner = anchor_component(h, k + h)
        assert partner.left == base.left + Fraction(1, 2)
        assert partner.right == base.right + Fraction(1, 2)

        low: list[Interval] = []
        high: list[Interval] = []
        for w in speeds:
            for tooth in danger_intervals_meeting(w, base):
                part = clipped(tooth, base, Fraction())
                if part is not None:
                    low.append(part)
            for tooth in danger_intervals_meeting(w, partner):
                part = clipped(tooth, base, Fraction(1, 2))
                if part is not None:
                    high.append(part)

        walls = {base.left, base.right}
        for interval in low + high:
            walls.add(interval.left)
            walls.add(interval.right)
        ordered = sorted(walls)

        for left, right in zip(ordered, ordered[1:]):
            if left == right:
                continue
            midpoint = (left + right) / 2
            a = sum(interval.left < midpoint < interval.right for interval in low)
            b = sum(interval.left < midpoint < interval.right for interval in high)
            # For odd w, danger at t+1/2 means ||wt||>3/7, which is disjoint
            # from danger at t.  The interval implementation audits this
            # aggregate consequence on every exact cell.
            assert a + b <= 12
            q = min(a, b)
            assert 0 <= q <= 6
            width = right - left
            total_mass += width
            q_masses[q] += width
            free_mass += width * ((a == 0) + (b == 0))
            both_zero_mass += width * (a == 0 and b == 0)
            tail_load += width * (a + b)
            absolute_current += width * abs(a - b)

    assert total_mass == Fraction(3, 7)
    assert free_mass == q_masses[0] + both_zero_mass
    assert tail_load <= Fraction(12, 7)

    moments = tuple(
        sum(
            (mass * comb(q, degree) for q, mass in enumerate(q_masses)),
            start=Fraction(),
        )
        for degree in range(7)
    )

    def integrated(coefficients: tuple[Fraction, ...]) -> Fraction:
        return sum(
            (
                coefficient * moments[degree]
                for degree, coefficient in enumerate(coefficients)
            ),
            start=Fraction(),
        )

    quadratic = integrated(QUADRATIC)
    cubic = integrated(CUBIC)
    exact_zero = integrated(EXACT_ZERO)
    current_chord_lower = (
        total_mass - Fraction(11, 12) * (tail_load - absolute_current)
    )
    assert quadratic <= free_mass
    assert cubic <= free_mass
    assert exact_zero == q_masses[0] <= free_mass
    assert current_chord_lower <= cubic

    missing_denominators = tuple(
        modulus
        for modulus in range(2, 15)
        if not any(speed % modulus == 0 for speed in (2 * h,) + speeds)
    )
    fixed_halfturn_blockers = tuple(
        w
        for w in speeds
        if 12 * h < (w % (28 * h)) < 16 * h
    )

    print(label)
    print(f" h={h} primitive={gcd(2 * h, *speeds)}")
    print(" speeds=" + ",".join(map(str, speeds)))
    print(" missing_denominators=" + (",".join(map(str, missing_denominators)) or "none"))
    print(
        " fixed_halfturn_blockers="
        + (",".join(map(str, fixed_halfturn_blockers)) or "none")
    )
    print(" q_masses=" + ",".join(f"{q}:{mass}" for q, mass in enumerate(q_masses)))
    print(f" q_max={max(q for q, mass in enumerate(q_masses) if mass)}")
    print(f" free_mass={free_mass}")
    print(f" paired_at_least_one_free_mass={q_masses[0]}")
    print(f" paired_both_free_mass={both_zero_mass}")
    print(f" tail_load={tail_load}")
    print(f" absolute_sheet_current={absolute_current}")
    print(f" current_chord_lower={current_chord_lower}")
    print(f" quadratic_dual={quadratic}")
    print(f" cubic_dual={cubic}")
    print(f" exact_degree6_dual={exact_zero}")

    return {
        "free_mass": free_mass,
        "q_zero": q_masses[0],
        "both_zero": both_zero_mass,
        "tail_load": tail_load,
        "absolute_current": absolute_current,
        "current_chord_lower": current_chord_lower,
        "quadratic": quadratic,
        "cubic": cubic,
        "q_max": max(q for q, mass in enumerate(q_masses) if mass),
    }


def canonical_ap_scope_guard() -> None:
    speeds = tuple(range(1, 14))
    tails = tuple(w for w in speeds if w != 2)
    assert len(tails) == 12 and any(w % 2 == 0 for w in tails)
    # The even tail 4 is dangerous at both members of the half-turn pair.
    first_phase = Fraction()
    second_phase = Fraction(1, 2)
    first_distance = min(
        (4 * first_phase) % 1, 1 - (4 * first_phase) % 1
    )
    second_distance = min(
        (4 * second_phase) % 1, 1 - (4 * second_phase) % 1
    )
    assert first_distance == second_distance == 0 < Fraction(1, 14)
    # The usual AP witness is equality-safe, so this scope guard is a genuine
    # tight-row warning rather than a claim that the AP is unsafe.
    witness = Fraction(1, 14)
    assert min(min((w * witness) % 1, 1 - (w * witness) % 1) for w in speeds) == Fraction(1, 14)
    print("canonical_AP_1_to_13_scope_guard")
    print(" odd_tail_hypothesis=FALSE even_tail_4_hits_both_halfturns_at_t0")
    print(" equality_witness=1/14 clearance=1/14")


def exact_small_bank_cubic_census() -> None:
    """Complete exact p3 census for 12-subsets of the first 15 odds."""
    pool = tuple(range(1, 30, 2))
    subsets = tuple(combinations(range(len(pool)), 12))
    assert len(pool) == 15 and len(subsets) == 455
    best: tuple[Fraction, int, tuple[int, ...]] | None = None
    checked = 0

    for h in range(1, 31):
        # Aggregate exact cell masses by the two labelled occupancy masks.
        state_mass: dict[tuple[int, int], Fraction] = {}
        for k in range(h):
            base = anchor_component(h, k)
            partner = anchor_component(h, k + h)
            low: list[tuple[int, Interval]] = []
            high: list[tuple[int, Interval]] = []
            for index, w in enumerate(pool):
                for tooth in danger_intervals_meeting(w, base):
                    part = clipped(tooth, base, Fraction())
                    if part is not None:
                        low.append((index, part))
                for tooth in danger_intervals_meeting(w, partner):
                    part = clipped(tooth, base, Fraction(1, 2))
                    if part is not None:
                        high.append((index, part))

            walls = {base.left, base.right}
            for _, interval in low + high:
                walls.add(interval.left)
                walls.add(interval.right)
            ordered = sorted(walls)
            for left, right in zip(ordered, ordered[1:]):
                if left == right:
                    continue
                midpoint = (left + right) / 2
                low_mask = sum(
                    1 << index
                    for index, interval in low
                    if interval.left < midpoint < interval.right
                )
                high_mask = sum(
                    1 << index
                    for index, interval in high
                    if interval.left < midpoint < interval.right
                )
                assert low_mask & high_mask == 0
                key = (low_mask, high_mask)
                state_mass[key] = state_mass.get(key, Fraction()) + right - left

        assert sum(state_mass.values(), start=Fraction()) == Fraction(3, 7)
        for subset in subsets:
            mask = sum(1 << index for index in subset)
            cubic = Fraction()
            for (low_mask, high_mask), mass in state_mass.items():
                a = (low_mask & mask).bit_count()
                b = (high_mask & mask).bit_count()
                q = min(a, b)
                cubic += mass * dual_value(q, CUBIC)
            row = tuple(pool[index] for index in subset)
            candidate = (cubic, h, row)
            if best is None or candidate < best:
                best = candidate
            checked += 1

    assert best is not None
    cubic, h, row = best
    assert cubic > 0
    print("small_bank_cubic_census")
    print(f" universe=12_subsets_first_15_odds h=1..30 checked={checked}")
    print(f" minimum={cubic} h={h} speeds=" + ",".join(map(str, row)))
    print(" cubic_negative=0")


def main() -> None:
    audit_pointwise_duals()
    abstract_marginal_extremizer()
    canonical_ap_scope_guard()
    exact_small_bank_cubic_census()

    odd_ap = tuple(range(1, 24, 2))
    small_ap = profile(2, odd_ap, "odd_AP_small_quadratic_hostile")
    assert small_ap["quadratic"] < 0 < small_ap["cubic"]
    assert small_ap["q_max"] == 3
    assert small_ap["cubic"] == Fraction(2646715, 56787276)

    residual_ap = profile(420, odd_ap, "odd_AP_h420_halfturn_positive_control")
    assert residual_ap["quadratic"] < 0 < residual_ap["cubic"]
    assert residual_ap["q_max"] == 3
    assert residual_ap["cubic"] == Fraction(16046063, 189290920)

    affine_ap = (9,) + tuple(1 + 840 * j for j in range(1, 12))
    affine = profile(420, affine_ap, "h420_affine_AP_quadratic_hostile")
    assert affine["quadratic"] < 0 < affine["cubic"]
    assert affine["q_max"] == 6
    assert affine["quadratic"] == Fraction(
        -164696925438924148369133272304940284525104,
        4325895791028726433155124208050350034149135,
    )
    assert affine["cubic"] == Fraction(
        13755395733596017612242204554383223168195,
        288393052735248428877008280536690002276609,
    )

    modular_nearzero = profile(
        4,
        (1, 7, 9, 17, 25, 33, 41, 49, 57, 65, 73, 81),
        "h4_modular_cubic_nearzero_scope_hostile",
    )
    assert modular_nearzero["quadratic"] < 0 < modular_nearzero["cubic"]
    assert modular_nearzero["q_max"] == 6
    assert modular_nearzero["cubic"] == Fraction(
        373566525247, 19204120390455
    )

    large_modular_nearzero = profile(
        420,
        (1, 839, 841, 1681, 2521, 3361, 4201, 5041, 5881, 6721, 7561, 8401),
        "h420_modular_cubic_nearzero_scope_hostile",
    )
    assert large_modular_nearzero["quadratic"] < 0 < large_modular_nearzero["cubic"]
    assert large_modular_nearzero["q_max"] == 6
    assert large_modular_nearzero["cubic"] == Fraction(
        28959147561368969682651461711510842816,
        844629198923707571893175739849000114881,
    )

    height_sharp = profile(
        420,
        (1, 735, 1155, 1365, 1679, 1681, 3255, 3359, 3361, 3609, 5039, 5041),
        "h420_height_sharp_THM366_complete_control",
    )
    assert height_sharp["quadratic"] > 0 and height_sharp["cubic"] > 0
    assert height_sharp["q_max"] == 5
    assert height_sharp["quadratic"] == Fraction(
        1877339701087410921940449985333,
        133247162596101650207195390175870,
    )
    assert height_sharp["cubic"] == Fraction(
        17339729627705789770092924238,
        180551710834826084291592669615,
    )

    common = (11, 1691, 3371, 5051, 6731, 8411, 10091, 525, 945, 1365, 1575)
    first = profile(420, common + (1287,), "joint_hostile_P1287")
    second = profile(420, common + (9009,), "joint_hostile_P9009")
    assert first["quadratic"] > 0 and first["cubic"] > 0
    assert second["quadratic"] > 0 and second["cubic"] > 0
    assert first["q_max"] == second["q_max"] == 4
    assert first["quadratic"] == Fraction(
        551054571711629821461569179, 23938321074154037836892907300
    )
    assert first["cubic"] == Fraction(
        14783540775565989336572164, 153450776116372037415980175
    )
    assert second["quadratic"] == Fraction(
        37049456465357700644747, 1579774373005611947264100
    )
    assert second["cubic"] == Fraction(
        16362496058066898986258, 169261539964886994349725
    )

    print("PASS")


if __name__ == "__main__":
    main()
