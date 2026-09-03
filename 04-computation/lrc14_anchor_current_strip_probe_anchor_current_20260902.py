#!/usr/bin/env python3
"""Exact anchor-strip/current probe for the LRC(14) minority branch.

For a row ``{2h} union W`` with twelve odd tails, write on the half-base

    s_w(t) = 1_{||wt||<1/14} - 1_{||w(t+1/2)||<1/14},
    C(t)   = sum_w s_w(t).

This script sweeps every rational wall once and records the current energy
on the anchor-safe core and on its complementary danger strip.  It is meant
to test how much information the full-circle pair kernel loses when the
anchor strip is deleted.  Every reported integral uses ordinary ``dt`` on
``[0,1/2)``; endpoints are null.
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
from math import ceil, floor, gcd


def clipped_arcs(
    speed: int,
    shift: Fraction,
    radius: Fraction,
    left: Fraction = Fraction(),
    right: Fraction = Fraction(1, 2),
) -> tuple[tuple[Fraction, Fraction], ...]:
    """Arcs for ``||speed*t+shift||<radius``, clipped to a base arc."""
    lo = floor(speed * left + shift - radius) - 1
    hi = ceil(speed * right + shift + radius) + 1
    answer: list[tuple[Fraction, Fraction]] = []
    for index in range(lo, hi + 1):
        tooth_left = (index - shift - radius) / speed
        tooth_right = (index - shift + radius) / speed
        a = max(left, tooth_left)
        b = min(right, tooth_right)
        if a < b:
            answer.append((a, b))
    return tuple(answer)


def clipped_teeth(
    speed: int, shift: Fraction, left: Fraction = Fraction(), right: Fraction = Fraction(1, 2)
) -> tuple[tuple[Fraction, Fraction], ...]:
    """Danger teeth for ``||speed*(t+shift)||<1/14``, clipped to a base arc."""
    # ``speed*(t+1/2)`` is the shifted phase, not ``speed*t+1/2``.  For the
    # only nonzero shift used here speed is odd, so these agree modulo one.
    return clipped_arcs(speed, speed * shift, Fraction(1, 14), left, right)


def d14(value: int) -> int:
    residue = value % 14
    return min(residue, 14 - residue)


def circle_distance(value: Fraction) -> Fraction:
    residue = value % 1
    return min(residue, 1 - residue)


def audit_nested_wall_transducer() -> int:
    """Exact finite audit of the general radius/phase Euclidean transducer."""
    checked = 0
    for degree in range(1, 40, 2):
        for numerator in range(7):
            quotient, remainder = divmod(numerator * degree, 7)
            for phase_bit in (0, 1):
                next_phase = (degree * phase_bit + quotient) % 2
                for sample in range(60):
                    y = Fraction(2 * sample + 1, 199)
                    radius = Fraction(numerator, 14)
                    residual_radius = Fraction(remainder, 14)
                    left_distances = tuple(
                        circle_distance(
                            y + Fraction(index, degree) + Fraction(phase_bit, 2)
                        )
                        for index in range(degree)
                    )
                    right_distance = circle_distance(
                        degree * y + Fraction(next_phase, 2)
                    )
                    # The identity is almost-everywhere; omit any sampled null
                    # boundary instead of imposing an endpoint convention.
                    if radius in left_distances or right_distance == residual_radius:
                        continue
                    left = sum(distance < radius for distance in left_distances)
                    right = quotient + int(right_distance < residual_radius)
                    assert left == right
                    checked += 1
    return checked


def covariance(u: int, v: int) -> Fraction:
    common = gcd(u, v)
    U, V = u // common, v // common
    return Fraction(d14(U + V) - d14(U - V), 7 * U * V)


def bonf3(depth: int) -> int:
    return 1 - depth + depth * (depth - 1) // 2 - depth * (depth - 1) * (depth - 2) // 6


def optimal_current_rebate(current_size: int) -> Fraction:
    """Sharp lower envelope of |bonf3(a)-bonf3(b)|/2 at |a-b| fixed."""
    assert current_size >= 0
    if current_size <= 2:
        return Fraction()
    if current_size == 3:
        return Fraction(1, 2)
    return Fraction(
        current_size * (current_size * current_size - 6 * current_size + 11), 12
    )


def audit_optimal_current_envelope() -> None:
    for current_size in range(13):
        values = [
            Fraction(abs(bonf3(a) - bonf3(b)), 2)
            for a in range(13)
            for b in range(13)
            if abs(a - b) == current_size
        ]
        assert min(values) == optimal_current_rebate(current_size)
    assert tuple(optimal_current_rebate(c) for c in range(13)) == (
        Fraction(0), Fraction(0), Fraction(0), Fraction(1, 2),
        Fraction(1), Fraction(5, 2), Fraction(11, 2), Fraction(21, 2),
        Fraction(18), Fraction(57, 2), Fraction(85, 2), Fraction(121, 2),
        Fraction(83),
    )


def tail_events(speeds: tuple[int, ...]) -> dict[Fraction, list[int]]:
    events: dict[Fraction, list[int]] = defaultdict(lambda: [0, 0])
    for speed in speeds:
        for a, b in clipped_teeth(speed, Fraction()):
            events[a][0] += 1
            events[b][0] -= 1
        for a, b in clipped_teeth(speed, Fraction(1, 2)):
            events[a][1] += 1
            events[b][1] -= 1
    events[Fraction()]
    events[Fraction(1, 2)]
    return events


def observable_and_wall_integral(
    speeds: tuple[int, ...],
    wall_speed: int,
    wall_shift: Fraction,
    wall_radius: Fraction,
    observable,
) -> tuple[Fraction, Fraction]:
    """Integrate a sheet-symmetric tail observable and its wall restriction."""
    events: dict[Fraction, list[int]] = defaultdict(lambda: [0, 0, 0])
    for wall, delta in tail_events(speeds).items():
        events[wall][0] += delta[0]
        events[wall][1] += delta[1]
    if wall_radius:
        for a, b in clipped_arcs(wall_speed, wall_shift, wall_radius):
            events[a][2] += 1
            events[b][2] -= 1
    events[Fraction()]
    events[Fraction(1, 2)]

    lower = upper = wall_depth = 0
    total = restricted = Fraction()
    walls = sorted(events)
    for wall, next_wall in zip(walls, walls[1:]):
        da, db, dw = events[wall]
        lower += da
        upper += db
        wall_depth += dw
        width = next_wall - wall
        value = observable(lower, upper)
        total += width * value
        restricted += width * value * wall_depth
        assert wall_depth in (0, 1)
    return total, restricted


def exact_profile(h: int, speeds: tuple[int, ...]) -> dict[str, Fraction | int]:
    assert len(speeds) == len(set(speeds)) == 12
    assert all(speed > 0 and speed % 2 for speed in speeds)

    # Event payload is (delta lower depth, delta upper depth, delta anchor-danger).
    events: dict[Fraction, list[int]] = defaultdict(lambda: [0, 0, 0])
    for wall, delta in tail_events(speeds).items():
        events[wall][0] += delta[0]
        events[wall][1] += delta[1]
    for a, b in clipped_teeth(2 * h, Fraction()):
        events[a][2] += 1
        events[b][2] -= 1

    zero, half = Fraction(), Fraction(1, 2)
    events[zero]
    events[half]
    walls = sorted(events)
    lower = upper = anchor = 0
    core_mass = strip_mass = Fraction()
    core_energy = strip_energy = Fraction()
    core_b3 = Fraction()
    core_abs_current = Fraction()
    core_positive_quadratic_current = Fraction()
    core_nonlinear_rebate = Fraction()
    core_pure_cubic_current = Fraction()
    core_sum_current_correction = Fraction()
    core_optimal_current_only_rebate = Fraction()
    core_q_cubic = Fraction()
    core_max_abs_current = strip_max_abs_current = 0
    core_state_mass: dict[tuple[int, int], Fraction] = defaultdict(Fraction)
    strip_state_mass: dict[tuple[int, int], Fraction] = defaultdict(Fraction)

    for wall, next_wall in zip(walls, walls[1:]):
        delta_lower, delta_upper, delta_anchor = events[wall]
        lower += delta_lower
        upper += delta_upper
        anchor += delta_anchor
        width = next_wall - wall
        if not width:
            continue
        assert 0 <= lower <= 12 and 0 <= upper <= 12 and lower + upper <= 12
        assert anchor in (0, 1)
        current = lower - upper
        if anchor:
            strip_mass += width
            strip_energy += width * current * current
            strip_max_abs_current = max(strip_max_abs_current, abs(current))
            strip_state_mass[(lower, upper)] += width
        else:
            core_mass += width
            core_energy += width * current * current
            core_b3 += width * bonf3(lower)
            core_abs_current += width * abs(current)
            core_positive_quadratic_current += width * max(current * current - 4, 0)
            core_nonlinear_rebate += width * abs(bonf3(lower) - bonf3(upper)) / 2
            core_pure_cubic_current += width * abs(current) * (current * current - 4) / 48
            core_sum_current_correction += (
                width * abs(current) * (lower + upper - 4) ** 2 / 16
            )
            current_size = abs(current)
            core_optimal_current_only_rebate += width * optimal_current_rebate(
                current_size
            )
            core_q_cubic += width * bonf3(min(lower, upper))
            core_max_abs_current = max(core_max_abs_current, abs(current))
            core_state_mass[(lower, upper)] += width

    full_half_formula = sum(
        (covariance(u, v) for u in speeds for v in speeds), start=Fraction()
    ) / 2
    assert core_mass == Fraction(3, 7)
    assert strip_mass == Fraction(1, 14)
    assert core_energy + strip_energy == full_half_formula
    assert core_q_cubic == core_b3 + core_nonlinear_rebate
    assert core_nonlinear_rebate == core_pure_cubic_current + core_sum_current_correction

    return {
        "full_half_energy": full_half_formula,
        "core_energy": core_energy,
        "strip_energy": strip_energy,
        "core_b3": core_b3,
        "core_nonlinear_rebate": core_nonlinear_rebate,
        "pure_cubic_current_rebate": core_pure_cubic_current,
        "sum_current_correction": core_sum_current_correction,
        "optimal_current_only_rebate": core_optimal_current_only_rebate,
        "core_q_cubic": core_q_cubic,
        "quadratic_current_rebate": core_positive_quadratic_current / 12,
        "positive_quadratic_certificate": core_b3 + core_positive_quadratic_current / 12,
        "optimal_current_only_certificate": core_b3 + core_optimal_current_only_rebate,
        "variance_certificate": core_b3 + (core_energy - Fraction(12, 7)) / 12,
        "core_abs_current": core_abs_current,
        "core_max_abs_current": core_max_abs_current,
        "strip_max_abs_current": strip_max_abs_current,
        "core_states": len(core_state_mass),
        "strip_states": len(strip_state_mass),
    }


def half_energy_on_residual_wall(
    h: int, base: tuple[int, ...], multiplier: int
) -> tuple[Fraction, Fraction, int, int, int]:
    """Return the exact Euclidean-remainder expression for strip energy.

    Put ``c=gcd(multiplier,2h)``, ``D=multiplier/c=7q+r``.  Averaging the
    anchor indicator over the D inverse images of multiplication by D gives

      (q + 1_{|| (2h/c)x + q/2 || < r/14}) / D.

    This routine computes the energy of ``C_base^2`` on that one residual
    wall and returns both the resulting strip value and the residual energy.
    """
    common = gcd(multiplier, 2 * h)
    reduced_multiplier = multiplier // common
    reduced_anchor = 2 * h // common
    quotient, remainder = divmod(reduced_multiplier, 7)

    events: dict[Fraction, list[int]] = defaultdict(lambda: [0, 0, 0])
    for speed in base:
        for a, b in clipped_teeth(speed, Fraction()):
            events[a][0] += 1
            events[b][0] -= 1
        for a, b in clipped_teeth(speed, Fraction(1, 2)):
            events[a][1] += 1
            events[b][1] -= 1
    if remainder:
        for a, b in clipped_arcs(
            reduced_anchor, Fraction(quotient, 2), Fraction(remainder, 14)
        ):
            events[a][2] += 1
            events[b][2] -= 1
    events[Fraction()]
    events[Fraction(1, 2)]

    lower = upper = residual = 0
    full_energy = residual_energy = Fraction()
    walls = sorted(events)
    for wall, next_wall in zip(walls, walls[1:]):
        da, db, dr = events[wall]
        lower += da
        upper += db
        residual += dr
        width = next_wall - wall
        current_squared = (lower - upper) ** 2
        full_energy += width * current_squared
        residual_energy += width * current_squared * residual
        assert residual in (0, 1)

    strip_formula = (
        Fraction(quotient, reduced_multiplier) * full_energy
        + Fraction(1, reduced_multiplier) * residual_energy
    )
    return strip_formula, residual_energy, reduced_multiplier, quotient, remainder


def audit_general_observable_identity(
    h: int, base: tuple[int, ...], multiplier: int
) -> tuple[int, int, int]:
    """Audit the identity on positive and signed symmetric observables."""
    common = gcd(multiplier, 2 * h)
    assert common == gcd(multiplier, h)  # multiplier is odd
    reduced_multiplier = multiplier // common
    reduced_anchor = 2 * h // common
    quotient, remainder = divmod(reduced_multiplier, 7)
    assert reduced_anchor % 2 == 0

    observables = {
        "constant": lambda _a, _b: 1,
        "current_energy": lambda a, b: (a - b) ** 2,
        "cross_current_sum": lambda a, b: Fraction((a - b) ** 2 - (a + b), 2),
        "q_cubic": lambda a, b: bonf3(min(a, b)),
        "free_sheet_count": lambda a, b: int(a == 0) + int(b == 0),
        "symmetric_bonf3": lambda a, b: Fraction(bonf3(a) + bonf3(b), 2),
    }
    scaled = tuple(multiplier * speed for speed in base)
    for name, observable in observables.items():
        base_total, residual = observable_and_wall_integral(
            base,
            reduced_anchor,
            Fraction(quotient, 2),
            Fraction(remainder, 14),
            observable,
        )
        scaled_total, direct_strip = observable_and_wall_integral(
            scaled, 2 * h, Fraction(), Fraction(1, 14), observable
        )
        predicted_strip = (
            Fraction(quotient, reduced_multiplier) * base_total
            + Fraction(1, reduced_multiplier) * residual
        )
        assert scaled_total == base_total, (name, multiplier)
        assert direct_strip == predicted_strip, (name, multiplier)

        discrepancy = direct_strip - base_total / 7
        absolute_bound = (
            Fraction(
                0 if remainder == 0 else max(remainder, 7 - remainder),
                7 * reduced_multiplier,
            )
            * abs(base_total)
        )
        # This L1 bound can use |integral f| only for the four nonnegative
        # observables.  The signed symmetric Bonferroni observable is audited
        # for the exact identity but deliberately excluded from this check.
        if name != "symmetric_bonf3" and all(
            observable(a, b) >= 0 for a in range(13) for b in range(13 - a)
        ):
            assert abs(discrepancy) <= absolute_bound

    # The residual wall itself has ordinary half-base mass r/14.
    _, residual_mass = observable_and_wall_integral(
        base,
        reduced_anchor,
        Fraction(quotient, 2),
        Fraction(remainder, 14),
        lambda _a, _b: 1,
    )
    assert residual_mass == Fraction(remainder, 14)
    return reduced_multiplier, quotient, remainder


def print_profile(label: str, h: int, speeds: tuple[int, ...]) -> dict[str, Fraction | int]:
    profile = exact_profile(h, speeds)
    print(label)
    print(f" h={h} primitive_gcd={gcd(2*h, *speeds)} speeds=" + ",".join(map(str, speeds)))
    for key, value in profile.items():
        print(f" {key}={value}")
    return profile


def main() -> None:
    transducer_checks = audit_nested_wall_transducer()
    print("NESTED_WALL_TRANSDUCER")
    print(f" exact_nonboundary_checks={transducer_checks}")
    print(" state_update=sD=7Q+R;phase_next=D*phase+Q_mod_2;residual_radius=R/14")
    audit_optimal_current_envelope()
    print("OPTIMAL_CURRENT_ENVELOPE")
    print(" g(0..12)=0,0,0,1/2,1,5/2,11/2,21/2,18,57/2,85/2,121/2,83")
    print(" pointwise_minimum_over_all_depth_pairs=PASS")
    base = tuple(range(1, 22, 2)) + (45,)
    profiles = []
    by_multiplier: dict[int, dict[str, Fraction | int]] = {}
    for multiplier in (1, 3, 9, 11, 13, 17, 19, 23, 49, 127):
        speeds = tuple(multiplier * speed for speed in base)
        profile = print_profile(f"AP_STEP_MULTIPLIER_{multiplier}", 420, speeds)
        formula, residual_energy, reduced, quotient, remainder = (
            half_energy_on_residual_wall(420, base, multiplier)
        )
        assert formula == profile["strip_energy"]
        audited_reduced, audited_quotient, audited_remainder = (
            audit_general_observable_identity(420, base, multiplier)
        )
        assert (audited_reduced, audited_quotient, audited_remainder) == (
            reduced,
            quotient,
            remainder,
        )
        # With exactly two tails, cross_current_sum is sigma_u*sigma_v, so
        # this separately audits the signed restricted-covariance formula.
        assert audit_general_observable_identity(420, (1, 45), multiplier) == (
            reduced,
            quotient,
            remainder,
        )
        print(
            " euclidean_strip_formula="
            f"D:{reduced}=7*{quotient}+{remainder};"
            f"residual_energy:{residual_energy};strip:{formula}"
        )
        profiles.append(profile)
        by_multiplier[multiplier] = profile

    full_values = {profile["full_half_energy"] for profile in profiles}
    assert len(full_values) == 1
    assert len({profile["strip_energy"] for profile in profiles}) > 1
    full_energy = Fraction(35666, 15015)
    base_strip = Fraction(15361, 45220)
    assert full_values == {full_energy}
    assert by_multiplier[1]["strip_energy"] == base_strip
    assert by_multiplier[1]["core_energy"] == Fraction(39490603, 19399380)
    assert by_multiplier[1]["variance_certificate"] == Fraction(
        -71225183, 162954792
    )
    assert by_multiplier[1]["positive_quadratic_certificate"] == Fraction(
        -3392354543, 9777287520
    )
    assert by_multiplier[1]["optimal_current_only_certificate"] == Fraction(
        2572403, 33948915
    )
    assert by_multiplier[1]["core_q_cubic"] == Fraction(526429, 5819814)
    assert by_multiplier[13]["strip_energy"] == (2 * full_energy - base_strip) / 13
    assert by_multiplier[49]["strip_energy"] == full_energy / 7
    assert by_multiplier[127]["strip_energy"] == (
        18 * full_energy + base_strip
    ) / 127
    assert all(
        profile["optimal_current_only_certificate"] <= profile["core_q_cubic"]
        for profile in profiles
    )
    print("COMMON_DILATION_AUDIT")
    print(" observables=constant,current_energy,cross_current_sum,q_cubic,free_sheet_count,symmetric_bonf3")
    print(" signed_pair_covariance_control=base_pair_1,45")
    print(" multipliers=1,3,9,11,13,17,19,23,49,127")
    print(" all_remainders_mod_7=0,1,2,3,4,5,6")
    print(" nontrivial_gcd_controls=3,9,49")
    print(" full_half_current_energy=INVARIANT")
    print(" anchor_strip_current_energy=NOT_INVARIANT")
    print(" D13_complement_relation=strip13=(2*full-strip1)/13")
    print(" D49_reduced_D7_exact_independence=strip49=full/7")
    print(" D127_nested_relation=strip127=(18*full+strip1)/127")
    print(" exact_variance_certificate_m1=-71225183/162954792<0")
    print(" exact_positive_quadratic_certificate_m1=-3392354543/9777287520<0")
    print(" exact_optimal_current_only_certificate_m1=2572403/33948915>0")
    print(" exact_q_cubic_m1=526429/5819814>0")
    print(" verdict=full_pair_kernel_does_not_determine_anchor_strip")
    print("PASS")


if __name__ == "__main__":
    main()
