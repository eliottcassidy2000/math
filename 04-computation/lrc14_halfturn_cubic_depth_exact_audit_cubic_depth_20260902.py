#!/usr/bin/env python3
"""Exact audit of half-turn cubic-depth duals at the h=420 residual.

This companion does four things, all over strict open danger teeth:

* proves the seven-state q=min(a,b) cubic is the maximum of two ordinary
  third Bonferroni polynomials and isolates its exact nonlinear-current term;
* proves a proposed symmetric bivariate cubic is a pointwise lower dual, then
  gives primitive denominator-complete h=420 rows where its integral is
  negative (including one failing the fixed half-turn clock gate);
* evaluates the original q-cubic exactly on the same rows; and
* audits the reduced-residue Brownian factorization of the signed pair-current
  covariance and checks it against a direct exact wall sweep.

The calculations are finite exact certificates.  They do not prove universal
positivity of the q-cubic and do not prove LRC(14).
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import ceil, comb, floor, gcd


@dataclass(frozen=True)
class Interval:
    left: Fraction
    right: Fraction


def choose(n: int, k: int) -> int:
    return comb(n, k) if n >= k else 0


def bonf3(n: int) -> int:
    return 1 - n + choose(n, 2) - choose(n, 3)


def q_cubic(a: int, b: int) -> int:
    return bonf3(min(a, b))


def bivariate_cubic(a: int, b: int) -> Fraction:
    return (
        Fraction(4, 5)
        - Fraction(a + b, 5)
        + Fraction(7, 15) * (choose(a, 2) + choose(b, 2))
        - Fraction(2, 5) * a * b
        - Fraction(3, 5) * (choose(a, 3) + choose(b, 3))
        + Fraction(2, 15) * (choose(a, 2) * b + a * choose(b, 2))
    )


def anchor_component(h: int, k: int) -> Interval:
    return Interval(
        Fraction(14 * k + 1, 28 * h),
        Fraction(14 * k + 13, 28 * h),
    )


def teeth_meeting(w: int, target: Interval, shift: Fraction) -> tuple[Interval, ...]:
    """Strict intervals where ||w(t+shift)||<1/14, clipped to target."""
    lo = floor(w * (target.left + shift) - Fraction(1, 14))
    hi = ceil(w * (target.right + shift) + Fraction(1, 14))
    answer = []
    for n in range(lo - 1, hi + 2):
        left = max(target.left, Fraction(14 * n - 1, 14 * w) - shift)
        right = min(target.right, Fraction(14 * n + 1, 14 * w) - shift)
        if left < right:
            answer.append(Interval(left, right))
    return tuple(answer)


def state_masses_on_intervals(
    speeds: tuple[int, ...], targets: tuple[Interval, ...]
) -> dict[tuple[int, int], Fraction]:
    masses: dict[tuple[int, int], Fraction] = {}
    for target in targets:
        low = []
        high = []
        for w in speeds:
            low.extend(teeth_meeting(w, target, Fraction()))
            high.extend(teeth_meeting(w, target, Fraction(1, 2)))

        walls = {target.left, target.right}
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
            assert a + b <= len(speeds)
            key = (a, b)
            masses[key] = masses.get(key, Fraction()) + right - left
    return masses


def integrate(
    masses: dict[tuple[int, int], Fraction], function
) -> Fraction:
    return sum(
        (mass * function(a, b) for (a, b), mass in masses.items()),
        start=Fraction(),
    )


def residue_distance_14(value: int) -> int:
    residue = value % 14
    return min(residue, 14 - residue)


def residue_level_sign(residue: int) -> tuple[int, int]:
    residue %= 14
    table = {
        1: (1, 1),
        5: (2, 1),
        3: (3, 1),
        13: (1, -1),
        9: (2, -1),
        11: (3, -1),
        7: (0, 0),
    }
    return table[residue]


def covariance_kernel(u: int, v: int) -> Fraction:
    common = gcd(u, v)
    U, V = u // common, v // common
    numerator = residue_distance_14(U + V) - residue_distance_14(U - V)
    return Fraction(numerator, 7 * U * V)


def audit_pointwise_identities() -> None:
    values = tuple(bonf3(q) for q in range(13))
    assert values[:7] == (1, 0, 0, 0, -1, -4, -10)
    assert all(values[q + 1] <= values[q] for q in range(12))

    bivariate_equalities = []
    current_equalities = []
    for a in range(13):
        for b in range(13 - a):
            free = int(a == 0) + int(b == 0)
            assert bivariate_cubic(a, b) <= free
            assert q_cubic(a, b) == max(bonf3(a), bonf3(b))

            n = a + b
            d = a - b
            bivar_factor = -Fraction(
                (n - 8) * (n - 6) * (n - 2) + (11 * n - 48) * d * d,
                120,
            )
            assert bivariate_cubic(a, b) == bivar_factor

            difference = abs(bonf3(a) - bonf3(b))
            bracket = 3 * (n - 4) ** 2 + d * d - 4
            assert Fraction(abs(d) * bracket, 24) == difference
            assert difference >= Fraction(max(d * d - 4, 0), 6)

            if bivariate_cubic(a, b) == free:
                bivariate_equalities.append((a, b))
            if difference == Fraction(max(d * d - 4, 0), 6):
                current_equalities.append((a, b))

    # The Brownian factorization is an exact seven-by-seven identity.
    brownian_order = (1, 5, 3, 13, 9, 11, 7)
    brownian_rows = []
    for r in brownian_order:
        row = []
        ell_r, eps_r = residue_level_sign(r)
        for s in brownian_order:
            ell_s, eps_s = residue_level_sign(s)
            numerator = residue_distance_14(r + s) - residue_distance_14(r - s)
            brownian = 2 * eps_r * eps_s * min(ell_r, ell_s)
            assert numerator == brownian
            row.append(numerator)
        brownian_rows.append(tuple(row))

    print("POINTWISE_AUDIT")
    print(" bonf3_0_to_12=" + ",".join(map(str, values)))
    print(" q_cubic=max(bonf3(a),bonf3(b))")
    print(" plateau=bonf3(1)=bonf3(2)=bonf3(3)=0")
    print(
        " nonlinear_current=abs(a-b)*(3*(a+b-4)^2+(a-b)^2-4)/24"
    )
    print(" current_lower=abs_bonf3_difference>=(max((a-b)^2-4,0))/6")
    print(" bivariate_factor=-((n-8)(n-6)(n-2)+(11n-48)d^2)/120")
    print(" bivariate_equality_states=" + repr(tuple(bivariate_equalities)))
    print(" current_lower_equality_states=" + repr(tuple(current_equalities)))
    print(" brownian_residue_order=" + ",".join(map(str, brownian_order)))
    print(" brownian_rows=" + repr(tuple(brownian_rows)))


def profile(h: int, speeds: tuple[int, ...], label: str) -> dict[str, Fraction | int]:
    assert len(speeds) == len(set(speeds)) == 12
    assert all(w > 0 and w % 2 for w in speeds)

    core_targets = tuple(anchor_component(h, k) for k in range(h))
    core = state_masses_on_intervals(speeds, core_targets)
    assert sum(core.values(), start=Fraction()) == Fraction(3, 7)
    assert all(core.get((a, b), Fraction()) == core.get((b, a), Fraction()) for a, b in core)

    q_masses = tuple(
        sum(
            (mass for (a, b), mass in core.items() if min(a, b) == q),
            start=Fraction(),
        )
        for q in range(7)
    )
    free_mass = integrate(core, lambda a, b: int(a == 0) + int(b == 0))
    q_value = integrate(core, q_cubic)
    bivariate_value = integrate(core, bivariate_cubic)
    one_sheet_bonf3 = integrate(core, lambda a, _b: bonf3(a))
    symmetric_check = integrate(core, lambda _a, b: bonf3(b))
    nonlinear_current = integrate(
        core, lambda a, b: abs(bonf3(a) - bonf3(b))
    )
    current_variance = integrate(core, lambda a, b: (a - b) ** 2)
    current_positive_part = integrate(
        core, lambda a, b: max((a - b) ** 2 - 4, 0)
    )
    assert one_sheet_bonf3 == symmetric_check
    assert q_value == one_sheet_bonf3 + nonlinear_current / 2
    assert q_value >= one_sheet_bonf3 + current_positive_part / 12
    variance_certificate = one_sheet_bonf3 + (current_variance - Fraction(12, 7)) / 12
    assert q_value >= variance_certificate

    # Pair-current formula on the full circle; C^2 is half-turn invariant,
    # so its first-half energy is one half of the exact pair sum.
    global_variance_formula = sum(
        (covariance_kernel(u, v) for u in speeds for v in speeds),
        start=Fraction(),
    )
    assert global_variance_formula >= 0
    full_half_variance = global_variance_formula / 2
    removed_strip_variance = full_half_variance - current_variance
    assert removed_strip_variance >= 0

    missing = tuple(
        modulus
        for modulus in range(2, 15)
        if not any(w % modulus == 0 for w in (2 * h,) + speeds)
    )
    blockers = tuple(
        w for w in speeds if 12 * h < (w % (28 * h)) < 16 * h
    )
    primitive = gcd(2 * h, *speeds)

    contributions = sorted(
        (
            (mass * bivariate_cubic(a, b), a, b, mass)
            for (a, b), mass in core.items()
        )
    )
    worst_contributions = tuple(contributions[:5])

    print(label)
    print(f" h={h} primitive_gcd={primitive}")
    print(" speeds=" + ",".join(map(str, speeds)))
    print(" missing_denominators=" + (",".join(map(str, missing)) or "none"))
    print(" fixed_halfturn_blockers=" + (",".join(map(str, blockers)) or "none"))
    print(" q_masses=" + ",".join(f"{q}:{mass}" for q, mass in enumerate(q_masses)))
    print(f" q_max={max(q for q, mass in enumerate(q_masses) if mass)}")
    print(f" free_mass={free_mass}")
    print(f" q_cubic={q_value}")
    print(f" bivariate_cubic={bivariate_value}")
    print(f" one_sheet_bonf3={one_sheet_bonf3}")
    print(f" nonlinear_current_rebate={nonlinear_current / 2}")
    print(f" current_variance_core={current_variance}")
    print(f" current_positive_part_core={current_positive_part}")
    print(f" current_variance_certificate={variance_certificate}")
    print(f" current_variance_full_circle={global_variance_formula}")
    print(f" current_variance_removed_anchor_strip={removed_strip_variance}")
    print(" worst_bivariate_state_contributions=" + repr(worst_contributions))

    return {
        "q_cubic": q_value,
        "bivariate": bivariate_value,
        "q_max": max(q for q, mass in enumerate(q_masses) if mass),
        "primitive": primitive,
        "missing": len(missing),
        "blockers": len(blockers),
    }


def main() -> None:
    audit_pointwise_identities()

    base = tuple(range(1, 22, 2)) + (45,)
    base_full = state_masses_on_intervals(
        base, (Interval(Fraction(), Fraction(1, 2)),)
    )
    assert max(min(a, b) for a, b in base_full) == 3
    print("AP_STEP_GLOBAL")
    print(" base=" + ",".join(map(str, base)))
    print(" global_q_max=3")
    print(" odd_multiplier_preserves_paired_state_via_t_mapsto_multiplier*t")
    family_rows = []
    for multiplier in (1, 19, 127):
        row = tuple(multiplier * w for w in base)
        result = profile(420, row, f"AP_STEP_FAMILY_multiplier_{multiplier}")
        assert result["primitive"] == 1
        assert result["missing"] == 0
        assert result["q_max"] == 3
        assert result["q_cubic"] > 0
        assert result["bivariate"] < 0
        family_rows.append(result)

    # At multiplier 127, tail 45*127=5715 lies in the failed fixed-clock
    # residue band (5040,6720); this keeps the clean common-multiplier family
    # inside that sharper necessary gate as well.
    assert family_rows[-1]["blockers"] > 0

    near_cubic = profile(
        420,
        (1, 837, 841, 1681, 2521, 3361, 4201, 5041, 5881, 6721, 7561, 8401),
        "DENOMINATOR_COMPLETE_MODULAR_CUBIC_HOSTILE",
    )
    assert near_cubic["primitive"] == 1
    assert near_cubic["missing"] == 0
    assert near_cubic["blockers"] > 0
    assert near_cubic["q_cubic"] > 0

    print("CONCLUSIONS")
    print(" bivariate_integral_uniform_positivity=REFUTED")
    print(" q_cubic_integral_uniform_positivity=OPEN")
    print(" q_cubic_AP_step_mechanism=q_max_at_most_3_so_q_cubic=H0")
    print(" brownian_pair_kernel_controls_full_current_variance_not_core_variance")
    print(" anchor_strip_variance_is_required_sidecar")
    print("PASS")


if __name__ == "__main__":
    main()
