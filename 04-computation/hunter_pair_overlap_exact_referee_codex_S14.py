#!/usr/bin/env python3
"""Exact referee for THM-856's global pair-overlap claim.

For ``D_q={t in R/Z: ||qt||<=1/13}``, write ``x=ga``, ``y=gb`` with
``g=gcd(x,y)`` and ``gcd(a,b)=1``.  Haar invariance under ``t -> gt`` first
reduces the global intersection to ``D_a cap D_b``.

The intersection of two line intervals of radii ``A,B`` whose centres are
distance ``d`` apart is

    (A+B-d)_+ - (|A-B|-d)_+.

The second, containment-plateau term was omitted in the original THM-856
formula.  Bezout identifies the ``ab`` pairs of teeth with all residues
``m mod ab``.  Since ``2(a+b)<13ab``, only their least integer
representatives contribute, and hence

    mu(D_x cap D_y)
      = [T((a+b)/13)-T(|a-b|/13)]/(ab),
    T(z)=sum_{m in Z}(z-|m|)_+.

Writing ``z=floor(z)+t`` gives ``T(z)=z^2+psi(t)`` with
``psi(t)=t(1-t)``.  Therefore

    mu(D_x cap D_y)
      = 4/169
        + [psi({(a+b)/13})-psi({|a-b|/13})]/(ab).

The correction has either sign.  In particular ``(a,b)=(6,7)`` and
``(5,9)`` disprove the claimed lower bound ``4/169``.

The second audit records the projective-ratio obstruction to THM-856's
proposed raw-speed deficit lemma.  On the legal scale-one ray

    g=1+13k, x=6g=6+13(6k), y=7g=7+13(7k),

the reduced pair remains ``(6,7)``, so its global density is always ``2/91``.
For a fixed finite interval union E, periodic averaging instead gives

    mu(E cap D_{ga} cap D_{gb})
      = mu(E) mu(D_a cap D_b) + O_E(1/g).

Thus the correct baseline is reduced-ratio dependent.  For ``E=E_[5]`` the
deficit from ``(4/169)mu(E)`` tends the nonzero value ``14/17745`` even while
the raw speeds tend to infinity.

All comparisons below use ``fractions.Fraction``.  The direct oracle forms
literal closed danger teeth and intersects interval unions; it is independent
of the triangular-sum formula being refereed.
"""

from __future__ import annotations

from fractions import Fraction as F
from math import gcd


BASELINE = F(4, 169)


def measure(intervals: tuple[tuple[F, F], ...]) -> F:
    return sum((hi - lo for lo, hi in intervals), F(0))


def meet(
    left: tuple[tuple[F, F], ...], right: tuple[tuple[F, F], ...]
) -> tuple[tuple[F, F], ...]:
    out: list[tuple[F, F]] = []
    i = j = 0
    while i < len(left) and j < len(right):
        lo = max(left[i][0], right[j][0])
        hi = min(left[i][1], right[j][1])
        if lo < hi:
            out.append((lo, hi))
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return tuple(out)


def danger(speed: int) -> tuple[tuple[F, F], ...]:
    """Literal closed danger teeth, split at the circle cut zero."""
    radius = F(1, 13 * speed)
    out = [(F(0), radius)]
    for centre_numerator in range(1, speed):
        centre = F(centre_numerator, speed)
        out.append((centre - radius, centre + radius))
    out.append((F(1) - radius, F(1)))
    return tuple(out)


def safe_set(speeds: tuple[int, ...]) -> tuple[tuple[F, F], ...]:
    current: tuple[tuple[F, F], ...] = ((F(0), F(1)),)
    for speed in speeds:
        safe = tuple(
            (F(13 * k + 1, 13 * speed), F(13 * (k + 1) - 1, 13 * speed))
            for k in range(speed)
        )
        current = meet(current, safe)
    return current


def direct_overlap(x: int, y: int) -> F:
    return measure(meet(danger(x), danger(y)))


def psi(t: F) -> F:
    assert 0 <= t < 1
    return t * (1 - t)


def triangular_sum(z: F) -> F:
    """T(z)=sum_m (z-|m|)_+, evaluated without truncation ambiguity."""
    assert z >= 0
    integer = z.numerator // z.denominator
    fractional = z - integer
    return F(integer * integer) + (2 * integer + 1) * fractional


def overlap_formula(x: int, y: int) -> F:
    g = gcd(x, y)
    a, b = x // g, y // g
    return (
        triangular_sum(F(a + b, 13))
        - triangular_sum(F(abs(a - b), 13))
    ) / (a * b)


def sawtooth_formula(x: int, y: int) -> F:
    g = gcd(x, y)
    a, b = x // g, y // g
    plus = F(a + b, 13) % 1
    minus = F(abs(a - b), 13) % 1
    return BASELINE + (psi(plus) - psi(minus)) / (a * b)


def global_formula_audit() -> None:
    # T(z)=z^2+psi({z}) is exact, including integral breakpoints.
    triangular_checks = 0
    for numerator in range(401):
        z = F(numerator, 13)
        assert triangular_sum(z) == z * z + psi(z % 1)
        triangular_checks += 1

    # Direct literal intervals check both the formula and gcd reduction.
    direct_checks = 0
    for x in range(1, 61):
        for y in range(x + 1, 61):
            exact = direct_overlap(x, y)
            assert exact == overlap_formula(x, y) == sawtooth_formula(x, y)
            direct_checks += 1
    assert direct_checks == 1_770

    # Freeze the two smallest transparent counterexamples supplied by review.
    assert overlap_formula(6, 7) == F(2, 91)
    assert overlap_formula(5, 9) == F(4, 195)
    assert F(2, 91) - BASELINE == -F(2, 1183)
    assert F(4, 195) - BASELINE == -F(8, 2535)

    below = equal = above = 0
    correction_bound_checks = 0
    equality_character_checks = 0
    for a in range(1, 81):
        for b in range(a + 1, 81):
            if gcd(a, b) != 1:
                continue
            value = overlap_formula(a, b)
            correction = value - BASELINE
            # Integer residues make the sharp sawtooth amplitude 42/169
            # (attained at residues 6 and 7), slightly sharper than 1/4.
            assert abs(correction) <= F(42, 169 * a * b)
            correction_bound_checks += 1
            # psi(r/13)=psi(s/13) iff r=+/-s mod 13.  Since 2 is a unit
            # mod 13, equality occurs exactly when 13 divides a or b.
            assert (value == BASELINE) == (a % 13 == 0 or b % 13 == 0)
            equality_character_checks += 1
            if value < BASELINE:
                below += 1
            elif value == BASELINE:
                equal += 1
            else:
                above += 1

    # Frozen census over the declared coprime domain.
    assert correction_bound_checks == equality_character_checks
    assert (below, equal, above) == (837, 282, 846)

    print("GLOBAL EXACT FORMULA")
    print(f"triangular_identity_checks={triangular_checks}")
    print(f"direct_closed_tooth_checks={direct_checks} domain=1<=x<y<=60")
    print(
        "formula=mu(D_x cap D_y)="
        "[T((a+b)/13)-T(|a-b|/13)]/(ab), gcd(a,b)=1"
    )
    print(
        "sawtooth=4/169+"
        "[psi({(a+b)/13})-psi({|a-b|/13})]/(ab)"
    )
    print("T(z)=z^2+psi({z}), psi(t)=t(1-t)")
    print(
        f"coprime_sign_census_a_lt_b_le_80="
        f"below:{below},equal:{equal},above:{above}"
    )
    print("equality_character=13_divides_a_or_b")
    print("uniform_error_bound=abs(mu-4/169)<=42/(169ab)")
    print(
        "counterexample_6_7="
        f"{overlap_formula(6,7)} baseline_gap={overlap_formula(6,7)-BASELINE}"
    )
    print(
        "counterexample_5_9="
        f"{overlap_formula(5,9)} baseline_gap={overlap_formula(5,9)-BASELINE}"
    )


def projective_ray_audit() -> None:
    core = safe_set((1, 2, 3, 4, 5))
    core_measure = measure(core)
    assert len(core) == 10 and core_measure == F(7, 15)

    reduced_density = overlap_formula(6, 7)
    reduced_target = core_measure * reduced_density
    ideal_target = core_measure * BASELINE
    persistent_deficit = ideal_target - reduced_target
    assert reduced_density == F(2, 91)
    assert persistent_deficit == F(14, 17745)

    rows: list[tuple[int, int, int, F, F, F]] = []
    for k in (1, 2, 5, 10, 30, 31, 32, 35, 40, 60):
        g = 1 + 13 * k
        x, y = 6 * g, 7 * g
        # These really are legal scale-one lifts of labels 6 and 7.
        assert x == 6 + 13 * (6 * k)
        assert y == 7 + 13 * (7 * k)
        assert gcd(x, y) == g
        assert direct_overlap(x, y) == reduced_density
        restricted = measure(meet(meet(core, danger(x)), danger(y)))
        error = restricted - reduced_target
        assert abs(error) <= F(len(core), g)
        rows.append((k, g, x, restricted, error, g * error))

    # g*(restricted-reduced_target) is endpoint-periodic.  Since the core
    # endpoint denominator lcm is 390 and g increases by 13 per k, k-period
    # 30 is exact.  Freeze five independently reconstructed pairs.
    periodic_checks = 0
    for k in range(1, 6):
        values = []
        for shift in (0, 30):
            kk = k + shift
            g = 1 + 13 * kk
            restricted = measure(
                meet(meet(core, danger(6 * g)), danger(7 * g))
            )
            values.append(g * (restricted - reduced_target))
        assert values[0] == values[1]
        periodic_checks += 1
    assert periodic_checks == 5

    print("\nPROJECTIVE-RATIO OBSTRUCTION")
    print(f"core=E_[5] components={len(core)} measure={core_measure}")
    print(
        "legal_ray=g=1+13k,x=6g=6+13(6k),y=7g=7+13(7k)"
    )
    print(
        f"fixed_reduced_density={reduced_density} "
        f"reduced_target_muE={reduced_target}"
    )
    print(
        f"ideal_4_over_169_target={ideal_target} "
        f"persistent_deficit={persistent_deficit}"
    )
    for k, g, x, restricted, error, scaled in rows:
        print(
            f"k={k:2d} g={g:3d} x={x:4d} "
            f"restricted={restricted} error_to_reduced={error} "
            f"g_times_error={scaled}"
        )
    print(f"scaled_error_period_30_checks={periodic_checks}")
    print(
        "correct_asymptotic=mu(E cap D_ga cap D_gb)="
        "mu(E)*rho(a,b)+O_E(1/g)"
    )
    print(
        "REFUTED: a uniform O_E(1/min(x,y)) deficit from "
        "(4/169)mu(E)"
    )


def scope_verdict() -> None:
    print("\nTHM-856 REFEREE VERDICT")
    print("RETAIN: Hunter/Kounias inequality and conditional coefficient arithmetic")
    print("RETAIN: first-moment seven-comb no-go and exact one-comb periodicity")
    print("RETAIN: the finite exact pilot packets as packet-specific computations")
    print("DELETE: global pair-overlap >=4/169 and (a+b)^2 leading term")
    print("DELETE: near-equal speeds are precisely the global minimum/failure locus")
    print("DELETE: raw-speed C(E)/x_min deficit lemma")
    print("REPLACE: stratify by reduced projective ratio (a,b) and common scale g")
    print("PASS: exact correction and projective obstruction independently verified")


if __name__ == "__main__":
    global_formula_audit()
    projective_ray_audit()
    scope_verdict()
