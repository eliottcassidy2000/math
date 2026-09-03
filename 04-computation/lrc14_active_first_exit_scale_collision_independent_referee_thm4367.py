#!/usr/bin/env python3
"""Independent exact referee for the THM-4367 arithmetic classification.

The only inherited mathematical input is THM-4365's proved first-exit law on
the common quotient fibre (odd P >= 11019): with

    1303 P = 47194 n + rho,       -23597 < rho <= 23597,

the strict active cells are |rho| < 3371 and have

    E_x(P) = 1303/47194 + (3371-rho)/(47194 P).

This program does not import the discovery/scout program.  It independently
checks the primitive-pair/scale round trip, strict boundaries, the least
collision, the sharp multiplicity bound and its complete equality case, and
physical binder/address controls.  The all-parameter implications themselves
are elementary divisibility arguments recorded in the theorem proof.

Nothing here advances LRC(14): every tested exit is a safe endpoint inherited
from THM-4365.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
from math import gcd, isqrt
import sys


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


CHECKS = 0


def require(condition: bool, message: str) -> None:
    """Optimization-safe assertion with an exact check count."""
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(message)


PRIME = 3371
MODULUS = 47194
NUMERATOR = 1303
THRESHOLD = 11019
ACTIVE_C_MIN = 2
ACTIVE_C_MAX = 6740
X1 = Fraction(NUMERATOR, MODULUS)

ROW_FIXED = (840, 3, 39, 11, 1691, 3371, 5051, 6731, 8411, 10091, 525, 945)
UNITS_MOD_14 = (1, 3, 5, 9, 11, 13)


@dataclass(frozen=True)
class ActiveRecord:
    parameter: int
    rho: int
    c: int
    gcd_scale: int
    a: int
    b: int
    kappa: int
    tooth_address: int
    exit_x: Fraction


def centered_rho(parameter: int) -> int:
    raw = (NUMERATOR * parameter) % MODULUS
    return raw if raw <= MODULUS // 2 else raw - MODULUS


def ceil_div(x: int, y: int) -> int:
    return -((-x) // y)


def candidate_a_values(b: int) -> tuple[int, ...]:
    """All even 2<=a<=6740 with 3371 | a+1303b."""
    residue = (-NUMERATOR * b) % PRIME
    answer = []
    for shift in (0, 1, 2):
        a = residue + shift * PRIME
        if ACTIVE_C_MIN <= a <= ACTIVE_C_MAX and a % 2 == 0:
            answer.append(a)
    return tuple(answer)


def admissible_scales(a: int, b: int) -> tuple[int, ...]:
    """All scales realizing primitive pair (a,b) in the proved fibre."""
    if (
        a < ACTIVE_C_MIN
        or a > ACTIVE_C_MAX
        or a % 2 != 0
        or b <= 0
        or b % 2 != 1
        or gcd(a, b) != 1
    ):
        return ()
    total = a + NUMERATOR * b
    if total % PRIME != 0:
        return ()
    kappa = total // PRIME
    if gcd(kappa, 14) != 1:
        return ()
    residue = pow(kappa, -1, 14)
    lower = ceil_div(THRESHOLD, b)
    upper = ACTIVE_C_MAX // a
    first = lower + ((residue - lower) % 14)
    if first > upper:
        return ()
    return tuple(range(first, upper + 1, 14))


def scale_count_formula(a: int, b: int) -> int:
    """Floor-function form of the admissible-scale count."""
    if (
        a < ACTIVE_C_MIN
        or a > ACTIVE_C_MAX
        or a % 2 != 0
        or b <= 0
        or b % 2 != 1
        or gcd(a, b) != 1
    ):
        return 0
    total = a + NUMERATOR * b
    if total % PRIME != 0:
        return 0
    kappa = total // PRIME
    if gcd(kappa, 14) != 1:
        return 0
    residue = pow(kappa, -1, 14)
    lower = ceil_div(THRESHOLD, b)
    upper = ACTIVE_C_MAX // a
    if upper < lower:
        return 0
    return max(
        0,
        (upper - residue) // 14 - ((lower - 1 - residue) // 14),
    )


def record_from_class(a: int, b: int, scale: int) -> ActiveRecord:
    scales = admissible_scales(a, b)
    require(scale in scales, "attempted reconstruction from an inadmissible scale")
    kappa = (a + NUMERATOR * b) // PRIME
    parameter = b * scale
    c = a * scale
    rho = PRIME - c
    require(scale * kappa % 14 == 1, "scale-address congruence failed")
    address = (scale * kappa - 1) // 14
    exit_x = Fraction(kappa, 14 * b)
    require(centered_rho(parameter) == rho, "reconstructed centered residue failed")
    require(
        NUMERATOR * parameter == MODULUS * address + rho,
        "reconstructed Euclidean address failed",
    )
    require(
        X1 + Fraction(c, MODULUS * parameter) == exit_x,
        "reconstructed metric exit failed",
    )
    require(
        Fraction(14 * address + 1, 14 * parameter) == exit_x,
        "reconstructed tooth wall failed",
    )
    return ActiveRecord(parameter, rho, c, scale, a, b, kappa, address, exit_x)


def record_from_parameter(parameter: int) -> ActiveRecord | None:
    require(parameter >= THRESHOLD and parameter % 2 == 1, "parameter outside fibre")
    rho = centered_rho(parameter)
    if abs(rho) >= PRIME:
        return None
    c = PRIME - rho
    scale = gcd(c, parameter)
    a = c // scale
    b = parameter // scale
    require(ACTIVE_C_MIN <= c <= ACTIVE_C_MAX and c % 2 == 0, "active c range failed")
    require(scale % 2 == 1 and a % 2 == 0 and b % 2 == 1, "parity split failed")
    require(gcd(a, b) == 1, "reduced pair is not primitive")
    require(scale < PRIME, "active scale reached the modulus prime")
    total = a + NUMERATOR * b
    require(total % PRIME == 0, "primitive divisibility failed")
    kappa = total // PRIME
    require(gcd(kappa, 14) == 1, "kappa is not a unit modulo fourteen")
    address = (NUMERATOR * parameter - rho) // MODULUS
    require(scale * kappa == 14 * address + 1, "address-scale identity failed")
    require(scale in admissible_scales(a, b), "forward scale missed converse set")
    exit_x = X1 + Fraction(c, MODULUS * parameter)
    require(exit_x == Fraction(kappa, 14 * b), "common-exit reduction failed")
    require(gcd(kappa, 14 * b) == 1, "common exit is not in lowest terms")
    return ActiveRecord(parameter, rho, c, scale, a, b, kappa, address, exit_x)


def circle_distance(speed: int, x: Fraction) -> Fraction:
    y = speed * x
    residue = y - (y.numerator // y.denominator)
    return min(residue, 1 - residue)


def binder_set(parameter: int, x: Fraction) -> tuple[Fraction, tuple[int, ...]]:
    values = tuple(
        (speed, circle_distance(speed, x))
        for speed in ROW_FIXED + (parameter,)
    )
    level = min(value for _, value in values)
    return level, tuple(speed for speed, value in values if value == level)


def residue_population(upper: int, residue: int, lower: int = 1) -> int:
    if upper < lower:
        return 0
    return max(
        0,
        (upper - residue) // 14 - ((lower - 1 - residue) // 14),
    )


def main() -> None:
    # Constant and primality firewall used when cancelling 3371 from a scale.
    require(MODULUS == 14 * PRIME, "modulus factorization changed")
    require(all(PRIME % p for p in range(2, isqrt(PRIME) + 1)), "3371 is not prime")
    require(gcd(NUMERATOR, MODULUS) == 1, "1303 ceased to be a modulus unit")
    require(pow(NUMERATOR, -1, PRIME) == 1485, "inverse modulo 3371 changed")

    # Strict boundary controls: neither c=0 nor c=6742 belongs to the active
    # classification.  Both physical rows exit at x1 with a two-binder tie.
    boundary_controls = []
    for parameter, expected_rho in ((43823, -3371), (50565, 3371)):
        require(centered_rho(parameter) == expected_rho, "boundary rho changed")
        require(record_from_parameter(parameter) is None, "boundary entered strict activity")
        level, binders = binder_set(parameter, X1)
        require(level == Fraction(1, 14), "boundary exit is not safe")
        require(set(binders) == {3371, parameter}, "boundary binder set changed")
        boundary_controls.append((parameter, expected_rho, binders))

    # Both strict active endpoints c=2 and c=6740 occur in the cofinite fibre.
    endpoint_controls = []
    for parameter, expected_rho, expected_c in (
        (47595, 3369, 2),
        (46793, -3369, 6740),
    ):
        rec = record_from_parameter(parameter)
        require(rec is not None, "strict endpoint unexpectedly inactive")
        require((rec.rho, rec.c) == (expected_rho, expected_c), "strict endpoint changed")
        level, binders = binder_set(parameter, rec.exit_x)
        require(level == Fraction(1, 14) and binders == (parameter,), "strict endpoint binder failed")
        endpoint_controls.append((parameter, rec.rho, rec.c, rec.tooth_address, rec.exit_x))

    # Three complete residue periods give a broad direct round-trip audit.
    direct_records = []
    first_active = None
    finite_limit = THRESHOLD + 3 * MODULUS - 2
    for parameter in range(THRESHOLD, finite_limit + 1, 2):
        rec = record_from_parameter(parameter)
        if rec is None:
            continue
        if first_active is None:
            first_active = rec
        rebuilt = record_from_class(rec.a, rec.b, rec.gcd_scale)
        require(rebuilt == rec, "parameter/class round trip failed")
        direct_records.append(
            (rec.parameter, rec.rho, rec.c, rec.gcd_scale, rec.a, rec.b, rec.kappa)
        )
    require(first_active is not None and first_active.parameter == 11045, "least active P changed")
    require(len(direct_records) == 10110, "three-period active count changed")

    # Independently sweep primitive pairs through b=50001 and reconstruct
    # every admissible scale.  This checks the converse without enumerating P
    # in its native order.
    reverse_pairs = 0
    reverse_scales = 0
    reverse_max = 0
    equality_pairs = []
    for b in range(1, 50002, 2):
        for a in candidate_a_values(b):
            if gcd(a, b) != 1:
                continue
            scales = admissible_scales(a, b)
            require(len(scales) == scale_count_formula(a, b), "count formula disagrees")
            if not scales:
                continue
            reverse_pairs += 1
            reverse_scales += len(scales)
            reverse_max = max(reverse_max, len(scales))
            if len(scales) == 241:
                equality_pairs.append((a, b))
            for scale in scales:
                rec = record_from_class(a, b, scale)
                require(
                    record_from_parameter(rec.parameter) == rec,
                    "class/parameter converse round trip failed",
                )
    require(reverse_pairs == 4634, "finite reverse-pair census changed")
    require(reverse_scales == 11373, "finite reverse-scale census changed")
    require(reverse_max == 241, "finite reverse sweep missed the sharp bound")

    # Least occupied denominator and least repeated parameter.  The b=75
    # class has four scales, and no b<75 has even one admissible scale.
    occupied_below_75 = []
    for b in range(1, 75, 2):
        for a in range(ACTIVE_C_MIN, ACTIVE_C_MAX + 1, 2):
            scales = admissible_scales(a, b)
            if scales:
                occupied_below_75.append((a, b, scales))
    require(not occupied_below_75, "an occupied primitive denominator below 75 appeared")

    least_scales = admissible_scales(34, 75)
    require(least_scales == (155, 169, 183, 197), "least collision scale set changed")
    least_records = tuple(record_from_class(34, 75, g) for g in least_scales)
    require(
        tuple((r.parameter, r.rho, r.c, r.tooth_address) for r in least_records)
        == (
            (11625, -1899, 5270, 321),
            (12675, -2375, 5746, 350),
            (13725, -2851, 6222, 379),
            (14775, -3327, 6698, 408),
        ),
        "least collision physical records changed",
    )
    require(all(r.exit_x == Fraction(29, 1050) for r in least_records), "least exits differ")
    for rec in least_records:
        level, binders = binder_set(rec.parameter, rec.exit_x)
        require(level == Fraction(1, 14), "least collision exit is not safe")
        require(binders == (rec.parameter,), "least collision lost its unique physical binder")

    # Least parameter which participates in any non-singleton metric class.
    # This is a different order statistic from the least occupied denominator
    # and from the least second parameter at which a repeated exit is visible.
    least_participating = None
    for parameter in range(THRESHOLD, 12676, 2):
        rec = record_from_parameter(parameter)
        if rec is None:
            continue
        scales = admissible_scales(rec.a, rec.b)
        if len(scales) >= 2:
            least_participating = (rec.parameter, rec.a, rec.b, rec.kappa, scales)
            break
    require(
        least_participating == (11085, 196, 2217, 857, (5, 19, 33)),
        "least participating parameter changed",
    )
    participant_records = tuple(record_from_class(196, 2217, g) for g in (5, 19, 33))
    require(
        tuple((r.parameter, r.rho, r.c, r.tooth_address) for r in participant_records)
        == (
            (11085, 2391, 980, 306),
            (42123, -353, 3724, 1163),
            (73161, -3097, 6468, 2020),
        ),
        "least participating class records changed",
    )
    require(
        all(r.exit_x == Fraction(857, 31038) for r in participant_records),
        "least participating class exits differ",
    )

    seen: dict[tuple[int, int], int] = {}
    first_repeat = None
    for parameter in range(THRESHOLD, 12676, 2):
        rec = record_from_parameter(parameter)
        if rec is None:
            continue
        key = (rec.a, rec.b)
        if key in seen:
            first_repeat = (seen[key], parameter, rec.a, rec.b, rec.kappa)
            break
        seen[key] = parameter
    require(
        first_repeat == (11625, 12675, 34, 75, 29),
        "least second-parameter collision changed",
    )

    # Dynamic hostile: translating every parameter in the least-denominator
    # quartet by two retains strict activity and the common Q fibre, but not
    # the metric collision.  Here the tooth addresses happen to stay fixed;
    # the four changed denominators already separate the right walls.
    shifted_expectations = (
        (707, Fraction(4495, 162778)),
        (231, Fraction(4901, 177478)),
        (-245, Fraction(5307, 192178)),
        (-721, Fraction(5713, 206878)),
    )
    shifted_records = []
    for old, (expected_rho, expected_exit) in zip(least_records, shifted_expectations):
        shifted = record_from_parameter(old.parameter + 2)
        require(shifted is not None, "translated quartet member became inactive")
        require(shifted.rho == expected_rho, "translated quartet rho changed")
        require(shifted.exit_x == expected_exit, "translated quartet exit changed")
        require(
            shifted.tooth_address == old.tooth_address,
            "translated quartet did not retain the hostile address",
        )
        level, binders = binder_set(shifted.parameter, shifted.exit_x)
        require(
            level == Fraction(1, 14) and binders == (shifted.parameter,),
            "translated quartet physical binder failed",
        )
        shifted_records.append((shifted.parameter, shifted.rho, shifted.tooth_address, shifted.exit_x))
    require(
        len({row[3] for row in shifted_records}) == 4,
        "translated quartet retained an exit collision",
    )

    # Residue phase alone is also too coarse for the metric consumer.
    phase_records = tuple(record_from_parameter(parameter) for parameter in (48679, 95873))
    require(all(rec is not None for rec in phase_records), "same-phase control became inactive")
    require(tuple(rec.rho for rec in phase_records) == (1, 1), "same-phase rhos changed")
    require(
        tuple(rec.exit_x for rec in phase_records)
        == (Fraction(18817, 681506), Fraction(37059, 1342222)),
        "same-phase exits changed",
    )
    require(phase_records[0].exit_x != phase_records[1].exit_x, "phase alone determined exit")

    # Global upper bound: a is positive even, hence g<=6740/a<=3370.
    # A unit residue modulo 14 has at most 241 representatives in that range;
    # for a>=4 it has at most 121, so equality forces a=2.
    for a in range(ACTIVE_C_MIN, ACTIVE_C_MAX + 1, 2):
        upper = ACTIVE_C_MAX // a
        maximum_residue_count = max(residue_population(upper, r) for r in UNITS_MOD_14)
        require(maximum_residue_count <= 241, "global scale bound failed")
        if a >= 4:
            require(maximum_residue_count <= 121, "a>=4 equality exclusion failed")

    # Complete a=2 equality analysis.  Odd integrality is exactly
    # b=401+6742t, kappa=155+2606t, t>=0.  The t=0 threshold exception and
    # every other residue class modulo seven are checked explicitly.
    require((-2 * pow(NUMERATOR, -1, PRIME)) % PRIME == 401, "a=2 base residue changed")
    a2_cycle = []
    for t in range(0, 701):
        b = 401 + 6742 * t
        kappa = 155 + 2606 * t
        require(2 + NUMERATOR * b == PRIME * kappa, "a=2 parametrization identity failed")
        scales = admissible_scales(2, b)
        if t == 0:
            expected = 239
        elif t % 7 in (0, 1, 2, 5):
            expected = 241
        elif t % 7 in (4, 6):
            expected = 240
        else:
            expected = 0
        require(len(scales) == expected, "a=2 multiplicity cycle failed")
        require(len(scales) == scale_count_formula(2, b), "a=2 floor formula failed")
        if t <= 6:
            residue = None if gcd(kappa, 14) != 1 else pow(kappa, -1, 14)
            a2_cycle.append((t, kappa % 14, residue, len(scales)))
    require(
        tuple(a2_cycle)
        == (
            (0, 1, 1, 239),
            (1, 3, 5, 241),
            (2, 5, 3, 241),
            (3, 7, None, 0),
            (4, 9, 11, 240),
            (5, 11, 9, 241),
            (6, 13, 13, 240),
        ),
        "a=2 seven-class table changed",
    )

    # Every equality pair in the independent finite sweep obeys the claimed
    # iff.  The preceding bound and a=2 parametrization prove the unbounded
    # direction once the same elementary identities are read symbolically.
    for a, b in equality_pairs:
        require(a == 2 and (b - 401) % 6742 == 0, "finite equality pair left a=2 family")
        t = (b - 401) // 6742
        require(t >= 1 and t % 7 in (0, 1, 2, 5), "finite equality residue failed")

    max_scales = admissible_scales(2, 7143)
    require(
        len(max_scales) == 241
        and max_scales[0] == 5
        and max_scales[-1] == 3365,
        "least maximum-multiplicity class changed",
    )
    max_first = record_from_class(2, 7143, max_scales[0])
    max_last = record_from_class(2, 7143, max_scales[-1])
    require(max_first.exit_x == max_last.exit_x == Fraction(2761, 100002), "max exits differ")
    require(
        (max_first.parameter, max_first.rho, max_first.tooth_address)
        == (35715, 3361, 986),
        "first maximum-class record changed",
    )
    require(
        (max_last.parameter, max_last.rho, max_last.tooth_address)
        == (24036195, -3359, 663626),
        "last maximum-class record changed",
    )
    for rec in (max_first, max_last):
        level, binders = binder_set(rec.parameter, rec.exit_x)
        require(level == Fraction(1, 14) and binders == (rec.parameter,), "max-class binder failed")

    # Infinite sharp family: t=1+7s gives b=7143+47194s,
    # kappa=2761+18242s and the same 241 scales 5+14j (0<=j<=240).
    infinite_controls = []
    previous_exit = None
    for s in (0, 1, 2, 17):
        t = 1 + 7 * s
        b = 7143 + MODULUS * s
        kappa = 2761 + 18242 * s
        require((b, kappa) == (401 + 6742 * t, 155 + 2606 * t), "infinite family formula failed")
        scales = admissible_scales(2, b)
        require(scales == tuple(5 + 14 * j for j in range(241)), "infinite scale set changed")
        exit_x = Fraction(kappa, 14 * b)
        require(exit_x == X1 + Fraction(2, MODULUS * b), "infinite exit identity failed")
        if previous_exit is not None:
            require(exit_x < previous_exit, "infinite maximum exits do not decrease")
        previous_exit = exit_x
        first = record_from_class(2, b, scales[0])
        last = record_from_class(2, b, scales[-1])
        infinite_controls.append((s, b, kappa, len(scales), first.parameter, last.parameter, exit_x))
    require(2606 * 401 - 6742 * 155 == -4, "a=2 rational classes ceased to be distinct")

    direct_hash = sha256(repr(tuple(direct_records)).encode("ascii")).hexdigest()
    reverse_summary = (reverse_pairs, reverse_scales, reverse_max, tuple(equality_pairs))
    reverse_hash = sha256(repr(reverse_summary).encode("ascii")).hexdigest()

    print("THM4367 ACTIVE FIRST-EXIT COLLISION CLEANROOM REFEREE")
    print("scope odd P>=11019 in THM-4365 common quotient fibre; strict |rho|<3371")
    print("classification c=3371-rho=ag; P=bg; gcd(a,b)=1; a even; b odd")
    print("kappa=(a+1303b)/3371; gcd(kappa,14)=1; g*kappa=1 mod 14")
    print("scale_window ceil(11019/b)<=g<=floor(6740/a)")
    print("metric_exit kappa/(14b); physical_binder bg; physical_address (g*kappa-1)/14")
    print("strict_boundary_controls", tuple(boundary_controls))
    print("strict_active_endpoint_controls", tuple(endpoint_controls))
    print("direct_three_period_records", len(direct_records), "first_active", first_active.parameter)
    print("direct_three_period_sha256", direct_hash)
    print("reverse_b_le_50001", reverse_pairs, reverse_scales, reverse_max)
    print("reverse_summary_sha256", reverse_hash)
    print("least_occupied_primitive_denominator", (34, 75, 29), "scales", least_scales)
    print(
        "least_collision_records",
        tuple((r.parameter, r.rho, r.c, r.tooth_address) for r in least_records),
        "exit",
        least_records[0].exit_x,
    )
    print("least_parameter_in_nonsingleton_class", least_participating)
    print(
        "least_participating_records",
        tuple((r.parameter, r.rho, r.c, r.tooth_address) for r in participant_records),
        "exit",
        participant_records[0].exit_x,
    )
    print("least_second_parameter_collision", first_repeat)
    print("plus_two_dynamic_hostile", tuple(shifted_records))
    print(
        "same_residue_phase_hostile",
        tuple((r.parameter, r.rho, r.exit_x) for r in phase_records),
    )
    print("sharp_global_multiplicity", 241)
    print("a2_cycle_t_0_to_6", tuple(a2_cycle))
    print(
        "equality_iff",
        "a=2; b=401+6742t; kappa=155+2606t; t>=1; t mod 7 in {0,1,2,5}",
    )
    print(
        "least_241_class",
        (2, 7143, 2761),
        "scale_count",
        len(max_scales),
        "scale_endpoints",
        (max_scales[0], max_scales[-1]),
        "exit",
        max_first.exit_x,
    )
    print("infinite_241_family_controls", tuple(infinite_controls))
    print("sidecar_separation", "same metric; distinct binder P=bg and address n=(g*kappa-1)/14")
    print("checks", CHECKS)
    print("PASS")
    print("NO_LRC14_DECREMENT; NO_COUNTEREXAMPLE; ACTIVE METRIC COINCIDENCE ONLY")


if __name__ == "__main__":
    main()
