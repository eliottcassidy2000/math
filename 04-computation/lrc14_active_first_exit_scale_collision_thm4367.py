#!/usr/bin/env python3
"""Exact audit for THM-4367's active first-exit collision classification.

The input is only THM-4365's proved arithmetic law

    rho = <1303 P>_(47194),  |rho| < 3371,
    E_x = 1303/47194 + (3371-rho)/(47194 P)

on odd P >= 11019.  This script does not rerun the 828-component geometry
and makes no LRC(14) claim.
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction as F
from hashlib import sha256
from math import gcd, isqrt


A = 3371
M = 14 * A
S = 1303
TAIL = 11019
ACTIVE_C_MAX = 2 * A - 2
CHECKS = 0


def require(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(message)


def ceil_div(x: int, y: int) -> int:
    return (x + y - 1) // y


def centered(value: int) -> int:
    residue = value % M
    if 2 * residue > M:
        residue -= M
    return residue


def active_record(P: int):
    """Return (rho,c,a,b,g,kappa,n,exit), or None off an active cell."""
    require(P >= TAIL and P % 2 == 1, "parameter outside odd cofinite tail")
    rho = centered(S * P)
    if abs(rho) >= A:
        return None
    c = A - rho
    g = gcd(c, P)
    a, b = c // g, P // g
    require(a > 0 and a % 2 == 0, "a is not positive even")
    require(b > 0 and b % 2 == 1 and g > 0 and g % 2 == 1, "bad odd scale")
    require(gcd(a, b) == 1, "fraction was not reduced")
    require(gcd(g, A) == 1, "active scale unexpectedly divisible by A")
    require((a + S * b) % A == 0, "A does not divide reduced wall numerator")
    kappa = (a + S * b) // A
    require(kappa > 0 and (g * kappa - 1) % 14 == 0, "scale congruence failed")
    n = (g * kappa - 1) // 14
    exit_point = F(S, M) + F(c, M * P)
    require(exit_point == F(kappa, 14 * b), "canonical exit reduction failed")
    require(gcd(kappa, 14 * b) == 1, "canonical exit is not in lowest terms")
    require(n == (S * P - rho) // M, "physical address changed")
    return rho, c, a, b, g, kappa, n, exit_point


def scale_class(a: int, b: int, kappa: int) -> tuple[int, ...]:
    """All tail-active scales in one reduced metric class."""
    require(a > 0 and a % 2 == 0, "a must be positive even")
    require(b > 0 and b % 2 == 1 and gcd(a, b) == 1, "bad reduced pair")
    require(a + S * b == A * kappa, "bad kappa")
    require(kappa > 0 and gcd(kappa, 14) == 1, "kappa is not a positive unit")
    g0 = pow(kappa, -1, 14)
    lower = ceil_div(TAIL, b)
    upper = ACTIVE_C_MAX // a
    if lower > upper:
        return ()
    first = g0
    if first < lower:
        first += 14 * ceil_div(lower - first, 14)
    return tuple(range(first, upper + 1, 14))


def main() -> None:
    require(M == 47194, "modulus changed")
    require(all(A % p for p in range(2, isqrt(A) + 1)), "3371 is not prime")
    require(pow(S, -1, M) == 1485, "modular inverse changed")
    require(gcd(S, M) == 1, "residue action not invertible")

    # If P0 is the least member of a nontrivial class, its next scale is
    # g+14, so P1/P0 <= 15.  This interval proves both minimality claims.
    minimality_bound = 15 * 11625
    by_correction: dict[F, list[int]] = defaultdict(list)
    by_exit: dict[F, list[int]] = defaultdict(list)
    records = {}
    active_count = 0
    for P in range(TAIL, minimality_bound + 1, 2):
        record = active_record(P)
        if record is None:
            continue
        active_count += 1
        records[P] = record
        rho, c, a, b, g, kappa, n, exit_point = record
        correction = F(c, P)
        require(correction == F(a, b), "reduced correction changed")
        require(A - a * g == rho and P == b * g, "source reconstruction failed")
        require(n == (g * kappa - 1) // 14, "address reconstruction failed")
        by_correction[correction].append(P)
        by_exit[exit_point].append(P)

    require(tuple(sorted(by_correction.values())) == tuple(sorted(by_exit.values())),
            "correction and metric partitions differ")
    for correction, parameters in by_correction.items():
        first_record = records[parameters[0]]
        a, b, kappa = first_record[2], first_record[3], first_record[5]
        predicted = tuple(
            b * g for g in scale_class(a, b, kappa)
            if b * g <= minimality_bound
        )
        require(tuple(parameters) == predicted, f"class mismatch at {correction}")

    collision_groups = tuple(
        (correction, tuple(parameters))
        for correction, parameters in sorted(by_correction.items())
        if len(parameters) >= 2
    )
    least_participating_pair = min(
        (parameters[0], parameters[1], correction)
        for correction, parameters in collision_groups
    )
    require(least_participating_pair == (11085, 42123, F(196, 2217)),
            "least participating collision changed")
    earliest_repeat_pair = min(
        (parameters[1], parameters[0], correction)
        for correction, parameters in collision_groups
    )
    require(earliest_repeat_pair == (12675, 11625, F(34, 75)),
            "earliest repeated exit changed")

    least_class = (11085, 42123, 73161)
    least_records = tuple(active_record(P) for P in least_class)
    require(tuple(record[0] for record in least_records) == (2391, -353, -3097),
            "least-class residues changed")
    require(tuple(record[4] for record in least_records) == (5, 19, 33),
            "least-class scales changed")
    require(tuple(record[6] for record in least_records) == (306, 1163, 2020),
            "least-class addresses changed")
    require({record[7] for record in least_records} == {F(857, 31038)},
            "least-class exits differ")
    require(scale_class(196, 2217, 857) == (5, 19, 33), "least class incomplete")

    quartet = (11625, 12675, 13725, 14775)
    quartet_records = tuple(active_record(P) for P in quartet)
    require(tuple(record[0] for record in quartet_records)
            == (-1899, -2375, -2851, -3327), "quartet residues changed")
    require(tuple(record[4] for record in quartet_records) == (155, 169, 183, 197),
            "quartet scales changed")
    require(tuple(record[6] for record in quartet_records) == (321, 350, 379, 408),
            "quartet addresses changed")
    require({record[7] for record in quartet_records} == {F(29, 1050)},
            "quartet exits differ")
    require(scale_class(34, 75, 29) == (155, 169, 183, 197),
            "quartet not complete")

    # Equal metric exits are not an operation congruence for the natural
    # odd-tail translation P -> P+2: one step separates the whole quartet.
    successors = tuple(P + 2 for P in quartet)
    successor_records = tuple(active_record(P) for P in successors)
    require(all(record is not None for record in successor_records),
            "a successor left the strict-active chart")
    require(tuple(record[0] for record in successor_records) == (707, 231, -245, -721),
            "successor residues changed")
    successor_exits = tuple(record[7] for record in successor_records)
    require(successor_exits == (F(4495, 162778), F(4901, 177478),
                                F(5307, 192178), F(5713, 206878)),
            "successor exits changed")
    require(len(set(successor_exits)) == 4, "translation preserved the metric collision")

    same_phase = (48679, 95873)
    same_phase_records = tuple(active_record(P) for P in same_phase)
    require(tuple(record[0] for record in same_phase_records) == (1, 1),
            "rho=1 hostile changed")
    require(tuple(record[7] for record in same_phase_records)
            == (F(18817, 681506), F(37059, 1342222)),
            "same-phase hostile exits changed")

    class_239 = scale_class(2, 401, 155)
    require(len(class_239) == 239 and class_239[0] == 29 and class_239[-1] == 3361,
            "239-cell class changed")
    require(F(155, 5614) == active_record(401 * class_239[0])[7]
            == active_record(401 * class_239[-1])[7], "239-cell exit changed")

    class_241 = scale_class(2, 7143, 2761)
    require(len(class_241) == 241 and class_241[0] == 5 and class_241[-1] == 3365,
            "maximum collision class changed")
    require(F(2761, 100002) == active_record(7143 * class_241[0])[7]
            == active_record(7143 * class_241[-1])[7], "maximum exit changed")
    require(max(sum(1 for g in range(r, 3371, 14))
                for r in (1, 3, 5, 9, 11, 13)) == 241,
            "universal multiplicity bound changed")

    # Equality forces a=2.  The divisibility and oddness conditions then give
    # b=401+6742t and kappa=155+2606t.  Residues of t modulo seven determine
    # whether the inverse scale class contains 241, 240, or no points.
    require(2 + S * 401 == A * 155, "equality-family seed changed")
    equality_residues = {0, 1, 2, 5}
    equality_controls = []
    for t in range(15):
        b = 401 + 2 * A * t
        kappa = 155 + 2 * S * t
        if t % 7 == 3:
            require(gcd(kappa, 14) != 1, "excluded equality residue became a unit")
            equality_controls.append((t, "inadmissible"))
            continue
        scales = scale_class(2, b, kappa)
        expected = 239 if t == 0 else (241 if t % 7 in equality_residues else 240)
        require(len(scales) == expected, "equality-class residue count changed")
        equality_controls.append((t, len(scales)))

    boundary = ((43823, -A, 1210), (50565, A, 1396))
    for P, expected_rho, expected_n in boundary:
        require(centered(S * P) == expected_rho, "boundary residue changed")
        require((S * P - expected_rho) // M == expected_n, "boundary address changed")
        require(active_record(P) is None, "open boundary became active")
        require(P - M < TAIL, "boundary is not least in the cofinite tail")

    infinite_controls = []
    for t in (0, 1, 2, 4, 5):
        a = 2
        b = 401 + 2 * A * t
        kappa = 155 + 2 * S * t
        scales = scale_class(a, b, kappa)
        require(len(scales) >= 2, f"family control t={t} is not a collision")
        infinite_controls.append((t, b, kappa, len(scales), scales[0], scales[-1]))

    serialized_collisions = repr(collision_groups).encode()
    print("THM-4367 PRIMARY: ACTIVE FIRST-EXIT SCALE COLLISIONS CLASSIFIED")
    print("SCOPE=THM-4365 odd active cells P>=11019; arithmetic only; LRC14_OPEN")
    print(f"BOUNDED=[{TAIL},{minimality_bound}]_odd;active={active_count};"
          f"collision_classes={len(collision_groups)}")
    print("NORMAL_FORM=c=ag,P=bg,a+1303b=3371kappa,gkappa=1_mod14")
    print("METRIC_EXIT=kappa/(14b); already reduced; equal iff (kappa,b) equal")
    print("SCALE_CLASS=g=inverse(kappa)_mod14 within ceil(11019/b)..floor(6740/a)")
    print("LEAST_PARTICIPANT=11085;class=(11085,42123,73161);exit=857/31038")
    print("LEAST_CLASS_BINDERS=(11085@306,42123@1163,73161@2020)")
    print("EARLIEST_REPEAT=(11625,12675);class_through=14775;exit=29/1050")
    print("EARLIEST_BINDERS=(11625@321,12675@350,13725@379,14775@408)")
    print("SHIFT_HOSTILE=P_PLUS_2:(11627,12677,13727,14777);four_distinct_exits")
    print("PHASE_HOSTILE=P:(48679,95873);rho=1;exits_distinct")
    print("MAX_MULTIPLICITY=241 attained_by=(a,b,kappa)=(2,7143,2761)")
    print("MAX_CLASS=g=5..3365_step14;P=35715..24036195;exit=2761/100002")
    print("EQUALITY_CLASSES=a=2,b=401+6742t,kappa=155+2606t;"
          "t>=1_and_t_mod7_in_{0,1,2,5}")
    print(f"BOUNDARY_CONTROLS={boundary};exit=1303/47194;binders={{3371,P}}")
    print(f"INFINITE_FAMILY_CONTROLS={tuple(infinite_controls)}")
    print(f"BOUNDED_COLLISION_SHA256={sha256(serialized_collisions).hexdigest()}")
    print(f"CHECKS={CHECKS}")
    print("NO_LRC14_DECREMENT;NO_828_GEOMETRY_REAUDIT")


if __name__ == "__main__":
    main()
