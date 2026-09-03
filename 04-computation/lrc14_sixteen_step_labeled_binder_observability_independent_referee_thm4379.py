#!/usr/bin/env python3
"""Import-free clean-room referee for the reserved THM-4379 candidate.

This script was designed before the candidate verifier was inspected.  It
uses normalized integer rational pairs, a fresh open-interval renewal, an
exhaustive phase-circle decoder, and a constructive audit of every active
metric-fibre shape.  It does not import candidate or canonical computation.
All claims remain relative to proved THM-4365/4367 geometry, and LRC(14)
remains open.
"""

from __future__ import annotations

from collections import Counter
from math import gcd
import sys


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")


A = 3371
S = 1303
M = 47194
N = 23597
TAIL = 11019
HEIGHT = 420
ANCHOR = 840
FIXED = (3, 39, 11, 1691, 3371, 5051, 6731, 8411, 10091, 525, 945)
K = 23
EPSILON = 0
UNITS_14 = (1, 3, 5, 9, 11, 13)

CHECKS = 0


def check(condition: bool, label: str) -> None:
    """Counted assertion which remains live under python -O."""
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(f"clean-room check failed: {label}")


def rational(numerator: int, denominator: int) -> tuple[int, int]:
    if denominator == 0:
        raise ZeroDivisionError
    if denominator < 0:
        numerator = -numerator
        denominator = -denominator
    divisor = gcd(abs(numerator), denominator)
    return numerator // divisor, denominator // divisor


def qeq(left: tuple[int, int], right: tuple[int, int]) -> bool:
    return left[0] * right[1] == right[0] * left[1]


def qlt(left: tuple[int, int], right: tuple[int, int]) -> bool:
    return left[0] * right[1] < right[0] * left[1]


def qle(left: tuple[int, int], right: tuple[int, int]) -> bool:
    return left[0] * right[1] <= right[0] * left[1]


def qsub(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
    return rational(left[0] * right[1] - right[0] * left[1], left[1] * right[1])


def qfloor(value: tuple[int, int]) -> int:
    return value[0] // value[1]


BASE_EXIT = rational(S, M)
CLEARANCE = rational(1, 14)
COMPONENT_LEFT = rational(14 * K + 1, 28 * HEIGHT)
COMPONENT_RIGHT = rational(14 * K + 13, 28 * HEIGHT)
SECOND_FRONTIER = rational(73, 2646)
NEXT_FIXED_LEFT = rational(2603, 94234)


def centered(value: int) -> int:
    residue = value % M
    return residue - M if residue > N else residue


def rho_of(parameter: int) -> int:
    return centered(S * parameter)


def address_of(parameter: int) -> int:
    rho = rho_of(parameter)
    quotient, remainder = divmod(S * parameter - rho, M)
    check(remainder == 0, "address integrality")
    return quotient


def phase_of(parameter: int) -> int:
    rho = rho_of(parameter)
    check((A - rho) % 2 == 0, "half-phase parity")
    return ((A - rho) // 2) % N


def tooth(speed: int, address: int) -> tuple[int, int, tuple[int, int], tuple[int, int]]:
    return (
        speed,
        address,
        rational(14 * address - 1, 14 * speed),
        rational(14 * address + 1, 14 * speed),
    )


def intersects_open(
    item: tuple[int, int, tuple[int, int], tuple[int, int]],
    left: tuple[int, int],
    right: tuple[int, int],
) -> bool:
    return qlt(item[2], right) and qlt(left, item[3])


def fixed_teeth() -> tuple[tuple[int, int, tuple[int, int], tuple[int, int]], ...]:
    answer = []
    for speed in FIXED:
        first = qfloor(rational(speed * COMPONENT_LEFT[0], COMPONENT_LEFT[1])) - 3
        last = qfloor(rational(speed * COMPONENT_RIGHT[0], COMPONENT_RIGHT[1])) + 3
        for address in range(first, last + 1):
            item = tooth(speed, address)
            if intersects_open(item, COMPONENT_LEFT, COMPONENT_RIGHT):
                answer.append(item)
    return tuple(answer)


LOCAL_FIXED_TEETH = fixed_teeth()


def active_at(
    item: tuple[int, int, tuple[int, int], tuple[int, int]],
    point: tuple[int, int],
) -> bool:
    return qlt(item[2], point) and qlt(point, item[3])


def parameter_tooth_at(
    parameter: int, point: tuple[int, int]
) -> tuple[int, int, tuple[int, int], tuple[int, int]] | None:
    center_floor = (parameter * point[0]) // point[1]
    live = []
    for address in range(center_floor - 2, center_floor + 3):
        item = tooth(parameter, address)
        if active_at(item, point):
            live.append(item)
    check(len(live) <= 1, "one open tooth per speed at a frontier")
    return live[0] if live else None


def better_tooth(
    challenger: tuple[int, int, tuple[int, int], tuple[int, int]],
    incumbent: tuple[int, int, tuple[int, int], tuple[int, int]],
) -> bool:
    """Farthest right, then wider, then smaller physical label/address."""
    if not qeq(challenger[3], incumbent[3]):
        return qlt(incumbent[3], challenger[3])
    if not qeq(challenger[2], incumbent[2]):
        return qlt(challenger[2], incumbent[2])
    if challenger[0] != incumbent[0]:
        return challenger[0] < incumbent[0]
    return challenger[1] < incumbent[1]


def owner(speed: int, address: int) -> int:
    return (2 * address - EPSILON * speed) % 2


def renew_direct(parameter: int) -> tuple[str, tuple[tuple[int, int, int], ...], tuple[tuple[int, int], ...], tuple[int, int]]:
    frontier = COMPONENT_LEFT
    chain = []
    frontiers = []
    while True:
        live = [item for item in LOCAL_FIXED_TEETH if active_at(item, frontier)]
        varying = parameter_tooth_at(parameter, frontier)
        if varying is not None:
            live.append(varying)
        if not live:
            return "missing", tuple(chain), tuple(frontiers), frontier
        winner = live[0]
        for item in live[1:]:
            if better_tooth(item, winner):
                winner = item
        frontiers.append(frontier)
        chain.append((winner[0], winner[1], owner(winner[0], winner[1])))
        check(qlt(frontier, winner[3]), "renewal advances")
        frontier = winner[3]
        if qlt(COMPONENT_RIGHT, frontier):
            return ("span" if len(chain) == 1 else "renew"), tuple(chain), tuple(frontiers), frontier


def distance_and_wall(
    speed: int, point: tuple[int, int]
) -> tuple[tuple[int, int], int, str | None]:
    numerator = speed * point[0]
    denominator = point[1]
    address_floor, residue = divmod(numerator, denominator)
    if residue <= denominator - residue:
        distance = rational(residue, denominator)
        address = address_floor
        side = "R" if 14 * residue == denominator else None
    else:
        distance = rational(denominator - residue, denominator)
        address = address_floor + 1
        side = "L" if 14 * residue == 13 * denominator else None
    return distance, address, side


def all_binders(parameter: int, point: tuple[int, int]) -> tuple[tuple[int, int], tuple[tuple[int, int, str, int], ...]]:
    rows = []
    for speed in (ANCHOR,) + FIXED + (parameter,):
        distance, address, side = distance_and_wall(speed, point)
        rows.append((speed, distance, address, side))
    least = rows[0][1]
    for _, distance, _, _ in rows[1:]:
        if qlt(distance, least):
            least = distance
    binders = []
    for speed, distance, address, side in rows:
        if qeq(distance, least):
            check(side in ("L", "R"), "minimum is an oriented 1/14 wall")
            binders.append((speed, address, side, owner(speed, address)))
    return least, tuple(sorted(binders))


def make_record(
    status: str,
    chain: tuple[tuple[int, int, int], ...],
    frontiers: tuple[tuple[int, int], ...],
    exit_point: tuple[int, int],
    clearance: tuple[int, int],
    binders: tuple[tuple[int, int, str, int], ...],
) -> tuple[object, ...]:
    return (K, K % HEIGHT, EPSILON, status, chain, frontiers, exit_point, clearance, binders)


def direct_record(parameter: int) -> tuple[object, ...]:
    status, chain, frontiers, exit_point = renew_direct(parameter)
    clearance, binders = all_binders(parameter, exit_point)
    return make_record(status, chain, frontiers, exit_point, clearance, binders)


FIXED_CHAIN = ((945, 26, 0), (A, 93, 0))
FIXED_FRONTIERS = (COMPONENT_LEFT, SECOND_FRONTIER)
FIXED_BINDER = (A, 93, "R", 0)


def formula_record(parameter: int) -> tuple[object, ...]:
    rho = rho_of(parameter)
    address = address_of(parameter)
    phase = phase_of(parameter)
    chain = FIXED_CHAIN
    frontiers = FIXED_FRONTIERS
    if 1 <= phase <= A - 1:
        chain = chain + ((parameter, address, owner(parameter, address)),)
        frontiers = frontiers + (BASE_EXIT,)
        exit_point = rational(14 * address + 1, 14 * parameter)
        binders = ((parameter, address, "R", owner(parameter, address)),)
    elif phase == 0:
        check(rho == A, "phase zero is the positive boundary")
        exit_point = BASE_EXIT
        binders = tuple(sorted((FIXED_BINDER, (parameter, address, "R", owner(parameter, address)))))
    elif phase == A:
        check(rho == -A, "phase A is the negative boundary")
        exit_point = BASE_EXIT
        binders = tuple(sorted((FIXED_BINDER, (parameter, address, "L", owner(parameter, address)))))
    else:
        check(abs(rho) > A, "remaining phase is strictly inactive")
        exit_point = BASE_EXIT
        binders = (FIXED_BINDER,)
    return make_record("missing", chain, frontiers, exit_point, CLEARANCE, binders)


def record_binders(record: tuple[object, ...]) -> tuple[tuple[int, int, str, int], ...]:
    return record[8]  # type: ignore[return-value]


def record_exit(record: tuple[object, ...]) -> tuple[int, int]:
    return record[6]  # type: ignore[return-value]


def visible(record: tuple[object, ...], parameter: int) -> bool:
    return any(atom[0] == parameter for atom in record_binders(record))


def marked(phase: int) -> bool:
    return 0 <= phase % N <= A


def output_word(parameter: int, horizon: int) -> tuple[tuple[object, ...], ...]:
    return tuple(formula_record(parameter + 2 * time) for time in range(horizon + 1))


def decode(word: tuple[tuple[object, ...], ...]) -> int | None:
    for time, record in enumerate(word):
        labels = tuple(atom[0] for atom in record_binders(record) if atom[0] != A)
        if labels:
            check(len(labels) == 1, "one varying binder label per marked record")
            return labels[0] - 2 * time
    return None


def extended_gcd(a: int, b: int) -> tuple[int, int, int]:
    old_r, r = a, b
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r:
        quotient = old_r // r
        old_r, r = r, old_r - quotient * r
        old_s, s = s, old_s - quotient * s
        old_t, t = t, old_t - quotient * t
    return old_r, old_s, old_t


def inverse_mod(a: int, modulus: int) -> int:
    divisor, coefficient, _ = extended_gcd(a, modulus)
    check(divisor == 1, "modular inverse exists")
    return coefficient % modulus


def role_normalized_record(parameter: int) -> tuple[object, ...]:
    """Erase the varying numeric speed/address everywhere in the full record."""
    record = formula_record(parameter)
    chain = tuple(
        (("P", None, tooth_owner) if speed == parameter else (speed, address, tooth_owner))
        for speed, address, tooth_owner in record[4]  # type: ignore[union-attr]
    )
    binders = tuple(
        (("P", None, side, tooth_owner) if speed == parameter else (speed, address, side, tooth_owner))
        for speed, address, side, tooth_owner in record_binders(record)
    )
    return record[:4] + (chain,) + record[5:8] + (binders,)


def main() -> None:
    check(M == 14 * A and N == M // 2 == 7 * A, "published moduli")
    check(gcd(S, M) == 1, "primitive phase multiplier")
    check(BASE_EXIT == rational(1303, 47194), "base exit")
    check(COMPONENT_LEFT == rational(323, 11760), "component left")
    check(COMPONENT_RIGHT == rational(67, 2352), "component right")
    check(tooth(945, 26)[3] == SECOND_FRONTIER, "first selected frontier")
    check(tooth(A, 93)[3] == BASE_EXIT, "fixed terminal right wall")
    check(tooth(6731, 186)[2] == NEXT_FIXED_LEFT, "next fixed left wall")
    check(qlt(rational(1, 7 * TAIL), qsub(NEXT_FIXED_LEFT, BASE_EXIT)), "whole-tail no-bridge bound")

    fixed_equalities = []
    for speed in (ANCHOR,) + FIXED:
        distance, address, side = distance_and_wall(speed, BASE_EXIT)
        check(qle(CLEARANCE, distance), "fixed bank is safe at the base exit")
        if qeq(distance, CLEARANCE):
            fixed_equalities.append((speed, address, side, owner(speed, address)))
    check(tuple(fixed_equalities) == (FIXED_BINDER,), "unique fixed binder at base exit")

    # Build the phase inverse by a full odd-residue enumeration, not by the
    # candidate's closed inverse formula.
    phase_to_residue: dict[int, int] = {}
    for residue in range(1, M, 2):
        z = phase_of(residue)
        check(z not in phase_to_residue, "phase map injective on odd residues")
        phase_to_residue[z] = residue
    check(set(phase_to_residue) == set(range(N)), "phase map is onto")

    def representative(z: int, lift: int = 0) -> int:
        parameter = phase_to_residue[z % N]
        if parameter < TAIL:
            parameter += ((TAIL - parameter + M - 1) // M) * M
        parameter += lift * M
        check(parameter >= TAIL and parameter % 2 == 1, "odd tail representative")
        check(phase_of(parameter) == z % N, "representative has requested phase")
        return parameter

    # Fresh direct renewal and all-thirteen-distance comparison at every
    # phase and at two M-separated scales.
    geometry_rows = 0
    owner_values = set()
    for z in range(N):
        for lift in (0, 1):
            parameter = representative(z, lift)
            raw = direct_record(parameter)
            closed = formula_record(parameter)
            check(raw == closed, f"direct/formula record P={parameter}")
            check(raw[3] == "missing", "component 23 status")
            check(raw[7] == CLEARANCE, "exact safe clearance")
            owner_values.update(item[2] for item in raw[4])  # type: ignore[union-attr]
            owner_values.update(item[3] for item in record_binders(raw))
            check(visible(raw, parameter) == marked(z), "closed binder visibility")
            geometry_rows += 1
    check(geometry_rows == 2 * N, "direct geometry row count")
    check(owner_values == {0}, "copy-zero owner parity is constant")

    # Closed interval cover and exact first-hit partition.
    first_hits: dict[int, int] = {}
    hit_counts: Counter[int] = Counter()
    for z in range(N):
        hits = tuple(time for time in range(17) if marked(z - time * S))
        check(bool(hits), "every phase marked through time 16")
        first_hits[z] = hits[0]
        hit_counts[hits[0]] += 1
        check(all(not marked(z - time * S) for time in range(hits[0])), "first marked time minimal")

    expected_hits = {0: A + 1, 16: 680}
    expected_hits.update({time: S for time in range(1, 16)})
    check(dict(sorted(hit_counts.items())) == expected_hits, "exact first-hit counts")
    check(tuple(z for z in range(N) if first_hits[z] == 0) == tuple(range(0, A + 1)), "time-zero closed phase interval")
    for time in range(1, 16):
        lower = A + (time - 1) * S + 1
        upper = A + time * S
        check(tuple(z for z in range(N) if first_hits[z] == time) == tuple(range(lower, upper + 1)), f"first-hit slab t={time}")
    check(tuple(z for z in range(N) if first_hits[z] == 16) == tuple(range(A + 15 * S + 1, N)), "time-16 terminal slab")

    running_cover: set[int] = set()
    unresolved_counts: dict[int, int] = {}
    unresolved_sets: dict[int, tuple[int, ...]] = {}
    for horizon in range(17):
        running_cover.update((point + horizon * S) % N for point in range(A + 1))
        if horizon <= 15:
            expected_cover = set(range(0, A + horizon * S + 1))
            check(running_cover == expected_cover, f"closed lifted interval H={horizon}")
        else:
            check(running_cover == set(range(N)), "closed intervals cover circle at H=16")
        unresolved = tuple(sorted(set(range(N)) - running_cover))
        unresolved_sets[horizon] = unresolved
        unresolved_counts[horizon] = len(unresolved)
        expected_count = max(0, N - (A + 1) - horizon * S)
        check(len(unresolved) == expected_count, f"unresolved phase count H={horizon}")
        if horizon <= 15:
            check(unresolved == tuple(range(A + horizon * S + 1, N)), f"unique blind interval H={horizon}")

    # Whole-tail word partition: two scales in each phase.  Before its first
    # mark a row is the common baseline word; afterward the first numeric
    # binder decodes its literal starting parameter.
    baseline_record = formula_record(11123)
    check(record_binders(baseline_record) == (FIXED_BINDER,), "fixed-only baseline binder")
    blind_sample_counts = Counter()
    decoded_sample_rows = 0
    for z in range(N):
        hit = first_hits[z]
        for lift in (0, 1):
            parameter = representative(z, lift)
            full_word = output_word(parameter, 16)
            check(decode(full_word) == parameter, "global W16 decoder")
            decoded_sample_rows += 1
            for horizon in range(17):
                prefix = full_word[: horizon + 1]
                if hit > horizon:
                    check(prefix == (baseline_record,) * (horizon + 1), f"blind word H={horizon}")
                    check(decode(prefix) is None, f"blind word has no label H={horizon}")
                    blind_sample_counts[horizon] += 1
                else:
                    check(decode(prefix) == parameter, f"visible word decodes H={horizon}")
    check(decoded_sample_rows == 2 * N, "W16 decoder sample count")
    for horizon in range(17):
        check(blind_sample_counts[horizon] == 2 * unresolved_counts[horizon], f"sample partition size H={horizon}")

    # Sharp H=15 hostile and its exact time-16 split.
    hostile_1, hostile_2 = 11123, 58317
    check(hostile_2 - hostile_1 == M, "hostile is one scale lift")
    check(phase_of(hostile_1) == phase_of(hostile_2) == 22927, "hostile phase")
    hostile_word_1 = output_word(hostile_1, 16)
    hostile_word_2 = output_word(hostile_2, 16)
    check(hostile_word_1[:16] == hostile_word_2[:16], "equal W15 hostile")
    check(hostile_word_1[16] != hostile_word_2[16], "W16 hostile split")
    check(record_binders(hostile_word_1[16]) == ((11155, 308, "R", 0),), "first hostile binder")
    check(record_binders(hostile_word_2[16]) == ((58349, 1611, "R", 0),), "second hostile binder")
    check(record_exit(hostile_word_1[16]) == rational(4313, 156170), "first hostile exit")
    check(record_exit(hostile_word_2[16]) == rational(22555, 816886), "second hostile exit")

    # The two boundary phases must retain both physical atoms and orientation.
    boundaries = []
    for target_rho in (-A, A):
        candidates = tuple(
            parameter
            for parameter in range(TAIL, TAIL + M, 2)
            if rho_of(parameter) == target_rho
        )
        check(len(candidates) == 1, "one boundary representative per odd period")
        parameter = candidates[0]
        record = formula_record(parameter)
        check(record == direct_record(parameter), "boundary direct geometry")
        check(record_exit(record) == BASE_EXIT, "boundary metric baseline")
        check(len(record_binders(record)) == 2, "boundary has two binders")
        check(visible(record, parameter), "boundary varying tooth is visible")
        boundaries.append((parameter, target_rho, address_of(parameter), record_binders(record)))
    expected_boundaries = [
        (43823, -3371, 1210, ((3371, 93, "R", 0), (43823, 1210, "L", 0))),
        (50565, 3371, 1396, ((3371, 93, "R", 0), (50565, 1396, "R", 0))),
    ]
    check(boundaries == expected_boundaries, "exact boundary conventions")

    # Construct one actual metric fibre for every possible (a,g0) scale
    # shape and verify horizon-zero full-record singleton separation.
    inv_s = inverse_mod(S, M)
    check(inv_s == 1485, "independent extended-Euclid inverse")
    structural_types = 0
    structural_points = 0
    structural_pairs = 0
    maximum_metric_fibre = 0
    for a in range(2, 6741, 2):
        cap = 6740 // a
        for g0 in UNITS_14:
            if g0 > cap:
                continue
            scales = tuple(range(g0, cap + 1, 14))
            kappa_residue = inverse_mod(g0, 14)
            b = (inv_s * (A * kappa_residue - a)) % M
            if b == 0:
                b = M
            attempts = 0
            while b < TAIL or gcd(a, b) != 1:
                b += M
                attempts += 1
                check(attempts < 10000, "primitive representative terminates")
            numerator = a + S * b
            check(numerator % A == 0, "fibre equation divisible by A")
            kappa = numerator // A
            check((g0 * kappa) % 14 == 1, "fibre scale congruence")
            check(b % 2 == 1 and gcd(a, b) == 1, "fibre primitive parity")
            records = set()
            metrics = set()
            for g in scales:
                parameter = b * g
                check(rho_of(parameter) == A - a * g, "constructed strict-active residue")
                record = formula_record(parameter)
                expected_address = (g * kappa - 1) // 14
                check((g * kappa - 1) % 14 == 0, "constructed address integrality")
                check(record_binders(record) == ((parameter, expected_address, "R", 0),), "numeric active binder")
                check(record_exit(record) == rational(kappa, 14 * b), "constant metric fibre")
                records.add(record)
                metrics.add(record_exit(record))
            check(len(metrics) == 1, "metric quotient collapses constructed fibre")
            check(len(records) == len(scales), "full record is singleton on metric fibre")
            structural_types += 1
            structural_points += len(scales)
            structural_pairs += len(scales) * (len(scales) - 1) // 2
            maximum_metric_fibre = max(maximum_metric_fibre, len(scales))
    check((structural_types, structural_points, structural_pairs) == (6106, 13427, 281073), "structural census")
    check(maximum_metric_fibre == 241, "sharp metric-fibre multiplicity")

    # A genuinely full role normalization (chain and binder occurrences) keeps
    # THM-4374's metric hostile through H=16 and splits at H=17.
    metric_hostile_1, metric_hostile_2 = 253031, 258645
    role_word_1 = tuple(role_normalized_record(metric_hostile_1 + 2 * time) for time in range(18))
    role_word_2 = tuple(role_normalized_record(metric_hostile_2 + 2 * time) for time in range(18))
    check(role_word_1[:17] == role_word_2[:17], "full role-normalized W16 hostile")
    check(role_word_1[17] != role_word_2[17], "full role-normalized split at H17")
    check(formula_record(metric_hostile_1) != formula_record(metric_hostile_2), "physical labels split metric hostile at H0")

    print("THM-4379 clean-room labeled-binder referee: PASS")
    print(f"checks={CHECKS}")
    print(f"constants=A:{A},S:{S},M:{M},N:{N},tail:{TAIL},component:{K}")
    print(f"direct_integer_geometry_rows={geometry_rows};all_13_distances=yes")
    print("closed_visible_phase=J:[0,3371];J_H:[0,3371+1303H]_for_H<=15;H16:all")
    print("first_hits=" + ",".join(f"{time}:{hit_counts[time]}" for time in range(17)))
    print("unresolved=" + ",".join(f"{horizon}:{unresolved_counts[horizon]}" for horizon in range(17)))
    print("whole_tail_partition=H0..15:one_infinite_blind_fibre_plus_singletons;H16:all_singletons")
    print("sharp_W15_pair=11123,58317;phase=22927;W16_binders=11155@308R,58349@1611R")
    print("boundaries=43823:-3371:1210:L+3371@93R;50565:+3371:1396:R+3371@93R")
    print(f"strict_active_shapes={structural_types};points={structural_points};pairs={structural_pairs};full_record_max_fibre=1")
    print("role_normalized_hostile=253031,258645;equal_through_W16;split_at_W17")
    print("scope=relative_to_THM-4365/4367_fixed_h420_odd_tail;all_safe;LRC(14)_OPEN")


if __name__ == "__main__":
    main()
