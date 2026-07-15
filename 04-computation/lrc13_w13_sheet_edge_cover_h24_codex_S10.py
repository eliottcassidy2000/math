#!/usr/bin/env python3
"""Exact moving-edge obstruction for the forced odd exception w=13.

Consider a two-sheet packet

    A = 2U union {13, y},       |U| = 10,

where no member of U is divisible by 13.  If the odd runner 13 is eligible
throughout the strict quotient loose set G_U, then its contrapositive can be
read on the thirteen lifts

    tau = (s+j)/13,       j in Z/13Z.

Whenever ||s|| > 2/13, runner 13 is ineligible on every lift, so the danger
teeth of U must cover all thirteen sheets.  Away from events u*s in Z, a
runner u kills exactly the two sheets

    u*j = -floor(u*s), -ceil(u*s)  (mod 13).

At an event u*s=k its edge changes, with the common endpoint retained,

    {-k+1,-k}/u  -->  {-k,-k-1}/u  (mod 13).

It is therefore enough, and exact, to initialize immediately to the right of
s=2/13 and sweep the grouped rational events k/u in (2/13,11/13).  An event
itself has the union of its left and right edges, so coverage of every open
chamber is equivalent to the desired closed-danger predicate.

This certificate tests every ten-subset U of [1,24] avoiding 13.  It also
reports the primitive divisor-complete subatlas forced by THM-772.  A
canonical decision digest records the initial missing-sheet mask and the
first failing chamber for every core.

Tournament Analysis / assumption challenge
-------------------------------------------
There is no clean runner-vertex tournament here.  The theorem-bearing object
is a time-ordered, runner-labelled edge cover on the thirteen sheet vertices.
In a covered generic chamber it has twenty incidences, hence exactly seven
degree-excess tokens.  At a simple event one excess token must move from the
departing sheet to the entering sheet.  A runner tournament would discard
sheet degrees, simultaneous events, and this token current.  The cyclic sheet
order 0,...,12 may be used as a tie Hamiltonian path only if the full labelled
incidence/event-word sidecar is retained; without that sidecar, scores, SCCs,
cycles, and Hamiltonian paths do not preserve the LRC predicate.

The script also certifies the exact degree-collision energy

    K = sum_j binom(d_j-1,2),        X_sheet = 8K.

For a covered simple slide from departure sheet a to entry sheet b,
Delta K=d_b-d_a+1.  Two explicit cores have the same labelled initial degree
vector and X_sheet but different first tears, proving that this energy is a
useful flux coordinate rather than a sufficient quotient.

Writing e^0=d^0-1 and C for cumulative entry-minus-departure sheet current,
the exact carrier is

    e=e^0+C,    cover iff C_j >= -e^0_j for every j,
    K=K_0+<e^0,C>+||C||^2/2.

The full height-24 atlas tears much earlier than the end of the ineligible
interval: every initial cover tears by the event at 3/8, and every
divisor-complete initial cover tears by 4/11.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import reduce
from hashlib import sha256
from itertools import combinations
from math import comb, gcd
from pathlib import Path


MODULUS = 13
CORE_SIZE = 10
MAX_U = 24
LEFT_BOUNDARY = F(2, 13)
RIGHT_BOUNDARY = F(11, 13)
ALL_SHEETS = (1 << MODULUS) - 1

# Filled after the first exact run and asserted on every canonical replay.
EXPECTED_DECISION_DIGEST = (
    "01dafb23a9562fc14f21dbe85bc36d6dcb4e93e4b346dde65a74e82fe2addf0c"
)
EXPECTED_EVENT_DIGEST = (
    "b7b9f5930d28c2dd7e34464851fe00af941561252585d8889b8737682de01cca"
)
ENERGY_LIAR_CORES = {
    (1, 2, 3, 4, 5, 6, 7, 8, 20, 23): ("4/23", 1),
    (1, 2, 3, 4, 5, 6, 8, 10, 11, 14): ("2/11", 7),
}
ENERGY_LIAR_DEGREES = (6, 1, 1, 2, 1, 2, 1, 1, 1, 1, 1, 1, 1)
PREFIX_TEAR_TIME = F(3, 8)
DIVISOR_PREFIX_TEAR_TIME = F(4, 11)
LATEST_TEAR_CORES = {
    (1, 2, 3, 5, 6, 7, 8, 9, 10, 11),
    (1, 2, 3, 5, 6, 7, 8, 9, 11, 20),
    (1, 2, 5, 6, 7, 8, 9, 10, 14, 22),
}
LATEST_DIVISOR_TEAR_CORES = {
    (1, 2, 6, 7, 8, 9, 10, 14, 22, 24),
    (2, 6, 8, 10, 14, 16, 18, 20, 22, 24),
}
EXPECTED_INITIAL_ENERGY_HISTOGRAM = {
    0: 14184, 1: 43927, 2: 26518, 3: 11036, 4: 4838, 5: 201,
    6: 922, 7: 193, 9: 1, 10: 28, 11: 2,
}
EXPECTED_PRETEAR_ENERGY_HISTOGRAM = {
    0: 15270, 1: 45098, 2: 26853, 3: 9518, 4: 4254, 5: 204,
    6: 511, 7: 127, 9: 1, 10: 12, 11: 2,
}


def set_gcd(values: tuple[int, ...]) -> int:
    return reduce(gcd, values)


def divisor_complete(values: tuple[int, ...]) -> bool:
    return all(any(value % modulus == 0 for value in values)
               for modulus in range(2, 13))


def edge_at_initial_chamber(speed: int) -> tuple[int, int]:
    """The two killed sheets immediately to the right of s=2/13."""
    assert speed % MODULUS
    inverse = pow(speed % MODULUS, -1, MODULUS)
    floor_value = (2 * speed) // 13
    return (
        (-floor_value * inverse) % MODULUS,
        (-(floor_value + 1) * inverse) % MODULUS,
    )


def event_update(speed: int, integer: int):
    """Return the left and right two-edges at speed*s=integer."""
    inverse = pow(speed % MODULUS, -1, MODULUS)
    old = (
        (-(integer - 1) * inverse) % MODULUS,
        (-integer * inverse) % MODULUS,
    )
    new = (
        (-integer * inverse) % MODULUS,
        (-(integer + 1) * inverse) % MODULUS,
    )
    return old, new


def build_atlas(candidates: tuple[int, ...]):
    initial_edges = {
        speed: edge_at_initial_chamber(speed) for speed in candidates
    }
    grouped: dict[F, list[tuple[int, tuple[int, int], tuple[int, int]]]] = (
        defaultdict(list)
    )
    for speed in candidates:
        for integer in range(speed + 1):
            if 13 * integer <= 2 * speed or 13 * integer >= 11 * speed:
                continue
            old, new = event_update(speed, integer)
            grouped[F(integer, speed)].append((speed, old, new))

    event_groups = tuple(
        (time, tuple(sorted(updates)))
        for time, updates in sorted(grouped.items())
    )
    digest = sha256()
    for speed in candidates:
        digest.update(
            f"I|{speed}|{initial_edges[speed][0]},{initial_edges[speed][1]}\n".encode()
        )
    for time, updates in event_groups:
        for speed, old, new in updates:
            digest.update(
                f"E|{time.numerator}/{time.denominator}|{speed}|"
                f"{old[0]},{old[1]}>{new[0]},{new[1]}\n".encode()
            )
    return initial_edges, event_groups, digest.hexdigest()


def missing_mask(degrees: list[int]) -> int:
    return sum(1 << sheet for sheet, degree in enumerate(degrees) if degree == 0)


def fmt_fraction(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def choose_two(value: int) -> int:
    return value * (value - 1) // 2


def collision_energy(degrees: tuple[int, ...]) -> int:
    """Collision count among degree-excess chips in a covered chamber."""
    assert len(degrees) == MODULUS and min(degrees) >= 1
    return sum(choose_two(degree - 1) for degree in degrees)


def slide_endpoints(
    old: tuple[int, int], new: tuple[int, int]
) -> tuple[int, int]:
    """Return the unique departure and entry sheets of one edge slide."""
    departure = tuple(set(old) - set(new))
    entry = tuple(set(new) - set(old))
    assert len(departure) == len(entry) == 1
    return departure[0], entry[0]


def audit_collision_flux_identity() -> None:
    """Audit every locally possible covered simple-slide degree pair."""
    for departure_degree in range(2, 21):
        for entry_degree in range(1, 21):
            before = (
                choose_two(departure_degree - 1)
                + choose_two(entry_degree - 1)
            )
            after = (
                choose_two(departure_degree - 2)
                + choose_two(entry_degree)
            )
            assert after - before == entry_degree - departure_degree + 1


def main() -> None:
    candidates = tuple(
        value for value in range(1, MAX_U + 1) if value % MODULUS
    )
    assert len(candidates) == 23
    initial_edges, event_groups, event_digest = build_atlas(candidates)

    decision_digest = sha256()
    total = 0
    initial_covers = 0
    survivors = 0
    divisor_rows = 0
    primitive_divisor_rows = 0
    divisor_initial_covers = 0
    primitive_divisor_initial_covers = 0
    divisor_survivors = 0
    primitive_divisor_survivors = 0
    failure_histogram: Counter[str] = Counter()
    initial_missing_cardinalities: Counter[int] = Counter()
    first_failure_missing_cardinalities: Counter[int] = Counter()
    initial_energy_histogram: Counter[int] = Counter()
    pretear_energy_histogram: Counter[int] = Counter()
    energy_liar_records = {}
    latest_tear_time = None
    latest_tear_cores: set[tuple[int, ...]] = set()
    latest_divisor_tear_time = None
    latest_divisor_tear_cores: set[tuple[int, ...]] = set()

    audit_collision_flux_identity()

    for core in combinations(candidates, CORE_SIZE):
        total += 1
        core_mask = sum(1 << (speed - 1) for speed in core)
        is_divisor = divisor_complete(core)
        is_primitive_divisor = is_divisor and set_gcd(core) == 1
        divisor_rows += int(is_divisor)
        primitive_divisor_rows += int(is_primitive_divisor)

        degrees = [0] * MODULUS
        for speed in core:
            for sheet in initial_edges[speed]:
                degrees[sheet] += 1
        initial_degree_vector = tuple(degrees)

        first_missing = missing_mask(degrees)
        initial_missing_cardinalities[first_missing.bit_count()] += 1
        if first_missing:
            failure_label = "initial"
            failure_missing = first_missing
        else:
            initial_covers += 1
            divisor_initial_covers += int(is_divisor)
            primitive_divisor_initial_covers += int(is_primitive_divisor)
            initial_energy_histogram[collision_energy(initial_degree_vector)] += 1
            failure_label = "survive"
            failure_missing = 0
            initial_excess = tuple(degree - 1 for degree in degrees)
            current = [0] * MODULUS

            for time, updates in event_groups:
                changed = False
                before_degrees = tuple(degrees)
                increment = [0] * MODULUS
                for speed, old, new in updates:
                    if not ((core_mask >> (speed - 1)) & 1):
                        continue
                    changed = True
                    departure, entry = slide_endpoints(old, new)
                    increment[departure] -= 1
                    increment[entry] += 1
                    for sheet in old:
                        degrees[sheet] -= 1
                    for sheet in new:
                        degrees[sheet] += 1
                if not changed:
                    continue
                for sheet in range(MODULUS):
                    current[sheet] += increment[sheet]
                    assert degrees[sheet] == initial_degree_vector[sheet] + current[sheet]
                assert sum(current) == sum(increment) == 0
                failure_missing = missing_mask(degrees)
                if failure_missing:
                    failure_label = fmt_fraction(time)
                    pretear_energy_histogram[collision_energy(before_degrees)] += 1
                    for sheet in range(MODULUS):
                        if failure_missing & (1 << sheet):
                            # The singleton cut has overdrawn its initial
                            # excess-chip capacity by exactly one.
                            assert current[sheet] == -initial_degree_vector[sheet]
                            assert initial_excess[sheet] + current[sheet] == -1
                    if latest_tear_time is None or time > latest_tear_time:
                        latest_tear_time = time
                        latest_tear_cores = {core}
                    elif time == latest_tear_time:
                        latest_tear_cores.add(core)
                    if is_divisor:
                        if (latest_divisor_tear_time is None
                                or time > latest_divisor_tear_time):
                            latest_divisor_tear_time = time
                            latest_divisor_tear_cores = {core}
                        elif time == latest_divisor_tear_time:
                            latest_divisor_tear_cores.add(core)
                    break

                # Grouped current/energy identity.  This includes simultaneous
                # events, for which increment can have entries larger than one.
                energy_before = collision_energy(before_degrees)
                energy_after = collision_energy(tuple(degrees))
                local_flux = sum(
                    (before_degrees[sheet] - 1) * increment[sheet]
                    for sheet in range(MODULUS)
                ) + sum(value * value for value in increment) // 2
                assert energy_after - energy_before == local_flux
                integrated_energy = (
                    collision_energy(initial_degree_vector)
                    + sum(initial_excess[sheet] * current[sheet]
                          for sheet in range(MODULUS))
                    + sum(value * value for value in current) // 2
                )
                assert energy_after == integrated_energy

            if failure_label == "survive":
                survivors += 1
                divisor_survivors += int(is_divisor)
                primitive_divisor_survivors += int(is_primitive_divisor)

        failure_histogram[failure_label] += 1
        first_failure_missing_cardinalities[failure_missing.bit_count()] += 1
        decision_digest.update(
            (",".join(map(str, core))
             + f"|{first_missing:04x}|{failure_label}|{failure_missing:04x}\n").encode()
        )
        if core in ENERGY_LIAR_CORES:
            energy_liar_records[core] = (
                initial_degree_vector,
                collision_energy(initial_degree_vector),
                failure_label,
                failure_missing,
            )

    decision_hexdigest = decision_digest.hexdigest()
    script_digest = sha256(Path(__file__).read_bytes()).hexdigest()

    print("w=13 moving sheet-edge cover certificate through quotient height 24")
    print("strict w-ineligibility / closed U-danger endpoint convention")
    print(f"script_sha256={script_digest}")
    print(f"event_atlas_sha256={event_digest}")
    print(f"decision_atlas_sha256={decision_hexdigest}")
    print()
    print("methodology:")
    print("  sheet vertices: j in Z/13Z, tau=(s+j)/13")
    print("  w=13 ineligible interval: 2/13 < s < 11/13")
    print("  generic runner edge: u*j=-floor(us),-ceil(us) mod 13")
    print("  exact grouped events: k/u with 2u < 13k < 11u")
    print("  simultaneous updates are applied before testing the next chamber")
    print("  event endpoints are safe to omit: the closed event edge is old union new")
    print()
    print(
        f"candidate_speeds={len(candidates)} range=[1,{MAX_U}] excluding=13 "
        f"event_groups={len(event_groups)}"
    )
    print(f"all_10_subsets={total} expected={comb(len(candidates), CORE_SIZE)}")
    print(f"initial_chamber_edge_covers={initial_covers}")
    print(f"full_moving_edge_cover_survivors={survivors}")
    print(
        f"prefix_event_groups_through_{fmt_fraction(PREFIX_TEAR_TIME)}="
        f"{sum(time <= PREFIX_TEAR_TIME for time, _ in event_groups)}"
    )
    print(
        f"latest_first_tear={fmt_fraction(latest_tear_time)} "
        f"cores_at_latest_tear={len(latest_tear_cores)}"
    )
    print("latest_tear_cores=" + str(sorted(latest_tear_cores)))
    print(
        f"divisor_complete_rows={divisor_rows} "
        f"initial_covers={divisor_initial_covers} survivors={divisor_survivors}"
    )
    print(
        f"primitive_divisor_complete_rows={primitive_divisor_rows} "
        f"initial_covers={primitive_divisor_initial_covers} "
        f"survivors={primitive_divisor_survivors}"
    )
    print(
        f"divisor_prefix_event_groups_through_"
        f"{fmt_fraction(DIVISOR_PREFIX_TEAR_TIME)}="
        f"{sum(time <= DIVISOR_PREFIX_TEAR_TIME for time, _ in event_groups)}"
    )
    print(
        f"latest_divisor_complete_first_tear="
        f"{fmt_fraction(latest_divisor_tear_time)} "
        f"cores_at_latest_tear={len(latest_divisor_tear_cores)}"
    )
    print(
        "latest_divisor_complete_tear_cores="
        + str(sorted(latest_divisor_tear_cores))
    )
    print(
        "initial_missing_sheet_cardinality_histogram="
        + str(dict(sorted(initial_missing_cardinalities.items())))
    )
    print(
        "first_failure_missing_sheet_cardinality_histogram="
        + str(dict(sorted(first_failure_missing_cardinalities.items())))
    )
    noninitial_failures = [
        (label, count) for label, count in failure_histogram.items()
        if label not in {"initial", "survive"}
    ]
    noninitial_failures.sort(
        key=lambda row: (-row[1], F(*map(int, row[0].split("/"))))
    )
    print(f"initial_failures={failure_histogram['initial']}")
    print(f"distinct_post_event_first_failure_times={len(noninitial_failures)}")
    print("top_post_event_first_failures=" + str(noninitial_failures[:12]))
    print()
    print("Tournament Analysis / assumption challenge:")
    print("  exact object: 13 sheet obligations + 10 labelled edges + event word")
    print("  covered generic chamber: 20 incidences, degree excess 7")
    print("  switch: a simple event transports one excess token departure->entry")
    print("  tie Hamiltonian path: cyclic sheet labels 0,1,...,12 (telemetry only)")
    print("  no clean runner tournament: pairwise orientations discard coverage")
    print("  required sidecar: sheet degrees, runner labels, simultaneous updates")
    print("  collision energy: K=sum binom(d_j-1,2), X_sheet=8K")
    print("  simple-slide flux: Delta K=d_entry-d_departure+1")
    print("  cumulative current: e=e_initial+C; cover iff C_j>=-e_initial_j")
    print("  quadratic cocycle: K=K_initial+<e_initial,C>+||C||^2/2")
    print(
        "  initial_cover_K_histogram="
        + str(dict(sorted(initial_energy_histogram.items())))
    )
    print(
        "  last_covered_chamber_K_histogram="
        + str(dict(sorted(pretear_energy_histogram.items())))
    )
    for core in ENERGY_LIAR_CORES:
        degrees, energy, failure_label, failure_mask = energy_liar_records[core]
        print(
            f"  energy liar {core}: degrees={degrees} X_sheet={8*energy} "
            f"first_tear={failure_label} missing_sheet={failure_mask.bit_length()-1}"
        )
    print()

    assert total == 1_144_066
    assert initial_covers == 101_850
    assert survivors == 0
    assert divisor_rows == 131_183
    assert primitive_divisor_rows == 131_149
    assert divisor_initial_covers == 20_612
    assert primitive_divisor_initial_covers == 20_604
    assert divisor_survivors == 0
    assert primitive_divisor_survivors == 0
    assert latest_tear_time == PREFIX_TEAR_TIME
    assert latest_tear_cores == LATEST_TEAR_CORES
    assert latest_divisor_tear_time == DIVISOR_PREFIX_TEAR_TIME
    assert latest_divisor_tear_cores == LATEST_DIVISOR_TEAR_CORES
    assert dict(initial_energy_histogram) == EXPECTED_INITIAL_ENERGY_HISTOGRAM
    assert dict(pretear_energy_histogram) == EXPECTED_PRETEAR_ENERGY_HISTOGRAM
    assert set(energy_liar_records) == set(ENERGY_LIAR_CORES)
    for core, (failure_label, missing_sheet) in ENERGY_LIAR_CORES.items():
        degrees, energy, actual_label, actual_mask = energy_liar_records[core]
        assert degrees == ENERGY_LIAR_DEGREES
        assert energy == 10 and 8 * energy == 80
        assert actual_label == failure_label
        assert actual_mask == 1 << missing_sheet
    if EXPECTED_EVENT_DIGEST != "TO_BE_FILLED":
        assert event_digest == EXPECTED_EVENT_DIGEST
    if EXPECTED_DECISION_DIGEST != "TO_BE_FILLED":
        assert decision_hexdigest == EXPECTED_DECISION_DIGEST
    print(
        "FINAL: PASS - every static edge-cover liar tears by 3/8 "
        "(divisor-complete by 4/11)"
    )


if __name__ == "__main__":
    main()
