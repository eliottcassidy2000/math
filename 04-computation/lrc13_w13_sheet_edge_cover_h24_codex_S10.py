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

        first_missing = missing_mask(degrees)
        initial_missing_cardinalities[first_missing.bit_count()] += 1
        if first_missing:
            failure_label = "initial"
            failure_missing = first_missing
        else:
            initial_covers += 1
            divisor_initial_covers += int(is_divisor)
            primitive_divisor_initial_covers += int(is_primitive_divisor)
            failure_label = "survive"
            failure_missing = 0

            for time, updates in event_groups:
                changed = False
                for speed, old, new in updates:
                    if not ((core_mask >> (speed - 1)) & 1):
                        continue
                    changed = True
                    for sheet in old:
                        degrees[sheet] -= 1
                    for sheet in new:
                        degrees[sheet] += 1
                if not changed:
                    continue
                failure_missing = missing_mask(degrees)
                if failure_missing:
                    failure_label = fmt_fraction(time)
                    break

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
        f"divisor_complete_rows={divisor_rows} "
        f"initial_covers={divisor_initial_covers} survivors={divisor_survivors}"
    )
    print(
        f"primitive_divisor_complete_rows={primitive_divisor_rows} "
        f"initial_covers={primitive_divisor_initial_covers} "
        f"survivors={primitive_divisor_survivors}"
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
    if EXPECTED_EVENT_DIGEST != "TO_BE_FILLED":
        assert event_digest == EXPECTED_EVENT_DIGEST
    if EXPECTED_DECISION_DIGEST != "TO_BE_FILLED":
        assert decision_hexdigest == EXPECTED_DECISION_DIGEST
    print("FINAL: PASS - every static edge-cover liar tears during the exact event word")


if __name__ == "__main__":
    main()
