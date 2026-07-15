#!/usr/bin/env python3
"""Uniform-in-c height-24 moving-edge obstruction for w=13c.

Let U be a ten-subset of [1,24] with no multiple of 13.  Persistent
two-sheet ownership by the odd exception w=13c would force the ten quotient
runners to cover all thirteen lifts whenever w is ineligible.  In the
coordinate

    tau_j=(s+j)/13,       j in Z/13Z,

the ineligible set is the union of the c windows

    I_(c,m)=((13m+2)/(13c),(13m+11)/(13c)),  0<=m<c.

Immediately right of a generic point s, runner u supplies the Cayley edge

    E_u(s)={-floor(us)u^-1, -(floor(us)+1)u^-1} mod 13.

At an event s=k/u its root-current increment is

    delta_(-(k+1)u^-1)-delta_(-(k-1)u^-1).

Thus the realizable root word in I_(c,m) consists exactly of the grouped
fractions k/u satisfying

    u(13m+2) < 13ck < u(13m+11).

For g=gcd(u,c), D=u/g, and C=c/g, all c windows together select the phase
coset residues k mod D for which ||Ck/D||>2/13, each with multiplicity g.
The number of selected u-events is therefore

    g*(D-2*floor(2D/13)-1).

It is enough to disprove coverage on the first window I_(c,0).  This exact
certificate sweeps every one of the C(23,10)=1,144,066 cores for every odd
scale c=1,3,...,19.  There are no survivors.  Two static observations close
all larger scales (and are audited here):

  * if 13c>2B, B=max(U), then immediately right of 2/(13c) every edge is
    {0,-u^-1}; their union has at most eleven sheets;
  * if 13c>=11B, the same description holds throughout the whole first
    window.

In particular B<=24 makes every odd c>=5 fail in the initial chamber; only
c=1 and c=3 require event sweeps.  The c=1 slice reproduces THM-792, while
c=3 has 2,528 initial covers and all tear by 1/7.

Root-current / tropical transfer audit
--------------------------------------
For an initially covered chamber put e0=d0-1.  If C_r is cumulative root
current after grouped prefix r and b_j=min_r C_(r,j), then

    d_r-1=e0+C_r,
    the whole prefix is covered iff e0_j+b_j>=0 for all j.

The block transfer T(W)=(C_final,b) composes tropically:

    T(UV)=(C_U+C_V, min(b_U,C_U+b_V)).

The scalar D(W)=sum_j max(0,-b_j) is only the minimum number of freely
placeable chips needed by the word.  The certificate records D at first tear;
most tears have D<=7 and therefore require the actual boundary allocation e0.

Tournament Analysis / assumption challenge
-------------------------------------------
The predicate-preserving vertices are the thirteen sheet cuts, not runners.
For telemetry, orient the pair (j,k) by the lexicographic deficit rank
(d_j,j)<(d_k,k), with time reversal as the current switch/gauge and sheet
order 0,...,12 as the tie Hamiltonian path.  This tournament is transitive:
score histogram 0,...,12, zero directed cycles, thirteen singleton SCCs, and
one Hamiltonian path.  Edge flips between the initial and tear tournaments
are recorded, but they do not determine failure.  Runner tournaments,
unlabelled Cayley edges, the scalar energy, and D(W) all discard actual cut
capacity.  The exact object is (window phase coset, e0, grouped A_12 root
word), equivalently (e0,T(prefix)) at a fixed cursor.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import comb, gcd
from pathlib import Path


MODULUS = 13
CORE_SIZE = 10
MAX_U = 24
SCALES = tuple(range(1, 20, 2))
ALL_SHEETS = (1 << MODULUS) - 1

# Filled after the first exact run and asserted on every canonical replay.
EXPECTED_EVENT_DIGEST = (
    "d51dbb090ee095368f1b240ade9c23168532d6f2db41479348f052c263f465d9"
)
EXPECTED_DECISION_DIGEST = (
    "39449f2c4d98e920ef5768b245fa66f84ae8f3278e6afd5fb9eb99611b599b49"
)


def fmt_fraction(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def missing_mask(degrees: list[int] | tuple[int, ...]) -> int:
    return sum(
        1 << sheet for sheet, degree in enumerate(degrees) if degree == 0
    )


def edge_right_of(speed: int, boundary: F) -> tuple[int, int]:
    """Two-edge in the open chamber immediately right of boundary."""
    inverse = pow(speed % MODULUS, -1, MODULUS)
    floor_value = (speed * boundary.numerator) // boundary.denominator
    return (
        (-floor_value * inverse) % MODULUS,
        (-(floor_value + 1) * inverse) % MODULUS,
    )


def event_update(
    speed: int, integer: int
) -> tuple[tuple[int, int], tuple[int, int], int, int]:
    """Left/right edges and the departure/entry root coordinates."""
    inverse = pow(speed % MODULUS, -1, MODULUS)
    old = (
        (-(integer - 1) * inverse) % MODULUS,
        (-integer * inverse) % MODULUS,
    )
    new = (
        (-integer * inverse) % MODULUS,
        (-(integer + 1) * inverse) % MODULUS,
    )
    departure = (-(integer - 1) * inverse) % MODULUS
    entry = (-(integer + 1) * inverse) % MODULUS
    assert set(old) - set(new) == {departure}
    assert set(new) - set(old) == {entry}
    return old, new, departure, entry


def build_first_window_atlas(
    scale: int, candidates: tuple[int, ...]
) -> tuple[
    F,
    F,
    dict[int, tuple[int, int]],
    tuple[
        tuple[
            F,
            tuple[
                tuple[int, tuple[int, int], tuple[int, int], int, int], ...
            ],
        ],
        ...,
    ],
]:
    left = F(2, 13 * scale)
    right = F(11, 13 * scale)
    initial_edges = {
        speed: edge_right_of(speed, left) for speed in candidates
    }
    grouped: dict[
        F, list[tuple[int, tuple[int, int], tuple[int, int], int, int]]
    ] = defaultdict(list)
    for speed in candidates:
        for integer in range(speed + 1):
            time = F(integer, speed)
            if left < time < right:
                old, new, departure, entry = event_update(speed, integer)
                grouped[time].append(
                    (speed, old, new, departure, entry)
                )
    event_groups = tuple(
        (time, tuple(sorted(updates)))
        for time, updates in sorted(grouped.items())
    )
    return left, right, initial_edges, event_groups


def divisor_complete(core: tuple[int, ...]) -> bool:
    return all(
        any(speed % modulus == 0 for speed in core)
        for modulus in range(2, 13)
    )


def selected_event_count_formula(speed: int, scale: int) -> int:
    common = gcd(speed, scale)
    reduced = speed // common
    return common * (reduced - 2 * ((2 * reduced) // 13) - 1)


def audit_phase_coset_count(candidates: tuple[int, ...]) -> None:
    """Audit the all-window arithmetic root-event count formula."""
    for scale in SCALES:
        for speed in candidates:
            actual = sum(
                13 * min(
                    (scale * integer) % speed,
                    (-(scale * integer)) % speed,
                )
                > 2 * speed
                for integer in range(speed)
            )
            assert actual == selected_event_count_formula(speed, scale)


def audit_static_lemmas() -> None:
    """Finite symbolic-range audit of the two elementary inequalities."""
    for bound in range(1, 65):
        speeds = tuple(
            speed for speed in range(1, bound + 1) if speed % MODULUS
        )
        for scale in range(1, 65):
            left = F(2, 13 * scale)
            right = F(11, 13 * scale)
            if 13 * scale > 2 * bound:
                for speed in speeds:
                    assert edge_right_of(speed, left)[0] == 0
            if 13 * scale >= 11 * bound:
                for speed in speeds:
                    assert speed * right <= 1
                    assert edge_right_of(speed, left)[0] == 0


def tournament_edge_flips(
    before: tuple[int, ...], after: tuple[int, ...]
) -> int:
    """Flips in the transitive lexicographic deficit tournament."""
    flips = 0
    for left in range(MODULUS):
        for right in range(left + 1, MODULUS):
            old_orientation = (before[left], left) < (before[right], right)
            new_orientation = (after[left], left) < (after[right], right)
            flips += old_orientation != new_orientation
    return flips


def main() -> None:
    candidates = tuple(
        speed for speed in range(1, MAX_U + 1)
        if speed % MODULUS
    )
    assert len(candidates) == 23
    assert len(SCALES) == 10
    audit_phase_coset_count(candidates)
    audit_static_lemmas()

    event_digest = sha256()
    decision_digest = sha256()
    scale_records = {}

    for scale in SCALES:
        left, right, initial_edges, event_groups = build_first_window_atlas(
            scale, candidates
        )
        event_digest.update(
            f"W|{scale}|{fmt_fraction(left)}|{fmt_fraction(right)}\n".encode()
        )
        for speed in candidates:
            edge = initial_edges[speed]
            event_digest.update(
                f"I|{scale}|{speed}|{edge[0]},{edge[1]}\n".encode()
            )
        for time, updates in event_groups:
            for speed, old, new, departure, entry in updates:
                event_digest.update(
                    f"E|{scale}|{fmt_fraction(time)}|{speed}|"
                    f"{old[0]},{old[1]}>{new[0]},{new[1]}|"
                    f"{departure}>{entry}\n".encode()
                )

        total = 0
        initial_covers = 0
        divisor_initial_covers = 0
        survivors = 0
        divisor_survivors = 0
        latest_tear = None
        latest_tear_cores: set[tuple[int, ...]] = set()
        first_failure_times: Counter[str] = Counter()
        tropical_demand: Counter[int] = Counter()
        tournament_flips: Counter[int] = Counter()

        for core in combinations(candidates, CORE_SIZE):
            total += 1
            core_mask = sum(1 << (speed - 1) for speed in core)
            is_divisor = divisor_complete(core)
            degrees = [0] * MODULUS
            for speed in core:
                for sheet in initial_edges[speed]:
                    degrees[sheet] += 1
            initial_degrees = tuple(degrees)
            initial_missing = missing_mask(degrees)

            if initial_missing:
                failure_label = "initial"
                failure_missing = initial_missing
            else:
                initial_covers += 1
                divisor_initial_covers += is_divisor
                initial_excess = tuple(degree - 1 for degree in degrees)
                current = [0] * MODULUS
                prefix_minimum = [0] * MODULUS
                failure_label = "survive"
                failure_missing = 0

                for time, updates in event_groups:
                    changed = False
                    increment = [0] * MODULUS
                    for speed, old, new, departure, entry in updates:
                        if not (core_mask & (1 << (speed - 1))):
                            continue
                        changed = True
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
                        prefix_minimum[sheet] = min(
                            prefix_minimum[sheet], current[sheet]
                        )
                        assert (
                            degrees[sheet] - 1
                            == initial_excess[sheet] + current[sheet]
                        )
                    assert sum(current) == 0
                    failure_missing = missing_mask(degrees)
                    if failure_missing:
                        failure_label = fmt_fraction(time)
                        demand = sum(
                            max(0, -value) for value in prefix_minimum
                        )
                        tropical_demand[demand] += 1
                        tournament_flips[
                            tournament_edge_flips(
                                initial_degrees, tuple(degrees)
                            )
                        ] += 1
                        for sheet in range(MODULUS):
                            if failure_missing & (1 << sheet):
                                assert (
                                    initial_excess[sheet]
                                    + prefix_minimum[sheet] <= -1
                                )
                                assert degrees[sheet] == 0
                        if latest_tear is None or time > latest_tear:
                            latest_tear = time
                            latest_tear_cores = {core}
                        elif time == latest_tear:
                            latest_tear_cores.add(core)
                        break

                if failure_label == "survive":
                    survivors += 1
                    divisor_survivors += is_divisor

            first_failure_times[failure_label] += 1
            decision_digest.update(
                (
                    f"{scale}|{','.join(map(str, core))}|"
                    f"{initial_missing:04x}|{failure_label}|"
                    f"{failure_missing:04x}\n"
                ).encode()
            )

        assert total == comb(len(candidates), CORE_SIZE)
        assert survivors == divisor_survivors == 0
        if scale >= 5:
            assert 13 * scale > 2 * MAX_U
            assert initial_covers == 0
            assert first_failure_times == Counter({"initial": total})

        scale_records[scale] = {
            "left": left,
            "right": right,
            "event_groups": len(event_groups),
            "total": total,
            "initial_covers": initial_covers,
            "divisor_initial_covers": divisor_initial_covers,
            "survivors": survivors,
            "divisor_survivors": divisor_survivors,
            "latest_tear": latest_tear,
            "latest_tear_cores": latest_tear_cores,
            "first_failure_times": first_failure_times,
            "tropical_demand": tropical_demand,
            "tournament_flips": tournament_flips,
        }

    event_hexdigest = event_digest.hexdigest()
    decision_hexdigest = decision_digest.hexdigest()
    script_digest = sha256(Path(__file__).read_bytes()).hexdigest()

    print("uniform w=13c moving sheet-edge cover certificate, quotient height 24")
    print("strict w-ineligibility / closed U-danger endpoint convention")
    print(f"script_sha256={script_digest}")
    print(f"event_atlas_sha256={event_hexdigest}")
    print(f"decision_atlas_sha256={decision_hexdigest}")
    print()
    print("methodology:")
    print("  sheet vertices: j in Z/13Z, tau=(s+j)/13")
    print("  scale-c first window: 2/(13c) < s < 11/(13c)")
    print("  generic edge: {-floor(us)/u, -ceil(us)/u} mod 13")
    print("  root event k/u: delta_(-(k+1)/u)-delta_(-(k-1)/u)")
    print("  grouped event criterion: 2u < 13ck < 11u")
    print("  all-window event count: g*(D-2*floor(2D/13)-1)")
    print("  transfer: T(UV)=(C_U+C_V,min(b_U,C_U+b_V))")
    print("  safety: e0+b>=0 coordinatewise; D=sum max(0,-b) is not enough")
    print()
    print("static lemmas:")
    print("  13c>2B: immediately-right edges all {0,-u^-1}, union<=11")
    print("  13c>=11B: those edges are constant on the whole first window")
    print("  B<=24 and odd c>=5: initial chamber cannot cover")
    print()
    print(
        f"candidate_speeds={len(candidates)} range=[1,{MAX_U}] excluding=13 "
        f"cores_per_scale={comb(len(candidates), CORE_SIZE)}"
    )
    for scale in SCALES:
        record = scale_records[scale]
        latest = (
            "none" if record["latest_tear"] is None
            else fmt_fraction(record["latest_tear"])
        )
        print(
            f"c={scale} window=({fmt_fraction(record['left'])},"
            f"{fmt_fraction(record['right'])}) "
            f"event_groups={record['event_groups']} "
            f"initial_covers={record['initial_covers']} "
            f"divisor_initial_covers={record['divisor_initial_covers']} "
            f"survivors={record['survivors']} "
            f"latest_tear={latest} "
            f"latest_tear_cores={len(record['latest_tear_cores'])}"
        )
        if record["tropical_demand"]:
            print(
                "  tropical_D_at_first_tear="
                + str(dict(sorted(record["tropical_demand"].items())))
            )
            print(
                "  deficit_tournament_edge_flips="
                + str(dict(sorted(record["tournament_flips"].items())))
            )
            print(
                "  latest_tear_core_list="
                + str(sorted(record["latest_tear_cores"]))
            )
    print()
    print("Tournament Analysis / assumption challenge:")
    print("  exact vertices: thirteen sheet cuts, not ten runners")
    print("  pair observable: lexicographic sheet deficit (d_j,j) vs (d_k,k)")
    print("  switch/gauge: reverse root current under time reversal")
    print("  tie Hamiltonian path: sheet order 0,1,...,12")
    print("  fingerprint: transitive, scores 0..12, cycles=0, SCCs=13, HP=1")
    print("  edge flips are telemetry; e0 plus root-prefix minima decide coverage")
    print("  challenged vertices: runners, sheets, cuts, events, windows, obligations")
    print()

    expected_initial = {1: 101_850, 3: 2_528}
    expected_latest = {1: F(3, 8), 3: F(1, 7)}
    expected_demand = {
        1: {1: 29_674, 2: 20_776, 3: 23_430, 4: 15_153, 5: 7_783,
            6: 3_530, 7: 1_210, 8: 279, 9: 14, 10: 1},
        3: {1: 17, 2: 124, 3: 378, 4: 643, 5: 678,
            6: 459, 7: 194, 8: 34, 9: 1},
    }
    for scale in (1, 3):
        assert scale_records[scale]["initial_covers"] == expected_initial[scale]
        assert scale_records[scale]["latest_tear"] == expected_latest[scale]
        assert dict(scale_records[scale]["tropical_demand"]) == expected_demand[scale]
    if EXPECTED_EVENT_DIGEST != "TO_BE_FILLED":
        assert event_hexdigest == EXPECTED_EVENT_DIGEST
    if EXPECTED_DECISION_DIGEST != "TO_BE_FILLED":
        assert decision_hexdigest == EXPECTED_DECISION_DIGEST
    print(
        "FINAL: PASS - every height-24 core tears for every odd c "
        "(event sweeps c=1,3; static initial obstruction c>=5)"
    )


if __name__ == "__main__":
    main()
