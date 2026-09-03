#!/usr/bin/env python3
"""Exact parameter-family extension of the THM-4363 LRC14 fibre.

This scratch audit is deliberately self-contained.  It rebuilds the h=420
anchor components and open danger teeth, proves a phase-free cofinite bound
for the completed-chain quotient, sharpens the eventual odd threshold by a
finite exact window, and derives the centered-residue law for the first exit.

Nothing here is an LRC(14) decrement.  Every row in the eventual family has
a literal safe first exit.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from hashlib import sha256
from math import ceil, floor, gcd
import sys


CHECKS = 0


def require(condition: bool, message: str) -> None:
    """Assertion which remains live under python -O."""
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(message)


H = 420
ANCHOR = 840
U = 3
MULTIPLE = 39
FIXED_NINE = (11, 1691, 3371, 5051, 6731, 8411, 10091, 525, 945)
BASE = (U, MULTIPLE) + FIXED_NINE
FOUR_CONTROLS = (241, 255, 761, 1015)
ROLE = {
    3: "U",
    39: "M",
    11: "D0",
    1691: "D1",
    3371: "D2",
    5051: "D3",
    6731: "D4",
    8411: "D5",
    10091: "D6",
    525: "C0",
    945: "C1",
}
WALL_RESIDUES = (1, 13, 15, 27, 29, 41)
COLLARS = tuple(
    sorted(
        {
            k
            for r in WALL_RESIDUES
            for k in ((20 * r - 1) % (2 * H), (20 * r) % (2 * H))
        }
    )
)
RESIDUAL = tuple(k for k in range(2 * H) if k not in set(COLLARS))


@dataclass(frozen=True)
class Tooth:
    speed: int
    address: int
    left: F
    right: F


@dataclass(frozen=True)
class Trace:
    status: str
    chain: tuple[Tooth, ...]
    frontiers: tuple[F, ...]
    exit_or_end: F


def tooth(speed: int, address: int) -> Tooth:
    return Tooth(
        speed,
        address,
        F(14 * address - 1, 14 * speed),
        F(14 * address + 1, 14 * speed),
    )


def component(k: int) -> tuple[F, F]:
    require(0 <= k < 2 * H, "component index outside anchor period")
    return F(14 * k + 1, 28 * H), F(14 * k + 13, 28 * H)


def meeting_teeth(speed: int, left: F, right: F) -> tuple[Tooth, ...]:
    lo = floor(speed * left - F(1, 14))
    hi = ceil(speed * right + F(1, 14))
    answer = []
    for address in range(lo - 1, hi + 2):
        item = tooth(speed, address)
        if item.left < right and left < item.right:
            answer.append(item)
    return tuple(answer)


def selection_key(item: Tooth) -> tuple[F, F, int, int]:
    return item.right, -item.left, -item.speed, -item.address


FIXED_TEETH = tuple(
    tuple(
        item
        for speed in BASE
        for item in meeting_teeth(speed, *component(k))
    )
    for k in range(2 * H)
)


def variable_tooth(parameter: int, frontier: F) -> Tooth | None:
    """The unique P-tooth active at a frontier, if one exists."""
    numerator = parameter * frontier.numerator
    denominator = frontier.denominator
    q = numerator // denominator
    active = []
    for address in (q, q + 1):
        if 14 * abs(numerator - address * denominator) < denominator:
            active.append(tooth(parameter, address))
    require(len(active) <= 1, "two teeth of one speed are simultaneously active")
    return active[0] if active else None


def trace(parameter: int | None, k: int) -> Trace:
    left, right = component(k)
    fixed = FIXED_TEETH[k]
    frontier = left
    chain = []
    frontiers = []
    while True:
        active = [item for item in fixed if item.left < frontier < item.right]
        if parameter is not None:
            extra = variable_tooth(parameter, frontier)
            if extra is not None:
                active.append(extra)
        if not active:
            return Trace("missing", tuple(chain), tuple(frontiers), frontier)
        winner = max(active, key=selection_key)
        frontiers.append(frontier)
        chain.append(winner)
        frontier = winner.right
        if frontier > right:
            status = "span" if len(chain) == 1 else "renew"
            return Trace(status, tuple(chain), tuple(frontiers), frontier)


BASE_TRACES = {k: trace(None, k) for k in RESIDUAL}
COMPLETED_K = tuple(k for k in RESIDUAL if BASE_TRACES[k].status != "missing")
MISSING_K = tuple(k for k in RESIDUAL if BASE_TRACES[k].status == "missing")


def physical_chain(row: Trace) -> tuple[tuple[int, int], ...]:
    return tuple((item.speed, item.address) for item in row.chain)


def same_quotient(parameter: int) -> bool:
    """THM-4363's status + completed-physical-chain quotient."""
    # Completed components reject most non-fibre parameters early.
    for k in COMPLETED_K:
        candidate = trace(parameter, k)
        reference = BASE_TRACES[k]
        if (
            candidate.status != reference.status
            or physical_chain(candidate) != physical_chain(reference)
        ):
            return False
    # A missing partial chain is deliberately erased; only its status matters.
    for k in MISSING_K:
        if trace(parameter, k).status != "missing":
            return False
    return True


def fixed_free_gaps(k: int) -> tuple[tuple[F, F, F], ...]:
    """Positive gaps (length,left,right) in the fixed-tooth union on I_k."""
    left, right = component(k)
    intervals = sorted(
        (max(left, item.left), min(right, item.right))
        for item in FIXED_TEETH[k]
        if item.right > left and item.left < right
    )
    cursor = left
    gaps = []
    for start, end in intervals:
        if end <= cursor:
            continue
        if start > cursor:
            gaps.append((start - cursor, cursor, start))
        cursor = max(cursor, end)
    if cursor < right:
        gaps.append((right - cursor, cursor, right))
    return tuple(gaps)


def distance(x: F) -> F:
    residue = x % 1
    return min(residue, 1 - residue)


def centered(value: int, modulus: int) -> int:
    """Representative in (-modulus/2,modulus/2]."""
    residue = value % modulus
    if 2 * residue > modulus:
        residue -= modulus
    return residue


def first_signature(parameter: int) -> tuple[object, ...]:
    for k in RESIDUAL:
        row = trace(parameter, k)
        if row.status == "missing":
            roles = tuple(
                "P" if item.speed == parameter else ROLE[item.speed]
                for item in row.chain
            )
            addressed = tuple(
                (
                    "P" if item.speed == parameter else ROLE[item.speed],
                    item.address,
                )
                for item in row.chain
            )
            return k, roles, addressed, row.exit_or_end
    raise AssertionError("row has no missing component")


def predicted_first(parameter: int) -> tuple[tuple[tuple[int, int], ...], F, int, int]:
    """Eventual k=23 chain, exit, centered residue and P-address."""
    x1 = F(1303, 47194)
    rho = centered(1303 * parameter, 47194)
    address = (1303 * parameter - rho) // 47194
    chain = ((945, 26), (3371, 93))
    if abs(rho) < 3371:
        chain += ((parameter, address),)
        exit_point = F(14 * address + 1, 14 * parameter)
    else:
        exit_point = x1
    return chain, exit_point, rho, address


def safe_grid_record(parameter: int) -> tuple[int, F]:
    speeds = (ANCHOR,) + BASE + (parameter,)
    safe = []
    for j in range(28 * H):
        t = F(14 * H + j, 28 * H)
        if min(distance(speed * t) for speed in speeds) >= F(1, 14):
            safe.append(t)
    require(bool(safe), f"control P={parameter} lost complete-grid safety")
    return len(safe), safe[0]


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(newline="\n")

    require(
        COLLARS == (19, 20, 259, 260, 299, 300, 539, 540, 579, 580, 819, 820),
        "collar set changed",
    )
    require(len(RESIDUAL) == 828, "residual body is not 828")
    census = Counter(row.status for row in BASE_TRACES.values())
    require(census == Counter({"missing": 546, "span": 276, "renew": 6}), "base census")
    require(len(COMPLETED_K) == 282 and len(MISSING_K) == 546, "base split")
    for parameter in FOUR_CONTROLS:
        require(same_quotient(parameter), f"THM-4363 control P={parameter} left fibre")

    # Every active P-tooth advances a frontier by strictly less than its full
    # width 1/(7P).  The minimum reach on a completed base trace is DELTA.
    completed_reaches = []
    for k in COMPLETED_K:
        row = BASE_TRACES[k]
        for frontier, winner in zip(row.frontiers, row.chain):
            completed_reaches.append((winner.right - frontier, k, frontier, winner))
    delta = min(item[0] for item in completed_reaches)
    delta_hits = tuple(
        (k, frontier, winner.speed, winner.address)
        for reach, k, frontier, winner in completed_reaches
        if reach == delta
    )
    require(delta == F(1, 28665), "minimum completed reach changed")
    require(
        delta_hits
        == (
            (216, F(47, 182), 945, 244),
            (496, F(323, 546), 945, 559),
            (776, F(505, 546), 945, 874),
        ),
        "minimum completed reach locations changed",
    )

    # A fixed-missing component stays missing if one of its fixed-free gaps is
    # at least as long as a whole P-tooth.  Take the maximum gap per component,
    # then the minimum of those hostile gaps over all 546 missing components.
    max_gaps = []
    for k in MISSING_K:
        gaps = fixed_free_gaps(k)
        require(bool(gaps), f"fixed-missing component k={k} has no positive gap")
        largest = max(gaps)
        max_gaps.append((largest[0], k, largest[1], largest[2]))
    gap = min(item[0] for item in max_gaps)
    gap_hits = tuple(item[1:] for item in max_gaps if item[0] == gap)
    require(gap == F(7003, 594127807), "minimum maximum-gap changed")
    require(
        gap_hits
        == (
            (410, F(57597, 117754), F(69103, 141274)),
            (429, F(72171, 141274), F(60157, 117754)),
        ),
        "minimum maximum-gap locations changed",
    )

    phase_free_threshold = max(ceil(F(1, 7) / delta), ceil(F(1, 7) / gap))
    require(phase_free_threshold == 12120, "phase-free quotient threshold changed")
    require(F(1, 7 * phase_free_threshold) <= gap, "threshold does not block a gap")
    require(F(1, 7 * (phase_free_threshold - 1)) > gap, "gap bound not boundary-sharp")
    require(F(1, 7 * phase_free_threshold) < delta, "threshold can alter completed chain")
    require(
        gap - F(1, 7 * phase_free_threshold) == F(137, 1028689860120),
        "phase-free gap margin changed",
    )
    require(
        F(1, 7 * (phase_free_threshold - 1)) - gap
        == F(6044, 7200234893033),
        "phase-free hostile margin changed",
    )

    # The width-gap certificate starts at 12120.  Direct exact enumeration of
    # the last odd window sharpens the actual eventual odd fibre onset.
    for parameter in range(11019, 12120, 2):
        require(same_quotient(parameter), f"odd onset window fails at P={parameter}")
    require(not same_quotient(11017), "hostile P=11017 unexpectedly lies in fibre")
    hostile_diffs = []
    for k in RESIDUAL:
        reference = BASE_TRACES[k]
        candidate = trace(11017, k)
        if reference.status != candidate.status or (
            candidate.status != "missing"
            and physical_chain(reference) != physical_chain(candidate)
        ):
            hostile_diffs.append(
                (k, reference.status, candidate.status, physical_chain(candidate))
            )
    require(
        tuple(hostile_diffs)
        == (
            (
                147,
                "missing",
                "renew",
                ((1691, 296), (11017, 1929), (10091, 1767), (525, 92), (11, 2)),
            ),
            (
                692,
                "missing",
                "renew",
                ((11, 9), (525, 433), (10091, 8324), (11017, 9088), (1691, 1395)),
            ),
        ),
        "P=11017 hostile quotient witness changed",
    )

    # First component 23.  Its base trace ends at the D2 wall x1.  Beyond x1,
    # the nearest fixed left wall is separated by FIRST_GAP.
    base_first = BASE_TRACES[23]
    require(
        physical_chain(base_first) == ((945, 26), (3371, 93)),
        "base first chain changed",
    )
    x0 = F(73, 2646)
    x1 = F(1303, 47194)
    first_gap = F(2110, 158831407)
    require(base_first.frontiers == (F(323, 11760), x0), "base frontiers changed")
    require(base_first.exit_or_end == x1, "base first exit changed")
    require(
        min(
            item.left - x1
            for item in FIXED_TEETH[23]
            if item.left > x1
        )
        == first_gap,
        "next fixed wall gap changed",
    )
    require(
        ceil(F(1, 7) / (x0 - F(323, 11760))) == 1164,
        "first-step stability threshold changed",
    )
    require(
        ceil(F(1, 7) / (x1 - x0)) == 6926,
        "second-step stability threshold changed",
    )

    # For odd P the active centered residues are odd and have |rho|<=3369;
    # hence the largest terminal extension is 6740/(47194 P).
    require(
        F(6740, 47194 * 10751) <= first_gap,
        "uniform odd terminal threshold fails",
    )
    require(
        F(6740, 47194 * 10749) > first_gap,
        "uniform odd terminal hostile bound no longer fails",
    )
    for parameter in range(10141, 10751, 2):
        candidate = trace(parameter, 23)
        expected_chain, expected_exit, _, _ = predicted_first(parameter)
        require(
            physical_chain(candidate) == expected_chain
            and candidate.exit_or_end == expected_exit,
            f"finite first-exit bridge window fails at P={parameter}",
        )
    hostile_first = trace(10139, 23)
    require(
        physical_chain(hostile_first)
        == (
            (945, 26),
            (3371, 93),
            (10139, 280),
            (6731, 186),
            (10091, 279),
        )
        and hostile_first.exit_or_end == F(3907, 141274),
        "P=10139 hostile first-exit chain changed",
    )
    _, predicted_hostile_exit, hostile_rho, _ = predicted_first(10139)
    require(hostile_rho == -3203, "P=10139 centered residue changed")
    require(
        predicted_hostile_exit == F(3921, 141946)
        and predicted_hostile_exit - x1 - first_gap == F(31, 68245609),
        "P=10139 bridge excess changed",
    )

    modulus = 47194
    numerator = 1303
    require(gcd(numerator, modulus) == 1, "first-exit residue map not invertible")
    inverse = pow(numerator, -1, modulus)
    require(inverse == 1485, "residue inverse changed")
    odd_residues = tuple(r for r in range(-23595, 23598, 2))
    require(len(odd_residues) == 23597, "odd centered residue universe changed")
    active_residues = tuple(r for r in odd_residues if abs(r) < 3371)
    require(len(active_residues) == 3370, "active odd residue count changed")
    require(
        len(odd_residues) - len(active_residues) == 20227,
        "inactive odd residue count changed",
    )
    require(
        (inverse * 3371) % modulus == 3371
        and (inverse * (-3371)) % modulus == 43823,
        "open-boundary residue cells changed",
    )

    # An explicit infinite active cell: P_j=1485+47194j, rho=1,
    # n_j=41+1303j.  j>=1 is already in the cofinite quotient fibre.
    family_records = []
    for j in (1, 2, 3):
        parameter = 1485 + modulus * j
        expected_chain, exit_point, rho, address = predicted_first(parameter)
        require(rho == 1 and address == 41 + 1303 * j, "explicit cell law changed")
        require(parameter >= 11019 and same_quotient(parameter), "explicit family left fibre")
        candidate = trace(parameter, 23)
        require(
            physical_chain(candidate) == expected_chain
            and candidate.exit_or_end == exit_point,
            "explicit family first exit changed",
        )
        values = tuple(
            (speed, distance(speed * exit_point))
            for speed in (ANCHOR,) + BASE + (parameter,)
        )
        clearance = min(value for _, value in values)
        binding = tuple(speed for speed, value in values if value == clearance)
        require(clearance == F(1, 14), "explicit exit is not safe")
        require(binding == (parameter,), "explicit active exit binding changed")
        family_records.append((j, parameter, address, exit_point))
    require(
        family_records[0][3] > family_records[1][3] > family_records[2][3] > x1,
        "same-cell metric exits are not strictly distinct",
    )

    # Preserve the original theorem's complete translated-grid firewall.  It
    # is asserted only for the four declared controls, not the whole family.
    grid_records = tuple((parameter, *safe_grid_record(parameter)) for parameter in FOUR_CONTROLS)
    require(
        grid_records
        == (
            (241, 214, F(121, 240)),
            (255, 222, F(121, 240)),
            (761, 218, F(281, 560)),
            (1015, 170, F(281, 560)),
        ),
        "complete translated-grid control changed",
    )

    # Reproduce the earlier bounded observation, now beside the proof that it
    # cannot be interpreted as a finite global signature law.
    prefix_fibre = []
    signatures = Counter()
    for parameter in range(5, 2001, 2):
        if parameter in BASE:
            continue
        if same_quotient(parameter):
            signature = first_signature(parameter)
            prefix_fibre.append(parameter)
            signatures[signature] += 1
    require(len(prefix_fibre) == 365, "bounded fibre count changed")
    require(len(signatures) == 51, "bounded first-signature count changed")
    require(Counter(signatures.values()) == Counter({1: 50, 315: 1}), "multiplicities")

    print("LRC14_828_PARAMETER_FIBRE=COFINITE_ODD_AND_RESIDUE_CLASSIFIED")
    print("scope=exact h=420 fixed bank; LRC14_OPEN; all eventual rows are safe controls")
    print(f"base_census={tuple(sorted(census.items()))};completed={len(COMPLETED_K)}")
    print(f"completed_min_reach={delta};hits={delta_hits}")
    print(f"missing_min_of_max_gaps={gap};hits={gap_hits}")
    print(
        "phase_free_threshold=12120;"
        f"margin_at_12120={gap - F(1, 7 * 12120)};"
        f"failure_at_12119={F(1, 7 * 12119) - gap}"
    )
    print("sharp_eventual_odd_Q_onset=11019;hostile_predecessor=11017")
    print(f"P11017_differences={tuple(hostile_diffs)}")
    print(
        "first_exit_sharp_odd_onset=10141;hostile_predecessor=10139;"
        f"hostile_chain={physical_chain(hostile_first)};hostile_rho={hostile_rho}"
    )
    print(
        "eventual_first_exit: x1=1303/47194; rho=<1303P>_47194; "
        "active iff |rho|<3371; n=(1303P-rho)/47194; "
        "exit=x1+(3371-rho)/(47194P), else x1"
    )
    print(
        f"odd_residue_cells={len(odd_residues)};active={len(active_residues)};"
        f"inactive={len(odd_residues)-len(active_residues)};inverse1303={inverse}"
    )
    print(f"explicit_rho1_family={tuple(family_records)}")
    print(f"complete_grid_controls={grid_records}")
    print(
        "bounded_P_lt_2001:tested=993;fibre=365;signatures=51;"
        "multiplicities=((1,50),(315,1))"
    )
    print(f"bounded_fibre_sha256={sha256(repr(tuple(prefix_fibre)).encode()).hexdigest()}")
    print(
        "bounded_signature_sha256="
        + sha256(repr(tuple(sorted(signatures.items(), key=repr))).encode()).hexdigest()
    )
    print(f"checks={CHECKS}")
    print(
        "CONCLUSION=FINITE_STATE_PERIODIC_ACTIVITY;AFFINE_ADDRESS;"
        "RECIPROCAL_AFFINE_METRIC_EXIT;INFINITE_CONSUMER_RANGE_ON_ONE_Q_FIBRE"
    )
    print("NO_LRC14_DECREMENT;NO_COUNTEREXAMPLE;NO_ARBITRARY_BANK_CLASSIFICATION")


if __name__ == "__main__":
    main()
