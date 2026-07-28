#!/usr/bin/env python3
"""Compact exact referee for the W-relative D^2 factorization defect.

This intentionally hard-codes the first nonconstant 4/17 <-> 13/17
macro-cycle found by the exhaustive scout.  The proof of universal failure of
integer one-step subdivision is algebraic; the script checks its two delayed
midpoints and every physical endpoint factor needed by the displayed witness.
"""

from fractions import Fraction as F
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
HALF_PATH = (ROOT / "04-computation" /
             "lrc14_central_half_odometer_full_local_cycle_thm2698.py")
SPEC = importlib.util.spec_from_file_location("half", HALF_PATH)
half = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(half)

P = 13
R = P ** 6
TARGET = P ** 3
GUARDS = (14, 27, 40, 53, 66, 2 * P ** 5)
ONE_FOURTEENTH = F(1, 14)
CONTROL_LABELS = ("H", "q1", "q2", "q3", "q4", "q5", "c1", "c2", "c3")
CONTROL_SPEEDS = (1, 14, 27, 40, 53, 66, 13, P ** 3, 2 * P ** 5)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def centered_distance(value):
    value = value - (value.numerator // value.denominator)
    return min(value, 1 - value)


def in_raw_w(y):
    return (centered_distance(TARGET * y) < ONE_FOURTEENTH
            and all(centered_distance(q * y) > ONE_FOURTEENTH
                    for q in GUARDS))


def danger_signature(y):
    return tuple(
        label for index, (label, speed) in enumerate(
            zip(CONTROL_LABELS, CONTROL_SPEEDS)
        )
        if centered_distance(speed * y) < (F(1, 7) if index == 0
                                             else ONE_FOURTEENTH)
    )


def raw_w_local_slack(y):
    """Symmetric y-radius before the first displayed W inequality binds."""
    margins = [(ONE_FOURTEENTH - centered_distance(TARGET * y)) / TARGET]
    margins.extend((centered_distance(q * y) - ONE_FOURTEENTH) / q
                   for q in GUARDS)
    require(all(value > 0 for value in margins),
            "raw-W slack requested outside the strict carrier")
    return min(margins)


def b0(y):
    return half.frac(P * y)


def macro(value, lift):
    return half.frac(P * P * value + F(lift, R))


def main():
    # Central-involution control: H^2=1, HT=TH and S=HT.
    y_half = F(11, 24)
    h = lambda y: half.frac(y + F(1, 2))
    s = lambda y: h(b0(y))
    require(h(h(y_half)) == y_half and h(b0(y_half)) == b0(h(y_half)),
            "central involution law changed")
    require(s(y_half) == y_half and s(s(y_half)) == b0(b0(y_half)),
            "S^2=T^2 or the half fixed point changed")
    require(in_raw_w(y_half) and not in_raw_w(b0(y_half)),
            "central W-relative factorization hostile changed")
    source_radius = raw_w_local_slack(y_half)
    endpoint_radius = source_radius / P ** 2
    midpoint_exclusion_radius = (
        centered_distance(TARGET * b0(y_half)) - ONE_FOURTEENTH
    ) / (TARGET * P)
    raw_radius = min(source_radius, endpoint_radius,
                     midpoint_exclusion_radius)
    require(raw_radius == F(1, 10541750856),
            "central raw-W endpoint radius changed")
    reflected = F(13, 24)
    require(s(reflected) == reflected
            and raw_w_local_slack(reflected) / P ** 2 == raw_radius
            and not in_raw_w(b0(reflected)),
            "reflected central defect interval changed")

    module, _, _, _, rails, present, present_starts = (
        half.core.build_carrier_data()
    )
    pair_prefixes = half.private.build_pair_prefixes(module)
    # Only rails j=3 and j=8 occur in the hard-coded witness.  Building the
    # intervening six-row shard is substantially cheaper than replaying the
    # exhaustive fourteen-row source bank.
    source_rows = half.private.shard((3, 9))
    rows_by_j = {j: source_rows[6][j - 3] for j in range(3, 9)}
    rail_starts = [[left for left, _, _ in rail[3]]
                   for rail in rails[:14]]

    first = {
        "x": F(39123022, 82055753), "N": 2301354, "y": F(4, 17),
        "j": 8, "shallow": 1, "owner": 4, "carry": 3,
        "h": 3, "kappa": 0, "edge": 0, "root": 7,
    }
    second = {
        "x": F(41305372, 82055753), "N": 2429727, "y": F(13, 17),
        "j": 3, "shallow": 4, "owner": 1, "carry": 1,
        "h": 9, "kappa": 1, "edge": 0, "root": 4,
    }
    k0, k1 = 4472391, 1956127
    require(macro(first["x"], k0) == second["x"]
            and macro(second["x"], k1) == first["x"],
            "the nonconstant D^2 macro-cycle changed")
    require(first["owner"] == second["shallow"]
            and second["owner"] == first["shallow"],
            "owner-to-next-shallow gluing changed")

    slacks = []
    for row in (first, second):
        digit = half.floor_fraction(2 * P * row["y"])
        require(half.floor_fraction(R * row["x"]) == row["N"]
                and row["N"] % P == row["carry"]
                and half.frac(R * row["x"]) == row["y"]
                and half.shallow(row["x"]) == row["shallow"]
                and half.owner(row["x"]) == row["owner"]
                and half.floor_fraction(P * row["y"]) == row["h"]
                and digit - 2 * row["h"] == row["kappa"],
                "hard-coded endpoint chart changed")
        require(rails[row["j"]][:2] == (1, row["owner"]),
                "source-one rail/owner typing changed")
        slack, _ = half.packet_slack(
            module, rails, rail_starts, present, present_starts,
            pair_prefixes, row,
        )
        values = rows_by_j[row["j"]][0][row["edge"]][row["carry"]][
            row["kappa"]][row["h"]]
        require(half.private.is_unit(values, row["root"], 26),
                "endpoint primitive-unit row changed")
        slacks.append(slack)
    require(tuple(slacks) == (F(11, 853068347561612),) * 2,
            "endpoint packet slack changed")

    roots0 = half.half_roots(module, half.frac(P * P * first["x"]), 0)
    roots1 = half.half_roots(module, half.frac(P * P * second["x"]), 0)
    step0, step1 = 2 * k0 % P, 2 * k1 % P
    require((roots0, roots1, step0, step1) == ((2,), (12,), 2, 8),
            "macro root input or increment changed")
    require((roots0[0] + step0) % P == second["root"]
            and (roots1[0] + step1) % P == first["root"],
            "macro root landing changed")

    endpoints = (first["y"], second["y"])
    midpoints = tuple(b0(y) for y in endpoints)
    require(tuple(b0(b0(y)) for y in endpoints) == endpoints[::-1],
            "13^2 stopped swapping the endpoint phases")
    require(all(in_raw_w(y) for y in endpoints)
            and midpoints == (F(1, 17), F(16, 17))
            and not any(in_raw_w(y) for y in midpoints),
            "W-relative endpoint/midpoint classification changed")
    require(tuple(centered_distance(TARGET * y) for y in midpoints)
            == (F(4, 17), F(4, 17)),
            "forced midpoint target failure changed")
    danger_signatures = tuple(danger_signature(y)
                              for y in (endpoints[0], midpoints[0],
                                        endpoints[1], midpoints[1]))
    require(danger_signatures == (("c1", "c2"), ("H",),
                                  ("c1", "c2"), ("H",)),
            "four-phase event/ghost danger quotient changed")

    physical_radius = min(slacks[0] / P ** 4, slacks[1] / P ** 2)
    require(physical_radius == F(11, 24364485074707200332),
            "three-macro-state physical radius changed")

    print("LRC14 COMPACT D2 MACRO / RELATIVE-SEGAL REFEREE")
    print(f"central_raw_radius={raw_radius}; central_defect_measure_floor={4 * raw_radius}; central_midpoints={(b0(y_half), b0(reflected))}")
    print(f"macro_endpoints={endpoints}; forced_midpoints={midpoints}")
    print(f"c4_danger_signatures={danger_signatures}")
    print(f"macro_lifts={(k0, k1)}; root_data={(roots0, roots1, step0, step1)}")
    print(f"packet_slacks={tuple(slacks)}")
    print(f"three_state_radius={physical_radius}; open_length={2 * physical_radius}")
    print("factorization_rule=M_K=F_b*F_a for 13*a+b=K mod R; every middle delayed phase is B0(y), independent of a,b")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
