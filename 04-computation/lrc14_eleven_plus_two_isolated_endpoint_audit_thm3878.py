#!/usr/bin/env python3
"""Exact endpoint-only audit for the 58 conditional THM-3878 survivors.

This is deliberately pointwise: finite isolated safe walls have Haar mass zero
and are invisible to every integral moment, but may still witness the non-strict
Lonely Runner inequality.  No canonical search is imported.
"""

from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
import json
import sys


sys.stdout.reconfigure(newline="\n")

SCALE1 = (
    (1, 3), (1, 4), (1, 9), (1, 10),
    (2, 3), (2, 9), (2, 15), (2, 21), (2, 23),
    (3, 7), (3, 8), (3, 14), (3, 17), (3, 19), (3, 20),
    (3, 22), (3, 26), (3, 31), (3, 38),
    (4, 7), (4, 13), (4, 19), (4, 21), (4, 25), (4, 37),
    (4, 43), (4, 49), (4, 51),
    (5, 6), (5, 12), (5, 17), (5, 18), (5, 24), (5, 29),
    (5, 36), (5, 39), (5, 41), (5, 42), (5, 48), (5, 53),
    (5, 54), (5, 63),
    (6, 11), (6, 17), (6, 19), (6, 23), (6, 41), (6, 47),
    (6, 53), (6, 65),
    (7, 10), (7, 13), (7, 15), (7, 22),
    (8, 9), (8, 21), (9, 11),
)

ISOLATED = tuple(Q(r, 14) for r in (3, 5, 9, 11))
DELTA = Q(1, 14)
CHECKS = 0


def require(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(label)


def dist(x: Q) -> Q:
    r = x % 1
    return min(r, 1 - r)


def ap11_safe(y: Q) -> bool:
    return all(dist(v * y) >= DELTA for v in range(1, 12))


def scale1_endpoint_safe(p: int, q: int, t: int, y: Q) -> bool:
    return dist(t * p * y) >= DELTA and dist(t * q * y) >= DELTA


def scale2_safe_lifts(t: int, y: Q) -> tuple[int, ...]:
    """Safe k in x=(y+k)/2 for speeds 2*AP11 plus t*(1,9)."""
    return tuple(
        k for k in (0, 1)
        if dist(t * Q(y + k, 2)) >= DELTA
        and dist(9 * t * Q(y + k, 2)) >= DELTA
    )


def danger_union_metrics(frequencies: tuple[int, ...]) -> tuple[Q, int]:
    """Exact measure and positive circle-component count of open danger."""
    pieces = []
    for w in frequencies:
        radius = Q(1, 14 * w)
        for k in range(w + 1):
            left = max(Q(0), Q(k, w) - radius)
            right = min(Q(1), Q(k, w) + radius)
            if left < right:
                pieces.append((left, right))
    merged = []
    for left, right in sorted(pieces):
        if not merged or left >= merged[-1][1]:
            merged.append([left, right])
        elif right > merged[-1][1]:
            merged[-1][1] = right
    measure = sum((right - left for left, right in merged), Q(0))
    count = len(merged)
    if count >= 2 and merged[0][0] == 0 and merged[-1][1] == 1:
        count -= 1
    return measure, count


def main() -> None:
    require(len(SCALE1) == 57 and len(set(SCALE1)) == 57, "scale-one universe")
    require(all(ap11_safe(y) for y in ISOLATED), "AP11 isolated walls")

    # Each isolated numerator is a unit modulo 14.  Hence all four endpoints
    # have the same scale-one verdict, namely 14 does not divide tp or tq.
    residue_counts = {}
    residue_sets = {}
    for p, q in SCALE1:
        good = tuple(
            t for t in range(14)
            if all(scale1_endpoint_safe(p, q, t, y) for y in ISOLATED)
        )
        predicted = tuple(
            t for t in range(14)
            if (t * p) % 14 != 0 and (t * q) % 14 != 0
        )
        require(good == predicted, f"endpoint residue rule {(p, q)}")
        # If one isolated endpoint works then all four work.
        for t in range(14):
            verdicts = {scale1_endpoint_safe(p, q, t, y) for y in ISOLATED}
            require(len(verdicts) == 1, f"four-wall coherence {(p, q, t)}")
        residue_sets[(p, q)] = good
        residue_counts[(p, q)] = len(good)

    histogram = Counter(residue_counts.values())
    require(histogram == Counter({0: 2, 6: 7, 7: 4, 12: 33, 13: 11}), "residue histogram")
    never = tuple(pair for pair, count in residue_counts.items() if count == 0)
    require(never == ((3, 14), (5, 42)), "never-help pairs")
    require(all(0 not in good for good in residue_sets.values()), "common residue-zero hostile")

    # A sharp liar for scalar measure + positive component count + inversion
    # symmetry.  The two atlas pairs even have the same sum.  At the same
    # legal t=15 they have identical pullback obstruction mass/count, but
    # opposite incidence on every AP11 isolated wall.
    require(danger_union_metrics((3, 14)) == (Q(13, 49), 14), "liar pair one")
    require(danger_union_metrics((7, 10)) == (Q(13, 49), 14), "liar pair two")
    require(danger_union_metrics((45, 210)) == (Q(13, 49), 210), "liar pullback one")
    require(danger_union_metrics((105, 150)) == (Q(13, 49), 210), "liar pullback two")
    require(not any(scale1_endpoint_safe(3, 14, 15, y) for y in ISOLATED), "liar killed walls")
    require(all(scale1_endpoint_safe(7, 10, 15, y) for y in ISOLATED), "liar safe walls")

    # In the scale-two row gcd(2,t)=1, so only odd residues occur.  Directly
    # check both coherent lifts for every endpoint and every odd residue modulo
    # 28 (the lift denominator).  Adding 28 to t changes neither phase.
    lift_packet = {}
    for t in range(1, 28, 2):
        lift_packet[t] = tuple(scale2_safe_lifts(t, y) for y in ISOLATED)
        require(all(lifts for lifts in lift_packet[t]), f"scale-two endpoint lifts t={t}")
    lift_histogram = Counter(len(lifts) for packet in lift_packet.values() for lifts in packet)
    require(lift_histogram == Counter({1: 32, 2: 24}), "scale-two lift histogram")

    semantic = {
        "universe": "AP11 isolated closed-safe walls against 58 conditional rows",
        "scale1_good_residue_sets": {
            f"{p},{q}": list(residue_sets[(p, q)]) for p, q in SCALE1
        },
        "scale2_safe_lifts": {
            str(t): [list(lifts) for lifts in lift_packet[t]] for t in lift_packet
        },
        "scalar_component_liar": {
            "t": 15,
            "killed_pair": [3, 14],
            "witness_pair": [7, 10],
            "pair_sum": 17,
            "measure": "13/49",
            "pullback_positive_components": 210,
        },
        "scope": "endpoint witnesses only; no universal u conclusion",
    }
    digest = sha256(json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")).hexdigest()

    print("THM3878_AP11_ISOLATED_ENDPOINT_WITNESS_AUDIT_20260823")
    print("scope=conditional_t>=U_hostile_control;AP11_only;LRC14=OPEN")
    print("closed_isolated_safe_walls=3/14,5/14,9/14,11/14;haar_mass=0")
    print("scale1_endpoint_rule=all_four_safe_iff_14_not_divide_tp_and_14_not_divide_tq")
    print("scale1_good_residue_count_histogram=0:2,6:7,7:4,12:33,13:11")
    print("scale1_never_help_pairs=(3,14),(5,42)")
    print("scale1_common_hostile=t_congruent_0_mod_14;all_57_rows;take_AP11_U=11,t=14")
    print("scale1_universal_endpoint_closures=0")
    print("endpoint_liar=t15,sum17:(3,14)_killed_vs_(7,10)_all4safe;measure=13/49;pullback_components=210;both_inversion_symmetric")
    print("scale2_pair=(1,9);allowed_t=odd;endpoint_residue_cases_mod28=56")
    print("scale2_safe_lift_histogram=one_lift:32,two_lifts:24;failures:0")
    print("scale2_AP11_endpoint_closure=all_odd_t;special_body_only")
    print("finite_points_change_Haar_moments=0;endpoint_witness_value=pointwise_only")
    print("semantic_sha256=" + digest)
    print(f"checks={CHECKS}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
