#!/usr/bin/env python3
"""Independent exact THM-3910 audit of the fixed-radius auxiliary-centre filter.

No THM-3878 interval code is imported.  Open pair-danger components are
rebuilt directly as unions of rational teeth on the circle.  For a candidate
centre w0 the conditions tested are

    closed B(w0,1/182) subset D_p union D_q,
    ||a w0|| >= 1/13.

Because the target danger union is open, the centre belongs to the *strict*
erosion (left+1/182,right-1/182).  Equality contacts with an erosion boundary
do not count.
"""

from __future__ import annotations

from fractions import Fraction as Q
from hashlib import sha256
import json
import math
import sys


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")

H = Q(1, 14)
DEEP = Q(1, 13)
RHO = Q(1, 182)

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

EXPECTED_REMAINING = (
    (1, 3), (1, 4), (1, 9), (1, 10),
    (2, 3), (2, 9), (3, 7), (3, 8),
    (4, 7), (5, 6), (5, 12), (6, 11),
    (7, 10), (8, 9), (8, 21), (9, 11),
)

CHECKS = 0


def require(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(label)


def raw_danger(w: int) -> list[tuple[Q, Q]]:
    """Open D_w teeth clipped at the chosen [0,1] circle cut."""
    radius = Q(1, 14 * w)
    pieces = [(Q(0), radius), (Q(1) - radius, Q(1))]
    for k in range(1, w):
        centre = Q(k, w)
        pieces.append((centre - radius, centre + radius))
    return pieces


def circle_open_components(frequencies: tuple[int, ...]) -> tuple[tuple[Q, Q], ...]:
    """Open components as lifted intervals; strict contacts are not merged."""
    merged: list[list[Q]] = []
    pieces = sorted(sum((raw_danger(w) for w in frequencies), []))
    for left, right in pieces:
        # Equality is a missing wall point for a union of open intervals.
        if not merged or left >= merged[-1][1]:
            merged.append([left, right])
        elif right > merged[-1][1]:
            merged[-1][1] = right
    if len(merged) >= 2 and merged[0][0] == 0 and merged[-1][1] == 1:
        wrapping = (merged[-1][0], merged[0][1] + 1)
        return (wrapping,) + tuple((left, right) for left, right in merged[1:-1])
    return tuple((left, right) for left, right in merged)


def eroded_cores(components: tuple[tuple[Q, Q], ...]) -> tuple[tuple[Q, Q], ...]:
    return tuple(
        (left + RHO, right - RHO)
        for left, right in components
        if left + RHO < right - RHO
    )


def deep_intersections(
    cores: tuple[tuple[Q, Q], ...], a: int
) -> tuple[tuple[tuple[Q, Q], ...], tuple[Q, ...]]:
    """Positive intersections and boundary-only contacts with ||a w||>=1/13."""
    hits: list[tuple[Q, Q]] = []
    contacts: set[Q] = set()
    for left, right in cores:
        k_lo = math.floor(a * left) - 2
        k_hi = math.ceil(a * right) + 2
        for k in range(k_lo, k_hi + 1):
            safe_left = (Q(k) + DEEP) / a
            safe_right = (Q(k) + 1 - DEEP) / a
            lo = max(left, safe_left)
            hi = min(right, safe_right)
            if lo < hi:
                hits.append((lo, hi))
            elif lo == hi:
                contacts.add(lo)
    return tuple(hits), tuple(sorted(contacts))


def exhaustive_limit(beta: Q) -> int:
    """First a for which the largest erosion beats an a-danger tooth."""
    gamma = beta - 2 * RHO
    require(gamma > 0, f"positive erosion beta={beta}")
    return math.floor((2 * DEEP) / gamma) + 1


def fmt(value: Q) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def main() -> None:
    require(len(SCALE1) == 57 and len(set(SCALE1)) == 57, "57 distinct pairs")
    require(all(math.gcd(p, q) == 1 and p < q for p, q in SCALE1), "primitive ordered pairs")

    closed = []
    remaining = []
    total_a = 0
    strict_contact_rows = []
    limits = []

    for p, q in SCALE1:
        components = circle_open_components((p, q))
        beta = max(right - left for left, right in components)
        cores = eroded_cores(components)
        gamma = max(right - left for left, right in cores)
        require(gamma == beta - 2 * RHO, f"erosion width {(p, q)}")

        limit = exhaustive_limit(beta)
        limits.append(limit)
        total_a += limit
        killers = []
        contacts_by_a = []
        for a in range(1, limit + 1):
            hits, contacts = deep_intersections(cores, a)
            if not hits:
                killers.append(a)
                if contacts:
                    contacts_by_a.append((a, contacts))

        # If a*gamma > 2/13, the largest open erosion core cannot fit in
        # one open component of {||a w||<1/13}; hence it meets the deep set.
        require(limit * gamma > 2 * DEEP, f"finite-a exhaustion {(p, q)}")
        require(deep_intersections(cores, limit)[0], f"limit positive control {(p, q)}")

        record = (p, q, beta, gamma, limit, tuple(killers))
        if killers:
            closed.append(record)
        else:
            remaining.append(record)
        if contacts_by_a:
            strict_contact_rows.append((p, q, tuple(contacts_by_a)))

    require(len(closed) == 41 and len(remaining) == 16, "41/16 split")
    require(tuple((p, q) for p, q, *_ in remaining) == EXPECTED_REMAINING, "remaining list")
    require(all(killers == (p,) for p, q, beta, gamma, limit, killers in closed), "unique killer a=p")
    require(total_a == 398 and min(limits) == 2 and max(limits) == 12, "finite scan size")

    # These are genuine endpoint cases: the a-deep set only touches an
    # erosion boundary.  The closed source arc cannot touch the boundary of
    # the open target, so these do not count as surviving centres.
    require(
        tuple((p, q, data[0][0], data[0][1]) for p, q, data in strict_contact_rows)
        == ((4, 13, 4, (Q(3, 13), Q(10, 13))),
            (7, 13, 7, (Q(2, 13), Q(11, 13)))),
        "strict erosion contact controls",
    )

    scale2_components = ((Q(2, 21), Q(8, 63)), (Q(55, 63), Q(19, 21)))
    scale2_beta = Q(2, 63)
    scale2_cores = eroded_cores(scale2_components)
    scale2_gamma = scale2_beta - 2 * RHO
    require(scale2_gamma == Q(17, 819), "scale-two erosion")
    scale2_limit = exhaustive_limit(scale2_beta)
    scale2_killers = tuple(
        a for a in range(1, scale2_limit + 1)
        if not deep_intersections(scale2_cores, a)[0]
    )
    require(scale2_limit == 8 and not scale2_killers, "scale-two survives")
    require(scale2_limit * scale2_gamma > 2 * DEEP, "scale-two exhaustion")

    closed_packet = tuple((p, q, killers) for p, q, _, _, _, killers in closed)
    remaining_packet = tuple((p, q) for p, q, *_ in remaining)
    semantic = {
        "scope": "THM3878 t>=U fixed-radius auxiliary-centre necessary condition",
        "radius": "1/182",
        "scale1_closed": [[p, q, list(k)] for p, q, k in closed_packet],
        "scale1_remaining": [list(pair) for pair in remaining_packet],
        "strict_contact_rows": [
            [p, q, [[a, [fmt(x) for x in contacts]] for a, contacts in data]]
            for p, q, data in strict_contact_rows
        ],
        "scale2_killers": list(scale2_killers),
        "finite_a_tests": total_a,
    }
    digest = sha256(
        json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()

    print("LRC14_AUXILIARY_CENTER_FIXED_RADIUS_INDEPENDENT_AUDIT_20260823")
    print("scope=THM3878_11+2;t>=U;necessary_filter_not_counterexample;LRC14=OPEN")
    print("source_arc=closed_radius_1/182;target_pair_danger=open;erosion=strict")
    print(f"scale1_input=57;center_depth_closed={len(closed)};remaining={len(remaining)}")
    print("scale1_closed=" + repr(closed_packet))
    print("scale1_remaining=" + repr(remaining_packet))
    print("killer_law=each_closed_pair_has_unique_killer_a=p")
    print("strict_boundary_contacts=(4,13,a4):3/13,10/13;(7,13,a7):2/13,11/13;contacts_excluded")
    print(f"finite_a_tests={total_a};limits={min(limits)}..{max(limits)}")
    print("large_a_exhaustion=a*(beta-1/91)>2/13_forces_a_deep_center")
    print(f"scale2_beta={fmt(scale2_beta)};erosion={fmt(scale2_gamma)};limit={scale2_limit};killers={scale2_killers}")
    print("semantic_sha256=" + digest)
    print(f"checks={CHECKS}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
