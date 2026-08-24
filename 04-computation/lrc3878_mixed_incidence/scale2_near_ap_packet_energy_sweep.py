#!/usr/bin/env python3
"""Scratch-only exact packet-energy sweep for THM-3878's (2,1,9) row.

Universe: every 11-subset E of {1,...,15}, every physical odd t>=max(E),
in the row 2E union {t,9t}.  The finite range is determined separately for
each E by the signed-endpoint tail bound; the exact t-grid energy closes every
cell below the tail.  This is conditional t>=U only and proves no LRC(14).
"""

from __future__ import annotations

from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations
import json
from math import isqrt
import sys


sys.stdout.reconfigure(newline="\n")
DELTA = Q(1, 14)
ALPHA = Q(4, 63)
CHECKS = 0
EXPECTED_SEMANTIC_SHA256 = "59ee93b99a53d20d7dd034a5d15d3d8f77199b328a3b6fd86dea15792fff317f"


def require(ok: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)


def fmt(x: Q) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def merge(pieces: list[tuple[Q, Q]]) -> list[tuple[Q, Q]]:
    out: list[list[Q]] = []
    for a0, b0 in sorted(pieces):
        a, b = Q(a0), Q(b0)
        if a >= b:
            continue
        if not out or a > out[-1][1]:
            out.append([a, b])
        elif b > out[-1][1]:
            out[-1][1] = b
    return [(a, b) for a, b in out]


def bad(speed: int) -> list[tuple[Q, Q]]:
    radius = Q(1, 14 * speed)
    out = []
    for k in range(speed):
        c = Q(k, speed)
        a, b = c - radius, c + radius
        if a < 0:
            out.extend(((Q(0), b), (a + 1, Q(1))))
        elif b > 1:
            out.extend(((a, Q(1)), (Q(0), b - 1)))
        else:
            out.append((a, b))
    return merge(out)


def safe_set(body: tuple[int, ...]) -> list[tuple[Q, Q]]:
    danger = merge(sum((bad(v) for v in body), []))
    out = []
    x = Q(0)
    for a, b in danger:
        if x < a:
            out.append((x, a))
        x = max(x, b)
    if x < 1:
        out.append((x, Q(1)))
    return out


def endpoint_word(pieces: list[tuple[Q, Q]]) -> tuple[tuple[Q, int], ...]:
    weights: dict[Q, int] = {}
    for a, b in pieces:
        aa, bb = a % 1, b % 1
        weights[aa] = weights.get(aa, 0) + 1
        weights[bb] = weights.get(bb, 0) - 1
    return tuple(sorted((x, s) for x, s in weights.items() if s))


def b2(x: Q) -> Q:
    x %= 1
    return x * x - x + Q(1, 6)


def grid_energy(edges: tuple[tuple[Q, int], ...], t: int) -> Q:
    return Q(1, 2 * t * t) * sum(
        se * sf * b2(t * (e - f))
        for e, se in edges for f, sf in edges
    )


def occupancy_energy(pieces: list[tuple[Q, Q]], t: int) -> Q:
    """Independently integrate the exact orbit-occupancy word N_t."""
    events: dict[Q, int] = {}
    initial = 0
    for a, b in pieces:
        aa, bb = t * a, t * b
        na, ra = divmod(aa.numerator, aa.denominator)
        nb, rb = divmod(bb.numerator, bb.denominator)
        left = na if ra == 0 else na + 1
        right = nb - 1 if rb == 0 else nb
        initial += max(0, right - left + 1)
        xa, xb = aa % 1, bb % 1
        if xa:
            events[xa] = events.get(xa, 0) + 1
        if xb:
            events[xb] = events.get(xb, 0) - 1

    n = initial
    previous = Q(0)
    first = second = Q(0)
    for x, jump in sorted(events.items()):
        width = x - previous
        first += width * n
        second += width * n * n
        n += jump
        previous = x
    width = Q(1) - previous
    first += width * n
    second += width * n * n
    mass = measure(pieces)
    require(first == t * mass, f"occupancy first moment t={t}")
    return second / (t * t) - mass * mass


def measure(pieces: list[tuple[Q, Q]]) -> Q:
    return sum((b - a for a, b in pieces), Q(0))


def strict_tail_start(mass: Q, arcs: int) -> int:
    # disc_t<=arcs^2/(3t^2), while containment would require
    # disc_t >= ((1-alpha)/alpha) mass^2.
    wall = Q(arcs * arcs) * ALPHA / (3 * mass * mass * (1 - ALPHA))
    t = isqrt(wall.numerator // wall.denominator)
    while Q(t * t) <= wall:
        t += 1
    return t


def main() -> None:
    # Independent reconstruction of the quotient obstruction recorded in
    # THM-3878: two open arcs, each of length 2/63.
    obstruction = ((Q(2, 21), Q(8, 63)), (Q(55, 63), Q(19, 21)))
    require(sum((b - a for a, b in obstruction), Q(0)) == ALPHA,
            "scale-two obstruction mass")
    require(all(b - a == Q(2, 63) for a, b in obstruction),
            "scale-two component widths")

    bodies = tuple(combinations(range(1, 16), 11))
    require(len(bodies) == 1365, "near-AP body universe")
    cells = []
    endpoint_audit_cells = 0
    max_tail = None
    closest = None
    for body in bodies:
        safe = safe_set(body)
        mass = measure(safe)
        edges = endpoint_word(safe)
        require(mass > 0, f"positive body mass {body}")
        require(len(edges) % 2 == 0, f"even endpoint word {body}")
        arcs = len(edges) // 2
        tail = strict_tail_start(mass, arcs)
        tail_record = (tail, body, arcs, mass)
        if max_tail is None or tail_record > max_tail:
            max_tail = tail_record

        first_t = max(body)
        if first_t % 2 == 0:
            first_t += 1
        # One independent Bernoulli-endpoint check for every body.  The full
        # finite sweep below uses the orbit word, so this is a separate path
        # rather than a second spelling of the certificate statistic.
        first_disc = occupancy_energy(safe, first_t)
        require(first_disc == grid_energy(edges, first_t),
                f"endpoint/occupancy energy agreement {(body,first_t)}")
        endpoint_audit_cells += 1
        for t in range(first_t, tail, 2):
            # If t is already a body speed then 2E union {t,9t} still has
            # distinct actual speeds (2E is even and t is odd); no exclusion
            # is needed.  Oddness is exactly gcd(2,t)=1.
            disc = first_disc if t == first_t else occupancy_energy(safe, t)
            floor = (1 - ALPHA) * mass * mass / ALPHA
            margin = floor - disc
            require(margin > 0, f"packet-energy closure {(body,t)}")
            record = (margin, body, t, mass, arcs, disc, tail)
            if closest is None or record < closest:
                closest = record
            cells.append((body, t, fmt(mass), arcs, fmt(disc), fmt(margin), tail))

    require(len(cells) == 5470, f"finite packet-energy cell count={len(cells)}")
    require(endpoint_audit_cells == len(bodies), "one endpoint audit per body")
    require(max_tail is not None and max_tail[0] == 64, "maximum tail threshold")
    require(closest is not None, "closest finite cell")

    semantic = {
        "scope": "THM3878 (2,1,9), t>=U, E any 11-subset of 1..15",
        "obstruction": [[fmt(a), fmt(b)] for a, b in obstruction],
        "cells": cells,
        "tail_rule": "disc_t<=r^2/(3t^2)<((1-alpha)/alpha)m^2",
    }
    digest = sha256(json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")).hexdigest()
    if digest != EXPECTED_SEMANTIC_SHA256:
        raise RuntimeError("frozen semantic transcript")

    margin, body, t, mass, arcs, disc, tail = closest
    print("THM3878_SCALE2_NEAR_AP_PACKET_ENERGY_SWEEP_20260823")
    print("scope=conditional_t>=U;scale2=(2,1,9);1365_bodies_E_subset_1..15;LRC14=OPEN")
    print("obstruction=(2/21,8/63)U(55/63,19/21);alpha=4/63;components=2")
    print(f"bodies={len(bodies)};finite_exact_cells={len(cells)};exact_energy_failures=0")
    print("tail_certificate=disc_t_le_r_squared_over_3t_squared;body_specific_strict_threshold")
    print("max_tail_start=" + repr((max_tail[0], max_tail[1], max_tail[2], fmt(max_tail[3]))))
    print("closest_finite_margin=" + repr((body, t, fmt(mass), arcs, fmt(disc), fmt(margin), tail)))
    print(f"independent_energy_path=exact_orbit_occupancy_word_N_t;endpoint_B2_controls={endpoint_audit_cells}")
    print("conclusion=every_2E_union_t{1,9}_with_E_subset_1..15_size11_and_odd_t>=maxE_has_positive_safe_measure")
    print("nonconsequence=unrestricted_t_below_U_and_general_eleven_speed_body_remain_open")
    print("semantic_sha256=" + digest)
    print(f"checks={CHECKS}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
