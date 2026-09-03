#!/usr/bin/env python3
"""Clean-room exact referee for the two THM-4391 nonunit l1=16 sectors.

This script is intentionally self-contained.  It imports no repository module
and does not read the producer packet.  Every logical check remains active
under ``python -O``.
"""

from fractions import Fraction as Q
from hashlib import sha256
from itertools import product
from math import gcd
import json


R = Q(3, 14)
BULK = Q(6, 343)
CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def qtext(x):
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def odd_three_unit(n):
    return n > 0 and n % 2 == 1 and n % 3 != 0


def gcd3(a, b, c):
    return gcd(gcd(a, b), c)


def admissible_partitions(total):
    """Primitive sorted full-support magnitude patterns of given l1 norm."""
    ans = []
    for a in range(1, total):
        for b in range(a, total):
            c = total - a - b
            if c < b:
                continue
            if c <= 0:
                continue
            if any(x % 3 == 0 for x in (a, b, c)):
                continue
            if gcd3(a, b, c) != 1:
                continue
            ans.append((a, b, c))
    return ans


def clip_halfplane(poly, A, B, C):
    """Clip a rational polygon by A*x+B*y <= C."""
    if not poly:
        return []
    out = []
    for P, S in zip(poly, poly[-1:] + poly[:-1]):
        # Traverse S -> P.
        fs = A * S[0] + B * S[1] - C
        fp = A * P[0] + B * P[1] - C
        sin = fs <= 0
        pin = fp <= 0
        if sin != pin:
            den = fs - fp
            require(den != 0, "crossing edge has zero affine difference")
            lam = fs / den
            X = (S[0] + lam * (P[0] - S[0]), S[1] + lam * (P[1] - S[1]))
            out.append(X)
        if pin:
            out.append(P)
    return out


def polygon_area(poly):
    if len(poly) < 3:
        return Q(0)
    twice = Q(0)
    for (x1, y1), (x2, y2) in zip(poly, poly[1:] + poly[:1]):
        twice += x1 * y2 - x2 * y1
    return abs(twice) / 2


def strip_area(h, m, ell, delta):
    # Reflection x -> -x removes the sign s, so use |delta-h*x+m*y| <= ell*r.
    square = [(-R, -R), (R, -R), (R, R), (-R, R)]
    # delta-h*x+m*y <= ell*r
    poly = clip_halfplane(square, -h, m, ell * R - delta)
    # -(delta-h*x+m*y) <= ell*r
    poly = clip_halfplane(poly, h, -m, ell * R + delta)
    return polygon_area(poly)


def ceil_q(x):
    return -((-x.numerator) // x.denominator)


def lift_length(p, b, q, s, h, m, ell, delta, k):
    """Six-term interval-intersection roof, derived from center distances."""
    terms = [
        2 * R / p,
        2 * R / b,
        2 * R / q,
        R / p + R / b - Q(abs(k), p * b),
        R / p + R / q - Q(abs(m * k + p * delta), ell * p * q),
        R / b + R / q - Q(abs(s * h * k + b * delta), ell * b * q),
    ]
    roof = min(terms)
    return roof if roof > 0 else Q(0)


def quotient_measure(p, b, q, s, h, m, ell):
    parts = []
    for delta in (-3, 0, 3):
        subtotal = Q(0)
        radius = ceil_q(R * (p + b))
        for k in range(-radius, radius + 1):
            if (m * k + p * delta) % ell:
                continue
            if k % 3 == 0:
                continue
            subtotal += lift_length(p, b, q, s, h, m, ell, delta, k)
        parts.append(subtotal)
    return sum(parts, Q(0)), tuple(parts)


def failure_intervals(w):
    """Definition-level nearest-integer failure intervals on y in [0,1]."""
    inv = pow(w % 3, -1, 3)
    ans = []
    # Only these centers can meet [0,1].  Endpoints have measure zero.
    for n in range(0, w + 1):
        lo = max(Q(0), Q(n, w) - R / w)
        hi = min(Q(1), Q(n, w) + R / w)
        if lo < hi:
            owner = (-inv * n) % 3
            ans.append((lo, hi, owner, n))
    return ans


def physical_measure(p, b, q, s, t, h, m, ell):
    """Definition-level physical comb, including independently labelled defects."""
    lists = [failure_intervals(w) for w in (p, b, q)]
    idx = [0, 0, 0]
    parts = {-3: Q(0), 0: Q(0), 3: Q(0)}
    while all(idx[j] < len(lists[j]) for j in range(3)):
        cur = [lists[j][idx[j]] for j in range(3)]
        lo = max(item[0] for item in cur)
        hi = min(item[1] for item in cur)
        if lo < hi and len({item[2] for item in cur}) == 3:
            n_p, n_b, n_q = (item[3] for item in cur)
            delta = m * n_b - s * h * n_p - t * ell * n_q
            require(delta in parts, "physical component has illegal defect")
            parts[delta] += hi - lo
        first_end = min(item[1] for item in cur)
        for j in range(3):
            if cur[j][1] == first_end:
                idx[j] += 1
    ordered = tuple(parts[d] for d in (-3, 0, 3))
    return sum(ordered, Q(0)), ordered


def presentations_below(h, m, ell, cutoff):
    units = [x for x in range(1, cutoff) if odd_three_unit(x)]
    found = set()
    for p, q, s, t in product(units, units, (-1, 1), (-1, 1)):
        numerator = s * h * p + t * ell * q
        if numerator <= 0 or numerator % m:
            continue
        b = numerator // m
        if b >= cutoff or not odd_three_unit(b):
            continue
        if len({p, b, q}) != 3 or gcd3(p, b, q) != 1:
            continue
        found.add((p, b, q, s, t))
    return sorted(found)


def primitive_relations_below(speeds, limit):
    """All full-support primitive coefficient rays with l1 < limit."""
    out = set()
    bound = limit - 2
    for a in range(-bound, bound + 1):
        if a == 0:
            continue
        for b in range(-bound, bound + 1):
            if b == 0:
                continue
            remaining = limit - abs(a) - abs(b)
            if remaining <= 1:
                continue
            # Solve c exactly rather than scan it.
            rhs = -(a * speeds[0] + b * speeds[1])
            if rhs % speeds[2]:
                continue
            c = rhs // speeds[2]
            if c == 0 or abs(c) >= remaining:
                continue
            if gcd3(abs(a), abs(b), abs(c)) != 1:
                continue
            v = (a, b, c)
            # Quotient global sign by making the first nonzero entry positive.
            if v[0] < 0:
                v = tuple(-x for x in v)
            out.add(v)
    return sorted(out)


def audit_sector(name, h, m, ell, cutoff, expected_presentations,
                 expected_triples, expected_max, expected_set, expected_parts):
    require(gcd(m, ell) == 1, f"{name}: m must invert modulo ell")
    require(ell in (2, 4), f"{name}: even coefficient sidecar not recognized")
    require(h + m + ell == 16, f"{name}: wrong l1 shell")

    # Exact polygon clipping independently certifies the continuous areas.
    areas = tuple(strip_area(h, m, ell, d) for d in (-3, 0, 3))
    area_sum = sum(areas, Q(0))
    require(Q(2, 3 * ell) * area_sum == BULK, f"{name}: wrong bulk area")

    # Six sampled residue classes (two per defect plane), each with roof <=3/(7W).
    envelope = lambda W: BULK + Q(18, 7 * W)
    require(envelope(cutoff) < expected_max, f"{name}: cutoff does not close tail")
    require(envelope(cutoff - 1) >= expected_max, f"{name}: cutoff not minimal")

    pres = presentations_below(h, m, ell, cutoff)
    require(len(pres) == expected_presentations, f"{name}: presentation count")
    by_set = {}
    for p, b, q, s, t in pres:
        relation = m * b - s * h * p - t * ell * q
        require(relation == 0, f"{name}: bad relation")
        require(gcd(p, b) == 1, f"{name}: cyclic chart gcd failure")

        q_measure, parts = quotient_measure(p, b, q, s, h, m, ell)
        direct, direct_parts = physical_measure(p, b, q, s, t, h, m, ell)
        require(q_measure == direct, f"{name}: quotient/physical disagreement at {(p,b,q)}")
        require(parts == direct_parts, f"{name}: defect-layer disagreement at {(p,b,q)}")
        key = tuple(sorted((p, b, q)))
        if key not in by_set:
            by_set[key] = (direct, q_measure, parts)
        else:
            old_direct, old_measure, old_parts = by_set[key]
            require(direct == old_direct, f"{name}: presentation-dependent physical measure")
            require(q_measure == old_measure, f"{name}: presentation-dependent measure")
            # Defect labels can reverse under a second presentation; compare as a multiset.
            require(sorted(parts) == sorted(old_parts), f"{name}: presentation-dependent defect masses")

    require(len(by_set) == expected_triples, f"{name}: unordered triple count")
    maximum = max(v[0] for v in by_set.values())
    winners = sorted(k for k, v in by_set.items() if v[0] == maximum)
    require(maximum == expected_max, f"{name}: maximum mismatch")
    require(winners == [expected_set], f"{name}: uniqueness mismatch")
    win_direct, win_formula, win_parts = by_set[expected_set]
    require(win_direct == win_formula, f"{name}: winning formulas disagree")
    require(win_parts == expected_parts, f"{name}: winning defect masses")

    shorter = primitive_relations_below(expected_set, 16)
    require(shorter, f"{name}: expected a shorter relation collision")
    require(all(any(c % 3 == 0 for c in v) for v in shorter),
            f"{name}: shorter ternary-unit relation found")

    return {
        "name": name,
        "shape": [h, m, ell],
        "areas": [qtext(x) for x in areas],
        "bulk": qtext(Q(2, 3 * ell) * area_sum),
        "cutoff": cutoff,
        "presentations": len(pres),
        "unordered_triples": len(by_set),
        "maximum": qtext(maximum),
        "winner": list(expected_set),
        "parts": [qtext(x) for x in win_parts],
        "winner_presentations": sum(1 for p, b, q, s, t in pres
                                    if tuple(sorted((p, b, q))) == expected_set),
        "shorter_relations": [list(v) for v in shorter],
    }


def main():
    patterns = admissible_partitions(16)
    expected_patterns = [
        (1, 1, 14), (1, 2, 13), (1, 4, 11), (1, 5, 10),
        (1, 7, 8), (2, 7, 7), (4, 5, 7),
    ]
    require(patterns == expected_patterns, "l1=16 partition classification")
    non_one = [x for x in patterns if 1 not in x]
    require(non_one == [(2, 7, 7), (4, 5, 7)], "non-coefficient-one classification")

    sectors = [
        audit_sector(
            "277", 7, 7, 2, 366, 1754, 877,
            Q(304, 12397), (1, 23, 77),
            (Q(31, 12397), Q(22, 1127), Q(31, 12397)),
        ),
        audit_sector(
            "457", 5, 7, 4, 416, 2389, 2388,
            Q(2178, 91945), (5, 37, 71),
            (Q(2, 2485), Q(58, 2627), Q(2, 2485)),
        ),
    ]

    require(sectors[0]["areas"] == ["9/4802", "117/2401", "9/4802"],
            "277 strip areas")
    require(sectors[1]["areas"] == ["9/3430", "171/1715", "9/3430"],
            "457 strip areas")
    require(sectors[0]["winner_presentations"] == 2, "277 winner presentation count")
    require(sectors[1]["winner_presentations"] == 1, "457 winner presentation count")

    # Confirm the two advertised shorter collisions in sorted-speed coordinates.
    require([8, 3, -1] in sectors[0]["shorter_relations"], "missing 8+3-1 collision")
    require([8, -3, 1] in sectors[1]["shorter_relations"], "missing 8-3+1 collision")

    payload = {
        "status": "VERIFIED-EXACT CLEANROOM; UNIVERSAL LRC14 REMAINS OPEN",
        "patterns_l1_16": [list(x) for x in patterns],
        "non_coefficient_one": [list(x) for x in non_one],
        "sectors": sectors,
    }
    semantic = sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    print("LRC14 L1=16 EXTRA-SECTOR INDEPENDENT CLEANROOM REFEREE")
    print(json.dumps(payload, sort_keys=True, separators=(",", ":")))
    print(f"explicit_checks={CHECKS}")
    print(f"semantic_sha256={semantic}")


if __name__ == "__main__":
    main()
