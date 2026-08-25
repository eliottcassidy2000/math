#!/usr/bin/env python3
"""Exact component geometry of the non-cover set for AP_m."""

from fractions import Fraction as Q
from math import floor, gcd


def covers_theta(theta: Q, m: int) -> bool:
    mask = 0
    for e in range(m):
        mask |= 1 << (floor(e * theta) % 7)
    return mask == 0x7F


def noncover_components(m: int):
    bps = {Q(0), Q(7)}
    for e in range(1, m):
        for r in range(7 * e + 1):
            bps.add(Q(r, e))
    bps = sorted(bps)
    components = []
    start = None
    end = None
    cells = 0
    for left, right in zip(bps, bps[1:]):
        bad = not covers_theta((left + right) / 2, m)
        if bad:
            if start is None:
                start = left
                cells = 0
            end = right
            cells += 1
        elif start is not None:
            components.append((start, end, cells))
            start = end = None
    if start is not None:
        components.append((start, end, cells))
    # Circle-glue endpoint components for semantic counting.
    if len(components) >= 2 and components[0][0] == 0 and components[-1][1] == 7:
        left = components[-1]
        right = components[0]
        wrapped = (left[0], right[1] + 7, left[2] + right[2])
        components = [wrapped] + components[1:-1]
    return components


def bad_rationals_theta():
    values = set()
    for q in range(1, 7):
        for p in range(q):
            if gcd(p, q) == 1:
                values.add(Q(7 * p, q))
    return sorted(values)


def circular_contains(left: Q, right: Q, x: Q) -> bool:
    if right <= 7:
        return left <= x <= right
    return left <= x + (7 if x < left else 0) <= right


def local_missing_sector_radii(p: int, q: int, m: int):
    """Guaranteed non-cover radii in theta around theta0=7p/q.

    For each sector missing at theta0, compute its first possible arrival under
    a positive/negative perturbation.  In each residue class e=s (mod q), the
    largest e<m is the earliest arrival, so this is a finite exact max-min.
    """
    n = m - 1
    data = []
    occupied = set()
    for s in range(q):
        A = Q(7 * p * s, q)
        integer = floor(A)
        frac = A - integer
        residue = integer % 7
        occupied.add(residue)
        E = n - ((n - s) % q)
        if E == 0:
            E += q
        data.append((residue, frac, E))
    missing = set(range(7)) - occupied
    plus_by_r = {}
    minus_by_r = {}
    for r in missing:
        plus_candidates = []
        minus_candidates = []
        for residue, frac, E in data:
            up = (r - residue) % 7
            down = (residue - r) % 7
            if not (1 <= up <= 6 and 1 <= down <= 6):
                raise RuntimeError((p, q, r, residue, "missing-sector direction"))
            plus_candidates.append((Q(up) - frac) / E)
            minus_candidates.append((Q(down - 1) + frac) / E)
        plus_by_r[r] = min(plus_candidates)
        minus_by_r[r] = min(minus_candidates)
    return max(minus_by_r.values()), max(plus_by_r.values()), minus_by_r, plus_by_r


def main():
    bad = bad_rationals_theta()
    print(f"persistent bad rational slopes ({len(bad)} on the circle): {bad}")
    for m in (14, 20, 27, 28, 40, 60, 80, 120):
        comps = noncover_components(m)
        total = sum((right - left for left, right, _ in comps), Q(0)) / 7
        print(f"\nm={m}: components={len(comps)}, deficit={total}")
        for left, right, cells in comps:
            owners = [x for x in bad if circular_contains(left, right, x)]
            owner = owners[0]
            owner_mod = owner % 7
            # theta owner 7p/q; recover the reduced x=p/q on the circle.
            x_owner = owner_mod / 7
            p, q = x_owner.numerator, x_owner.denominator
            minus, plus, _, _ = local_missing_sector_radii(p, q, m)
            actual_minus = (owner - left) if owner >= left else (owner + 7 - left)
            actual_plus = right - (owner if owner >= left else owner + 7)
            print(
                f"  theta=({left},{right}) width_x={(right-left)/7} "
                f"cells={cells} persistent_owners={owners} "
                f"radii=(-{actual_minus},+{actual_plus}) "
                f"maxmin=(-{minus},+{plus}) equal={actual_minus == minus and actual_plus == plus}"
            )


if __name__ == "__main__":
    main()
