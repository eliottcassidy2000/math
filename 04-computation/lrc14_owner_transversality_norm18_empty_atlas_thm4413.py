#!/usr/bin/env python3
"""Clean-room exact audit of the norm-18 carrier boundary and empty combs.

This program deliberately does not import repository code and does not use
the affine-fibre enumerator from the producer report.  It uses two unrelated
physical constructions:

* direct enumeration of the rank-two relation lattice C.w=0 inside the exact
  roof bounds; and
* a rational event sweep on y in [0,1), recomputing nearest integers and owner
  labels on every open cell.

All theorem gates use explicit exceptions and remain live under ``python -O``.
Only Python's standard library is used.
"""

from collections import defaultdict
from fractions import Fraction
from itertools import combinations, product
from math import gcd
import sys


R = Fraction(3, 14)
CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def gcd3(a, b, c):
    return gcd(gcd(abs(a), abs(b)), abs(c))


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def pattern(c):
    return tuple(sorted(abs(x) for x in c))


def pattern_code(c):
    # Hyphens avoid the traditional but ambiguous abbreviation 116, which can
    # be misread as (1,16) rather than (1,1,16).
    return "-".join(str(x) for x in pattern(c))


def owner_admissible(C):
    return all(x % 3 != 0 for x in C)


def margins(w, C):
    return (
        3 * (w[1] + w[2]) - 14 * abs(C[0]),
        3 * (w[0] + w[2]) - 14 * abs(C[1]),
        3 * (w[0] + w[1]) - 14 * abs(C[2]),
    )


def gap(w, C):
    return -min(margins(w, C))


def carrier_length(w, C):
    p = margins(w, C)
    if min(p) <= 0:
        return Fraction(0)
    candidates = [Fraction(3, 7 * x) for x in w]
    candidates.extend(
        (
            Fraction(p[0], 14 * w[1] * w[2]),
            Fraction(p[1], 14 * w[0] * w[2]),
            Fraction(p[2], 14 * w[0] * w[1]),
        )
    )
    return min(candidates)


def coefficient_rays(max_norm):
    """Primitive full-support ternary-unit coefficient rays, modulo +/-1."""
    by_norm = defaultdict(list)
    for norm in range(3, max_norm + 1):
        for a in range(1, norm - 1):
            if a % 3 == 0:
                continue
            for b in range(1, norm - a):
                d = norm - a - b
                if b % 3 == 0 or d <= 0 or d % 3 == 0:
                    continue
                if gcd3(a, b, d) != 1:
                    continue
                # Fix the global sign by requiring the first coordinate > 0.
                for sb, sd in product((-1, 1), repeat=2):
                    by_norm[norm].append((a, sb * b, sd * d))
    return dict(by_norm)


def speed_universe(height):
    values = [x for x in range(1, height + 1, 2) if x % 3 != 0]
    return [w for w in combinations(values, 3) if gcd3(*w) == 1]


def direct_carriers(w, keep_owner=True):
    """Enumerate C.w=0 in the strict roof, without choosing a relation chart."""
    bounds = (
        (3 * (w[1] + w[2]) - 1) // 14,
        (3 * (w[0] + w[2]) - 1) // 14,
        (3 * (w[0] + w[1]) - 1) // 14,
    )
    answer = {}
    for C0 in range(-bounds[0], bounds[0] + 1):
        for C1 in range(-bounds[1], bounds[1] + 1):
            numerator = w[0] * C0 + w[1] * C1
            if numerator % w[2] != 0:
                continue
            C2 = -numerator // w[2]
            if abs(C2) > bounds[2]:
                continue
            C = (C0, C1, C2)
            if keep_owner and not owner_admissible(C):
                continue
            length = carrier_length(w, C)
            require(length > 0, f"strict-bound carrier was not live: {w=} {C=}")
            answer[C] = length
    return answer


def nearest_integer_nonnegative(x):
    # Event midpoints never lie at a half-integer on an eligible sheet.
    return (x.numerator * 2 + x.denominator) // (2 * x.denominator)


def literal_event_sweep(w):
    """Direct physical construction on y in [0,1), independent of C.w=0."""
    events = {Fraction(0), Fraction(1)}
    for speed in w:
        for n in range(speed + 1):
            for sign in (-1, 1):
                endpoint = Fraction(14 * n + sign * 3, 14 * speed)
                if 0 < endpoint < 1:
                    events.add(endpoint)
    events = sorted(events)
    pieces = defaultdict(Fraction)
    for left, right in zip(events, events[1:]):
        middle = (left + right) / 2
        lifts = tuple(nearest_integer_nonnegative(speed * middle) for speed in w)
        errors = tuple(speed * middle - n for speed, n in zip(w, lifts))
        if not all(abs(e) < R for e in errors):
            continue
        owners = tuple((-pow(speed, -1, 3) * n) % 3 for speed, n in zip(w, lifts))
        if len(set(owners)) != 3:
            continue
        C = cross(w, lifts)
        pieces[C] += right - left
    return dict(pieces)


def egcd_nonnegative(a, b):
    old_r, r = a, b
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t
    return old_r, old_s, old_t


def signed_bezout_pair(a, b):
    g, x, y = egcd_nonnegative(abs(a), abs(b))
    return g, x if a >= 0 else -x, y if b >= 0 else -y


def bezout3(c):
    g12, x, y = signed_bezout_pair(c[0], c[1])
    g, u, z = signed_bezout_pair(g12, c[2])
    require(g == 1, f"coefficient row not primitive: {c}")
    answer = (u * x, u * y, z)
    require(dot(c, answer) == 1, f"bad Bezout vector: {c=} {answer=}")
    return answer


def ceil_fraction(x):
    return -((-x.numerator) // x.denominator)


def floor_fraction(x):
    return x.numerator // x.denominator


def fibre_min_gap(w, c, delta):
    """Exact min D on the owner-admissible defect fibre c x C=delta*w."""
    u = bezout3(c)
    C0 = cross(w, tuple(delta * x for x in u))
    require(cross(c, C0) == tuple(delta * x for x in w), "bad defect basepoint")
    upper = gap(w, C0)
    # Any minimizer has D<=D(C0).  The first coordinate alone then gives a
    # rigorous finite interval for k in C=C0+k*c.  Here c[0]>0 by convention.
    radius = Fraction(upper + 3 * (w[1] + w[2]), 14)
    lo = ceil_fraction((Fraction(-C0[0]) - radius) / c[0])
    hi = floor_fraction((Fraction(-C0[0]) + radius) / c[0])
    values = []
    for k in range(lo, hi + 1):
        C = tuple(C0[i] + k * c[i] for i in range(3))
        if owner_admissible(C):
            values.append((gap(w, C), C, k))
    require(values, f"no owner samples in defect fibre: {w=} {c=} {delta=}")
    return min(values)


def defect_of(w, c, C):
    z = cross(c, C)
    require(all(z[i] % w[i] == 0 for i in range(3)), "nonintegral defect")
    qs = tuple(z[i] // w[i] for i in range(3))
    require(len(set(qs)) == 1, f"cross product not parallel to w: {w=} {c=} {C=}")
    return qs[0]


def admissible_cutoff(M):
    bound = (14 * M) // 3
    while bound % 2 == 0 or bound % 3 == 0:
        bound -= 1
    return bound


def format_fraction(x):
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def main():
    sys.stdout.reconfigure(newline="\n")
    rays = coefficient_rays(18)
    # Odd speed triples can only support even l1 norms: modulo two,
    # c.w and sum(c_i) have the same parity, as does sum(abs(c_i)).
    short = [c for n, rows in rays.items() if n < 18 and n % 2 == 0 for c in rows]
    norm18 = rays[18]

    patterns18 = sorted({pattern(c) for c in norm18})
    expected_patterns = [
        (1, 1, 16),
        (1, 4, 13),
        (1, 7, 10),
        (2, 5, 11),
        (4, 7, 7),
        (5, 5, 8),
    ]
    require(patterns18 == expected_patterns, f"wrong norm-18 patterns: {patterns18}")
    short_patterns = sorted({pattern(c) for c in short})
    require(len(short_patterns) == 20, f"expected 20 shorter patterns, got {short_patterns}")

    cutoffs = {pat: admissible_cutoff(max(pat)) for pat in patterns18}
    expected_cutoffs = {
        (1, 1, 16): 73,
        (1, 4, 13): 59,
        (1, 7, 10): 43,
        (2, 5, 11): 49,
        (4, 7, 7): 31,
        (5, 5, 8): 37,
    }
    require(cutoffs == expected_cutoffs, f"wrong analytic cutoffs: {cutoffs}")

    universe = speed_universe(max(cutoffs.values()))
    require(len(universe) == 2289, f"wrong height-73 universe size: {len(universe)}")
    all_relations = {}
    minimal_relations = {}
    raw_incidence_count = 0
    for w in universe:
        rel18 = tuple(c for c in norm18 if dot(c, w) == 0)
        if not rel18:
            continue
        all_relations[w] = rel18
        raw_incidence_count += len(rel18)
        has_shorter = any(dot(c, w) == 0 for c in short)
        if not has_shorter:
            minimal_relations[w] = rel18
    require(raw_incidence_count == 232, f"raw incidence mismatch: {raw_incidence_count}")
    require(len(all_relations) == 209, f"raw triple mismatch: {len(all_relations)}")
    minimal_incidence_count = sum(map(len, minimal_relations.values()))
    require(minimal_incidence_count == 190, f"minimal incidence mismatch: {minimal_incidence_count}")
    require(len(minimal_relations) == 180, f"minimal triple mismatch: {len(minimal_relations)}")

    # Check the algebraic cutoff implication on every relation incidence in
    # the complete coarse core, before the minimality filter is applied.
    for w, rels in all_relations.items():
        for c in rels:
            if gap(w, c) < 2:
                continue
            pat = pattern(c)
            require(w[2] <= cutoffs[pat], f"central failure escaped its cutoff: {w=} {c=}")
            failed_coordinates = [i for i, p in enumerate(margins(w, c)) if p < 0]
            require(failed_coordinates, f"dead carrier lacks a failed coordinate: {w=} {c=}")
            for i in failed_coordinates:
                j, k = [q for q in range(3) if q != i]
                require(3 * (w[j] + w[k]) <= 14 * abs(c[i]), "pair-sum cutoff algebra failed")
                require(3 * w[i] <= 14 * max(abs(c[j]), abs(c[k])), "opposite-coordinate cutoff algebra failed")

    # A second physical implementation agrees carrier-by-carrier on the whole
    # finite core forced by the analytic cutoff.
    physical = {}
    grouped_component_count = 0
    for w in minimal_relations:
        lattice = direct_carriers(w)
        literal = literal_event_sweep(w)
        require(lattice == literal, f"literal/lattice mismatch: {w=} {lattice=} {literal=}")
        physical[w] = lattice
        grouped_component_count += len(lattice)
        for C, length in lattice.items():
            require(owner_admissible(C), f"owner failure among live carriers: {w=} {C=}")
            require(dot(C, w) == 0, f"nonrelation live carrier: {w=} {C=}")
            p = margins(w, C)
            require(all(x % 2 == 0 and x % 3 != 0 for x in p), f"margin alphabet failed: {w=} {C=} {p=}")
            require(gap(w, C) <= -2, f"live gap breached: {w=} {C=} D={gap(w,C)}")
            require(length >= Fraction(1, 7 * w[1] * w[2]), f"component floor breached: {w=} {C=}")
            adaptive_floor = Fraction(min(min(p), 6 * w[1]), 14 * w[1] * w[2])
            require(length >= adaptive_floor, f"adaptive component floor breached: {w=} {C=}")

    empty = sorted(w for w, pieces in physical.items() if not pieces)
    expected_empty = [(7, 11, 43), (7, 11, 47), (7, 25, 29)]
    require(empty == expected_empty, f"wrong empty physical triples: {empty}")

    central_dead = []
    by_pattern = {}
    for pat in patterns18:
        rows = [
            (w, c)
            for w, rels in minimal_relations.items()
            if w[2] <= cutoffs[pat]
            for c in rels
            if pattern(c) == pat
        ]
        triples = {w for w, _ in rows}
        dead_rows = [(w, c) for w, c in rows if gap(w, c) >= 2]
        by_pattern[pat] = (len(triples), len(rows), len(dead_rows))
        central_dead.extend(dead_rows)

    expected_pattern_census = {
        (1, 1, 16): (14, 14, 2),
        (1, 4, 13): (16, 17, 4),
        (1, 7, 10): (8, 9, 1),
        (2, 5, 11): (12, 13, 2),
        (4, 7, 7): (1, 1, 1),
        (5, 5, 8): (5, 5, 0),
    }
    require(by_pattern == expected_pattern_census, f"pattern census mismatch: {by_pattern}")
    require(len(central_dead) == 10, f"central-dead count mismatch: {len(central_dead)}")

    dead_rows_summary = []
    rescued = 0
    for w, c in sorted(central_dead):
        zero_gap = gap(w, c)
        require(zero_gap >= 2 and zero_gap % 2 == 0, "central dead gap parity failed")
        plus_gap, plus_C, _ = fibre_min_gap(w, c, 3)
        pieces = physical[w]
        counts = defaultdict(int)
        masses = defaultdict(Fraction)
        for C, length in pieces.items():
            delta = defect_of(w, c, C)
            require(delta in (-3, 0, 3), f"unexpected norm-18 defect: {w=} {c=} {C=} {delta=}")
            counts[delta] += 1
            masses[delta] += length
        profile = (counts[-3], counts[0], counts[3])
        total = sum(pieces.values(), Fraction())
        if total > 0:
            rescued += 1
            require(plus_gap <= -2, f"positive comb not rescued in defect 3: {w=} {c=}")
        else:
            require(plus_gap >= 2, f"empty comb has a live defect-3 sample: {w=} {c=} {plus_C=}")
        dead_rows_summary.append((w, c, pattern_code(c), zero_gap, plus_gap, profile, total))
    require(rescued == 6, f"wrong rescued presentation count: {rescued}")

    expected_dead_rows = [
        ((1, 17, 37), (2, -11, 5), "2-5-11", 40, -8, (1, 0, 1), Fraction(8, 4403)),
        ((1, 17, 55), (1, -13, 4), "1-4-13", 14, -34, (1, 0, 1), Fraction(2, 385)),
        ((1, 19, 35), (16, 1, -1), "1-1-16", 62, -4, (1, 0, 1), Fraction(64, 4655)),
        ((1, 25, 41), (16, 1, -1), "1-1-16", 26, -22, (2, 0, 2), Fraction(172, 7175)),
        ((5, 13, 41), (1, -13, 4), "1-4-13", 44, -22, (1, 0, 1), Fraction(22, 3731)),
        ((7, 11, 43), (5, -11, 2), "2-5-11", 4, 34, (0, 0, 0), Fraction(0)),
        ((7, 11, 47), (13, -4, -1), "1-4-13", 8, 20, (0, 0, 0), Fraction(0)),
        ((7, 25, 29), (4, 7, -7), "4-7-7", 2, 20, (0, 0, 0), Fraction(0)),
        ((7, 25, 29), (13, 1, -4), "1-4-13", 20, 2, (0, 0, 0), Fraction(0)),
        ((13, 23, 31), (1, -10, 7), "1-7-10", 8, -22, (1, 0, 1), Fraction(22, 4991)),
    ]
    require(dead_rows_summary == expected_dead_rows, f"central-dead atlas mismatch: {dead_rows_summary}")

    # Sharpening: after the analytic 73 cutoff makes the search complete, the
    # exact minimal shell has a substantially smaller attained central-dead
    # ceiling in every populated pattern.
    actual_dead_heights = {}
    for pat in patterns18:
        heights = [w[2] for w, c in central_dead if pattern(c) == pat]
        actual_dead_heights[pat] = max(heights) if heights else None
    expected_actual_heights = {
        (1, 1, 16): 41,
        (1, 4, 13): 55,
        (1, 7, 10): 31,
        (2, 5, 11): 43,
        (4, 7, 7): 29,
        (5, 5, 8): None,
    }
    require(actual_dead_heights == expected_actual_heights, f"actual dead ceilings mismatch: {actual_dead_heights}")

    # Sharp live/dead boundary witnesses inside the minimal norm-18 shell.
    live_w = (29, 73, 77)
    live_c = (7, -7, 4)
    live_C = (32, 1, -13)
    require(dot(live_c, live_w) == 0 and sum(map(abs, live_c)) == 18, "bad live witness chart")
    require(not any(dot(c, live_w) == 0 for c in short), "live witness has a shorter unit relation")
    require(owner_admissible(live_C) and dot(live_C, live_w) == 0, "bad live witness carrier")
    require(margins(live_w, live_C) == (2, 304, 124), "unexpected live witness margins")
    require(gap(live_w, live_C) == -2, "live witness does not attain D=-2")
    require(carrier_length(live_w, live_C) == Fraction(1, 7 * 73 * 77), "floor witness length mismatch")

    dead_w = (7, 25, 29)
    dead_C = (4, 7, -7)
    require(owner_admissible(dead_C) and dot(dead_C, dead_w) == 0, "bad dead witness")
    require(margins(dead_w, dead_C) == (106, 10, -2), "unexpected dead witness margins")
    require(gap(dead_w, dead_C) == 2 and carrier_length(dead_w, dead_C) == 0, "dead witness not sharp")

    # Hostile controls isolate the two hypotheses in the arithmetic gap.
    tangent_w = (1, 5, 13)
    tangent_C = (-2, 3, -1)
    require(dot(tangent_C, tangent_w) == 0, "bad tangent hostile")
    require(margins(tangent_w, tangent_C) == (26, 0, 4), "bad tangent margins")
    require(not owner_admissible(tangent_C) and gap(tangent_w, tangent_C) == 0, "owner hostile failed")

    even_w = (1, 5, 14)
    even_C = (-4, -2, 1)
    require(dot(even_C, even_w) == 0 and owner_admissible(even_C), "bad parity hostile")
    require(margins(even_w, even_C) == (1, 17, 4), "bad parity-hostile margins")
    require(carrier_length(even_w, even_C) == Fraction(1, 980), "parity hostile length mismatch")
    require(Fraction(1, 980) < Fraction(1, 7 * 5 * 14), "parity hostile did not break odd floor")

    # Failure anatomy: the three empty roofs have positive geometric lattice
    # points, but every such point is removed by the owner gate.  A distinguished
    # defect-three sample is the index-three one-zero complement 3q.
    complements = {
        (7, 11, 43): ((3, 2, -1), Fraction(18, 3311), ((5, -11, 2),)),
        (7, 11, 47): ((2, 3, -1), Fraction(18, 2303), ((13, -4, -1),)),
        (7, 25, 29): ((-3, 2, -1), Fraction(18, 5075), ((4, 7, -7), (13, 1, -4))),
    }
    for w, (q, expected_length, charts) in complements.items():
        C = tuple(3 * x for x in q)
        geometry = direct_carriers(w, keep_owner=False)
        owners = direct_carriers(w, keep_owner=True)
        require(all(cross(c, q) == w for c in charts), f"bad complementary row: {w=} {q=}")
        require(len(geometry) == 11, f"unexpected geometric-roof count: {w=} count={len(geometry)}")
        require(C in geometry and geometry[C] == expected_length, f"missing deleted geometric carrier: {w=} {C=}")
        require(not owner_admissible(C) and not owners, f"empty failure anatomy mismatch: {w=}")

    # The exact all-height shell-specific sharpening follows because every
    # central failure lies under the proved analytic cutoff, and the complete
    # finite census has now evaluated every such presentation.
    require(max(w[2] for w, _ in central_dead) == 55, "global central-dead ceiling is not 55")
    require(max(w[2] for w in empty) == 47, "global empty ceiling is not 47")

    print("LRC14 NORM-18 GAP/EMPTY-COMB CLEAN-ROOM REFEREE")
    print("status=PASS exact_all_height_local_shell; LRC14=OPEN")
    print("method=coefficient_box + direct_rank2_carriers + rational_circle_event_sweep")
    print("patterns=" + ",".join("-".join(map(str, p)) for p in patterns18))
    print("owner_boundary=margins_even_and_not_divisible_by_3; live_D<=-2; dead_D>=2; D=0_tangent_requires_owner_failure")
    print("component_floor=min(p_star,6*w2)/(14*w2*w3)>=1/(7*w2*w3); sharp=(29,73,77),C=(32,1,-13)")
    print("parity_hostile=(1,5,14),C=(-4,-2,1),L=1/980")
    print("analytic_central_cutoffs=" + ",".join(f"{'-'.join(map(str,p))}:{cutoffs[p]}" for p in patterns18))
    print(f"finite_core={len(universe)}_triples,{raw_incidence_count}_raw_incidences,{minimal_incidence_count}_minimal_incidences,{len(minimal_relations)}_minimal_triples")
    print(f"literal_lattice_components={grouped_component_count}")
    print("pattern_census=" + ";".join(f"{'-'.join(map(str,p))}:{by_pattern[p]}" for p in patterns18))
    print("central_dead_actual_height_ceiling=" + ",".join(f"{'-'.join(map(str,p))}:{actual_dead_heights[p]}" for p in patterns18))
    print(f"central_dead_presentations={len(central_dead)}; defect3_rescued={rescued}")
    for w, c, pat, d0, d3, profile, mass in dead_rows_summary:
        print(f"central_dead={w};c={c};pattern={pat};D0={d0};D3={d3};counts={profile};mass={format_fraction(mass)}")
    print("physical_empty=" + ",".join(map(str, empty)))
    print("sharp_complete_ceiling=central_dead_height_55; empty_height_47; first_empty_height_29")
    print("empty_geometric_roofs=11_carriers_each; deleted_defect3=3*(3,2,-1),3*(2,3,-1),3*(-3,2,-1)")
    print(f"checks={CHECKS}")


if __name__ == "__main__":
    main()
