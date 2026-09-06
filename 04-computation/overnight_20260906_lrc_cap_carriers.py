#!/usr/bin/env python3
"""Exact finite base for the all-height multi-direction cap-carrier theorem.

Universe: 1<=a<b<c<=100, odd ternary units, gcd(a,b,c)=1.
Independent carriers: rectangular integer box versus Bezout progression rows.
Cap means the FULL dictionary has no three collinear, never a selected subset.
"""
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations
from math import gcd
import json
import sys

CHECKS = 0
SHARP = Q(204, 5957)


def need(test, detail):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(detail)


def limits(w):
    return tuple((3*(sum(w)-x)-1)//14 for x in w)


def box(w):
    a, b, c = w
    bx, by, bz = limits(w)
    out = set()
    for x in range(-bx, bx+1):
        if not x % 3:
            continue
        for y in range(-by, by+1):
            if not y % 3:
                continue
            znum = -a*x-b*y
            if znum % c:
                continue
            z = znum//c
            if abs(z) <= bz and z % 3:
                out.add((x, y, z))
    return out


def egcd(a, b):
    if not b:
        return a, 1, 0
    g, x, y = egcd(b, a % b)
    return g, y, x-(a//b)*y


def ceildiv(a, b):
    return -((-a)//b)


def rows(w):
    # This does not invoke the box enumerator or its divisibility test.
    a, b, c = w
    d, u, v = egcd(b, c)
    need(gcd(a, d) == 1 and u*b+v*c == d, ("row primitive Bezout", w))
    bb, cc = b//d, c//d
    bx, by, bz = limits(w)
    out = set()
    for n in range(-(bx//d), bx//d+1):
        y0, z0 = -a*n*u, -a*n*v
        lo = max(ceildiv(-by-y0, cc), ceildiv(z0-bz, bb))
        hi = min((by-y0)//cc, (z0+bz)//bb)
        for t in range(lo, hi+1):
            point = (d*n, y0+cc*t, z0-bb*t)
            if all(x % 3 for x in point):
                out.add(point)
    return out


def cross(u, v):
    return (u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2],
            u[0]*v[1]-u[1]*v[0])


def sub(u, v):
    return tuple(x-y for x, y in zip(u, v))


def cap_triples(points):
    for p, q, r in combinations(sorted(points), 3):
        if cross(sub(q, p), sub(r, p)) == (0, 0, 0):
            return False, (p, q, r)
    return True, None


def cap_directions(points):
    for p in points:
        slopes = set()
        for q in points:
            if p == q:
                continue
            d = sub(q, p)
            gg = gcd(*d)
            d = tuple(x//gg for x in d)
            if next(x for x in d if x) < 0:
                d = tuple(-x for x in d)
            if d in slopes:
                return False
            slopes.add(d)
    return True


def phase(w, p):
    colors = tuple(a*b % 3 for a, b in zip(w, p))
    need(len(set(colors)) == 1 and colors[0] != 0, ("one owner color", w, p))
    return colors[0]


def midpoint_certificate(w, points):
    buckets = {}
    pair_count = 0
    witnesses = []
    for p in sorted(points):
        parity = tuple(x % 2 for x in p)
        need(sum(parity) % 2 == 0, ("rank two parity plane", w, p))
        key = phase(w, p), parity
        for q in buckets.get(key, []):
            midpoint = tuple((x+y)//2 for x, y in zip(p, q))
            need(midpoint in points and midpoint not in (p, q),
                 ("convex owner midpoint", w, p, q, midpoint))
            pair_count += 1
            witnesses.append(tuple(sorted((p, midpoint, q))))
        buckets.setdefault(key, []).append(p)
    need(len(set(witnesses)) == pair_count, ("endpoint pairs inject into AP triples", w))
    n = len(points)//2
    k, r = divmod(n, 4)
    lower = 2*(4*k*(k-1)//2+r*k)
    need(pair_count >= lower, ("balanced parity pair lower bound", w))
    if len(points) > 8:
        need(pair_count > 0, ("large dictionary forced collinear", w))
    return pair_count, lower


def projection(w, points):
    q = Q(3, 7*w[2])
    sums = []
    for i in range(3):
        j, k = [j for j in range(3) if j != i]
        sums.append(sum((min(q, Q(3*(w[j]+w[k])-14*abs(p[i]),
                                 14*w[j]*w[k])) for p in points), Q(0)))
    return tuple(sums)


def main():
    sys.stdout.reconfigure(newline="\n")
    values = [x for x in range(1, 101) if x % 2 and x % 3]
    universe = [w for w in combinations(values, 3) if gcd(*w) == 1]
    counts = Counter()
    cap_sizes = Counter()
    maxima = {}
    sharp_rows = []
    all_records = []
    forced_midpoints = 0
    for w in universe:
        points = box(w)
        alternate = rows(w)
        need(points == alternate, ("independent full carrier dictionary", w, points ^ alternate))
        need({tuple(-x for x in p) for p in points} == points, ("antipodal dictionary", w))
        pairs, lower = midpoint_certificate(w, points)
        forced_midpoints += pairs
        # Proven midpoint rejection is followed by an independent direction
        # test, so the cap test does not silently inherit a cardinality filter.
        is_cap = not pairs and cap_triples(points)[0]
        need(is_cap == cap_directions(points), ("independent cap test", w))
        multi = any(cross(p, q) != (0, 0, 0) for p, q in combinations(points, 2))
        if is_cap:
            need(len(points) <= 8, ("analytic cap bound", w))
            cap_sizes[len(points)] += 1
        counts["all"] += 1
        counts["cap"] += is_cap
        counts["multi_cap"] += is_cap and multi
        sums = None
        if is_cap and multi:
            sums = projection(w, points)
            value = min(sums)
            need(value <= SHARP, ("sharp finite cap ceiling", w, sums))
            if value == SHARP:
                sharp_rows.append((w, sums, sorted(points)))
            n = len(points)
            candidate = (value, w)
            if n not in maxima or candidate > maxima[n]:
                maxima[n] = candidate
        all_records.append((w, len(points), is_cap, multi,
                            None if sums is None else tuple(map(str, sums)), pairs, lower))
    need(len(sharp_rows) == 1 and sharp_rows[0][0] == (23, 29, 37), "unique cap equality")
    need(Q(24, 7*101) < SHARP < Q(6, 77), "analytic cap tail and target separation")
    # Selected subsets cannot be used as full cap certificates.
    hostile_w = (1, 11, 23)
    hostile = box(hostile_w)
    subset = {min(hostile), max(hostile)}
    need(not cap_directions(hostile) and cap_directions(subset), "full dictionary versus subset hostile")
    # The original equality row is a one-ray control outside our new scope.
    eq = box((1, 5, 11))
    need(len(eq) == 2 and projection((1, 5, 11), eq) == (Q(6, 77),)*3, "old one-ray equality control")
    # Physical interval construction is an extra independent control at the
    # new sharp row; map pair-name order back to coordinate-name order.
    from synthesis_20260905_lrc_sparse_transport import row
    network, mass, _ = row((23, 29, 37), audit_reference=True)
    need((network[2], network[1], network[0]) == sharp_rows[0][1], "physical network sharp control")
    print("status=PROVED_ANALYTICALLY cap bound; FINITE_EXACT sharp finite base")
    print("universe=all primitive sorted distinct positive odd ternary-unit speeds c<=100")
    print("counts="+str(dict(sorted(counts.items()))))
    print("cap_cardinalities="+str(dict(sorted(cap_sizes.items()))))
    print("cardinalitywise_multi_cap_maxima="+str(maxima))
    print("sharp_row="+str(sharp_rows[0]))
    print("sharp_physical_mass="+str(mass))
    print("universal_cap_bound=N<=8; universal_multi_cap_ceiling=204/5957")
    print("strict_target_gap="+str(Q(6, 77)-SHARP))
    print("tail_from_c=101; tail_bound_at_101="+str(Q(24, 707)))
    print("same_owner_midpoint_triples_checked="+str(forced_midpoints))
    print("subset_hostile="+str((hostile_w, sorted(hostile), sorted(subset))))
    print("checks="+str(CHECKS))
    print("semantic_sha256="+sha256(json.dumps(all_records, sort_keys=True).encode()).hexdigest())
    print("LRC14=OPEN; full carrier dictionaries with collinear triples remain")


if __name__ == "__main__":
    main()
