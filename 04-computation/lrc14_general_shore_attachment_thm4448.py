#!/usr/bin/env python3
"""Exact self-contained audit for THM-4448 general shore attachment.

The general r-shore theorem is proved topologically in the theorem/result
note.  This script supplies its exact signed-(1,1,2) tail cells, filtered and
unfiltered two-shore gap atlases, bounded body census, literal controls, and
cofinal prescribed-component hostiles.  Every interval endpoint is rational.
Danger intervals are open; touching endpoints remain separate.
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction as Q
from itertools import combinations
from math import gcd


TAILS = ((1, 5, 11), (2, 11, 20))


def dist(x: Q) -> Q:
    x %= 1
    return min(x, 1 - x)


def safe_at(speeds: tuple[int, ...], x: Q, threshold: Q = Q(1, 14)) -> bool:
    return all(dist(Q(v) * x) >= threshold for v in speeds)


def danger_intervals(speed: int, radius_numerator: int = 1):
    radius = Q(radius_numerator, 14 * speed)
    out = []
    for k in range(speed + 1):
        lo = max(Q(0), Q(k, speed) - radius)
        hi = min(Q(1), Q(k, speed) + radius)
        if lo < hi:
            out.append((lo, hi, speed, k))
    return out


def merge(intervals):
    out = []
    for lo, hi, *payload in sorted(intervals):
        # The intervals are open.  Equality leaves a safe singleton and must
        # not be merged.
        if not out or lo >= out[-1][1]:
            out.append([lo, hi])
        elif hi > out[-1][1]:
            out[-1][1] = hi
    return tuple((lo, hi) for lo, hi in out)


def safe_components(speeds: tuple[int, ...]):
    bad = merge([z for v in speeds for z in danger_intervals(v)])
    out = []
    cursor = Q(0)
    for lo, hi in bad:
        if cursor <= lo and cursor != 0:
            out.append((cursor, lo))
        cursor = max(cursor, hi)
    if cursor < 1:
        out.append((cursor, Q(1)))
    return tuple(out)


def owner(t: int, y: Q):
    """Return (active, nearest integer, killed sheet), away from a wall."""
    z = Q(t) * y
    n = (2 * z.numerator + z.denominator) // (2 * z.denominator)
    delta = abs(z - n)
    if delta >= Q(3, 14):
        return False, n, None
    return True, n, (-n * pow(t, -1, 3)) % 3


def fully_spoiled(tails: tuple[int, int, int], y: Q) -> bool:
    data = [owner(t, y) for t in tails]
    return all(x[0] for x in data) and len({x[2] for x in data}) == 3


def tail_walls(tails: tuple[int, int, int]):
    walls = {Q(0), Q(1)}
    events = defaultdict(list)
    for t in tails:
        for k in range(t + 1):
            for sign, kind in ((-1, "enter"), (1, "exit")):
                y = Q(k, t) + sign * Q(3, 14 * t)
                if 0 <= y <= 1:
                    walls.add(y)
                    events[y].append((t, k, kind))
    return tuple(sorted(walls)), events


def failure_components(tails: tuple[int, int, int]):
    walls, events = tail_walls(tails)
    cells = []
    for lo, hi in zip(walls, walls[1:]):
        mid = (lo + hi) / 2
        if fully_spoiled(tails, mid):
            signature = tuple(owner(t, mid)[2] for t in tails)
            # A wall point is equality-safe for its changing tail, so two
            # adjacent open cells are distinct F_T components even if their
            # midpoint signatures happen to agree.
            cells.append((lo, hi, signature))
    return tuple(cells), events


def component_address(component, tails):
    lo, hi = component
    mid = (lo + hi) / 2
    rows = []
    for t in tails:
        active, n, sheet = owner(t, mid)
        left = Q(n, t) - Q(3, 14 * t)
        right = Q(n, t) + Q(3, 14 * t)
        rows.append((t, active, n, sheet, left, right, lo-left, right-hi))
    return tuple(rows)


def containment_certificate(
    body: tuple[int, ...], tails: tuple[int, int, int], comps=None, fcomps=None
):
    """Exact G_C subset F_T decision and first escaping component/cell."""
    if comps is None:
        comps = safe_components(body)
    if fcomps is None:
        fcomps, _ = failure_components(tails)
    assignments = []
    for i, (lo, hi) in enumerate(comps):
        containing = [
            (j, flo, fhi, sig)
            for j, (flo, fhi, sig) in enumerate(fcomps)
            if flo < lo and hi < fhi
        ]
        if not containing:
            # Return an exact point in this component outside the open F_T.
            cuts = sorted(
                {lo, hi}
                | {x for flo, fhi, _ in fcomps for x in (flo, fhi) if lo < x < hi}
            )
            candidates = [lo, hi] + [(a+b)/2 for a, b in zip(cuts, cuts[1:])]
            y = next(y for y in candidates if lo <= y <= hi and not fully_spoiled(tails, y))
            safe_sheets = tuple(j for j in range(3) if safe_at(tails, (y+j)/3))
            return False, comps, assignments, (i, lo, hi, y, safe_sheets)
        assignments.append(containing[0])
    return True, comps, assignments, None


def cross_height(body: tuple[int, ...], pair: tuple[int, int]) -> int:
    shore = set(pair)
    return min(
        max(x // gcd(x, z), z // gcd(x, z))
        for x in pair for z in body if z not in shore
    )


def max_distinguished_cross_height(body: tuple[int, ...]):
    vals = [(cross_height(body, p), p) for p in combinations(body, 2)]
    return max(vals)


def factor(n: int):
    out = []
    p = 2
    while p*p <= n:
        if n % p:
            p += 1
            continue
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        out.append((p, e))
        p += 1
    if n > 1:
        out.append((n, 1))
    return out


def decoder_pair(pair):
    x, y = pair
    g = gcd(x, y)
    p, q = x//g, y//g
    return p < q and p+q <= 356 and all(r % 3 == 2 and e <= 2 for r, e in factor(p+q))


def best_decoder_cross_height(body):
    rows = []
    for pair in combinations(body, 2):
        if decoder_pair(pair):
            rows.append((cross_height(body, pair), pair))
    return max(rows, default=(-1, None))


def interval_intersection_measure(left, right):
    total = Q(0)
    for a, b in left:
        for c, d, *_ in right:
            total += max(Q(0), min(b, d)-max(a, c))
    return total


def longest_circle_danger_component(speeds: tuple[int, ...]):
    pieces = list(merge([z for v in speeds for z in danger_intervals(v)]))
    if pieces and pieces[0][0] == 0 and pieces[-1][1] == 1 and len(pieces) > 1:
        wrap = pieces[0][1] + (1-pieces[-1][0])
        middle = [hi-lo for lo,hi in pieces[1:-1]]
        return max([wrap] + middle)
    return max((hi-lo for lo,hi in pieces), default=Q(0))


def decoder_pair_gap_atlas():
    rows = []
    for p in range(1, 356):
        for q in range(p+1, 357-p):
            if gcd(p,q)==1 and decoder_pair((p,q)):
                rows.append((longest_circle_danger_component((p,q)),p,q))
    return sorted(rows, reverse=True)


def all_bounded_pair_gap_atlas():
    rows=[]
    for p in range(1,178):
        for q in range(p+1,357-p):
            if gcd(p,q)==1:
                rows.append((longest_circle_danger_component((p,q)),p,q))
    return sorted(rows,reverse=True)


EXPECTED_FAILURE = {
    (1, 5, 11): (
        (Q(25,154), Q(31,154), (0,1,2)),
        (Q(123,154), Q(129,154), (2,1,0)),
    ),
    (2, 11, 20): (
        (Q(5,56), Q(3,28), (0,1,2)),
        (Q(123,280), Q(129,280), (1,2,0)),
        (Q(151,280), Q(157,280), (1,0,2)),
        (Q(25,28), Q(51,56), (2,1,0)),
    ),
}
EXPECTED_MEASURE = {(1,5,11): Q(6,77), (2,11,20): Q(11,140)}
HEIGHTS = (13,14,16,18)


def req(x, msg):
    if not x:
        raise AssertionError(msg)


def integer_grid_gap(p: int, q: int) -> Q:
    """Independent integer-endpoint union on denominator 14*p*q."""
    n = 14*p*q
    intervals = []
    for speed, scale in ((p,q),(q,p)):
        for k in range(speed+1):
            lo = max(0, scale*(14*k-1))
            hi = min(n, scale*(14*k+1))
            if lo < hi:
                intervals.append((lo,hi))
    merged = []
    for lo,hi in sorted(intervals):
        if not merged or lo >= merged[-1][1]:
            merged.append([lo,hi])
        elif hi > merged[-1][1]:
            merged[-1][1] = hi
    if len(merged)>1 and merged[0][0]==0 and merged[-1][1]==n:
        lengths = [merged[0][1] + n-merged[-1][0]]
        lengths.extend(hi-lo for lo,hi in merged[1:-1])
    else:
        lengths = [hi-lo for lo,hi in merged]
    return Q(max(lengths), n)


def scan_head():
    records = {
        T: {h: defaultdict(int) for h in HEIGHTS} for T in TAILS
    }
    extremal = {T: {h: {"max_low_kappa":(-1,None,None), "min_escape":None,
                         "max_overlap":None} for h in HEIGHTS} for T in TAILS}
    failures = {T:failure_components(T)[0] for T in TAILS}
    for body in combinations(range(1,19),10):
        req(gcd(*body)==1, "all H18 ten-subsets are primitive")
        H=max(body)
        comps=safe_components(body)
        measure=sum((b-a for a,b in comps),Q(0))
        for T in TAILS:
            contained,_,_,escape_record=containment_certificate(body,T,comps,failures[T])
            if H <= 13:
                physical=tuple(3*c for c in body)+T
                literal=bool(safe_components(tuple(sorted(physical))))
                req(literal == (not contained), f"literal containment equivalence {body} {T}")
            overlap=interval_intersection_measure(comps,failures[T])
            escape_measure=measure-overlap
            for h in HEIGHTS:
                if H>h:
                    continue
                R=records[T][h]
                E=extremal[T][h]
                R["bodies"]+=1
                R["containments"]+=int(contained)
                R["components"]+=len(comps)
                if H <= 13:
                    R["literal_rows"]+=1
                if measure < EXPECTED_MEASURE[T]:
                    R["scalar_low"]+=1
                    kappa,pair=best_decoder_cross_height(body)
                    E["max_low_kappa"]=max(E["max_low_kappa"],(kappa,body,pair))
                elif measure == EXPECTED_MEASURE[T]:
                    R["scalar_equal"]+=1
                e=(escape_measure,body,escape_record)
                o=(overlap,body)
                E["min_escape"]=e if E["min_escape"] is None else min(E["min_escape"],e)
                E["max_overlap"]=o if E["max_overlap"] is None else max(E["max_overlap"],o)
    return records,extremal



def equality_overlap_checks():
    controls = (
        ((1,5,11), (1,2,3,4,7,8,9,13,14,19), Q(6,77)),
        ((2,11,20), (1,2,3,4,5,6,7,8,12,14), Q(11,140)),
    )
    rows = []
    for tails, body, expected in controls:
        comps = safe_components(body)
        fcomps = failure_components(tails)[0]
        overlap = interval_intersection_measure(comps, fcomps)
        contained = containment_certificate(body, tails, comps, fcomps)[0]
        req(overlap == expected, ("full failure-set overlap", tails, body, overlap))
        req(not contained, ("overlap is not containment", tails, body))
        rows.append((tails, body, overlap, contained))
    return tuple(rows)


def hostile_family_checks():
    A=tuple(range(1,9))
    rows=[]
    controls=(
        ((1,5,11),53,Q(2,11),Q(3,154)),
        ((2,11,20),121,Q(1,10),Q(1,140)),
    )
    for T,N,y0,boundary_distance in controls:
        body=A+(N,4*N)
        req(gcd(*body)==1 and len(set(body))==10,"primitive distinct hostile family")
        req(cross_height(body,(N,4*N))==N,"cross height")
        req(all(dist(Q(c)*y0)>Q(1,14) for c in body),"strict body safety")
        req(fully_spoiled(T,y0),"strict full spoilage")
        req(Q(6,7*N)<boundary_distance,"component diameter bound")
        req(14*N>=87*max(A) and 14*N>=29*max(T),"simultaneous cofinal closure")
        rows.append((T,N,y0,boundary_distance,cross_height(body,(N,4*N))))
    return rows


def main():
    print("LRC14_GENERAL_SHORE_ATTACHMENT_THM4448_PRIMARY")
    print("STATUS=FINITE_EXACT_SUPPORT_FOR_PROVED_ATTACHMENT;LRC14_OPEN")
    for T in TAILS:
        F,events=failure_components(T)
        req(F==EXPECTED_FAILURE[T],f"failure components {T}")
        measure=sum((b-a for a,b,_ in F),Q(0))
        req(measure==EXPECTED_MEASURE[T],f"failure measure {T}")
        print(f"tail={T} failure_components={F} measure={measure}")
        for lo,hi,sig in F:
            print(f"  interval=({lo},{hi}) owners={sig} left={tuple(sorted(events[lo]))} right={tuple(sorted(events[hi]))}")

    atlas=decoder_pair_gap_atlas()
    req(len(atlas)==5855,"decoder atlas count")
    for gap,p,q in atlas:
        req(gap==integer_grid_gap(p,q),f"gap engines {(p,q)}")
    leaders=[row for row in atlas if row[0]==atlas[0][0]]
    req(atlas[0]==(Q(29,196),1,28) and leaders==[atlas[0]],"sharp gap leader")
    print(f"decoder_pairs={len(atlas)} sharp_pair_danger_gap={atlas[0]} unique=True")
    all_atlas=all_bounded_pair_gap_atlas()
    for gap,p,q in all_atlas:
        req(gap==integer_grid_gap(p,q),f"all-pair gap engines {(p,q)}")
    all_leaders=[row for row in all_atlas if row[0]==all_atlas[0][0]]
    req(all_atlas[0]==(Q(15,98),1,14) and all_leaders==[all_atlas[0]],"all-pair leader")
    print(f"all_coprime_pairs={len(all_atlas)} sharp_unfiltered_gap={all_atlas[0]} unique=True")
    print("scope_guard=29/196_requires_THM3818_inert_sum_decoder_filter")
    print("unfiltered_cone_integer=7*h>=45*M_A and 7*h>=15*maxT")
    print("uniform_cone=2*h*rho>=29/196;rho=min(1/(84*M_A),1/(28*maxT))")
    print("uniform_cone_integer=14*h>=87*M_A and 14*h>=29*maxT")

    records,extremal=scan_head()
    for T in TAILS:
        print(f"HEAD tail={T}")
        for h in HEIGHTS:
            R=dict(records[T][h]);E=extremal[T][h]
            req(R["containments"]==0,f"head containment {T} H{h}")
            print(f"  H={h} counts={R} max_low_kappa={E['max_low_kappa']} min_escape={E['min_escape']} max_overlap={E['max_overlap']}")
    print(f"equality_overlap_controls={equality_overlap_checks()}")
    print(f"hostile_component_controls={hostile_family_checks()}")
    print("CONCLUSION=EXACT_ADDRESS_CRITERION;SHARP_DECODER_GAP_CONE;PREDETERMINED_COMPONENT_REFUTED_COFINALLY")


if __name__=="__main__":
    main()
