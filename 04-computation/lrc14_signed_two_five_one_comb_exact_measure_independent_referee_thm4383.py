#!/usr/bin/env python3
"""Clean-room exact audit of the reserved signed-(2,5,1) comb claim.

This file was written without inspecting the prior lrc251 probe.  It uses a
definition-level full-x-circle interval union, not a y-circle determinant
cell decomposition.  All arithmetic is exact (fractions and integers).
"""

from fractions import Fraction as Q
from functools import reduce
from hashlib import sha256
from itertools import permutations, product
from math import gcd


R = Q(3, 14)
TARGET = Q(12, 371)
CHECKS = 0


def check(condition, *context):
    """Optimization-safe assertion: remains active under python -O."""
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(context)


def gcd_many(values):
    return reduce(gcd, values)


def unit6(n):
    return n > 0 and n % 2 == 1 and n % 3 != 0


def pospart(x):
    return max(Q(0), x)


def E(t):
    """Period-three quadrature error, evaluated exactly."""
    rho = t - 3 * (t // 3)
    if rho <= 1:
        return -(rho * rho) / 3
    if rho <= 2:
        return rho - 1 - (rho * rho) / 3
    return 2 * rho - 3 - (rho * rho) / 3


def series_measure(p, q):
    A = Q(3 * abs(q - p), 70)
    B = Q(3 * (p + q), 70)
    total = Q(0)
    for k in range(1, int(B) + 2):
        if k % 3:
            f = Q(5, p * q) * (pospart(B - k) - pospart(A - k))
            total += 2 * f
    return total


def quadrature_measure(p, q):
    A = Q(3 * abs(q - p), 70)
    B = Q(3 * (p + q), 70)
    return Q(6, 245) + Q(10, p * q) * (E(B) - E(A))


def merge(intervals):
    """Union closed representatives; endpoints have measure zero."""
    if not intervals:
        return []
    intervals = sorted(intervals)
    out = [list(intervals[0])]
    for lo, hi in intervals[1:]:
        if lo <= out[-1][1]:
            out[-1][1] = max(out[-1][1], hi)
        else:
            out.append([lo, hi])
    return [(lo, hi) for lo, hi in out]


def circle_piece(center, radius):
    center %= 1
    lo, hi = center - radius, center + radius
    if lo < 0:
        return [(Q(0), hi), (lo + 1, Q(1))]
    if hi > 1:
        return [(Q(0), hi - 1), (lo, Q(1))]
    return [(lo, hi)]


def danger_intervals(w, j):
    radius = Q(1, 14 * w)
    pieces = []
    for n in range(w):
        pieces.extend(circle_piece(Q(n, w) - Q(j, 3), radius))
    return merge(pieces)


def intersect(left, right):
    out = []
    i = j = 0
    while i < len(left) and j < len(right):
        lo = max(left[i][0], right[j][0])
        hi = min(left[i][1], right[j][1])
        if lo < hi:
            out.append((lo, hi))
        if left[i][1] < right[j][1]:
            i += 1
        else:
            j += 1
    return out


def full_x_measure(p, b, q):
    speeds = (p, b, q)
    sheets = {(w, j): danger_intervals(w, j) for w in speeds for j in range(3)}
    pieces = []
    for js in permutations(range(3)):
        cur = [(Q(0), Q(1))]
        for w, j in zip(speeds, js):
            cur = intersect(cur, sheets[(w, j)])
        pieces.extend(cur)
    union = merge(pieces)
    return sum((hi - lo for lo, hi in union), Q(0)), len(union)


def in_intervals(x, intervals):
    left, right = 0, len(intervals)
    while left < right:
        mid = (left + right) // 2
        if intervals[mid][0] < x:
            left = mid + 1
        else:
            right = mid
    i = left - 1
    return i >= 0 and intervals[i][0] < x < intervals[i][1]


def nearest_info(w, y):
    z = w * y
    floor_z = z.numerator // z.denominator
    frac = z - floor_z
    n = floor_z if frac < Q(1, 2) else floor_z + 1
    error = z - n
    if abs(error) >= R:
        return None
    owner = (-pow(w, -1, 3) * n) % 3
    return n, error, owner


def full_x_cell_audit(p, b, q, sigma, tau):
    """Check the owner/defect/determinant iff on every definition-level x-cell."""
    speeds = (p, b, q)
    sheets = {(w, j): danger_intervals(w, j) for w in speeds for j in range(3)}
    walls = {Q(0), Q(1)}
    for intervals in sheets.values():
        for lo, hi in intervals:
            walls.add(lo)
            walls.add(hi)
    walls = sorted(walls)
    measure = Q(0)
    f_cells = 0
    for left, right in zip(walls, walls[1:]):
        if left == right:
            continue
        x = (left + right) / 2
        active = {}
        for w in speeds:
            js = [j for j in range(3) if in_intervals(x, sheets[(w, j)])]
            check(len(js) <= 1, "overlapping physical sheets", w, x, js)
            active[w] = js[0] if js else None

        y = (3 * x) % 1
        lift_shift = (3 * x).numerator // (3 * x).denominator
        info = {w: nearest_info(w, y) for w in speeds}
        for w in speeds:
            check((active[w] is not None) == (info[w] is not None), "eligibility mismatch", w, x)
            if info[w] is not None:
                # If 3x=y+m, the physical sheet is owner-m (mod 3).
                check(active[w] == (info[w][2] - lift_shift) % 3, "owner/lift mismatch", w, x)

        physical = all(active[w] is not None for w in speeds) and len(set(active.values())) == 3
        gate = False
        if info[p] is not None and info[q] is not None:
            np, _, _ = info[p]
            nq, _, _ = info[q]
            numerator = 2 * sigma * np + tau * nq
            # This integrality test is essential when gcd(p,q)=5; K is then
            # automatically divisible by five and is not by itself sufficient.
            if numerator % 5 == 0:
                nb_from_endpoints = numerator // 5
                k_from_endpoints = b * np - p * nb_from_endpoints
                gate = k_from_endpoints % 3 != 0
        check(physical == gate, "physical/gate mismatch", p, b, q, sigma, tau, x, y, active, info, physical, gate)
        if physical:
            np, ep, op = info[p]
            nb, eb, ob = info[b]
            nq, eq, oq = info[q]
            delta = 2 * sigma * np + tau * nq - 5 * nb
            k = b * np - p * nb
            K = q * np - p * nq
            check(delta == 0, "nonzero owner defect", p, b, q, x, delta)
            check(K == 5 * tau * k, "determinant quotient mismatch", p, b, q, x)
            check(k % 3 != 0, "owner gate mismatch", p, b, q, x, k)
            check({op, ob, oq} == {0, 1, 2}, "owner collision", p, b, q, x)
            check(5 * eb == 2 * sigma * ep + tau * eq, "error relation mismatch", p, b, q, x)
            measure += right - left
            f_cells += 1
    return measure, len(walls) - 1, f_cells


def presentations(product_bound):
    out = []
    # p*q < product_bound makes both p and q strictly smaller than it.
    for p in range(1, product_bound):
        if not unit6(p):
            continue
        for q in range(1, product_bound):
            if p * q >= product_bound or not unit6(q):
                continue
            for sigma, tau in product((1, -1), repeat=2):
                numerator = 2 * sigma * p + tau * q
                if numerator <= 0 or numerator % 5:
                    continue
                b = numerator // 5
                if not unit6(b) or len({p, b, q}) != 3:
                    continue
                if gcd_many((p, b, q)) != 1:
                    continue
                out.append((p, b, q, sigma, tau))
    return out


def canonical_vectors(abs_coeffs):
    """All coefficient vectors modulo overall sign, normalized first nonzero >0."""
    ans = set()
    for magnitudes in set(permutations(abs_coeffs)):
        for signs in product((1, -1), repeat=3):
            v = tuple(a * s for a, s in zip(magnitudes, signs))
            if v[0] < 0:
                v = tuple(-x for x in v)
            ans.add(v)
    return sorted(ans)


def cross(u, v):
    return (
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    )


def positive_primitive_ray(u, v):
    z = cross(u, v)
    if 0 in z or not (all(x > 0 for x in z) or all(x < 0 for x in z)):
        return None
    z = tuple(abs(x) for x in z)
    d = gcd_many(z)
    return tuple(x // d for x in z)


def ray_is_target_eligible(z):
    return len(set(z)) == 3 and all(unit6(x) for x in z)


def cross_audit(left_abs, right_abs=None):
    left = canonical_vectors(left_abs)
    right = left if right_abs is None else canonical_vectors(right_abs)
    examined = positive = eligible = 0
    rays = set()
    if right_abs is None:
        pairs = ((left[i], left[j]) for i in range(len(left)) for j in range(i + 1, len(left)))
    else:
        pairs = product(left, right)
    for u, v in pairs:
        z0 = cross(u, v)
        if z0 == (0, 0, 0):
            continue
        examined += 1
        z = positive_primitive_ray(u, v)
        if z is not None:
            positive += 1
            rays.add(z)
            if ray_is_target_eligible(z):
                eligible += 1
    return len(left), len(right), examined, positive, eligible, sorted(rays)


def all_primitive_triples(height):
    vals = [n for n in range(1, height + 1) if unit6(n)]
    for a, b, c in permutations(vals, 3):
        if a < b < c and gcd_many((a, b, c)) == 1:
            yield (a, b, c)


def relation_presentations_of_triple(triple, abs_coeffs):
    vecs = canonical_vectors(abs_coeffs)
    return [v for v in vecs if sum(x * y for x, y in zip(v, triple)) == 0]


def main():
    rows = presentations(425)
    check(len(rows) == len(set(rows)) == 75, "small-product count", len(rows), len(set(rows)))

    sign_counts = {}
    gcd_counts = {}
    component_total = 0
    wall_cell_total = 0
    failure_cell_total = 0
    measured = []
    for p, b, q, sigma, tau in rows:
        check(5 * b == 2 * sigma * p + tau * q, "relation", p, b, q, sigma, tau)
        check(gcd(p, b) == gcd(b, q) == 1, "carrier gcd", p, b, q)
        check(gcd(p, q) in (1, 5), "endpoint gcd", p, b, q)
        sign_counts[(sigma, tau)] = sign_counts.get((sigma, tau), 0) + 1
        gcd_counts[gcd(p, q)] = gcd_counts.get(gcd(p, q), 0) + 1

        ms = series_measure(p, q)
        mq = quadrature_measure(p, q)
        mx, components = full_x_measure(p, b, q)
        mc, wall_cells, failure_cells = full_x_cell_audit(p, b, q, sigma, tau)
        check(ms == mq == mx == mc, "measure mismatch", p, b, q, sigma, tau, ms, mq, mx, mc)
        component_total += components
        wall_cell_total += wall_cells
        failure_cell_total += failure_cells
        measured.append((mx, p, b, q, sigma, tau, components))

    maxima = [row for row in measured if row[0] == max(x[0] for x in measured)]
    check(max(x[0] for x in measured) == TARGET, "wrong maximum")
    check(maxima == [(TARGET, 1, 11, 53, 1, 1, maxima[0][-1])], "nonunique maximum", maxima)

    # The analytic ceiling is already strict at every integer product >=425.
    ceiling_425 = Q(6, 245) + Q(10, 3 * 425)
    check(ceiling_425 < TARGET, "large-product ceiling")

    # Required hostile for the exceptional endpoint-gcd seam.  Endpoint
    # determinant divisibility alone falsely admits this cell: K=-5, but the
    # middle nearest integer numerator is 1 and is not divisible by five.
    hp, hb, hq, hs, ht, hx = 5, 7, 25, 1, 1, Q(13, 1050)
    hy = (3 * hx) % 1
    hip, hib, hiq = nearest_info(hp, hy), nearest_info(hb, hy), nearest_info(hq, hy)
    check(hip is not None and hib is None and hiq is not None, "gcd5 hostile eligibility")
    hnp, hnq = hip[0], hiq[0]
    hK = hq * hnp - hp * hnq
    hnum = 2 * hs * hnp + ht * hnq
    check(hK == -5 and hK % 5 == 0 and hnum == 1 and hnum % 5 != 0, "gcd5 hostile arithmetic")

    # Algebraic residue and exceptional gcd checks for every finite row.
    exceptional = []
    for p, b, q, sigma, tau in rows:
        check(q % 3 == (-sigma * tau * p) % 3, "mod3 q", p, b, q, sigma, tau)
        check(b % 3 == (-sigma * p) % 3, "mod3 b", p, b, q, sigma, tau)
        check(q % 5 == (-2 * sigma * tau * p) % 5, "mod5 q", p, b, q, sigma, tau)
        if gcd(p, q) == 5:
            exceptional.append((p, b, q, sigma, tau, quadrature_measure(p, q)))
            check(p % 5 == q % 5 == 0 and b % 5 != 0, "gcd5 structure", p, b, q)

    # Finite coefficient-vector cross products give scale-independent audits:
    # duplicate 251 presentations, and overlaps with signed 121 / signed 141.
    dup = cross_audit((1, 2, 5))
    ov121 = cross_audit((1, 2, 5), (1, 2, 1))
    ov141 = cross_audit((1, 2, 5), (1, 4, 1))
    check(dup[4] == ov121[4] == ov141[4] == 0, "cross-product overlap", dup[4], ov121[4], ov141[4])

    # A separate bounded hostile census checks that the vector audit agrees
    # with literal relation evaluation, without being used for the proof.
    census_height = 199
    census_triples = census_251 = census_dup = census_121 = census_141 = 0
    for triple in all_primitive_triples(census_height):
        census_triples += 1
        p251 = relation_presentations_of_triple(triple, (1, 2, 5))
        p121 = relation_presentations_of_triple(triple, (1, 2, 1))
        p141 = relation_presentations_of_triple(triple, (1, 4, 1))
        if p251:
            census_251 += 1
            census_dup += len(p251) > 1
            census_121 += bool(p121)
            census_141 += bool(p141)
        check(len(p251) <= 1, "literal duplicate presentation", triple, p251)
        check(not (p251 and p121), "literal 251/121 overlap", triple)
        check(not (p251 and p141), "literal 251/141 overlap", triple)
    check(census_dup == census_121 == census_141 == 0, "bounded hostile overlap")

    print("THM-4383 CLEAN-ROOM EXACT AUDIT")
    print("status=PASS")
    print("scope=primitive distinct positive odd 3-units admitting signed abs-coefficients (2,5,1)")
    print("relation=5*b=2*sigma*p+tau*q, sigma,tau in {+1,-1}")
    print("owner_defect=delta=2*sigma*n_p+tau*n_q-5*n_b=0")
    print("quotient=k=b*n_p-p*n_b; K=q*n_p-p*n_q=5*tau*k; owner_gate=3_not_dividing_k")
    print("A=3*abs(q-p)/70; B=3*(q+p)/70")
    print("measure=6/245+10/(p*q)*(E(B)-E(A))")
    print(f"large_product_ceiling_at_425={ceiling_425}; target={TARGET}; strict={ceiling_425 < TARGET}")
    print(f"small_product_presentations={len(rows)}")
    print("sign_counts=" + ",".join(f"({s:+d},{t:+d}):{sign_counts.get((s,t),0)}" for s,t in product((1,-1), repeat=2)))
    print("endpoint_gcd_counts=" + ",".join(f"{d}:{gcd_counts[d]}" for d in sorted(gcd_counts)))
    print(f"gcd5_presentations={len(exceptional)}")
    for p, b, q, sigma, tau, mu in exceptional:
        print(f"gcd5_row=({p},{b},{q},{sigma:+d},{tau:+d}),mu={mu}")
    print("gcd5_hostile=(p,b,q,x)=(5,7,25,13/1050);K=-5;middle_numerator=1;endpoint_K_gate=FalsePositive")
    print(f"full_x_total_components={component_total}")
    print(f"full_x_wall_cells={wall_cell_total}; failure_cells={failure_cell_total}")
    print(f"small_product_max={TARGET}; unique=(p,b,q,sigma,tau)=(1,11,53,+1,+1)")
    print(f"cross_duplicate_251=vectors:{dup[0]},pairs:{dup[2]},positive_rays:{dup[3]},unique_positive_rays:{len(dup[5])},eligible_rays:{dup[4]}")
    print(f"cross_251x121=vectors:{ov121[0]}x{ov121[1]},pairs:{ov121[2]},positive_rays:{ov121[3]},unique_positive_rays:{len(ov121[5])},eligible_rays:{ov121[4]}")
    print(f"cross_251x141=vectors:{ov141[0]}x{ov141[1]},pairs:{ov141[2]},positive_rays:{ov141[3]},unique_positive_rays:{len(ov141[5])},eligible_rays:{ov141[4]}")
    print(f"bounded_hostile_height={census_height}; triples={census_triples}; signed251={census_251}; duplicate={census_dup}; overlap121={census_121}; overlap141={census_141}")
    print(f"optimization_safe_runtime_checks={CHECKS}")


if __name__ == "__main__":
    main()
