#!/usr/bin/env python3
"""Symbolic finite-state tail for all small primitive H two-star channels.

For every ordered label edge, every coprime P<Q<=26 with Q<=2P and
P+Q>=8, and every dilation g making both levels at least six, this proves
the cell-90 high-channel floor.  The sole weaker lane is edge labels (1,6)
on primitive channel 3:5, where the floor is 2030/280393.

The proof is infinite, not a large finite sample.  On a ray p=gP,q=gQ,
write k=r+tP and l=s+tQ.  Then the center determinant is A(g)+tC with
C=Qe-Pf fixed.  On each residue class of g modulo |C|, every floor endpoint
is affine in the new ray parameter.  Once the active max/min and absolute-
value branches stabilize, the cleared floor margin is an exact quadratic.
This referee proves branch stability from slopes, checks the finite head,
constructs every quadratic, and proves its integer tail is nonnegative.
"""

from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as F
from hashlib import sha256
from math import gcd
from pathlib import Path

L = 168
CELL = 90
H = (1, 2, 3, 4, 6, 12)
EDGES = ((0, 1), (0, 2), (0, 3), (0, 4), (0, 5),
         (1, 2), (1, 3), (1, 4), (1, 5))
MAX_Q = 26
WORKERS = 6


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def ceildiv(a, b):
    return -((-a) // b)


def floor_moments(n, m, a, b):
    if n == 0:
        return 0, 0, 0
    s1 = n * (n - 1) // 2
    s2 = n * (n - 1) * (2 * n - 1) // 6
    qa, a0 = divmod(a, m)
    qb, b0 = divmod(b, m)
    base0 = qa * s1 + qb * n
    base1 = qa * s2 + qb * s1
    base2 = qa * qa * s2 + 2 * qa * qb * s1 + qb * qb * n
    if a0 == 0:
        return base0, base1, base2
    height = (a0 * (n - 1) + b0) // m
    if height == 0:
        return base0, base1, base2
    u0, u1, u2 = floor_moments(height, a0, m, m - b0 + a0 - 1)
    r0 = n * height - u0
    r1 = height * s1 - (u2 - u0) // 2
    r2 = n * height * height - 2 * u1 - u0
    return (base0 + r0, base1 + r1,
            base2 + 2 * qa * r1 + 2 * qb * r0 + r2)


def residue_prefix(n, m, a, b, threshold, base):
    shifted = floor_moments(n, m, a, b + m - threshold)
    d0 = shifted[0] - base[0]
    d1 = shifted[1] - base[1]
    y0d = (shifted[2] - base[2] - d0) // 2
    high_sum = a * d1 + b * d0 - m * y0d
    total = a * n * (n - 1) // 2 + b * n - m * base[0]
    return n - d0, total - high_sum


def triangle_sum(n, m, a, b, peak, base, total):
    if peak <= 0:
        return 0
    radius = (peak - 1) // L
    require(2 * radius < m, ("overlapping triangle tails", n, m, peak))
    low_count, low_sum = residue_prefix(n, m, a, b, radius + 1, base)
    before_count, before_sum = residue_prefix(n, m, a, b, m - radius, base)
    high_count, high_sum = n - before_count, total - before_sum
    return (peak * low_count - L * low_sum
            + (peak - L * m) * high_count + L * high_sum)


def fast_mass(e, p, f, q):
    if p > q:
        return fast_mass(f, q, e, p)
    z, w = L * p - e, L * q - f
    r, s = e * CELL % L, f * CELL % L
    determinant = r * w - s * z
    require(determinant % L == 0, ("nonintegral phase", e, p, f, q))
    b = (determinant // L) % z
    a = w % z
    base = floor_moments(p, z, a, b)
    total = a * p * (p - 1) // 2 + b * p - z * base[0]
    numerator = (
        triangle_sum(p, z, a, b, 12 * (z + w), base, total)
        - triangle_sum(p, z, a, b, 12 * (w - z), base, total)
    )
    return numerator, z * w


def term_candidates(P, Q, e, f, residue_index, shift_index, kernel, g):
    z, w = L * g * P - e, L * g * Q - f
    r0, s0 = e * CELL % L, f * CELL % L
    base = (r0 * w - s0 * z) // L
    cross = P * w - Q * z
    lo = max(0, ceildiv(-shift_index, Q))
    hi = min(g - 1, (g * Q - 1 - shift_index) // Q)
    affine = base + residue_index * w - shift_index * z
    peak = 12 * ((z + w) if kernel == 0 else (w - z))
    radius = (peak - 1) // L
    if cross < 0:
        affine, cross = -affine, -cross
    if cross:
        positive_lo = ceildiv(-radius - affine, cross)
        positive_hi = (radius - affine) // cross
        negative_hi = ceildiv(-affine, cross) - 1
        return (lo, positive_lo, hi, positive_hi,
                negative_hi, negative_hi + 1), cross
    return (lo, hi, affine, peak - L * abs(affine)), 0


def stable_term(P, Q, e, f, residue_index, shift_index, kernel,
                residue, period, h):
    here, cross = term_candidates(
        P, Q, e, f, residue_index, shift_index, kernel, residue + period * h
    )
    nxt, cross2 = term_candidates(
        P, Q, e, f, residue_index, shift_index, kernel,
        residue + period * (h + 1)
    )
    require(cross == cross2, ("cross instability", P, Q, e, f, residue,
                              period, h))
    slopes = tuple(y - x for x, y in zip(here, nxt))
    # Every endpoint is genuinely affine on the fixed residue class.
    nxt2, _ = term_candidates(
        P, Q, e, f, residue_index, shift_index, kernel,
        residue + period * (h + 2)
    )
    require(tuple(y - x for x, y in zip(nxt, nxt2)) == slopes,
            ("endpoint nonaffinity", P, Q, e, f, residue, period, h))

    if cross:
        lo, positive_lo, hi, positive_hi, negative_hi, _ = here
        slo, splo, shi, sphi, sneg, _ = slopes
        if positive_lo > lo or (positive_lo == lo and splo > slo):
            left, sleft = positive_lo, splo
        else:
            left, sleft = lo, slo
        if not (sleft >= slo and sleft >= splo):
            return False
        if positive_hi < hi or (positive_hi == hi and sphi < shi):
            right, sright = positive_hi, sphi
        else:
            right, sright = hi, shi
        if not (sright <= shi and sright <= sphi):
            return False
        if left > right:
            return sleft >= sright
        if negative_hi < left:
            return sneg <= sleft
        if negative_hi >= right:
            return sneg >= sright
        return sleft <= sneg <= sright

    lo, hi, affine, value = here
    slo, shi, saffine, svalue = slopes
    if hi < lo or shi < slo:
        return False
    if affine > 0 and saffine < 0:
        return False
    if affine < 0 and saffine > 0:
        return False
    if affine == 0 and saffine != 0:
        return False
    if value > 0:
        return svalue >= 0
    if value < 0:
        return svalue <= 0
    return svalue == 0


def ray_is_stable(P, Q, e, f, residue, period, h):
    for residue_index in range(P):
        # The relative index differs from rQ/P by less than two; this
        # deliberately oversized bank retains hostile boundary candidates.
        for shift_index in range(-2, Q + 3):
            for kernel in (0, 1):
                if not stable_term(P, Q, e, f, residue_index, shift_index,
                                   kernel, residue, period, h):
                    return False
    return True


def cleared_margin(edge, P, Q, orientation, g):
    i, j = EDGES[edge]
    e, f = ((H[i], H[j]), (H[j], H[i]))[orientation]
    numerator, denominator = fast_mass(e, g * P, f, g * Q)
    if edge == 3 and (P, Q) == (3, 5):
        return 280393 * numerator - 2030 * denominator
    return 105 * numerator - denominator


def quadratic_minimum(v0, v1, v2):
    second = v2 - 2 * v1 + v0
    require(second % 2 == 0, ("nonintegral quadratic", v0, v1, v2))
    a = second // 2
    b = v1 - v0 - a
    require(a >= 0, ("concave tail", v0, v1, v2))
    candidates = {0}
    if a == 0:
        require(b >= 0, ("decreasing linear tail", v0, v1, v2))
    elif b < 0:
        center = (-b) // (2 * a)
        candidates.update((max(0, center - 1), center, center + 1, center + 2))
    minimum = min(a * n * n + b * n + v0 for n in candidates)
    return minimum, (a, b, v0)


def prove_ray(task):
    edge, P, Q, orientation = task
    i, j = EDGES[edge]
    e, f = ((H[i], H[j]), (H[j], H[i]))[orientation]
    cross = abs(Q * e - P * f)
    period = cross or 1
    g_min = max(1, ceildiv(6, P))
    head_checks = tail_classes = 0
    maximum_stable_g = 0
    minimum_cleared = None
    polynomial_digest = sha256()

    for residue in range(1, period + 1):
        h_min = max(0, ceildiv(g_min - residue, period))
        h = max(1, h_min)
        while not ray_is_stable(P, Q, e, f, residue, period, h):
            h *= 2
            require(h <= 1 << 20,
                    ("stabilization ceiling", edge, P, Q, orientation,
                     residue, period))
        stable_g = residue + period * h
        maximum_stable_g = max(maximum_stable_g, stable_g)

        for head_h in range(h_min, h):
            value = cleared_margin(edge, P, Q, orientation,
                                   residue + period * head_h)
            require(value >= 0,
                    ("negative finite head", edge, P, Q, orientation,
                     residue, period, head_h, value))
            minimum_cleared = value if minimum_cleared is None else min(minimum_cleared, value)
            head_checks += 1

        values = tuple(cleared_margin(edge, P, Q, orientation,
                                      residue + period * (h + k))
                       for k in range(4))
        minimum, polynomial = quadratic_minimum(*values[:3])
        a, b, c = polynomial
        require(values[3] == 9 * a + 3 * b + c,
                ("fourth-point failure", edge, P, Q, orientation,
                 residue, period, h, values, polynomial))
        require(minimum >= 0,
                ("negative quadratic tail", edge, P, Q, orientation,
                 residue, period, h, polynomial, minimum))
        minimum_cleared = minimum if minimum_cleared is None else min(minimum_cleared, minimum)
        polynomial_digest.update(repr((residue, h, polynomial, minimum)).encode())
        tail_classes += 1

    return (task, period, head_checks, tail_classes, maximum_stable_g,
            minimum_cleared, polynomial_digest.hexdigest())


def main():
    channels = tuple(
        (P, Q)
        for Q in range(2, MAX_Q + 1)
        for P in range((Q + 1) // 2, Q)
        if gcd(P, Q) == 1 and P + Q >= 8
    )
    tasks = tuple(
        (edge, P, Q, orientation)
        for edge in range(len(EDGES))
        for P, Q in channels
        for orientation in (0, 1)
    )
    with ProcessPoolExecutor(max_workers=WORKERS) as pool:
        rows = tuple(pool.map(prove_ray, tasks, chunksize=1))

    head_checks = sum(row[2] for row in rows)
    tail_classes = sum(row[3] for row in rows)
    maximum_stable_g = max(row[4] for row in rows)
    require(all(row[5] >= 0 for row in rows), "negative ray summary")
    sharp = tuple(row for row in rows if row[5] == 0)
    require(sharp and all(row[0][0] == 3 and row[0][1:3] == (3, 5)
                          for row in sharp),
            ("sharp-class mismatch", sharp))
    semantic = repr(rows).encode()
    print("LRC H SMALL HIGH-CHANNEL RAYS -- SYMBOLIC FINITE-STATE TAIL")
    print(f"primitive_channels={len(channels)};ordered_edge_rays={len(tasks)};Q<=26")
    print(f"finite_head_checks={head_checks};affine_residue_tail_classes={tail_classes};"
          f"maximum_stable_g={maximum_stable_g}")
    print(f"sharp_classes={len(sharp)};sharp_floor=2030/280393;"
          f"sharp_channel=edge(1,6),primitive(3,5)")
    print(f"semantic_sha256={sha256(semantic).hexdigest()}")


if __name__ == "__main__":
    main()
