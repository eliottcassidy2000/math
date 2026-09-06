#!/usr/bin/env python3
"""Exact generalized two-phase packet measure and unbounded residue suppliers.

No carrier/network producer is used. Fraction interval arithmetic verifies
54 literal measure controls; six-point arithmetic checks the unbounded
family symbolically and on named very large rows without interval expansion.
"""

from fractions import Fraction as Q
from hashlib import sha256
from math import gcd

CHECKS = 0


def need(test, payload):
    global CHECKS
    CHECKS += 1
    if not test:
        raise AssertionError(payload)


def distance(x):
    return min(x % 1, (-x) % 1)


def meet(A, B):
    out = []
    i = j = 0
    while i < len(A) and j < len(B):
        lo, hi = max(A[i][0], B[j][0]), min(A[i][1], B[j][1])
        if lo <= hi:
            out.append((lo, hi))
        ar, br = A[i][1], B[j][1]
        i += ar <= br
        j += br <= ar
    return out


def safe_parts(S, radius):
    parts = [(Q(0), Q(1))]
    for s in S:
        parts = meet(parts, [((k + radius) / s, (k + 1 - radius) / s) for k in range(s)])
    return parts


def shifted(parts, t):
    out = []
    for lo, hi in parts:
        lo += t
        hi += t
        if lo >= 1:
            lo -= 1
            hi -= 1
        if hi <= 1:
            out.append((lo, hi))
        else:
            out.extend(((lo, Q(1)), (Q(0), hi - 1)))
    return sorted(out)


def measure(parts):
    return sum((hi - lo for lo, hi in parts), Q(0))


def inside(parts, x):
    return any(lo <= x <= hi for lo, hi in parts)


def dyadic_parameters(T):
    common = gcd(gcd(T[0], T[1]), T[2])
    g = common & -common
    e = sum((w // g) % 2 == 0 for w in T)
    return g, e, Q(1, 2 * g)


def packet(T, y, delta, radius):
    points = tuple(((y - eps * delta + j) / 3) % 1 for eps in range(2) for j in range(3))
    masks = [{k for k, x in enumerate(points) if distance(w * x) < radius} for w in T]
    good = set(range(6)) - set().union(*masks)
    return points, masks, good


def integrate_survivors(A, T, delta, radius):
    ends = {Q(0), Q(1)} | {x for I in A for x in I}
    ends |= {((k + sign * 3 * radius) / w + eps * delta) % 1
             for w in T for k in range(w) for sign in (-1, 1) for eps in range(2)}
    ends = sorted(ends)
    value = Q(0)
    for lo, hi in zip(ends, ends[1:]):
        y = (lo + hi) / 2
        if inside(A, y):
            value += (hi - lo) * len(packet(T, y, delta, radius)[2])
    return value


def main():
    bodies = (
        (1, 2, 3, 4, 5, 7, 8, 9, 11, 13),
        (1, 2, 3, 5, 7, 8, 9, 11, 12, 13),
        (1, 3, 4, 10, 11, 13, 14, 16, 17, 18),
    )
    digest = sha256()
    rows = positive = 0
    for C in bodies:
        for base in ((1, 5, 7), (1, 4, 5), (1, 4, 10)):
            for scale in (1, 2, 4):
                T = tuple(scale * w for w in base)
                g, e, delta = dyadic_parameters(T)
                need(g == scale and e == sum(w % 2 == 0 for w in base), "actual dyadic normalization")
                for radius in (Q(1, 14), Q(1, 12)):
                    G = safe_parts(C, radius)
                    A = meet(G, shifted(G, delta))
                    actual = measure(safe_parts(tuple(3 * c for c in C) + T, radius))
                    lower = Q(3 - e, 6) * measure(A)
                    integral = integrate_survivors(A, T, delta, radius)
                    need((3 - e) * measure(A) <= integral <= 6 * actual,
                         ("independent survivor integral and multiplicity", C, T, radius))
                    need(actual >= lower, ("literal physical measure floor", C, T, radius))
                    positive += lower > 0
                    rows += 1
                    digest.update((repr((C, T, radius, measure(A), lower, integral, actual)) + "\n").encode())
    need(rows == 54, "complete stated measure bank")
    print("LITERAL MEASURE BANK", rows, "POSITIVE CERTIFICATE FLOORS", positive)
    print("MEASURE SEMANTIC SHA256", digest.hexdigest())

    for radius, data in (
        (Q(1, 14), (((1, 5, 7), Q(5, 84)), ((1, 4, 5), Q(83, 420)), ((1, 4, 10), Q(83, 420)))),
        (Q(1, 12), (((1, 5, 7), Q(13, 240)), ((1, 4, 5), Q(23, 120)), ((1, 4, 10), Q(23, 120)))),
    ):
        for base, base_y in data:
            for scale in (1, 2, 4, 8):
                T = tuple(scale * w for w in base)
                g, e, delta = dyadic_parameters(T)
                points, masks, good = packet(T, base_y / scale, delta, radius)
                need([len(mask) for mask in masks] == [gcd(w // g, 6) for w in T], "sharp actual multiplicities")
                need(len(good) == 3 - e, "sharp pointwise survivor count")
            print("SHARP COUNT", radius, base, base_y, "SURVIVORS", 3 - e)

    midpoint_data = ((Q(167, 2002), Q(15, 182)), (Q(695, 2912), Q(33, 364)), (Q(515, 6188), Q(39, 476)))
    for C, (y, alpha) in zip(bodies, midpoint_data):
        need(min(distance(c * z) for c in C for z in (y, y + Q(1, 2))) == alpha, "exact common midpoint margin")
        print("INTERIOR PACKET", y, "COMMON BODY MARGIN", alpha, "UNIFORM FULL-ROW MARGIN", min(alpha, Q(1, 12)))

    step = 480480
    base_body = (1, 3, 5, 7, 8, 9, 11, 12, 13, 11650)
    y = Q(695, 2912)
    alpha = Q(33, 364)
    need(step % 2912 == 0 and (11650 - 2) % 2912 == 0, "packet-preserving residue lift")
    need(step % 14 == 0 and step % 165 == 0, "small-clock divisor pins survive")
    need(alpha - Q(1, 14) == Q(1, 52), "body slack at LRC14")
    need(alpha - Q(1, 12) == Q(2, 273), "body slack at stronger margin")
    need(Q(4, 6) * Q(1, 52) == Q(1, 78), "uniform L14 coefficient")
    need(Q(4, 6) * Q(2, 273) == Q(4, 819), "uniform L12 coefficient")
    lifts = ((0,) * 10, tuple(range(10)), tuple(i * i for i in range(10)), tuple(10 ** i for i in range(10)))
    tails = ((1, 5, 14), (14, 1001, 1004), (14, 1000001, 1000012), (14, 5000003, 7000004))
    for ns, T in zip(lifts, tails):
        C = tuple(b + step * n for b, n in zip(base_body, ns))
        S = tuple(3 * c for c in C) + T
        need(len(set(S)) == 13 and all(w % 3 for w in T) and any(w % 2 for w in T), "typed unbounded row")
        need(gcd(3 * C[0], 14) == 1, "literal primitive full row")
        need(all(any(s % q == 0 for s in S) for q in range(2, 15)), "every clock denominator through14 killed")
        need(min(distance(c * z) for c in C for z in (y, y + Q(1, 2))) == alpha, "exact margin after independent body lifts")
        points, masks, good = packet(T, y, Q(1, 2), Q(1, 12))
        need(len(good) >= 3 - sum(w % 2 == 0 for w in T), "unbounded packet survivor count")
        margins = [(points[k], min(distance(s * points[k]) for s in S)) for k in sorted(good)]
        need(all(margin >= Q(1, 12) for x, margin in margins), "actual large-row stronger witnesses")
        H = max(C)
        need(Q(1, 52 * H) < Q(1, 4) and Q(2, 273 * H) < Q(1, 4), "disjoint Lipschitz body neighborhoods")
        print("UNBOUNDED FAMILY CONTROL", "MAX BODY", H, "T", T,
              "FIRST PHYSICAL", margins[0], "SMALL CLOCKS ALL KILLED", True)

    hostile = (4, 26, 34)
    for radius in (Q(1, 14), Q(1, 12)):
        need(not packet(hostile, y, Q(1, 2), radius)[2], "wrong all-even half-shift kills entire retained packet")
    print("ACTUAL-VALUATION HOSTILE", hostile, "AT", y, "WRONG HALF-SHIFT ALL SIX DEAD")
    print("UNBOUNDED FAMILY FLOORS: M>=1/12; L_1/14>=(3-e)/(78H); L_1/12>=4(3-e)/(819H)")
    print("CHECKS", CHECKS)
    print("PROVED inherited gcd-count transport with actual dyadic phase; FINITE-EXACT controls")
    print("OPEN unrestricted body supplier and chart entry; constructed residue families are not universal")


if __name__ == "__main__":
    main()
