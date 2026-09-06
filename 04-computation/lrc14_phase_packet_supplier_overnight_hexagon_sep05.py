#!/usr/bin/env python3
"""Exact eighth-wave phase-packet controls; no carrier or network producer.

The frozen quantitative-pair module supplies only rational interval utilities.
Every new packet, divisor pin, discrepancy and literal-row test is below.
Universe: declared finite controls, not a census of arbitrary ten-bodies.
"""

from fractions import Fraction as Q
from hashlib import sha256
from math import gcd

from lrc14_quantitative_pair_packet_overnight_hexagon_sep05 import (
    distance, inside, integrate_survivors, measure, meet, packet,
    safe_parts, shifted,
)

CHECKS = 0
DIGEST = sha256()


def need(test, payload):
    global CHECKS
    CHECKS += 1
    if not test:
        raise AssertionError(payload)


def record(*row):
    DIGEST.update((repr(row) + "\n").encode())


def third_safe(C, lam):
    G = safe_parts(C, lam)
    return meet(meet(G, shifted(G, Q(1, 3))), shifted(G, Q(2, 3)))


def nine(T, y, lam):
    X = tuple((y / 3 + Q(k, 9)) % 1 for k in range(9))
    masks = [{k for k, x in enumerate(X) if distance(w*x) < lam} for w in T]
    return X, masks, set(range(9)) - set().union(*masks)


def integral_nine(A, T, lam):
    ends = {Q(0), Q(1)} | {x for I in A for x in I}
    # For 3-unit w the nine physical families yield every residue modulo3w;
    # using only w walls would miss two thirds of the packet phase changes.
    ends |= {((k + sign*9*lam)/(3*w)) % 1
             for w in T for k in range(3*w) for sign in (-1, 1)}
    ends = sorted(ends)
    return sum(((hi-lo)*len(nine(T, (lo+hi)/2, lam)[2])
                for lo, hi in zip(ends, ends[1:]) if inside(A, (lo+hi)/2)), Q(0))


def omega(v):
    return 1 if v % 3 == 0 else 3


def ap8_interval(lam):
    return (2+3*lam)/21, (1-lam)/8


def outlier_budget(v, w, lam=Q(1, 14)):
    lo, hi = ap8_interval(lam)
    W = omega(v)+omega(w)
    return ((hi-lo)*(1-2*lam*W)
            - 2*lam*(1-2*lam*omega(v))/v
            - 2*lam*(1-2*lam*omega(w))/w)


def main():
    C = (1, 2, 3, 5, 7, 8, 9, 11, 12, 13)
    gap_rows = 0
    for T in ((1, 5, 11), (2, 7, 14), (4, 26, 34)):
        for delta in (Q(1, 2), Q(1, 3), Q(4, 7), Q(3, 14)):
            for lam in (Q(1, 14), Q(1, 12)):
                G = safe_parts(C, lam)
                A = meet(G, shifted(G, delta))
                h = sum(distance(w*delta) >= 6*lam for w in T)
                actual = measure(safe_parts(tuple(3*c for c in C)+T, lam))
                integ = integrate_survivors(A, T, delta, lam)
                need(h*measure(A) <= integ <= 6*actual, ("gap integral", T, delta, lam))
                for y in (Q(0), Q(1, 5), Q(13, 72), Q(9, 14)):
                    _, masks, good = packet(T, y, delta, lam)
                    need(all(len(mask) <= (1 if distance(w*delta) >= 6*lam else 2)
                             for w, mask in zip(T, masks)), "gap individual multiplicity")
                    need(len(good) >= h, "gap survivor count")
                gap_rows += 1
                record("gap", T, delta, lam, h, measure(A), integ, actual)
    need(gap_rows == 24, "complete declared gap bank")
    long_body = tuple(range(1, 14))
    y, delta = Q(9, 14), Q(4, 7)
    need(all(distance(c*z) >= Q(1, 14) for c in long_body for z in (y, y-delta)), "singleton gap endpoint")
    need(distance(delta) == Q(3, 7) and packet((1, 5, 11), y, delta, Q(1, 14))[2], "closed band escape")
    for delta in (Q(2, 5), Q(3, 7)-Q(1, 10000)):
        X, masks, _ = packet((1,), delta/2, delta, Q(1, 14))
        need({0, 3} <= masks[0], "below-band individual cap hostile")
    print("GAP BANK", gap_rows, "CLOSED BAND", Q(3, 7), "HOSTILE below-band cap fails")

    tri_rows = 0
    for C in (tuple(range(1, 11)), tuple(range(1, 20, 2)), (1, 3, 4, 7, 9, 12, 15, 30, 33, 39)):
        for T in ((1, 5, 11), (2, 8, 14), (4, 26, 34)):
            for lam in (Q(1, 14), Q(1, 9)):
                A = third_safe(C, lam)
                actual = measure(safe_parts(tuple(3*c for c in C)+T, lam))
                integ = integral_nine(A, T, lam)
                need(3*measure(A) <= integ <= 9*actual, ("nine-grid integral", C, T, lam))
                need(measure(shifted(A, Q(1, 3))) == measure(A), "third periodic measure")
                for y in (Q(0), Q(1, 8), Q(1, 17)):
                    _, masks, good = nine(T, y, lam)
                    need(all(len(mask) <= 2 for mask in masks) and len(good) >= 3, "nine-grid strict occupancy")
                tri_rows += 1
                record("tri", C, T, lam, measure(A), integ, actual)
    need(not third_safe(tuple(range(1, 11)), Q(1, 14)), "AP10 empty three-packet hostile")
    need(measure(third_safe(tuple(range(1, 20, 2)), Q(1, 14))) == Q(86, 2261), "positive but clock6-trivial control")
    need(all(distance(w*Q(r, 3)) == Q(1, 3) < Q(3, 7)
             for w in (1, 5, 11) for r in (1, 2)), "selected three-packet has no gap-band pair edge")
    print("THREE-PHASE BANK", tri_rows, "AP10 EMPTY; ODD BODY positive but clock6-trivial")

    step = 120120
    base = (1, 3, 4, 7, 9, 12, 15, 30, 33, 39)
    phases = (Q(1, 8), Q(11, 24), Q(19, 24))
    need(step % 24 == 0, "phase-preserving modulus")
    need(tuple(min(distance(c*y) for y in phases) for c in base)
         == (Q(1, 8), Q(3, 8), Q(1, 6), Q(1, 8), Q(1, 8), Q(1, 2), Q(1, 8), Q(1, 4), Q(1, 8), Q(1, 8)), "base body margin vector")
    nonprimitive_ns = (9, 10, 2, 12, 13, 6, 16, 15, 8, 11)
    nonprimitive_C = tuple(b+step*n for b, n in zip(base, nonprimitive_ns))
    need(all(c % 17 == 0 for c in nonprimitive_C), "body-gcd17 hostile to full-row/body primitivity conflation")
    need(gcd(3*nonprimitive_C[0], 8) == 1, "same hostile still has primitive full row")
    record("body-gcd hostile", nonprimitive_C)
    X, masks, retained = nine((8, 14), Q(1, 8), Q(1, 8))
    need(masks == [{2, 3, 4}, {6, 8}] and retained == {0, 1, 5, 7}, "literal retained quartet")
    quartet = tuple(X[k] for k in sorted(retained))
    need(quartet == (Q(3, 72), Q(11, 72), Q(43, 72), Q(59, 72)), "quartet physical coordinates")
    accepted = []
    for u in range(72):
        good = [x for x in quartet if distance(u*x) >= Q(1, 8)]
        need(bool(good) == (u % 72 != 0), ("complete72 criterion", u))
        if good:
            accepted.append(u)
        if u % 9:
            need(len(nine((u,), Q(1, 8), Q(1, 8))[1][0]) <= 3, "non9 appender cap")
        else:
            need(len({(u*x) % 1 for x in quartet}) == 1, "9-multiple common phase")
        record("72", u, tuple(good))
    need(len(accepted) == 71, "only zero residue fails")
    need(all(max(distance(u*x) for x in quartet) == Q(1, 8) for u in (27, 45)), "exact endpoint appender residues")
    print("QUARTET", quartet, "ACCEPTED MOD72", len(accepted), "ONLY HOSTILE 0; equality controls27,45")
    for ns, u in (((0,)*10, 2), (tuple(range(10)), 1000001), (tuple(10**i for i in range(10)), 1000012)):
        C = tuple(b+step*n for b, n in zip(base, ns))
        S = tuple(3*c for c in C)+(8, 14, u)
        need(len(set(S)) == 13 and u % 72, "typed distinct constructed row")
        need(gcd(gcd(3*C[0], 8), 14) == 1, "literal full-row primitivity")
        need(all(any(s % q == 0 for s in S) for q in range(2, 15)), "all small clocks killed")
        need(min(distance(c*y) for c in C for y in phases) == Q(1, 8), "independent body lifts preserve margin")
        good = [x for x in quartet if all(distance(s*x) >= Q(1, 8) for s in S)]
        need(bool(good), "literal very-large thirteen-row witness")
        print("UNBOUNDED QUARTET CONTROL", "H", max(C), "u", u, "FIRST", good[0])
        record("large", C, u, tuple(good))

    need(Q(3, 28)/Q(39, 8) == Q(2, 91), "point-bound crossover")
    need(Q(6, 49)-Q(24, 49)/Q(39, 8) == Q(2, 91), "comb-bound crossover")
    need(Q(6)* (Q(1, 8)-Q(1, 14))/3 == Q(3, 28), "unit-nine L14 floor")
    need(Q(6)* (Q(1, 8)-Q(1, 9))/3 == Q(1, 36), "unit-nine L19 floor")
    need(Q(2, 9)*6*(Q(1, 8)-Q(1, 14)) == Q(1, 14), "valuation-one appender L14 floor")
    need(Q(2, 9)*6*(Q(1, 8)-Q(1, 9)) == Q(1, 54), "valuation-one appender L19 floor")
    for u in (1, 2, 5, 13, 27, 40, 71, 73, 81, 117, 190, 191, 200, 1000):
        S = tuple(3*c for c in base)+(8, 14, u)
        if len(set(S)) != 13:
            continue
        actual = measure(safe_parts(S, Q(1, 14)))
        H = max(base)
        point = Q(3, 28*max(3*H, u))
        comb = Q(6, 49*H)-Q(24, 49*u)
        need(actual >= max(point, comb) >= Q(2, 91*H), ("literal uniform quartet floor", u))
        if u % 9:
            k = 3 if u % 3 else 2
            need(actual >= Q(k, 9)*measure(third_safe(base, Q(1, 14))), "typed stronger nine-grid measure")
        record("quartet measure", u, actual, point, comb)
    print("QUARTET UNIFORM ACTUAL L14 >=2/(91H); accepted appender height unrestricted")

    reflo, refhi = ap8_interval(Q(1, 14))
    need((reflo, refhi, refhi-reflo) == (Q(31, 294), Q(13, 112), Q(25, 2352)), "AP8 rational component")
    # Each addressed safe tooth is affine in lambda. Endpoint checks prove
    # the entire closed lambda interval [0,1/9], not a sample-only assertion.
    for c in range(1, 9):
        for r in range(3):
            j = (c*((reflo+refhi)/2+Q(r, 3))).__floor__()
            for lam in (Q(0), Q(1, 9)):
                lo, hi = ap8_interval(lam)
                need(j+lam <= c*(lo+Q(r, 3)) <= c*(hi+Q(r, 3)) <= j+1-lam, "affine AP8 tooth endpoints")
    comb_rows = 0
    for v in (3, 4, 5, 6, 10, 11, 12, 13, 33):
        for lam in (Q(1, 14), Q(1, 13)):
            safe = third_safe((v,), lam)
            for I in ((Q(0), Q(1)), (Q(1, 73), Q(17, 73)), ap8_interval(lam), (Q(0), 2*lam/v)):
                length = I[1]-I[0]
                bad = length-measure(meet([I], safe))
                bound = 2*lam*omega(v)*length + 2*lam*(1-2*lam*omega(v))/v
                need(bad <= bound, ("literal third-orbit comb discrepancy", v, lam, I))
                comb_rows += 1
                record("comb", v, lam, I, bad, bound)
    controls = ((110, 130, Q(1747, 11771760)), (109, 112, Q(73, 1794576)),
                (45, 46, Q(113, 1893360)), (33, 36, Q(29, 60368)))
    for v, w, expected in controls:
        b = outlier_budget(v, w)
        need(b == expected > 0, ("exact outlier budget", v, w))
        C = tuple(range(1, 9))+(v, w)
        residual = measure(meet([ap8_interval(Q(1, 14))], third_safe(C, Q(1, 14))))
        A = third_safe(C, Q(1, 14))
        need(residual >= b and measure(A) >= 3*residual, "three disjoint surviving translates")
        for T in ((1, 5, 14), (2, 8, 14)):
            actual = measure(safe_parts(tuple(3*c for c in C)+T, Q(1, 14)))
            need(actual >= measure(A)/3 >= b, "native full-row outlier floor")
            record("outlier", v, w, T, b, residual, measure(A), actual)
        print("OUTLIERS", v, w, "BUDGET", b)
    L = Q(25, 2352)
    need(Q(8, 49)/(L*Q(1, 7)) == Q(2688, 25), "unit-pair raw height threshold")
    need(Q(10, 49)/(L*Q(3, 7)) == Q(224, 5), "mixed-pair raw height threshold")
    need(Q(12, 49)/(L*Q(5, 7)) == Q(4032, 125), "divisible-pair raw height threshold")
    for m, n in ((1, 1), (2, 5), (1000, 1001), (1000000, 1000001)):
        v, w = 110*m, 130*n
        C = tuple(range(1, 9))+(v, w)
        S = tuple(3*c for c in C)+(1, 5, 14)
        need(v != w and m % 3 and n % 3, "typed unrestricted outlier lift")
        need(all(any(s % q == 0 for s in S) for q in range(2, 15)), "AP8 outlier family kills all small clocks")
        need(outlier_budget(v, w) >= controls[0][2], "absolute tail/body-height independent measure floor")
    lam = Q(1, 13)
    need(ap8_interval(lam)[1]-ap8_interval(lam)[0] == Q(5, 546), "stronger AP8 interval length")
    need(Q(28, 169)/Q(5, 7098) == Q(1176, 5), "stronger unit-outlier height threshold")
    for v, w in ((236, 238), (260, 440)):
        b = outlier_budget(v, w, lam)
        need(b == Q(5, 7098)-Q(14, 169)*(Q(1, v)+Q(1, w)) > 0, "L13 exact positive budget")
        C = tuple(range(1, 9))+(v, w)
        for T in ((1, 5, 14), (2, 8, 14)):
            actual = measure(safe_parts(tuple(3*c for c in C)+T, lam))
            need(actual >= b, "literal stronger full-row measure")
            record("outlier13", v, w, T, b, actual)
    print("COMB BANK", comb_rows, "MIXED HEIGHT BUDGETS109/45/33; UNIT L13 HEIGHT236")
    print("CHECKS", CHECKS)
    print("SEMANTIC SHA256", DIGEST.hexdigest())
    print("PROVED inherited occupancy/discrepancy transport; FINITE-EXACT declared native controls")
    print("OPEN arbitrary primitive clock3 body entry; packet rejection is not an LRC counterexample")


if __name__ == "__main__":
    main()
