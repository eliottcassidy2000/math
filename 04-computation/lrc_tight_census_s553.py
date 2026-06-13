#!/usr/bin/env python3
"""
Pinning the LRC tight family at n=14, and stress-testing the S552 spectral gap.
opus-2026-06-01-S553 (remote-control), complementary to oracle-S552.

oracle-S552 proved (exhaustively n<=8): the max-collar M(S)=max_t min_i ||s_i t||
is either 1/n (the "tight family", the minimax extremizers) or >= 2/(2n-1) -- a
SPECTRAL GAP with nothing in the open interval (1/n, 2/(2n-1)).  The gap edge is
achieved by the apex-doubled AP  A_n = {1,...,n-2, 2(n-1)}.

A set V is TIGHT  <=>  M(V) = 1/n exactly  <=>  SAFE_{1/n}(V) = {t: min ||v_i t||
>= 1/n} is nonempty but has measure ZERO (the peak just touches the floor).

This script, all exact (Fraction):
 (1) ENUMERATES the tight family at n=14 (and calibrators 5,6,7; neighbours
     12,15,16) over primitive (n-1)-subsets of [1,B] -- pins what oracle-S552
     could only do for n<=8;
 (2) characterises each tight set (AP? apex-doubled? (q,q)-necklace symmetry?);
 (3) STRESS-TESTS the spectral gap AT n=14: verifies no config in the search has
     M(V) in the open interval (1/14, 2/27);
 (4) confirms M(A_14) = 2/27 for the apex-doubled AP {1,...,12,26}.
"""

from fractions import Fraction
from math import gcd
from itertools import combinations


def safe_set(V, thr):
    """Exact (measure, nonempty) of {t in R/Z : min_i ||v_i t|| >= thr}.
    thr a Fraction. Endpoints where ||v t|| = thr: t = (k*den +/- num)/(v*den)."""
    num, den = thr.numerator, thr.denominator
    endpoints = set()
    for v in V:
        for k in range(v + 1):
            for s in (-1, 1):
                endpoints.add(Fraction(k * den + s * num, v * den) % 1)
    pts = sorted(endpoints)
    m = len(pts)
    meas = Fraction(0)
    nonempty = False
    for i in range(m):
        a = pts[i]
        b = pts[(i + 1) % m]
        length = (b - a) if b > a else (b - a + 1)
        mid = (a + length / 2) % 1
        if all(min((v * mid) % 1, 1 - (v * mid) % 1) >= thr for v in V):
            meas += length
            nonempty = True
    if not nonempty:                       # tight points sit exactly on endpoints
        for t in pts:
            if all(min((v * t) % 1, 1 - (v * t) % 1) >= thr for v in V):
                nonempty = True
                break
    return meas, nonempty


def primitive(V):
    g = 0
    for v in V:
        g = gcd(g, v)
    return tuple(sorted(v // g for v in V))


def classify(V, n):
    """Structural labels for a tight set."""
    m = n - 1
    red = primitive(V)
    labels = []
    if red == tuple(range(1, m + 1)):
        labels.append("AP/regular-polygon")
    # apex-doubled AP: {1,...,m-1, 2m}? (oracle gap edge is a SECOND value, not
    # tight, but check AP-with-one-doubled forms anyway)
    # (q,q) necklace symmetry: rotation-by-2 of the observer order has 2 q-cycles
    # crude proxy: speeds symmetric under v -> (something). Report the reduced set.
    if not labels:
        labels.append("non-AP")
    return red, labels


def census(n, B, label=""):
    m = n - 1
    thr = Fraction(1, n)
    tight = []
    counterex = []
    count = 0
    for combo in combinations(range(1, B + 1), m):
        if primitive(combo) != tuple(combo):
            continue
        count += 1
        meas, ne = safe_set(combo, thr)
        if meas == 0:
            (tight if ne else counterex).append(combo)
    aps = [t for t in tight if primitive(t) == tuple(range(1, m + 1))]
    non_ap = [t for t in tight if t not in aps]
    print(f"== n={n} {label}: primitive {m}-subsets of [1,{B}] ({count} sets) ==")
    print(f"   TIGHT family |M=1/n| size: {len(tight)}   "
          f"(AP: {len(aps)}, non-AP: {len(non_ap)})")
    if counterex:
        print(f"   *** COUNTEREXAMPLES (M<1/n): {counterex}")
    for t in non_ap[:40]:
        red, lab = classify(t, n)
        print(f"       non-AP tight {t}  reduced={red}")
    print()
    return tight, non_ap


def gap_check(n, B, label=""):
    """Stress-test the S552 spectral gap: no config has M in (1/n, 2/(2n-1))."""
    m = n - 1
    floor = Fraction(1, n)
    edge = Fraction(2, 2 * n - 1)
    violations = []
    edge_achievers = 0
    count = 0
    for combo in combinations(range(1, B + 1), m):
        if primitive(combo) != tuple(combo):
            continue
        count += 1
        meas1, ne1 = safe_set(combo, floor)       # M>=1/n iff ne1
        measg, neg = safe_set(combo, edge)         # M>=edge iff neg
        # M in (1/n, 2/(2n-1)) iff  SAFE_{1/n} positive measure (M>1/n)
        #                           AND not neg (M<edge)
        if meas1 > 0 and not neg:
            violations.append((combo, float(meas1)))
        if (not neg) is False and measg == 0:      # M == edge exactly
            edge_achievers += 1
    print(f"== n={n} {label}: SPECTRAL-GAP stress test over {count} configs ==")
    print(f"   floor 1/n = {floor},  gap edge 2/(2n-1) = {edge}")
    print(f"   configs with M in OPEN (1/n, 2/(2n-1)): {len(violations)}  "
          f"(gap holds iff 0)")
    for v in violations[:20]:
        print(f"        VIOLATION: {v}")
    print(f"   configs achieving M = 2/(2n-1) exactly (edge): {edge_achievers}")
    print()
    return violations


def verify_apex_doubled(n):
    """A_n = {1,...,n-2, 2(n-1)} should have M(A_n) = 2/(2n-1) (S552 theorem)."""
    m = n - 1
    A = tuple(sorted(list(range(1, m)) + [2 * m]))   # {1,..,m-1, 2m} = {1..n-2,2(n-1)}
    edge = Fraction(2, 2 * n - 1)
    # M = edge  <=>  nonempty at edge, empty just above
    me, ne = safe_set(A, edge)
    above = edge + Fraction(1, 100 * n * n)
    ma, na = safe_set(A, above)
    print(f"== apex-doubled A_{n} = {A} ==")
    print(f"   SAFE at 2/(2n-1)={edge}: nonempty={ne}, measure={me}")
    print(f"   SAFE just above edge: nonempty={na}  (should be False => M=edge)")
    print(f"   => M(A_{n}) = 2/(2n-1) confirmed: {ne and not na}\n")


if __name__ == "__main__":
    print("######## CALIBRATION (n=5,6 known NON-unique tight; n=7 proven) ########\n", flush=True)
    census(5, 12, "calib")
    census(6, 13, "calib")
    census(7, 13, "calib")

    print("######## TARGET: n=14 = 2*7 doubled prime ########\n", flush=True)
    t14, na14 = census(14, 18, "DOUBLED PRIME")
    verify_apex_doubled(14)
    gap_check(14, 16, "n=14 (bounded stress test)")

    print("######## NEIGHBOURS (doubled-prime comparison) ########\n", flush=True)
    census(12, 15, "4*3")
    census(16, 18, "2^4")
