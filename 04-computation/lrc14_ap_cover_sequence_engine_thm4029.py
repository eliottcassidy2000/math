#!/usr/bin/env python3
"""Independent exact audit of the THM-536 AP-cover threshold sequence.

Two semantically different engines are used:

1. direct x-space endpoint refinement at r/(7e);
2. Farey/Sturmian alpha-space cells for theta=j+alpha.

The second engine also retains the cover time and an explicit Farey interval
witness, rather than only outputting the scalar measure prefix.
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction as Q
from itertools import combinations
from math import floor


def covered_at_x(x: Q, m: int) -> bool:
    mask = 0
    for e in range(m):
        mask |= 1 << (floor(7 * e * x) % 7)
    return mask == 0x7F


def a_direct(m: int) -> Q:
    """Exact x-space measure of S_7({0,...,m-1})."""
    bps = {Q(0), Q(1)}
    for e in range(1, m):
        for r in range(7 * e + 1):
            bps.add(Q(r, 7 * e))
    bps = sorted(bps)
    total = Q(0)
    for left, right in zip(bps, bps[1:]):
        mid = (left + right) / 2
        if covered_at_x(mid, m):
            total += right - left
    return total


def farey_neighbors(n: int):
    """Consecutive pairs in the Farey sequence F_n on [0,1]."""
    a, b, c, d = 0, 1, 1, n
    previous = (a, b)
    while True:
        current = (c, d)
        yield previous, current
        if c == 1 and d == 1:
            break
        k = (n + b) // d
        a, b, c, d = c, d, k * c - a, k * d - b
        previous = current


def farey_cover_distribution(max_m: int):
    """Return exact first-cover-time masses and one Farey-cell witness each.

    On an open Farey cell of order max_m-1, floor(e*alpha) is constant for
    every 0 <= e < max_m. For theta=j+alpha, the visited residue is
    e*j+floor(e*alpha) (mod 7).
    """
    mass = defaultdict(Q)
    witness = {}
    unresolved = Q(0)
    n = max_m - 1
    for (a, b), (c, d) in farey_neighbors(n):
        alpha = Q(a + c, b + d)  # strict interior mediant
        width = Q(1, b * d)
        for j in range(7):
            mask = 0
            first = None
            trace = []
            for e in range(max_m):
                residue = (e * j + floor(e * alpha)) % 7
                mask |= 1 << residue
                trace.append(residue)
                if mask == 0x7F:
                    first = e + 1
                    break
            contribution = width / 7
            if first is None:
                unresolved += contribution
            else:
                mass[first] += contribution
                witness.setdefault(first, (j, Q(a, b), Q(c, d), tuple(trace)))
    if sum(mass.values(), Q(0)) + unresolved != 1:
        raise RuntimeError("Farey-cell mass partition")
    return mass, unresolved, witness


def good_measure(P) -> Q:
    """Exact measure of {x: ||p x|| >= 1/14 for all p in P}."""
    P = tuple(sorted(P))
    if not P:
        return Q(1)
    bps = {Q(0), Q(1)}
    for p in P:
        for r in range(p + 1):
            for sign in (-1, 1):
                x = (Q(r) + sign * Q(1, 14)) / p
                if 0 <= x <= 1:
                    bps.add(x)
    bps = sorted(bps)
    total = Q(0)
    for left, right in zip(bps, bps[1:]):
        mid = (left + right) / 2
        if all(min((p * mid) % 1, 1 - ((p * mid) % 1)) >= Q(1, 14) for p in P):
            total += right - left
    return total


def exact_caps():
    caps = {}
    minimizers = {}
    speeds = tuple(range(1, 14))
    for k in range(8, 14):
        values = [(good_measure(P), P) for P in combinations(speeds, 13 - k)]
        cap = min(v for v, _ in values)
        caps[k] = cap
        minimizers[k] = tuple(P for v, P in values if v == cap)
    return caps, minimizers


def main():
    max_m = 40
    mass, unresolved, witness = farey_cover_distribution(max_m)
    cumulative = Q(0)
    farey_a = {}
    for m in range(1, max_m + 1):
        cumulative += mass[m]
        farey_a[m] = cumulative

    print("EXACT AP COVER SEQUENCE a(m)=meas S7(AP_m)")
    print("m  a(m)                    increment")
    for m in range(1, max_m + 1):
        print(f"{m:2d} {str(farey_a[m]):>24s} {str(mass[m]):>24s}")
    print(f"unresolved through m={max_m}: {unresolved}")

    print("\nINDEPENDENT DIRECT-ENDPOINT CROSS-CHECK m=1..30")
    for m in range(1, 31):
        direct = a_direct(m)
        print(f"m={m:2d} direct={direct} farey={farey_a[m]} equal={direct == farey_a[m]}")
        if direct != farey_a[m]:
            raise RuntimeError((m, "direct/Farey mismatch"))

    caps, minimizers = exact_caps()
    print("\nEXHAUSTIVE TRUE CAPS")
    for k in range(8, 14):
        print(f"k={k}: cap={caps[k]}, minimizers={minimizers[k]}")

    print("\nFIRST CROSSINGS AND SEMANTIC WITNESSES")
    thresholds = {}
    for k in range(8, 14):
        cap = caps[k]
        if cap == 1:
            thresholds[k] = None
            print(f"k={k}: cap=1, a(m)<1 for every finite m (alpha in (0,1/(m-1)) misses residues), N*=infinity")
            continue
        first_above = next(m for m in range(1, max_m + 1) if farey_a[m] > cap)
        Nstar = first_above - 2  # N uses AP_{N+1}
        thresholds[k] = Nstar
        j, left, right, trace = witness[first_above]
        print(
            f"k={k}: a({first_above-1})={farey_a[first_above-1]} <= cap={cap} "
            f"< a({first_above})={farey_a[first_above]}; N*={Nstar}"
        )
        print(
            f"  positive jump witness: theta={j}+alpha, alpha in ({left},{right}), "
            f"first-cover trace={trace}"
        )

    if tuple(thresholds[k] for k in range(8, 13)) != (7, 8, 10, 13, 26):
        raise RuntimeError("finite threshold row")
    if thresholds[13] is not None:
        raise RuntimeError("cap-one threshold")
    print("\nAUDITED N*(k) = (7,8,10,13,26,infinity).")


if __name__ == "__main__":
    main()
