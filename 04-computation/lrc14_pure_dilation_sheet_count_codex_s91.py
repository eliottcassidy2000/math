#!/usr/bin/env python3
"""
lrc14_pure_dilation_sheet_count_codex_s91.py

Exact scout for the pure-dilation hard core

    S(b,V) = {b,2b,...,12b,V},   gcd(b,V)=1,   V == 0 mod 14.

This complements the mod-13 certificate in HYP-2858.  The point here is not
only to prove the AP core again, but to isolate a reusable sheet quotient:
after writing tau=(n+u)/b, the dilated block depends only on u, while each
extra coprime runner can spoil only about 1/7 of the b sheets.
"""

from __future__ import annotations

from fractions import Fraction as F
from math import ceil, gcd

HALF14 = F(1, 14)
DANGER_LEN = F(1, 7)


def fmt(x: F) -> str:
    return f"{x.numerator}/{x.denominator}" if x.denominator != 1 else str(x.numerator)


def norm(x: F) -> F:
    x %= 1
    return min(x, 1 - x)


def in_g12(u: F) -> bool:
    return all(norm(j * u) >= HALF14 for j in range(1, 13))


def g12_intervals() -> list[tuple[F, F]]:
    bps = {F(0), F(1)}
    for j in range(1, 13):
        for k in range(j):
            for sign in (-1, 1):
                bps.add(((F(k) + sign * HALF14) / j) % 1)
    cuts = sorted(bps)
    out: list[tuple[F, F]] = []
    for lo, hi in zip(cuts, cuts[1:]):
        if lo < hi and in_g12((lo + hi) / 2):
            out.append((lo, hi))
    return out


def safe_sheet_count(b: int, V: int, u: F) -> int:
    return sum(norm(F(V) * (F(n) + u) / b) >= HALF14 for n in range(b))


def sheet_breakpoints(b: int, V: int, intervals: list[tuple[F, F]]) -> list[F]:
    bps = {F(0), F(1)}
    for lo, hi in intervals:
        bps.add(lo)
        bps.add(hi)
    # Boundaries where one V-sheet enters/leaves the danger interval.
    for n in range(b):
        for k in range(V + 1):
            for sign in (-1, 1):
                u = F(b) * (F(k) + sign * HALF14) / V - n
                if F(0) < u < F(1):
                    bps.add(u)
    return sorted(bps)


def actual_lonely_measure(b: int, V: int, intervals: list[tuple[F, F]]) -> tuple[F, int, int]:
    cuts = sheet_breakpoints(b, V, intervals)
    total = F(0)
    min_safe: int | None = None
    max_killed = 0
    for lo, hi in zip(cuts, cuts[1:]):
        if lo >= hi:
            continue
        mid = (lo + hi) / 2
        if not in_g12(mid):
            continue
        safe = safe_sheet_count(b, V, mid)
        killed = b - safe
        total += F(safe, b) * (hi - lo)
        min_safe = safe if min_safe is None else min(min_safe, safe)
        max_killed = max(max_killed, killed)
    return total, (min_safe if min_safe is not None else 0), max_killed


def main() -> None:
    intervals = g12_intervals()
    mu = sum(hi - lo for lo, hi in intervals)
    widest = max(intervals, key=lambda ab: ab[1] - ab[0])
    r12 = len(intervals)
    kappa = F(r12, 1) / (6 * mu)

    print("LRC14 pure-dilation sheet-count scout (Codex S91)")
    print("Family S(b,V)={b,2b,...,12b,V}, gcd(b,V)=1, V == 0 mod 14")
    print()
    print("Exact G_12 = {u : ||j u|| >= 1/14 for j=1..12}")
    print(f"  arcs r12       = {r12}")
    print(f"  measure mu12   = {fmt(mu)} = {float(mu):.12f}")
    print(
        "  widest arc     = "
        f"[{fmt(widest[0])}, {fmt(widest[1])}] "
        f"length {fmt(widest[1] - widest[0])}"
    )
    print("  all intervals:")
    for lo, hi in intervals:
        print(f"    [{fmt(lo)}, {fmt(hi)}] length {fmt(hi - lo)}")
    print()

    print("Comb floor")
    print("  L(S) >= (6/7) mu12 - b*r12/(7V)")
    print(f"  positivity threshold V/b > r12/(6 mu12) = {fmt(kappa)} = {float(kappa):.6f}")
    print()

    print("Sheet-count floor")
    print("  For fixed u in G_12, the offsets frac(V*n/b) are a b-point grid.")
    print("  A danger interval of length 1/7 kills at most b/7+1 sheets.")
    print("  Hence safe sheets >= 6b/7 - 1, and")
    print("  L(S) >= mu12 * (6/7 - 1/b).")
    print()

    examples = [(3, 14), (5, 14), (13, 14), (13, 28), (13, 70), (5, 294)]
    print("Representative exact integrations")
    print("    b     V   V/b       actual_L        sheet_LB        comb_LB   min_safe max_killed")
    for b, V in examples:
        assert V % 14 == 0 and gcd(b, V) == 1
        actual, min_safe, max_killed = actual_lonely_measure(b, V, intervals)
        sheet_lb = mu * (F(6, 7) - F(1, b))
        comb_lb = F(6, 7) * mu - F(b * r12, 7 * V)
        print(
            f"  {b:3d} {V:5d} {float(F(V,b)):7.2f} "
            f"{fmt(actual):>14s} {fmt(sheet_lb):>14s} {fmt(comb_lb):>14s} "
            f"{min_safe:9d} {max_killed:10d}"
        )
        assert min_safe >= ceil(F(6 * b, 7) - 1)
        assert actual >= sheet_lb
    print()

    print("Proof connection")
    print("  The mod-13 certificate proves every AP boundary core directly.")
    print("  The sheet quotient is the reusable part: for bE plus h coprime parked")
    print("  runners, the same union bound gives a first floor")
    print("      meas(G_E) * (1 - h/7 - h/b),")
    print("  before any slow-fast V->infinity argument.  This is a finite-V bridge")
    print("  for non-AP coordinated-growth clusters with bounded h.")
    print()

    print("Tournament Analysis")
    print("  vertices: sheet quotient > mod-13 residue > comb floor > raw V/b growth > runner vertices")
    print("  pairwise observable: certified positive lonely measure for S(b,V)")
    print("  challenged assumption: the hard core needs V/b large; sheet counting closes small V/b.")


if __name__ == "__main__":
    main()
