#!/usr/bin/env python3
"""
lrc_antipodal_summand_units_s571.py

codex-2026-06-03 S571

Audit the bridge between:
  - summand-graph node C, with incoming pairs {a, C-a};
  - LRC antipodal witnesses at modulus C = 2n - 1;
  - multiplication by units modulo C.

Methodology / Tournament Analysis note:
The natural vertices here are antipodal summand pairs P_a={a,C-a}, not runners.
A pairwise tournament on P_a is not the clean object: multiplication by units is a
group action, so its information is orbit/gcd-stratum data rather than asymmetric
pair preference. A forced orientation by unit visibility or gcd would be a
transitive ledger and would destroy the nonunit-hole predicate. This script
therefore reports the orbit fingerprints that preserve the LRC predicate:
missed unit pair -> S553 witness; missed nonunit pair -> composite-modulus hole.
"""

from __future__ import annotations

from collections import defaultdict
from math import gcd


def factor(n: int) -> str:
    m = n
    out: list[str] = []
    p = 2
    while p * p <= m:
        if m % p == 0:
            e = 0
            while m % p == 0:
                m //= p
                e += 1
            out.append(f"{p}^{e}" if e > 1 else str(p))
        p += 1 if p == 2 else 2
    if m > 1:
        out.append(str(m))
    return " * ".join(out) if out else "1"


def phi(n: int) -> int:
    result = n
    m = n
    p = 2
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            result -= result // p
        p += 1 if p == 2 else 2
    if m > 1:
        result -= result // m
    return result


def pair_reps(C: int) -> list[int]:
    if C % 2:
        return list(range(1, (C + 1) // 2))
    return list(range(1, C // 2))


def pair_label(a: int, C: int) -> str:
    return f"{{{a},{C-a}}}"


def pair_count_summand_node(C: int) -> int:
    return (C - 1) // 2 if C % 2 else C // 2 - 1


def units(C: int) -> list[int]:
    return [u for u in range(1, C) if gcd(u, C) == 1]


def canonical_pair(x: int, C: int) -> int | None:
    r = x % C
    if r == 0:
        return None
    if C % 2 == 0 and r == C // 2:
        return None
    return min(r, (-r) % C)


def unit_orbits(C: int) -> list[set[int]]:
    reps = set(pair_reps(C))
    U = units(C)
    seen: set[int] = set()
    orbits: list[set[int]] = []
    for a in sorted(reps):
        if a in seen:
            continue
        orbit = {canonical_pair(u * a, C) for u in U}
        clean = {x for x in orbit if x is not None}
        seen |= clean
        orbits.append(clean)
    return orbits


def summarize_modulus(C: int) -> None:
    reps = pair_reps(C)
    strata: dict[int, list[int]] = defaultdict(list)
    for a in reps:
        strata[gcd(a, C)].append(a)

    fixed = "yes" if C % 2 == 0 else "no"
    print(f"C={C:<2} factors={factor(C):<8} summand_pairs={pair_count_summand_node(C):<2} "
          f"unit_pairs={phi(C)//2:<2} fixed_midpoint={fixed}")
    for d, vals in sorted(strata.items()):
        kind = "unit-visible" if d == 1 else f"nonunit gcd={d}"
        labels = ", ".join(pair_label(a, C) for a in vals[:8])
        more = "" if len(vals) <= 8 else f", ... (+{len(vals)-8})"
        print(f"  {kind:<18}: {len(vals):<2} {labels}{more}")

    orbit_bits = []
    for orbit in unit_orbits(C):
        ds = sorted({gcd(a, C) for a in orbit})
        orbit_bits.append(f"size={len(orbit)}, gcds={ds}, reps={[pair_label(a, C) for a in sorted(orbit)[:5]]}")
    print("  unit-action orbits:")
    for bit in orbit_bits:
        print(f"    {bit}")
    print()


def classify_speed_set(label: str, speeds: tuple[int, ...], C: int) -> None:
    R = [v % C for v in speeds]
    counts: dict[int, int] = defaultdict(int)
    zeros = 0
    midpoint = 0
    for r in R:
        if r == 0:
            zeros += 1
            continue
        if C % 2 == 0 and r == C // 2:
            midpoint += 1
            continue
        key = canonical_pair(r, C)
        assert key is not None
        counts[key] += 1

    missed = [a for a in pair_reps(C) if counts.get(a, 0) == 0]
    doubled = [a for a in pair_reps(C) if counts.get(a, 0) > 1]
    missed_unit = [a for a in missed if gcd(a, C) == 1]
    missed_nonunit = [a for a in missed if gcd(a, C) != 1]
    print(f"{label}: C={C}, speeds={speeds}")
    print(f"  zeros={zeros} midpoint={midpoint} missed={len(missed)} doubled={len(doubled)}")
    print(f"  missed_unit={', '.join(pair_label(a,C) for a in missed_unit) or '-'}")
    print(f"  missed_nonunit={', '.join(pair_label(a,C) for a in missed_nonunit) or '-'}")
    print(f"  doubled={', '.join(pair_label(a,C) for a in doubled) or '-'}")
    if missed_unit:
        print("  reading: a missed unit pair is invertible, so S553 gives a 2/C witness.")
    elif missed_nonunit:
        print("  reading: only nonunit pairs are missed; this is the composite-modulus hole.")
    else:
        print("  reading: perfect antipodal transversal; residual is flip-set structure.")
    print()


def main() -> None:
    print("Antipodal summand pairs and multiplicative unit visibility (S571)")
    print("=" * 78)
    print("Odd C=2n-1 nodes: summand pairs {a,C-a} are exactly antipodal residue pairs.")
    print("Units modulo C are the multiplicative clocks that can move a missed pair to {+1,-1}.")
    print("Even C nodes have a fixed midpoint C/2, the missing distinct-summand/apex pair.\n")

    print("Modulus fingerprints")
    print("-" * 78)
    for C in (9, 11, 13, 15, 27, 14, 28):
        summarize_modulus(C)

    print("Speed-set fingerprints")
    print("-" * 78)
    classify_speed_set("AP n=14", tuple(range(1, 14)), 27)
    classify_speed_set("V* n=14", (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 24), 27)
    classify_speed_set("n=8 sporadic A", (1, 2, 3, 4, 5, 7, 12), 15)
    classify_speed_set("n=8 sporadic B", (1, 4, 5, 6, 7, 11, 13), 15)

    print("Synthesis")
    print("-" * 78)
    print("Addition supplies the pair shell P_a={a,C-a}.")
    print("Multiplication by units permutes the visible shells and supplies inverse witness times.")
    print("Odd C has no midpoint shell; even C has a fixed midpoint/apex that must be treated separately.")
    print("Composite odd C splits shells by gcd; missed nonunit shells are invisible to the S553 inverse clock.")


if __name__ == "__main__":
    main()
