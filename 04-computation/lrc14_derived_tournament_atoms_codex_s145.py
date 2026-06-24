#!/usr/bin/env python3
"""Derived tournament fingerprints for AP/Goddyn-Wong LRC14 atoms.

This is a small, reproducible scout for the poke-forum post
20260624-084020Z-derived-tournament-ap-gw-conditions.

Tournament Analysis declarations
--------------------------------
Runner-pressure tournament:
    vertices: runners s in S
    observable: endpoint pressure profile over U_14
    switch: lexicographic profile comparison
    tie Hamiltonian path: 1,13,3,11,5,9,7,2,12,4,10,6,8,0 then speed

Residue-gap tournament:
    vertices: residues in Z/14
    observable: multiplicity, neighbor multiplicity, unit/apex flags,
                distance from observer hole
    switch/tie: lexicographic comparison and the same Hamiltonian path

Acceleration-gate tournament:
    vertices: acceleration sites v = 1..13
    observable: Jacobsthal admissibility of W(v)=[14-v,27-2v],
                prime coverage of the window, negative window length
    switch: lexicographic comparison

Node-squared local tournament:
    vertices: ordered pairs (outer runner, inner runner)
    observable: for each outer node, the inner distance order at the first
                unit time where the outer runner is most boundary-critical
    switch/tie: distance comparison and runner-pressure tie path

The script does not attempt a full graph-isomorphism canonization.  It emits
stable fingerprints sufficient for the finite family tested here: AP, GW, all
single AP accelerations, all double AP accelerations, and the AP-residue false
look-alike that loses the multiple-of-12 covering debt.  Exploratory
accelerations that duplicate an existing speed are retained but marked
distinct_speeds=False.
"""

from __future__ import annotations

import hashlib
from collections import Counter, defaultdict
from itertools import combinations
from math import gcd


N = 14
U = (1, 3, 5, 9, 11, 13)
RESIDUE_TIE_PATH = (1, 13, 3, 11, 5, 9, 7, 2, 12, 4, 10, 6, 8, 0)
RESIDUE_RANK = {r: i for i, r in enumerate(RESIDUE_TIE_PATH)}
AP = tuple(range(1, 14))
GW = tuple(list(range(1, 12)) + [13, 24])
FALSE_AP_RESIDUE = tuple(list(range(1, 12)) + [13, 26])


def clock_distance_num(speed: int, unit: int) -> int:
    """Return 14 * ||speed * unit / 14|| as an integer in 0..7."""

    r = (speed * unit) % N
    return min(r, N - r)


def pressure_core(speed: int) -> tuple[int, ...]:
    counts = Counter(clock_distance_num(speed, u) for u in U)
    return tuple(counts[j] for j in (1, 3, 5, 7, 2, 4, 6))


def pressure_profile(speed: int) -> tuple[int, ...]:
    return pressure_core(speed) + (gcd(speed, N), speed)


def runner_tie_key(speed: int) -> tuple[int, int]:
    return (RESIDUE_RANK[speed % N], speed)


def ordered_by_pressure(row: tuple[int, ...]) -> list[int]:
    return sorted(row, key=lambda s: (pressure_profile(s), tuple(-x for x in runner_tie_key(s))), reverse=True)


def residue_multiplicities(row: tuple[int, ...]) -> tuple[int, ...]:
    c = Counter(s % N for s in row)
    return tuple(c[r] for r in range(N))


def max_empty_gap(row: tuple[int, ...]) -> int:
    occupied = sorted({s % N for s in row})
    if not occupied:
        return N
    gaps = []
    for a, b in zip(occupied, occupied[1:] + [occupied[0] + N]):
        gaps.append(b - a)
    return max(gaps)


def covering_ok(row: tuple[int, ...]) -> bool:
    return all(any(s % d == 0 for s in row) for d in range(1, 14)) and all(s % N != 0 for s in row)


def distinct_speeds(row: tuple[int, ...]) -> bool:
    return len(set(row)) == len(row)


def full_unit_shell(row: tuple[int, ...]) -> bool:
    residues = {s % N for s in row}
    return set(U).issubset(residues)


def has_apex(row: tuple[int, ...]) -> bool:
    return any(s % N == 7 for s in row)


def endpoint_units_ok(row: tuple[int, ...]) -> bool:
    return all(min(clock_distance_num(s, u) for s in row) == 1 for u in U)


def jacobsthal_window(v: int) -> tuple[int, int]:
    return (N - v, 2 * N - 1 - 2 * v)


def gate_admissible(v: int) -> bool:
    lo, hi = jacobsthal_window(v)
    if lo > hi:
        return False
    return all(gcd(k, v) > 1 for k in range(lo, hi + 1))


def prime_divisors(n: int) -> set[int]:
    out = set()
    d = 2
    x = n
    while d * d <= x:
        while x % d == 0:
            out.add(d)
            x //= d
        d += 1
    if x > 1:
        out.add(x)
    return out


def gate_profile(v: int) -> tuple[int, int, int, int]:
    lo, hi = jacobsthal_window(v)
    if lo > hi:
        window = []
    else:
        window = list(range(lo, hi + 1))
    primes = prime_divisors(v)
    hit = {p for p in primes if any(k % p == 0 for k in window)}
    return (int(gate_admissible(v)), len(hit), -len(window), v)


def gate_order() -> list[int]:
    return sorted(range(1, 14), key=gate_profile, reverse=True)


def residue_gap_profiles(row: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    m = residue_multiplicities(row)
    profiles = []
    for r in range(N):
        dist0 = min(r, N - r)
        profiles.append((m[r], m[(r - 1) % N] + m[(r + 1) % N], int(r in U), int(r == 7), dist0))
    return tuple(profiles)


def residue_gap_order(row: tuple[int, ...]) -> list[int]:
    profiles = residue_gap_profiles(row)
    return sorted(range(N), key=lambda r: (profiles[r], -RESIDUE_RANK[r]), reverse=True)


def first_best_unit(speed: int) -> int:
    values = [(clock_distance_num(speed, u), RESIDUE_RANK[(speed * u) % N], u) for u in U]
    return min(values)[2]


def inner_score_sequence(row: tuple[int, ...], outer: int) -> tuple[int, ...]:
    u = first_best_unit(outer)
    scores = [0 for _ in row]
    for i, j in combinations(range(len(row)), 2):
        a = row[i]
        b = row[j]
        da = clock_distance_num(a, u)
        db = clock_distance_num(b, u)
        if da < db:
            scores[i] += 1
        elif db < da:
            scores[j] += 1
        else:
            winner = min((i, j), key=lambda k: runner_tie_key(row[k]) + (k,))
            scores[winner] += 1
    return tuple(sorted(scores, reverse=True))


def relative_profile(outer: int, inner: int) -> tuple[tuple[tuple[int, int, int], int], ...]:
    """Profile of an inner runner as seen from an outer runner over unit times."""

    counts = Counter(
        (
            clock_distance_num(outer, u),
            clock_distance_num(inner, u),
            clock_distance_num(inner - outer, u),
        )
        for u in U
    )
    return tuple(sorted(counts.items()))


def inner_profile_class(row: tuple[int, ...], outer: int) -> tuple[tuple[tuple[tuple[int, int, int], int], int], ...]:
    """Full node-squared local profile class at one outer runner."""

    return tuple(sorted(Counter(relative_profile(outer, inner) for inner in row).items()))


def short_hash(obj: object) -> str:
    return hashlib.sha1(repr(obj).encode("utf-8")).hexdigest()[:12]


def square_inner_histogram(row: tuple[int, ...]) -> tuple[tuple[tuple[int, ...], int], ...]:
    hist = Counter(inner_score_sequence(row, outer) for outer in row)
    return tuple(sorted(hist.items(), key=lambda item: (item[0], item[1])))


def square_profile_hash_histogram(row: tuple[int, ...]) -> tuple[tuple[str, int], ...]:
    hist = Counter(short_hash(inner_profile_class(row, outer)) for outer in row)
    return tuple(sorted(hist.items()))


def accelerate(row: tuple[int, ...], *sites: int) -> tuple[int, ...]:
    values = list(row)
    for v in sites:
        values.remove(v)
    values.extend(2 * v for v in sites)
    return tuple(sorted(values))


def row_summary(name: str, row: tuple[int, ...]) -> dict[str, object]:
    shell_hist = Counter(pressure_core(s) for s in row)
    return {
        "name": name,
        "row": row,
        "distinct_speeds": distinct_speeds(row),
        "covering_ok": covering_ok(row),
        "endpoint_units_ok": endpoint_units_ok(row),
        "full_unit_shell": full_unit_shell(row),
        "has_7_apex": has_apex(row),
        "max_empty_gap": max_empty_gap(row),
        "pressure_shell_hist": tuple(sorted(shell_hist.items())),
        "residue_multiplicities": residue_multiplicities(row),
        "pressure_order": tuple(ordered_by_pressure(row)),
        "gap_order": tuple(residue_gap_order(row)),
        "square_inner_histogram": square_inner_histogram(row),
        "square_profile_hash_histogram": square_profile_hash_histogram(row),
    }


def print_named_summary(summary: dict[str, object]) -> None:
    print(f"\n[{summary['name']}]")
    print(f"row={summary['row']}")
    for key in (
        "distinct_speeds",
        "covering_ok",
        "endpoint_units_ok",
        "full_unit_shell",
        "has_7_apex",
        "max_empty_gap",
    ):
        print(f"{key}={summary[key]}")
    print(f"residue_multiplicities={summary['residue_multiplicities']}")
    print(f"pressure_shell_hist={summary['pressure_shell_hist']}")
    print(f"square_inner_histogram={summary['square_inner_histogram']}")
    print(f"square_profile_hash_histogram={summary['square_profile_hash_histogram']}")


def main() -> None:
    print("LRC14 derived tournament atom scout")
    print("=" * 40)
    print("Gate profiles (v, W(v), profile, admissible):")
    for v in range(1, 14):
        print(f"  v={v:2d} W={jacobsthal_window(v)} profile={gate_profile(v)} admissible={gate_admissible(v)}")
    print(f"Gate tournament order={tuple(gate_order())}")
    print(f"Admissible sites={[v for v in range(1, 14) if gate_admissible(v)]}")

    named_rows = {
        "AP": AP,
        "GW_12_to_24": GW,
        "false_AP_residue_12_to_26": FALSE_AP_RESIDUE,
    }
    named_rows["near_12_to_36"] = tuple(sorted([s for s in AP if s != 12] + [36]))

    for name, row in named_rows.items():
        print_named_summary(row_summary(name, row))

    single_rows = {v: accelerate(AP, v) for v in range(1, 14)}
    double_rows = {sites: accelerate(AP, *sites) for sites in combinations(range(1, 14), 2)}

    def passes_basic(row: tuple[int, ...]) -> bool:
        return (
            distinct_speeds(row)
            and covering_ok(row)
            and endpoint_units_ok(row)
            and full_unit_shell(row)
            and has_apex(row)
            and max_empty_gap(row) <= 3
        )

    single_pass = [v for v, row in single_rows.items() if passes_basic(row)]
    double_pass = [sites for sites, row in double_rows.items() if passes_basic(row)]
    admissible_single_pass = [v for v in single_pass if gate_admissible(v)]

    print("\nFinite-family audit")
    print(f"single_accelerations={len(single_rows)}")
    print(f"double_accelerations={len(double_rows)}")
    print(f"single_basic_pass_sites={single_pass}")
    print(f"single_basic_and_gate_pass_sites={admissible_single_pass}")
    print(f"double_basic_pass_count={len(double_pass)}")
    print(f"double_basic_pass_sites={double_pass}")

    pressure_classes = defaultdict(list)
    square_classes = defaultdict(list)
    square_profile_classes = defaultdict(list)
    for label, row in [("AP", AP), ("GW", GW)] + [(f"single_{v}", row) for v, row in single_rows.items()]:
        pressure_classes[row_summary(label, row)["pressure_shell_hist"]].append(label)
        square_classes[row_summary(label, row)["square_inner_histogram"]].append(label)
        square_profile_classes[row_summary(label, row)["square_profile_hash_histogram"]].append(label)

    print("\nPressure-shell class buckets among AP/GW/singles:")
    for fingerprint, labels in sorted(pressure_classes.items(), key=lambda item: (len(item[1]), item[1])):
        print(f"  labels={labels} fingerprint={fingerprint}")

    print("\nNode-squared inner-class buckets among AP/GW/singles:")
    for fingerprint, labels in sorted(square_classes.items(), key=lambda item: (len(item[1]), item[1])):
        print(f"  labels={labels} fingerprint={fingerprint}")

    print("\nNode-squared relative-profile buckets among AP/GW/singles:")
    for fingerprint, labels in sorted(square_profile_classes.items(), key=lambda item: (len(item[1]), item[1])):
        print(f"  labels={labels} fingerprint={fingerprint}")

    print("\nReadout")
    print("- Only site 12 is Jacobsthal-admissible under W(v)=[14-v,27-2v].")
    print("- AP and GW both pass the basic endpoint/divisibility/unit/apex/gap filters.")
    print("- The AP-residue false look-alike passes endpoint residues but fails covering_ok.")
    print("- Some formal accelerations duplicate speeds; distinct_speeds filters them before summit use.")
    print("- Coarse pressure shells and inner score sequences are too weak; the gate tournament and")
    print("  node-squared relative-profile hashes are the discriminators.")


if __name__ == "__main__":
    main()
