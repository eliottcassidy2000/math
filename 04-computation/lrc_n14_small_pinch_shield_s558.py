#!/usr/bin/env python3
"""Small-pinch shield analysis for HYP-2060 / HYP-2059.

This script records the S558 proof-search result:

* THM-396: for n=14, a universal blocker of every safe pinch residue of a
  reduced-sum <=14 pair must be a multiple of the pair sum.
* N3 from HYP-2059 is too strong: some HYP-2060-like sets have no small
  reduced-sum clearing pinch, while remaining lonely.
* The residual is a finite blocker-cover problem.  Tournament Analysis is run
  on small pair-cells for one such set.
"""

from __future__ import annotations

from collections import Counter, deque
from fractions import Fraction
from itertools import combinations, permutations
from math import gcd
import random


N = 14
THR = Fraction(1, N)


def dist_num(x: int, den: int) -> int:
    r = x % den
    return min(r, den - r)


def dist_frac(x: Fraction) -> Fraction:
    x %= 1
    return min(x, 1 - x)


def collar_at(speeds: tuple[int, ...], t: Fraction) -> Fraction:
    return min(dist_frac(v * t) for v in speeds)


def reduced_sum(a: int, b: int) -> int:
    return (a + b) // gcd(a, b)


def small_pairs(speeds: tuple[int, ...], max_s: int = N) -> list[tuple[int, int, int, int]]:
    out = []
    for a, b in combinations(speeds, 2):
        s = reduced_sum(a, b)
        if s <= max_s:
            out.append((a, b, s, a + b))
    return out


def pair_safe_residues(a: int, b: int) -> list[int]:
    den = a + b
    return [
        m
        for m in range(1, den)
        if dist_num(a * m, den) * N >= den
        and dist_num(b * m, den) * N >= den
    ]


def dangerous_residues(c: int, den: int, residues: list[int]) -> set[int]:
    return {m for m in residues if dist_num(c * m, den) * N < den}


def universal_blockers(speeds: tuple[int, ...], a: int, b: int) -> list[int]:
    den = a + b
    residues = pair_safe_residues(a, b)
    return [
        c
        for c in speeds
        if c not in (a, b) and dangerous_residues(c, den, residues) == set(residues)
    ]


def sum_multiple_shields(speeds: tuple[int, ...], a: int, b: int) -> list[int]:
    den = a + b
    return [c for c in speeds if c not in (a, b) and c % den == 0]


def best_pinch_for_pair(speeds: tuple[int, ...], a: int, b: int) -> tuple[Fraction, Fraction, list[int]]:
    den = a + b
    best = Fraction(-1)
    best_t = Fraction(0)
    blockers: list[int] = []
    for m in range(1, den):
        t = Fraction(m, den)
        value = collar_at(speeds, t)
        if value > best:
            best = value
            best_t = t
            blockers = [v for v in speeds if dist_frac(v * t) < THR]
    return best, best_t, blockers


def small_pinch_witness(speeds: tuple[int, ...]) -> tuple[int, int, int, Fraction] | None:
    for a, b, s, den in small_pairs(speeds):
        for m in range(1, den):
            t = Fraction(m, den)
            if collar_at(speeds, t) >= THR:
                return a, b, s, t
    return None


def exact_M(speeds: tuple[int, ...]) -> tuple[Fraction, list[Fraction]]:
    candidates: set[Fraction] = set()
    for v in speeds:
        for k in range(2 * v):
            candidates.add(Fraction(2 * k + 1, 2 * v) % 1)
    for a, b in combinations(speeds, 2):
        for den in (a + b, abs(a - b)):
            if den:
                for m in range(den + 1):
                    candidates.add(Fraction(m, den) % 1)
    best = Fraction(-1)
    arg: list[Fraction] = []
    for t in candidates:
        value = collar_at(speeds, t)
        if value > best:
            best = value
            arg = [t]
        elif value == best:
            arg.append(t)
    return best, sorted(arg)


def verify_thm396(limit_D: int = 500) -> tuple[int, int]:
    checked = 0
    failures = 0
    for den in range(3, limit_D + 1):
        for a in range(1, den):
            b = den - a
            if a >= b:
                continue
            s = reduced_sum(a, b)
            if s > N:
                continue
            residues = pair_safe_residues(a, b)
            for c in range(1, den):
                if dangerous_residues(c, den, residues) == set(residues):
                    failures += 1
            checked += 1
    return checked, failures


def collective_cover_counterexample() -> tuple[tuple[int, int], int, list[int], list[int]]:
    a, b = 3, 12
    den = a + b
    blockers = [14, 8, 10, 11, 13]
    residues = pair_safe_residues(a, b)
    covered: set[int] = set()
    for c in blockers:
        covered |= dangerous_residues(c, den, residues)
    assert covered == set(residues)
    assert not any(c % den == 0 for c in blockers)
    return (a, b), den, residues, blockers


def primitive(speeds: tuple[int, ...]) -> bool:
    g = 0
    for v in speeds:
        g = gcd(g, v)
    return g == 1


def sieve_covered(speeds: tuple[int, ...]) -> bool:
    return all(any(v % q == 0 for v in speeds) for q in range(2, N + 1))


def prime_power_covered(speeds: tuple[int, ...]) -> bool:
    return all(any(v % q == 0 for v in speeds) for q in (8, 9, 5, 7, 11, 13))


def not_ap_scaled(speeds: tuple[int, ...]) -> bool:
    g = 0
    for v in speeds:
        g = gcd(g, v)
    return tuple(sorted(v // g for v in speeds)) != tuple(range(1, N))


def checkable_hyp2060(speeds: tuple[int, ...]) -> bool:
    return (
        primitive(speeds)
        and len(set(speeds)) == len(speeds)
        and sieve_covered(speeds)
        and any(v % N == 0 for v in speeds)
        and any(v % 7 == 0 for v in speeds)
        and prime_power_covered(speeds)
        and bool(small_pairs(speeds, max_s=6))
        and not_ap_scaled(speeds)
    )


def sample_checkable_failures(trials: int = 200, hi: int = 220, seed: int = 77) -> tuple[int, int, list[tuple[int, ...]]]:
    rng = random.Random(seed)
    checked = 0
    tries = 0
    failures: list[tuple[int, ...]] = []
    while checked < trials and tries < 500_000:
        tries += 1
        speeds = tuple(sorted(rng.sample(range(1, hi + 1), 13)))
        if not checkable_hyp2060(speeds):
            continue
        checked += 1
        if small_pinch_witness(speeds) is None:
            failures.append(speeds)
    return checked, tries, failures


def tournament_fingerprint(rows: list[tuple[tuple[int, int], tuple[Fraction, int, int, int, int]]]) -> dict[str, object]:
    # Higher best collar is better; then fewer universal blockers, fewer shields,
    # smaller reduced sum, smaller pair sum.
    def key(row: tuple[tuple[int, int], tuple[Fraction, int, int, int, int]]) -> tuple[Fraction, int, int, int, int, int, int]:
        (a, b), (best, universals, shields, s, den) = row
        return (best, -universals, -shields, -s, -den, -a, -b)

    n = len(rows)
    adj = [[False] * n for _ in range(n)]
    flips_vs_s = 0
    for i, left in enumerate(rows):
        for j, right in enumerate(rows):
            if i == j:
                continue
            left_wins = key(left) > key(right)
            adj[i][j] = left_wins
            # reduced-sum order would prefer smaller s.
            if i < j:
                ls = left[1][3]
                rs = right[1][3]
                if ls != rs:
                    s_order_left_wins = ls < rs
                    if left_wins != s_order_left_wins:
                        flips_vs_s += 1

    scores = [sum(adj[i]) for i in range(n)]
    c3 = 0
    for i, j, k in combinations(range(n), 3):
        cyc = (
            adj[i][j] and adj[j][k] and adj[k][i]
        ) or (
            adj[i][k] and adj[k][j] and adj[j][i]
        )
        c3 += int(cyc)

    def reaches(start: int) -> set[int]:
        seen = {start}
        todo = deque([start])
        while todo:
            u = todo.popleft()
            for v in range(n):
                if adj[u][v] and v not in seen:
                    seen.add(v)
                    todo.append(v)
        return seen

    remaining = set(range(n))
    sccs = []
    while remaining:
        u = next(iter(remaining))
        ru = reaches(u)
        comp = {v for v in remaining if u in reaches(v) and v in ru}
        sccs.append(len(comp))
        remaining -= comp

    hp: int | str = 0
    if n <= 8:
        for perm in permutations(range(n)):
            if all(adj[perm[i]][perm[i + 1]] for i in range(n - 1)):
                hp += 1
    else:
        hp = "skipped(n>8)"

    return {
        "score_hist": dict(sorted(Counter(scores).items())),
        "c3": c3,
        "sccs": sorted(sccs, reverse=True),
        "hamiltonian_paths": hp,
        "edge_flips_vs_reduced_sum_order": flips_vs_s,
    }


def analyze_pair_cells(speeds: tuple[int, ...]) -> None:
    rows = []
    print(f"\nPair-cell analysis for {speeds}:")
    for a, b, s, den in small_pairs(speeds):
        best, best_t, blockers = best_pinch_for_pair(speeds, a, b)
        universals = universal_blockers(speeds, a, b)
        shields = sum_multiple_shields(speeds, a, b)
        rows.append(((a, b), (best, len(universals), len(shields), s, den)))
        print(
            f"  pair=({a},{b}) s={s} D={den}: best={best} at {best_t}; "
            f"best-blockers={blockers}; universal={universals}; shields={shields}"
        )
    if rows:
        print("  Tournament Analysis:")
        print("    vertices: small reduced-sum pair-cells")
        print("    pairwise observable: (best collar, universal blockers, shields, reduced sum, pair sum)")
        print("    switch/gauge: higher best collar; ties lower blocker debt; tie path lexicographic")
        print(f"    fingerprints: {tournament_fingerprint(rows)}")


def main() -> None:
    print("S558 small-pinch shield analysis for HYP-2060")
    print("=" * 72)
    checked, failures = verify_thm396()
    print(f"THM-396 finite residue verification: checked {checked} pair types with D<=500; non-shield universal blockers={failures}")

    pair, den, residues, blockers = collective_cover_counterexample()
    print("\nCollective-cover obstruction to a stronger theorem:")
    print(f"  pair={pair}, D={den}, pair-safe residues={residues}")
    print(f"  non-shield blockers={blockers} cover every pair-safe residue")

    examples = {
        "AP tight wall": tuple(range(1, 14)),
        "V* tight sporadic": (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 24),
        "sieve-minimal lonely": (2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
        "checkable no-small-pinch": (1, 2, 9, 26, 110, 153, 166, 170, 178, 190, 192, 196, 201),
    }
    for name, speeds in examples.items():
        M, arg = exact_M(speeds)
        witness = small_pinch_witness(speeds)
        print(f"\n{name}:")
        print(f"  M={M}; argmax sample={[str(t) for t in arg[:4]]}; small-pinch witness={witness}")
        print(f"  checkable-HYP2060-proxy={checkable_hyp2060(speeds)}")
        analyze_pair_cells(speeds)

    checked, tries, failures = sample_checkable_failures()
    print("\nRandom checkable-HYP2060 proxy probe:")
    print(f"  checked={checked}, tries={tries}, no-small-pinch failures={len(failures)}")
    for speeds in failures[:5]:
        print(f"  failure={speeds}")
        analyze_pair_cells(speeds)


if __name__ == "__main__":
    main()
