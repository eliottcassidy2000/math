#!/usr/bin/env python3
"""
lrc_n14_n17_gate_skip_transfer_s560.py

codex-2026-06-02 S560

Carry the n=17 prime-gate / skip-8 idea back to n=14.

The n=17 attempt isolated a family

    {r} union {17*q : 1 <= q <= 16, q != 8}

as the closest structured prime-gate repair.  This script asks whether the
same "skip the half-gate" packet already exists in the n=14 work.

It does, with a composite twist: for n=14 the best row-parent, gate, and
double-gate ladders all skip q=6, the predecessor of the apex q=7.  The apex is
kept as a shield/bridge, so the missing half-gate shifts one step left.
"""

from __future__ import annotations

from collections import Counter, deque
from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from itertools import combinations, permutations
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
S356 = SourceFileLoader(
    "lonely_runner_residue_probe_s356",
    str(ROOT / "04-computation" / "lonely_runner_residue_probe_s356.py"),
).load_module()


@dataclass(frozen=True)
class LadderRow:
    n: int
    scale: int
    breaker: int
    skip: int
    speeds: tuple[int, ...]
    primitive: bool
    missing_moduli: tuple[int, ...]
    forbidden_length: Fraction
    max_gap: Fraction
    gap_ratio: Fraction
    witness: Fraction | None
    boundary_witnesses: int
    classification: str


def fmt(value: Fraction | None) -> str:
    return S356.fmt_frac(value)


def ffloat(value: Fraction) -> str:
    return f"{float(value):.6f}"


def gcd_all(values: tuple[int, ...]) -> int:
    g = 0
    for value in values:
        g = gcd(g, value)
    return g


def missing_moduli(n: int, speeds: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(m for m in range(2, n + 1) if all(speed % m for speed in speeds))


def classify(report: object) -> str:
    if report.max_gap > 0:
        return "positive_gap"
    if report.boundary_witness_count:
        return "boundary_only"
    return "open_cover_candidate"


def ladder_speeds(n: int, scale: int, breaker: int, skip: int) -> tuple[int, ...]:
    return tuple(sorted((breaker,) + tuple(scale * q for q in range(1, n) if q != skip)))


def ladder_row(n: int, scale: int, breaker: int, skip: int) -> LadderRow:
    raw = ladder_speeds(n, scale, breaker, skip)
    primitive = gcd_all(raw) == 1
    speeds = S356.normalize_speed_set(list(raw))
    report = S356.report(f"n{n}_scale{scale}_skip{skip}", list(speeds))
    return LadderRow(
        n=n,
        scale=scale,
        breaker=breaker,
        skip=skip,
        speeds=speeds,
        primitive=primitive,
        missing_moduli=missing_moduli(n, speeds),
        forbidden_length=report.forbidden_length,
        max_gap=report.max_gap,
        gap_ratio=report.max_gap / report.threshold,
        witness=report.witness or report.boundary_witness,
        boundary_witnesses=report.boundary_witness_count,
        classification=classify(report),
    )


def unit_wall_count(n: int, speeds: tuple[int, ...]) -> int:
    if any(speed % n == 0 for speed in speeds):
        return 0
    return sum(
        1
        for a in range(1, n)
        if gcd(a, n) == 1 and S356.is_lonely_witness(speeds, Fraction(a, n))
    )


def scan_ladders(n: int, scale: int) -> list[LadderRow]:
    rows: list[LadderRow] = []
    for breaker in range(1, n):
        for skip in range(1, n):
            row = ladder_row(n, scale, breaker, skip)
            if row.primitive:
                rows.append(row)
    rows.sort(
        key=lambda row: (
            row.gap_ratio,
            len(row.missing_moduli),
            -row.forbidden_length,
            row.skip,
            row.breaker,
        )
    )
    return rows


def print_best_table(title: str, rows: list[LadderRow], limit: int = 10) -> None:
    print(title)
    print("  label                 gap/th       max_gap       length boundary missing witness")
    for row in rows[:limit]:
        label = f"n{row.n} scale{row.scale} r{row.breaker} skip{row.skip}"
        miss = ",".join(str(m) for m in row.missing_moduli) or "-"
        print(
            f"  {label:<21} {fmt(row.gap_ratio):>10} {fmt(row.max_gap):>12} "
            f"{fmt(row.forbidden_length):>12} {row.boundary_witnesses:>8} "
            f"{miss:<8} {fmt(row.witness):>10}"
        )
    print()


def skip_hist(rows: list[LadderRow], top: int = 24) -> dict[int, int]:
    return dict(sorted(Counter(row.skip for row in rows[:top]).items()))


def tournament_fingerprint(rows: list[LadderRow]) -> dict[str, object]:
    # Vertices are whole gate-packet rows, not runners.
    # More counterexample-like rows win: smaller gap, fewer missing moduli,
    # longer forbidden length, then smaller scale.
    def key(row: LadderRow) -> tuple[Fraction, int, Fraction, int, int]:
        return (-row.gap_ratio, -len(row.missing_moduli), row.forbidden_length, -row.scale, -row.skip)

    n = len(rows)
    adj = [[False] * n for _ in range(n)]
    for i, left in enumerate(rows):
        for j, right in enumerate(rows):
            if i != j:
                adj[i][j] = key(left) > key(right)

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
    sccs: list[int] = []
    while remaining:
        u = next(iter(remaining))
        ru = reaches(u)
        comp = {v for v in remaining if v in ru and u in reaches(v)}
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
        "vertices": [f"n{row.n}:scale{row.scale}:skip{row.skip}" for row in rows],
        "score_hist": dict(sorted(Counter(scores).items())),
        "c3": c3,
        "sccs": sorted(sccs, reverse=True),
        "hamiltonian_paths": hp,
    }


def main() -> None:
    print("S560 n=14/n=17 gate-skip transfer")
    print("=" * 72)

    print("1. Top-gate unit wall transfer")
    n14_initial = tuple(range(1, 14))
    n17_initial = tuple(range(1, 17))
    print(
        f"  n=14 no-14-gate unit witnesses={unit_wall_count(14, n14_initial)} "
        f"(phi(14)=6)"
    )
    print(
        f"  n=17 no-17-gate unit witnesses={unit_wall_count(17, n17_initial)} "
        f"(phi(17)=16)"
    )
    print(
        "  Carryover: the n=17 prime-gate lemma is the top-q case of the "
        "denominator sieve, not a prime-only phenomenon."
    )
    print()

    n14_scale7 = scan_ladders(14, 7)
    n14_scale14 = scan_ladders(14, 14)
    n14_scale28 = scan_ladders(14, 28)
    n17_scale17 = scan_ladders(17, 17)

    print_best_table("2. n=14 row-parent ladder (scale 7)", n14_scale7)
    print_best_table("3. n=14 gate ladder (scale 14)", n14_scale14)
    print_best_table("4. n=14 double-gate ladder (scale 28)", n14_scale28)
    print_best_table("5. n=17 prime-gate ladder (scale 17)", n17_scale17)

    print("6. Skip histograms among the top 24 primitive ladder rows")
    print(f"  n=14 scale 7:  {skip_hist(n14_scale7)}")
    print(f"  n=14 scale 14: {skip_hist(n14_scale14)}")
    print(f"  n=14 scale 28: {skip_hist(n14_scale28)}")
    print(f"  n=17 scale 17: {skip_hist(n17_scale17)}")
    print(
        "  Surprise: n=17 skips the literal half-gate q=8, while n=14 skips "
        "q=6 and keeps the apex/bridge q=7."
    )
    print()

    selected = [
        n14_scale7[0],
        n14_scale14[0],
        n14_scale28[0],
        n17_scale17[0],
    ]
    print("7. Cross-row ratios")
    for row in selected:
        apex = row.n // 2 if row.n % 2 == 0 else (row.n - 1) // 2
        print(
            f"  n={row.n} scale={row.scale} skip={row.skip} apex_or_half={apex} "
            f"gap/th={fmt(row.gap_ratio)} witness={fmt(row.witness)}"
        )
    print(
        "  n=14 scale 7 -> 14 -> 28 preserves skip=6 and halves the exact "
        "gap/th at each gate-depth step."
    )
    print()

    print("8. Tournament Analysis")
    print("  vertices: best gate-packet rows (n14 scale 7/14/28, n17 scale17)")
    print("  pairwise observable: (gap/th, missing moduli, forbidden length, scale, skip)")
    print("  switch/gauge: smaller gap wins; ties prefer sieve-complete and longer coverage")
    print(f"  fingerprints: {tournament_fingerprint(selected)}")
    print()

    print("Synthesis")
    print(
        "  The n=17 skip-8 packet carries back to n=14 as a predecessor-of-apex "
        "packet: skip q=6, keep q=7.  In the 2*7 row the apex is a necessary "
        "shield/bridge, so the missing middle gate shifts left."
    )


if __name__ == "__main__":
    main()
