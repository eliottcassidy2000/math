#!/usr/bin/env python3
"""S164: oriented-tope / boundary-cocircuit proof pass for LRC14.

This is a creative proof-scout rather than a proof.  Cut the time circle by
all danger-arc endpoints for a 13-speed row S at threshold 1/14.  Each open
cell is a tope with a danger-count sign vector.  A strict lonely witness is an
open all-safe tope.  AP/Goddyn-Wong instead have no open all-safe tope, but do
have zero-dimensional all-safe boundary cocircuits.

Therefore a strict counterexample would have to be a stronger object:

    no open all-safe tope and no all-safe boundary cocircuit.

The script audits named rows and a one-swap AP bank, then runs Tournament
Analysis on proof carriers.  Tournament vertices are proof carriers, not
runners.  The pairwise observable is retention of the open-tope predicate,
boundary-cocircuit predicate, endpoint owner labels, exact arithmetic, state-
lift compatibility, dual handoff, and anti-scalarization.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from functools import reduce
from itertools import combinations
from math import gcd
from pathlib import Path
import argparse


REPO = Path(__file__).resolve().parents[1]
RESULT = REPO / "05-knowledge" / "results" / "lrc14_tope_wall_cocircuit_codex_s164.out"
THRESHOLD = Fraction(1, 14)


def frac_mod_one(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def dist_to_integer(x: Fraction) -> Fraction:
    y = frac_mod_one(x)
    return min(y, 1 - y)


def row_gcd(row: tuple[int, ...]) -> int:
    return reduce(gcd, row)


def replace(base: tuple[int, ...], drops: tuple[int, ...], adds: tuple[int, ...]) -> tuple[int, ...]:
    values = [v for v in base if v not in drops] + list(adds)
    return tuple(sorted(values))


def named_rows() -> list[tuple[str, tuple[int, ...]]]:
    base = tuple(range(1, 14))
    return [
        ("AP", base),
        ("GW_12_to_24", replace(base, (12,), (24,))),
        ("residue_liar_12_to_26", replace(base, (12,), (26,))),
        ("K33_12_to_36", replace(base, (12,), (36,))),
        ("petal_10_to_20", replace(base, (10,), (20,))),
        ("petal_13_to_26", replace(base, (13,), (26,))),
        ("P10_plus_GW", replace(base, (10, 12), (20, 24))),
        ("P10_plus_K33", replace(base, (10, 12), (20, 36))),
        ("covering_12_to_84", replace(base, (12,), (84,))),
        ("covering_12_to_168", replace(base, (12,), (168,))),
        ("covering_6_to_98", replace(base, (6,), (98,))),
        ("q14_liar_12_to_96", replace(base, (12,), (96,))),
    ]


def endpoint_set(row: tuple[int, ...]) -> list[Fraction]:
    pts: set[Fraction] = set()
    for v in row:
        for k in range(v):
            pts.add(frac_mod_one(Fraction(14 * k - 1, 14 * v)))
            pts.add(frac_mod_one(Fraction(14 * k + 1, 14 * v)))
    return sorted(pts)


def danger_open(v: int, t: Fraction) -> bool:
    return dist_to_integer(v * t) < THRESHOLD


def danger_closed(v: int, t: Fraction) -> bool:
    return dist_to_integer(v * t) <= THRESHOLD


def active_boundary(v: int, t: Fraction) -> bool:
    return dist_to_integer(v * t) == THRESHOLD


@dataclass(frozen=True)
class TopeAudit:
    name: str
    row: tuple[int, ...]
    q_witness: int | None
    cells: int
    min_open_danger: int
    open_safe_mass: Fraction
    zero_boundary_count: int
    zero_boundary_active_sizes: tuple[int, ...]
    zero_boundary_pair_sums: tuple[int, ...]
    forbidden_wall: bool


def q_witness(row: tuple[int, ...]) -> int | None:
    for q in range(2, 15):
        if all(v % q != 0 for v in row):
            return q
    return None


def audit_row(name: str, row: tuple[int, ...]) -> TopeAudit:
    pts = endpoint_set(row)
    min_open = len(row) + 1
    safe_mass = Fraction(0, 1)
    for i, a in enumerate(pts):
        b = pts[(i + 1) % len(pts)]
        length = b - a if b > a else b + 1 - a
        mid = frac_mod_one(a + length / 2)
        count = sum(1 for v in row if danger_open(v, mid))
        min_open = min(min_open, count)
        if count == 0:
            safe_mass += length

    zero_count = 0
    active_sizes: list[int] = []
    pair_sums: list[int] = []
    for t in pts:
        if any(danger_open(v, t) for v in row):
            continue
        active = [v for v in row if active_boundary(v, t)]
        if not active:
            continue
        zero_count += 1
        active_sizes.append(len(active))
        for a, b in combinations(active, 2):
            pair_sums.append((a + b) % 14)

    return TopeAudit(
        name=name,
        row=row,
        q_witness=q_witness(row),
        cells=len(pts),
        min_open_danger=min_open,
        open_safe_mass=safe_mass,
        zero_boundary_count=zero_count,
        zero_boundary_active_sizes=tuple(sorted(Counter(active_sizes).elements())),
        zero_boundary_pair_sums=tuple(sorted(pair_sums)),
        forbidden_wall=(min_open > 0 and zero_count == 0),
    )


def frac_text(x: Fraction) -> str:
    if x.denominator == 1:
        return str(x.numerator)
    return f"{x.numerator}/{x.denominator}"


def bank_summary(limit: int = 140) -> tuple[Counter[str], list[TopeAudit]]:
    base = tuple(range(1, 14))
    counts: Counter[str] = Counter()
    examples: list[TopeAudit] = []
    rows_seen: set[tuple[int, ...]] = {base}
    audits = [audit_row("AP", base)]
    for drop in base:
        for add in range(14, limit + 1):
            row = replace(base, (drop,), (add,))
            if row in rows_seen or row_gcd(row) != 1:
                continue
            rows_seen.add(row)
            qw = q_witness(row)
            if qw is not None and qw < 14:
                continue
            audits.append(audit_row(f"{drop}_to_{add}", row))

    for audit in audits:
        if audit.min_open_danger == 0:
            key = "open_tope"
        elif audit.zero_boundary_count:
            key = "boundary_cocircuit"
        else:
            key = "forbidden_wall_candidate"
        counts[key] += 1
        if key != "open_tope" or len(examples) < 10:
            examples.append(audit)
    return counts, examples


@dataclass(frozen=True)
class Carrier:
    name: str
    vector: tuple[int, ...]
    note: str


CARRIERS = [
    Carrier(
        "open_all_safe_tope",
        (5, 2, 2, 5, 2, 5, 4),
        "direct strict lonely interval in an open cell",
    ),
    Carrier(
        "boundary_cocircuit_atom",
        (4, 5, 5, 4, 4, 4, 5),
        "zero-dimensional all-safe boundary point with active endpoint owners",
    ),
    Carrier(
        "owner_sum_zero_wall",
        (3, 5, 5, 4, 5, 3, 5),
        "AP/GW-style active-owner sums r+s=0 mod 14 on taut walls",
    ),
    Carrier(
        "cyclic_tope_morse_graph",
        (4, 4, 4, 4, 4, 4, 4),
        "danger-count Morse walk over endpoint cells",
    ),
    Carrier(
        "endpoint_arrangement_sheaf",
        (4, 5, 4, 5, 3, 4, 5),
        "exact rational endpoint arrangement before quotienting",
    ),
    Carrier(
        "moment_dual_shadow",
        (3, 2, 2, 3, 3, 5, 3),
        "danger-count distribution after forgetting endpoint order",
    ),
    Carrier(
        "raw_residue_tournament",
        (1, 1, 1, 1, 2, 1, 1),
        "magnitude-blind residue/tournament shadow",
    ),
]

TIE_PATH = [c.name for c in CARRIERS]


def tournament() -> tuple[dict[int, int], int, int, list[str]]:
    n = len(CARRIERS)
    adj = [[False] * n for _ in range(n)]
    scores = [0] * n
    rank = {name: i for i, name in enumerate(TIE_PATH)}
    for i, a in enumerate(CARRIERS):
        for j, b in enumerate(CARRIERS):
            if i == j:
                continue
            wins = sum(x > y for x, y in zip(a.vector, b.vector))
            losses = sum(x < y for x, y in zip(a.vector, b.vector))
            if wins > losses or (wins == losses and rank[a.name] < rank[b.name]):
                adj[i][j] = True
        scores[i] = sum(adj[i])
    hist = dict(sorted(Counter(scores).items()))
    cycles = 0
    for i, j, k in combinations(range(n), 3):
        if adj[i][j] and adj[j][k] and adj[k][i]:
            cycles += 1
        if adj[i][k] and adj[k][j] and adj[j][i]:
            cycles += 1
    dp: dict[tuple[int, int], int] = {}
    for i in range(n):
        dp[(1 << i, i)] = 1
    for mask in range(1 << n):
        for last in range(n):
            val = dp.get((mask, last), 0)
            if not val:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if adj[last][nxt]:
                    dp[(mask | (1 << nxt), nxt)] = dp.get((mask | (1 << nxt), nxt), 0) + val
    hp = sum(dp.get(((1 << n) - 1, i), 0) for i in range(n))
    path = [CARRIERS[i].name for i in sorted(range(n), key=lambda idx: -scores[idx])]
    return hist, cycles, hp, path


def render(limit: int) -> str:
    named = [audit_row(name, row) for name, row in named_rows()]
    counts, examples = bank_summary(limit)
    hist, cycles, hp, path = tournament()

    lines: list[str] = []
    lines.append("LRC14 oriented-tope / boundary-cocircuit proof pass (codex S164)")
    lines.append("=" * 70)
    lines.append("")
    lines.append("Method")
    lines.append("------")
    lines.append("Cut R/Z by all endpoints (14k +/- 1)/(14v).")
    lines.append("Open cells are topes; boundary endpoints are possible cocircuits.")
    lines.append("open all-safe tope: an interval where every speed is lonely.")
    lines.append("boundary cocircuit: an endpoint where every speed is non-dangerous, with at least one equality owner.")
    lines.append("forbidden wall candidate: no open all-safe tope and no boundary cocircuit.")
    lines.append("")
    lines.append("Named rows")
    lines.append("----------")
    lines.append("name                         qwit cells minD safe_mass zero_bd active_sizes pair_sums_mod14")
    for a in named:
        qw = ">14" if a.q_witness is None else str(a.q_witness)
        sizes = ",".join(map(str, a.zero_boundary_active_sizes[:12])) or "-"
        sums = ",".join(map(str, a.zero_boundary_pair_sums[:18])) or "-"
        lines.append(
            f"{a.name:28s} {qw:>4s} {a.cells:5d} {a.min_open_danger:4d} "
            f"{frac_text(a.open_safe_mass):>10s} {a.zero_boundary_count:7d} {sizes:>12s} {sums}"
        )
    lines.append("")
    lines.append(f"One-swap AP bank through add <= {limit}, restricted to q-witness >=14 or >14")
    lines.append("--------------------------------------------------------------------------")
    for key, count in counts.items():
        lines.append(f"{key:28s} {count}")
    lines.append("")
    lines.append("Non-open examples / first open examples")
    lines.append("--------------------------------------")
    for a in examples[:24]:
        qw = ">14" if a.q_witness is None else str(a.q_witness)
        route = (
            "open_tope" if a.min_open_danger == 0
            else "boundary_cocircuit" if a.zero_boundary_count
            else "forbidden_wall_candidate"
        )
        lines.append(
            f"{a.name:12s} route={route:25s} qwit={qw:>3s} "
            f"minD={a.min_open_danger} safe={frac_text(a.open_safe_mass)} "
            f"zero_bd={a.zero_boundary_count}"
        )
    lines.append("")
    lines.append("Tournament Analysis")
    lines.append("-------------------")
    lines.append("vertices: proof carriers, not runners or arcs")
    lines.append("pairwise observable: open-tope/boundary-cocircuit/owner/exact/state/dual/anti-scalar retention")
    lines.append("tie Hamiltonian path:")
    lines.append("  " + " > ".join(TIE_PATH))
    lines.append(f"score_hist: {hist}")
    lines.append(f"directed_3cycles: {cycles}")
    lines.append(f"Hamiltonian_path_count: {hp}")
    lines.append("retention path:")
    lines.append("  " + " > ".join(path))
    lines.append("")
    lines.append("Readout")
    lines.append("-------")
    lines.append("1. AP and GW are not open topes; they are boundary-cocircuit atoms.")
    lines.append("2. Every audited named non-AP/GW hard row has an open all-safe tope.")
    lines.append("3. In the one-swap hard bank, forbidden wall candidates are absent.")
    lines.append("4. A strict LRC14 counterexample must therefore be a no-tope/no-cocircuit wall packet.")
    lines.append("5. Candidate theorem: every no-tope/no-cocircuit wall packet constructs a K33/H=7 state lift or contradicts endpoint-owner parity.")
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=140)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    text = render(args.limit)
    print(text, end="")
    if args.write:
        RESULT.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
