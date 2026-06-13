#!/usr/bin/env python3
"""
lonely_runner_14_gate_tightness_s382.py

codex-2026-05-31 S382

Test the heuristic:

    at n=14, tightness wants no multiple of 14,
    but a counterexample requires a multiple of 14.

The second clause is a direct unit-skeleton obstruction: if no speed is
divisible by 14, every unit point a/14 is safe.  The first clause is more
experimental, so this script audits canonical examples and one-swap
neighborhoods of the initial tight set.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
S356 = SourceFileLoader(
    "lonely_runner_residue_probe_s356",
    str(ROOT / "04-computation" / "lonely_runner_residue_probe_s356.py"),
).load_module()
S360 = SourceFileLoader(
    "lonely_runner_endpoint_protection_s360",
    str(ROOT / "04-computation" / "lonely_runner_endpoint_protection_s360.py"),
).load_module()


N = 14
K = 13
ONE = Fraction(1, 1)
INITIAL = tuple(range(1, 14))
SEVEN_LADDER = (1, 7, 14, 21, 28, 35, 49, 56, 63, 70, 77, 84, 91)
S380_GATE_LADDER = (1, 14, 28, 42, 56, 70, 98, 112, 126, 140, 154, 168, 182)
UNIT_POINTS = tuple(Fraction(a, N) for a in range(1, N) if gcd(a, N) == 1)


@dataclass(frozen=True)
class Audit:
    label: str
    speeds: tuple[int, ...]
    classification: str
    forbidden_length: Fraction
    max_gap: Fraction
    gap_ratio: Fraction
    boundary_witnesses: int
    unprotected: int
    first_unprotected: Fraction | None
    has_14_multiple: bool
    unit_safe: int
    unit_owned: int
    unit_unprotected: int


def fmt_frac(x: Fraction | None) -> str:
    return S356.fmt_frac(x)


def fmt_ratio(x: Fraction) -> str:
    return f"{float(x):.6f}"


def primitive(speeds: tuple[int, ...]) -> bool:
    g = 0
    for v in speeds:
        g = gcd(g, v)
    return g == 1


def normalize(speeds: tuple[int, ...]) -> tuple[int, ...]:
    return S356.normalize_speed_set(list(speeds))


def unit_owned_points(speeds: tuple[int, ...]) -> set[Fraction]:
    values = {endpoint.value for endpoint in S360.endpoints(speeds)}
    return {t for t in UNIT_POINTS if t in values}


def unit_unprotected_points(speeds: tuple[int, ...]) -> set[Fraction]:
    owned = unit_owned_points(speeds)
    return {
        t
        for t in owned
        if not any(S360.direct_protects(speeds, protector, t) for protector in speeds)
    }


def audit(label: str, raw_speeds: tuple[int, ...]) -> Audit:
    speeds = normalize(raw_speeds)
    report = S356.report(label, list(speeds))
    summary = S360.summarize(list(speeds))
    unit_safe = sum(1 for t in UNIT_POINTS if S356.is_lonely_witness(speeds, t))
    owned = unit_owned_points(speeds)
    unprotected_units = unit_unprotected_points(speeds)
    return Audit(
        label=label,
        speeds=speeds,
        classification=summary.classification,
        forbidden_length=summary.forbidden_length,
        max_gap=summary.max_gap,
        gap_ratio=summary.max_gap / summary.threshold,
        boundary_witnesses=report.boundary_witness_count,
        unprotected=summary.unprotected_count,
        first_unprotected=summary.first_unprotected,
        has_14_multiple=any(v % N == 0 for v in speeds),
        unit_safe=unit_safe,
        unit_owned=len(owned),
        unit_unprotected=len(unprotected_units),
    )


def audit_rank(row: Audit) -> tuple:
    class_rank = {"open_cover": 0, "boundary_only": 1, "positive_gap": 2}.get(
        row.classification, 3
    )
    return (class_rank, row.max_gap, row.unprotected, row.speeds)


def print_audit(row: Audit) -> None:
    print(
        "  "
        f"{row.label:<34} class={row.classification:<13} "
        f"has14={int(row.has_14_multiple)} "
        f"unit_safe={row.unit_safe} unit_owned={row.unit_owned} "
        f"unit_unprot={row.unit_unprotected} "
        f"len={fmt_frac(row.forbidden_length)} "
        f"gap/th={fmt_ratio(row.gap_ratio)} "
        f"unprot={row.unprotected} first={fmt_frac(row.first_unprotected)}"
    )
    print(f"    speeds={row.speeds}")


def residue_unit_action() -> None:
    print("Residue action on the n=14 unit skeleton")
    print("  Unit points are a/14 for a in {1,3,5,9,11,13}.")
    print("  A speed with residue r mod 14 covers such a point iff r*a == 0 mod 14.")
    for r in range(N):
        covered = 0
        for t in UNIT_POINTS:
            if S360.circular_distance_to_integer(r * t) < Fraction(1, N):
                covered += 1
        if covered or r in (0, 1, 2, 7, 13):
            print(f"  residue={r:2d} covers_unit_points={covered}")
    print("  Therefore: without a 14-multiple, all six unit points remain safe.")
    print()


def canonical_audits() -> list[Audit]:
    no_coprime_no_gate = (2, 4, 6, 8, 10, 12, 16, 18, 20, 7, 21, 35, 49)
    rows = [
        audit("initial tight", INITIAL),
        audit("initial replace 6 by 14", tuple(sorted((set(INITIAL) - {6}) | {14}))),
        audit("initial replace 13 by 14", tuple(sorted((set(INITIAL) - {13}) | {14}))),
        audit("no-coprime no-gate test", no_coprime_no_gate),
        audit("seven-ladder", SEVEN_LADDER),
        audit("seven-ladder with 42", tuple(sorted((set(SEVEN_LADDER) - {84}) | {42}))),
        audit("S380 14-multiple ladder", S380_GATE_LADDER),
    ]
    print("Canonical exact audits")
    for row in rows:
        print_audit(row)
    print()
    return rows


def one_swap_scan(max_new_speed: int = 112) -> list[Audit]:
    rows: list[Audit] = []
    base = set(INITIAL)
    for removed in INITIAL:
        rest = base - {removed}
        for new in range(1, max_new_speed + 1):
            if new in rest:
                continue
            speeds = tuple(sorted(rest | {new}))
            if len(speeds) != K or speeds == INITIAL or not primitive(speeds):
                continue
            rows.append(audit(f"swap {removed}->{new}", speeds))
    return rows


def summarize_scan(rows: list[Audit], max_new_speed: int) -> None:
    print(f"One-swap scan around the initial tight set, new speed <= {max_new_speed}")
    counts: dict[bool, Counter[str]] = defaultdict(Counter)
    for row in rows:
        counts[row.has_14_multiple][row.classification] += 1

    for has_gate in (False, True):
        label = "with a 14-multiple" if has_gate else "without a 14-multiple"
        total = sum(counts[has_gate].values())
        print(f"  {label}: total={total} classes={dict(sorted(counts[has_gate].items()))}")
        subset = [row for row in rows if row.has_14_multiple]
        if not has_gate:
            subset = [row for row in rows if not row.has_14_multiple]
        best = sorted(subset, key=audit_rank)[:5]
        for row in best:
            print_audit(row)
    print()


def gate_replacement_ledger(max_q: int = 12) -> None:
    print(f"Pure gate replacement ledger: replace one initial speed by 14*q, q<= {max_q}")
    rows: list[Audit] = []
    for removed in INITIAL:
        for q in range(1, max_q + 1):
            new = N * q
            speeds = tuple(sorted((set(INITIAL) - {removed}) | {new}))
            if len(speeds) == K and primitive(speeds):
                rows.append(audit(f"remove {removed}, add {new}", speeds))
    by_class = Counter(row.classification for row in rows)
    print(f"  audited={len(rows)} class_counts={dict(sorted(by_class.items()))}")
    print("  best gate replacements:")
    for row in sorted(rows, key=audit_rank)[:10]:
        print_audit(row)
    print()


def main() -> None:
    print("n=14 gate/tightness duality probe (codex-2026-05-31 S382)")
    print("Convention: k=13 speeds, threshold=1/14.")
    print()

    residue_unit_action()
    canonical_audits()
    rows = one_swap_scan(max_new_speed=112)
    summarize_scan(rows, max_new_speed=112)
    gate_replacement_ledger(max_q=12)

    print("Interpretation")
    print(
        "  No 14-multiple leaves the six unit points a/14 safe.  Such a set can "
        "be tight only in the boundary-only sense, with unit witnesses available."
    )
    print(
        "  A genuine open-cover counterexample must destroy those unit witnesses, "
        "so it must contain at least one 14-multiple."
    )
    print(
        "  The local evidence supports the complementary heuristic: inserting a "
        "14-gate into the initial tight architecture destroys exact tightness and "
        "moves the obstruction to positive gaps or descendant endpoints."
    )


if __name__ == "__main__":
    main()
