#!/usr/bin/env python3
"""
lrc_zeckendorf_bridge_s451.py

codex-2026-05-31 S451

Merge the repo's Zeckendorf thread into the LRC anti-Bohr analogy thread.

Zeckendorf is not an exact LRC formulation.  It is the canonical normal form
for independent sets in a path graph, and therefore a useful model for any LRC
proof route that turns endpoint debt into a no-adjacent-carry recursion.
"""

from __future__ import annotations

from dataclasses import dataclass


def fibs_upto(n: int) -> list[tuple[int, int]]:
    """Zeckendorf convention: F1=1, F2=2, then Fk=F{k-1}+F{k-2}."""
    fibs = [(1, 1), (2, 2)]
    while fibs[-1][1] < n:
        idx = fibs[-1][0] + 1
        fibs.append((idx, fibs[-1][1] + fibs[-2][1]))
    return fibs


def zeckendorf(n: int) -> list[tuple[int, int]]:
    rem = n
    out: list[tuple[int, int]] = []
    for idx, value in reversed(fibs_upto(n)):
        if value <= rem:
            out.append((idx, value))
            rem -= value
    return list(reversed(out))


def fmt_z(n: int) -> str:
    parts = zeckendorf(n)
    return " + ".join(f"F{idx}={value}" for idx, value in parts)


def gaps(parts: list[tuple[int, int]]) -> str:
    if len(parts) < 2:
        return "-"
    return ",".join(str(b[0] - a[0]) for a, b in zip(parts, parts[1:]))


def path_independence(m: int, x: int) -> int:
    """Independence polynomial I(P_m,x)."""
    if m == 0:
        return 1
    prev2 = 1
    prev1 = 1 + x
    if m == 1:
        return prev1
    for _ in range(2, m + 1):
        prev2, prev1 = prev1, prev1 + x * prev2
    return prev1


@dataclass(frozen=True)
class Shadow:
    label: str
    value: int
    reading: str


SHADOWS = [
    Shadow("initial n=8 debt", 4, "unit endpoints left by {1,...,7}"),
    Shadow("initial n=14 debt", 6, "unit endpoints left by {1,...,13}"),
    Shadow("initial n=16 debt", 8, "unit endpoints left by {1,...,15}"),
    Shadow("n=14 runners", 14, "first post-current-frontier denominator"),
    Shadow("n=14 lower speeds", 13, "columns before a forced 14-gate"),
    Shadow("seven-ladder debt", 84, "unprotected endpoints in S450 atlas"),
    Shadow("S380 gate-ladder debt", 168, "unprotected endpoints after gate export"),
    Shadow("seven-ladder gap denom", 924, "denominator of max_gap/th=5/924"),
    Shadow("S380 gap denom", 1848, "denominator of max_gap/th=5/1848"),
]


def print_path_fugacity_table() -> None:
    print("Path independence: Zeckendorf to tournament fugacity")
    print("=" * 88)
    print(f"{'m':>2} {'I(P_m,1)':>10} {'I(P_m,2)':>10} {'Zeckendorf of I(P_m,2)':<42} note")
    for m in range(1, 9):
        z1 = path_independence(m, 1)
        z2 = path_independence(m, 2)
        note = "H=21 path obstruction" if z2 == 21 else ""
        print(f"{m:>2} {z1:>10} {z2:>10} {fmt_z(z2):<42} {note}")


def print_lrc_shadow_table() -> None:
    print()
    print("LRC boundary quantities in Zeckendorf normal form")
    print("=" * 88)
    print(f"{'quantity':<26} {'value':>8} {'Zeckendorf':<44} gaps")
    for item in SHADOWS:
        parts = zeckendorf(item.value)
        print(f"{item.label:<26} {item.value:>8} {fmt_z(item.value):<44} {gaps(parts)}")
    print()
    print("Readings")
    print("=" * 88)
    for item in SHADOWS:
        print(f"- {item.label}: {item.reading}")


def print_synthesis() -> None:
    print()
    print("S451 synthesis")
    print("=" * 88)
    print(
        "Zeckendorf is the x=1 path-graph independence normal form; "
        "tournament OCF is the x=2 independence-polynomial regime."
    )
    print(
        "For LRC, the relevant import is not Fibonacci numerology.  It is "
        "the no-adjacent-carry rule: a proof can try to make endpoint debt "
        "move through a recurrence graph where adjacent repairs cannot both "
        "be used without exporting new debt."
    )
    print(
        "In irrational rotation language this is Ostrowski numeration; for "
        "the golden slope all continued-fraction digits are 1, so the "
        "Ostrowski expansion is exactly Zeckendorf.  That makes the golden "
        "case a clean model of anti-Bohr boundary recursion."
    )


def main() -> None:
    print_path_fugacity_table()
    print_lrc_shadow_table()
    print_synthesis()


if __name__ == "__main__":
    main()
