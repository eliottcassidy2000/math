#!/usr/bin/env python3
"""
lrc_zeckendorf_bridge_s460.py

codex-2026-05-31 S460

Connect Zeckendorf decomposition to the Lonely Runner endpoint program and
the tournament independence-polynomial thread.

Old repo threads give the dictionary:

* Zeckendorf representations are independent sets in the Fibonacci path,
  evaluated at fugacity x=1.
* Tournament OCF is H(T)=I(Omega(T),2), an independence polynomial at fugacity
  x=2 on the odd-cycle conflict graph.
* LRC endpoint protection is a finite row/column incidence problem.

This script adds the missing bridge: use Zeckendorf supports as coordinates on
the LRC speed columns, then re-audit the existing n=14/15/16 gate-cover rows.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from itertools import combinations
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
S360 = SourceFileLoader(
    "lonely_runner_endpoint_protection_s360",
    str(ROOT / "04-computation" / "lonely_runner_endpoint_protection_s360.py"),
).load_module()
S420 = SourceFileLoader(
    "lrc_integer_programming_modes_s420",
    str(ROOT / "04-computation" / "lrc_integer_programming_modes_s420.py"),
).load_module()


FRONTIER_NS = (14, 15, 16, 18, 21, 24)


@dataclass(frozen=True)
class CoverFamily:
    n: int
    exact_size: int | None
    forced: tuple[int, ...]
    covers: tuple[tuple[int, ...], ...]


def fibs_up_to(value: int) -> list[int]:
    fibs = [1, 2]
    while fibs[-1] <= value:
        fibs.append(fibs[-1] + fibs[-2])
    return fibs[:-1]


def zeck(value: int) -> tuple[tuple[int, int], ...]:
    """Return Zeckendorf terms as descending (index, Fibonacci value)."""
    if value <= 0:
        return tuple()
    remaining = value
    out: list[tuple[int, int]] = []
    fibs = fibs_up_to(value)
    for offset in range(len(fibs) - 1, -1, -1):
        fib = fibs[offset]
        if fib <= remaining:
            out.append((offset + 1, fib))
            remaining -= fib
        if remaining == 0:
            break
    return tuple(out)


def zeck_indices(value: int) -> tuple[int, ...]:
    return tuple(index for index, _fib in zeck(value))


def zeck_text(value: int) -> str:
    terms = zeck(value)
    if not terms:
        return "0"
    symbolic = "+".join(f"F{index}" for index, _fib in terms)
    numeric = "+".join(str(fib) for _index, fib in terms)
    return f"{symbolic}={numeric}"


def v2(value: int) -> int:
    out = 0
    while value % 2 == 0:
        out += 1
        value //= 2
    return out


def phi(value: int) -> int:
    return sum(1 for item in range(1, value + 1) if gcd(item, value) == 1)


def mode_text(value: int) -> str:
    h = v2(value)
    return f"2^{h}*{value >> h}"


def path_independence(m: int, x: int) -> int:
    if m == 0:
        return 1
    if m == 1:
        return 1 + x
    a, b = 1, 1 + x
    for _ in range(2, m + 1):
        a, b = b, b + x * a
    return b


def cycle_independence_at_two(m: int) -> int:
    return 2**m + (-1) ** m


def digit_load(values: tuple[int, ...]) -> Counter[int]:
    load: Counter[int] = Counter()
    for value in values:
        for index in zeck_indices(value):
            load[index] += 1
    return load


def load_text(load: Counter[int]) -> str:
    if not load:
        return "-"
    return " ".join(f"F{index}:{load[index]}" for index in sorted(load))


def exact_gate_cover_family(n: int) -> CoverFamily:
    targets = S420.endpoint_values_for_owner(n, n)
    candidates = tuple(range(1, n))
    columns = S420.build_cover_columns(n, targets, candidates)
    result = S420.set_cover_result(f"owner {n} lower cover", n, targets, candidates)
    speed_to_mask = {column.speed: column.mask for column in columns}
    full_mask = (1 << len(targets)) - 1
    covers: list[tuple[int, ...]] = []
    if result.exact_size is not None:
        for combo in combinations(candidates, result.exact_size):
            mask = 0
            for speed in combo:
                mask |= speed_to_mask.get(speed, 0)
            if mask == full_mask:
                covers.append(combo)
    return CoverFamily(
        n=n,
        exact_size=result.exact_size,
        forced=result.forced_columns,
        covers=tuple(covers),
    )


def cover_column_stats(n: int) -> dict[int, tuple[int, int]]:
    targets = S420.endpoint_values_for_owner(n, n)
    candidates = tuple(range(1, n))
    columns = S420.build_cover_columns(n, targets, candidates)
    union_without: dict[int, int] = {}
    for column in columns:
        mask = 0
        for other in columns:
            if other.speed != column.speed:
                mask |= other.mask
        union_without[column.speed] = mask
    return {
        column.speed: (column.size, (column.mask & ~union_without[column.speed]).bit_count())
        for column in columns
    }


def seven_ladder() -> tuple[int, ...]:
    return (1, 7, 14, 21, 28, 35, 49, 56, 63, 70, 77, 84, 91)


def s380_gate_ladder() -> tuple[int, ...]:
    return (1, 14, 28, 42, 56, 70, 98, 112, 126, 140, 154, 168, 182)


def unprotected_values(speeds: tuple[int, ...]) -> set[Fraction]:
    endpoints = {endpoint.value for endpoint in S360.endpoints(speeds)}
    return {
        value
        for value in endpoints
        if not any(S360.direct_protects(speeds, speed, value) for speed in speeds)
    }


def owner_debt_rows(speeds: tuple[int, ...], limit: int = 8) -> list[tuple[int, int, int]]:
    bad = unprotected_values(speeds)
    by_owner_labels: Counter[int] = Counter()
    by_owner_unique: dict[int, set[Fraction]] = defaultdict(set)
    for endpoint in S360.endpoints(speeds):
        if endpoint.value in bad:
            by_owner_labels[endpoint.speed] += 1
            by_owner_unique[endpoint.speed].add(endpoint.value)
    return [
        (owner, label_count, len(by_owner_unique[owner]))
        for owner, label_count in by_owner_labels.most_common(limit)
    ]


def print_header(title: str) -> None:
    print("=" * 96)
    print(title)
    print("=" * 96)


def print_old_thread_dictionary() -> None:
    print_header("OLD THREAD DICTIONARY")
    print(
        "Zeckendorf: integer columns are independent supports in the Fibonacci path P_infty."
    )
    print(
        "Tournaments: H(T)=I(Omega(T),2), with Omega the odd-cycle conflict graph."
    )
    print(
        "LRC: endpoint rows are covered by speed columns; this script labels each column by its Zeckendorf support."
    )
    print()
    print("path/cycle fugacity table")
    print("m  I(P_m,1)  I(P_m,2)  Zeck(I(P_m,2))         I(C_m,2) odd-cycle note")
    print("-" * 96)
    for m in range(1, 9):
        path_one = path_independence(m, 1)
        path_two = path_independence(m, 2)
        cycle_text = "-"
        if m >= 3:
            cycle_text = str(cycle_independence_at_two(m))
        note = "odd cycle branch" if m >= 3 and m % 2 == 1 else "even/not a tournament Omega branch"
        print(
            f"{m:<2} {path_one:>9} {path_two:>11} "
            f"{zeck_text(path_two):<24} {cycle_text:>9} {note}"
        )
    print()
    print(
        "Key anchors: 7 = I(C_3,2) = F4+F2, while 21 = I(P_4,2) = F7."
    )
    print(
        "So the old forbidden-H story and the LRC endpoint-cover story both pass through path independence."
    )
    print()


def print_frontier_shell() -> None:
    print_header("ZECKENDORF SHELL AROUND THE LRC FRONTIER")
    print("n   mode      phi(n)  Zeckendorf         top anchor  low payload  note")
    print("-" * 96)
    for n in range(10, 25):
        terms = zeck(n)
        top_index, top_value = terms[0]
        payload = n - top_value
        gaps = [terms[i][0] - terms[i + 1][0] for i in range(len(terms) - 1)]
        note = ""
        if n in (14, 15, 16):
            note = "frontier shell F6 + F1/F2/F3"
        elif len(terms) == 1:
            note = "pure Fibonacci reset"
        elif len(terms) == 2 and gaps == [2]:
            note = "Lucas/min-gap pair"
        elif n in FRONTIER_NS:
            note = "active LRC comparison point"
        print(
            f"{n:<3} {mode_text(n):<9} {phi(n):>6}  "
            f"{zeck_text(n):<18} F{top_index}={top_value:<4} {payload:>5}       {note}"
        )
    print()
    print(
        "Reframe: n=14,15,16 are consecutive low-payload moves over the same F6=13 anchor."
    )
    print(
        "n=18 is the Lucas/min-gap pair F6+F4, and n=21 is the next pure Fibonacci reset."
    )
    print()


def print_gate_cover_atlas() -> None:
    print_header("LOCAL n-GATE COVER FAMILIES IN ZECKENDORF COORDINATES")
    print(
        "For each n, cover the endpoints owned by the mandatory n-gate using lower columns 1..n-1."
    )
    print("n   exact  forced columns                 #minimum covers  free part(s)")
    print("-" * 96)
    families = {n: exact_gate_cover_family(n) for n in (14, 15, 16)}
    for n, family in families.items():
        forced = set(family.forced)
        free_parts = sorted({tuple(speed for speed in cover if speed not in forced) for cover in family.covers})
        free_text = ", ".join(str(part) for part in free_parts)
        print(
            f"{n:<3} {str(family.exact_size):>5}  "
            f"{str(family.forced):<30} {len(family.covers):>7}         {free_text}"
        )
    print()
    print(
        "Exact new fact: n=14 factors as all seven odd lower columns plus exactly one arbitrary even bridge."
    )
    print(
        "The n=16 row is rigid and unique; n=15 is intermediate, with two free columns after its forced core."
    )
    print()


def print_n14_fan_digit_load() -> None:
    print_header("THE n=14 GATE FAN AS A FIBONACCI-CUBE TRANSVERSAL")
    family = exact_gate_cover_family(14)
    stats = cover_column_stats(14)
    forced = family.forced
    even_bridges = tuple(speed for speed in range(1, 14) if speed % 2 == 0)

    print("forced odd fan columns")
    print("p   mode      covers  private  Zeckendorf")
    print("-" * 96)
    for speed in forced:
        covers, private = stats[speed]
        print(
            f"{speed:<3} {mode_text(speed):<9} {covers:>6} {private:>8}  {zeck_text(speed)}"
        )
    print(f"\nforced digit load: {load_text(digit_load(forced))}")
    print()

    print("six minimum-cover choices: forced odd fan plus one even bridge")
    print("bridge  Zeckendorf       full cover digit load")
    print("-" * 96)
    for bridge in even_bridges:
        cover = tuple(sorted(forced + (bridge,)))
        print(f"{bridge:<7} {zeck_text(bridge):<16} {load_text(digit_load(cover))}")
    print()
    print(
        "Interpretation: private endpoints force the whole odd residue fan. The even bridge is a free Zeckendorf-cube coordinate locally, but the global LRC rows should break that freedom."
    )
    print()


def print_owner_debt_zeckendorf() -> None:
    print_header("OWNER DEBT LEDGER WITH ZECKENDORF SUPPORTS")
    for label, speeds in (
        ("seven-ladder", seven_ladder()),
        ("S380 gate ladder", s380_gate_ladder()),
    ):
        print(label)
        print("owner  labels  unique  mode       Zeckendorf")
        print("-" * 96)
        for owner, labels, unique in owner_debt_rows(speeds):
            print(
                f"{owner:<6} {labels:>6} {unique:>7}  {mode_text(owner):<9} {zeck_text(owner)}"
            )
        print()
    print(
        "Reframe: gate-heavy repairs export debt to owners with high Zeckendorf top indices, not just to larger magnitudes."
    )
    print()


def print_hypothesis() -> None:
    print_header("WORKING HYPOTHESIS")
    print(
        "HYP-1920: The n=14 LRC proof should be a Zeckendorf-shell row-cover certificate."
    )
    print()
    print(
        "Speed columns live as vertices of the Fibonacci cube. The mandatory 14-gate is F6+F1, and its local endpoint invoice forces every odd lower column; the only local freedom is one even bridge. This is the LRC analogue of the tournament path-fugacity lift: Zeckendorf gives the x=1 path coordinate, while the tournament OCF gives the x=2 independence evaluation."
    )
    print()
    print(
        "Proof target: branch on the forced odd fan, quotient the six even bridges as a local Fibonacci-cube fiber, then prove coarse denominator rows and exported owner-debt rows give a positive dual weight on every bridge fiber."
    )
    print()


def main() -> None:
    print_header("LRC ZECKENDORF BRIDGE")
    print("codex-2026-05-31-S460")
    print()
    print_old_thread_dictionary()
    print_frontier_shell()
    print_gate_cover_atlas()
    print_n14_fan_digit_load()
    print_owner_debt_zeckendorf()
    print_hypothesis()


if __name__ == "__main__":
    main()
