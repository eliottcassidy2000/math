#!/usr/bin/env python3
"""
lrc_egyptian_q_bridge_s397.py

codex-2026-05-31 S397

Bridge the S394 hyperbola/q-counterterm picture to Egyptian fractions and the
pure n=16 Lonely Runner row.

The guiding dictionary:

* Egyptian fractions split reciprocal mass by divisor congruences.
* q(x)=x-1/x is the anti-reciprocal finite-part counterterm.
* LRC endpoint protection at n=16 is also a reciprocal/congruence problem:
  a lower speed protects an endpoint when a residue is small modulo 16u.

This script does not claim a proof.  It builds tables that make the shared
obstruction visible.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from itertools import combinations
from math import gcd, log
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
S372 = SourceFileLoader(
    "lonely_runner_creative_multiroute_s372",
    str(ROOT / "04-computation" / "lonely_runner_creative_multiroute_s372.py"),
).load_module()
S391 = SourceFileLoader(
    "lrc_n16_dyadic_endpoint_formula_s391",
    str(ROOT / "04-computation" / "lrc_n16_dyadic_endpoint_formula_s391.py"),
).load_module()


N = 16
K = 15
OWNER = 16
LOWER_NINE_COVER = (1, 3, 5, 7, 8, 9, 11, 13, 15)
EVEN_BREAKERS = (2, 4, 6, 10, 12, 14)


@dataclass(frozen=True)
class EgyptianSplit:
    numerator: int
    denominator: int
    divisor: int
    b: int
    c: int

    @property
    def parents(self) -> tuple[int, int]:
        return (self.b, self.c)


@dataclass(frozen=True)
class CoordinateRow:
    coord: int
    split: EgyptianSplit | None
    product_sum_arities: tuple[int, ...]
    owner16_endpoint_count: int
    half_turn_missed: int
    in_nine_cover: bool
    is_even_breaker: bool
    reciprocal_partner: int | None


def divisors(value: int) -> list[int]:
    out: list[int] = []
    for d in range(1, int(value**0.5) + 1):
        if value % d == 0:
            out.append(d)
            if d * d != value:
                out.append(value // d)
    return sorted(out)


def egyptian_split(k: int, denominator: int) -> EgyptianSplit | None:
    """Return a split k/N = 1/b + 1/c via (kb-N)(kc-N)=N^2."""

    square = denominator * denominator
    target_residue = (-denominator) % k
    for d in divisors(square):
        if d % k != target_residue:
            continue
        other = square // d
        if (denominator + d) % k or (denominator + other) % k:
            continue
        b = (denominator + d) // k
        c = (denominator + other) // k
        if b > 0 and c > 0:
            return EgyptianSplit(k, denominator, d, min(b, c), max(b, c))
    return None


def collision_seed_arities(max_arity: int, max_target: int) -> dict[int, tuple[int, ...]]:
    """Targets of product-sum cores padded by ones: product = sum + ones."""

    by_target: dict[int, set[int]] = defaultdict(set)

    def rec(start: int, core: tuple[int, ...], product: int, total: int) -> None:
        if len(core) >= 2:
            ones = product - total
            arity = len(core) + ones
            if ones >= 0 and 2 <= arity <= max_arity and product <= max_target:
                by_target[product].add(arity)
        for value in range(start, max_target // product + 1):
            rec(value, core + (value,), product * value, total + value)

    rec(2, (), 1, 0)
    return {target: tuple(sorted(arities)) for target, arities in by_target.items()}


def reciprocal_partner_mod16(coord: int) -> int | None:
    if gcd(coord, N) != 1:
        return None
    for candidate in range(1, N):
        if (coord * candidate) % N == 1:
            return candidate
    return None


def half_turn_misses() -> dict[int, int]:
    system = S372.build_pattern_system(N)
    half = N // 2
    out: dict[int, int] = {}
    for coord in range(1, K + 1):
        vector = [0] * K
        vector[coord - 1] = half
        out[coord] = S372.score_vector(system, tuple(vector)).missed
    return out


def half_turn_pair_misses() -> tuple[tuple[int, tuple[int, int]], ...]:
    system = S372.build_pattern_system(N)
    half = N // 2
    rows: list[tuple[int, tuple[int, int]]] = []
    for a, b in combinations(range(1, K + 1), 2):
        vector = [0] * K
        vector[a - 1] = half
        vector[b - 1] = half
        rows.append((S372.score_vector(system, tuple(vector)).missed, (a, b)))
    return tuple(sorted(rows))


def coordinate_table() -> tuple[CoordinateRow, ...]:
    target_arities = collision_seed_arities(max_arity=20, max_target=32)
    misses = half_turn_misses()
    rows: list[CoordinateRow] = []
    for coord in range(1, N):
        rows.append(
            CoordinateRow(
                coord=coord,
                split=egyptian_split(N, coord),
                product_sum_arities=target_arities.get(coord, ()),
                owner16_endpoint_count=len(S391.protected_endpoints(OWNER, coord)),
                half_turn_missed=misses[coord],
                in_nine_cover=coord in LOWER_NINE_COVER,
                is_even_breaker=coord in EVEN_BREAKERS,
                reciprocal_partner=reciprocal_partner_mod16(coord),
            )
        )
    return tuple(rows)


def harmonic_shells() -> tuple[tuple[int, Fraction, float], ...]:
    rows: list[tuple[int, Fraction, float]] = []
    for m in (1, 2, 4, 8, 16, 32, 64, 128):
        shell = sum(Fraction(1, j) for j in range(m + 1, 2 * m + 1))
        rows.append((m, shell, float(shell) - log(2)))
    return tuple(rows)


def classify_speeds(label: str, speeds: tuple[int, ...]) -> tuple[str, Fraction, int, int, int]:
    report = S356.report(label, list(speeds))
    endpoints = S360.summarize(list(speeds))
    if report.forbidden_length < 1:
        classification = "positive_gap"
    elif report.boundary_witness_count:
        classification = "boundary_only"
    else:
        classification = "open_cover"
    return (
        classification,
        report.max_gap / report.threshold,
        endpoints.unprotected_count,
        endpoints.unique_endpoint_count,
        endpoints.boundary_modulus,
    )


def breaker_speed_sets() -> tuple[tuple[str, tuple[int, ...]], ...]:
    seed = set(LOWER_NINE_COVER + (16,))
    rows: list[tuple[str, tuple[int, ...]]] = []
    for omitted in EVEN_BREAKERS:
        speeds = tuple(sorted(seed | (set(EVEN_BREAKERS) - {omitted})))
        rows.append((f"nine+16 omit {omitted}", speeds))

    # The coordinates 10 and 15 were the best support-2 half-turn defect in S393.
    # Make two deliberately weird probes that keep the nine-cover but replace a
    # low even breaker by a 16-multiple tied to that Egyptian/half-turn pair.
    rows.append(("nine+16 + {2,4,6,10,160}", tuple(sorted(seed | {2, 4, 6, 10, 160}))))
    rows.append(("nine+16 + {2,4,6,15,160}", tuple(sorted(seed | {2, 4, 6, 15, 160}))))
    rows.append(("nine+16 + {2,4,10,15,160}", tuple(sorted(seed | {2, 4, 10, 15, 160}))))
    return tuple(rows)


def fmt_split(split: EgyptianSplit | None) -> str:
    if split is None:
        return "-"
    return f"1/{split.b}+1/{split.c}"


def yesno(value: bool) -> str:
    return "Y" if value else "-"


def print_coordinate_table() -> None:
    print("1. Coordinate ledger: Egyptian split vs n=16 LRC pressure")
    rows = coordinate_table()
    print(
        "   c  16/c split      ps-arity inv16 endpoints16 halfmiss nine evenbreak"
    )
    for row in rows:
        arities = ",".join(str(a) for a in row.product_sum_arities) or "-"
        inv = "-" if row.reciprocal_partner is None else str(row.reciprocal_partner)
        print(
            f"   {row.coord:>2} {fmt_split(row.split):>13} {arities:>9} "
            f"{inv:>5} {row.owner16_endpoint_count:>11} {row.half_turn_missed:>8} "
            f"{yesno(row.in_nine_cover):>4} {yesno(row.is_even_breaker):>9}"
        )
    print()
    print("   Read: '16/c split' uses the master Egyptian criterion")
    print("         16/c = 1/b + 1/c'.  The same coordinates are also LRC")
    print("         half-turn defect coordinates and endpoint-cover residues.")
    print("   The best S393 half-turn coordinate is c=15: it is the")
    print("   anti-residue -1 mod 16 and has the trivial split 16/15=1+1/15.")
    print("   The best support-2 pair (10,15) mixes that anti-residue with")
    print("   c=10, which is not two-unit-fraction splitable at numerator 16.")
    print()


def print_harmonic_shells() -> None:
    print("2. Best support-2 half-turn pairs")
    target_arities = collision_seed_arities(max_arity=20, max_target=32)
    print("   rank pair    missed  split-flags product-sum-arities")
    for rank, (missed, pair) in enumerate(half_turn_pair_misses()[:12], start=1):
        split_flags = "".join("S" if egyptian_split(N, coord) else "-" for coord in pair)
        arities = "/".join(
            ",".join(str(a) for a in target_arities.get(coord, ())) or "-"
            for coord in pair
        )
        print(f"   {rank:>4} {str(pair):>8} {missed:>7} {split_flags:>11} {arities:>20}")
    print()
    print("   The best pair is (10,15): an unsplit even residue plus the")
    print("   anti-residue -1 mod 16.  This looks like a q-pattern: one")
    print("   coordinate carries the unsplit debt, the other supplies the")
    print("   reciprocal counterterm.")
    print()

    print("3. Egyptian shells approximating ln 2")
    print("   m     sum_{m<j<=2m} 1/j        error over ln2")
    for m, shell, error in harmonic_shells():
        print(f"   {m:>3}   {str(shell):>24}   {error: .9f}")
    print()
    print("   These finite Egyptian fractions are the rectangle-sum version of")
    print("   the S394 hyperbola area.  At m=16 the shell [17,32] is the first")
    print("   denominator octave beyond the n=16 gate; the error is endpoint-debt")
    print("   in area language.")
    print()


def print_breaker_probes() -> None:
    print("4. Exact LRC probes from Egyptian/dyadic breaker sets")
    print("   label                         class          gap/th unprot endpoints modulus")
    for label, speeds in breaker_speed_sets():
        classification, gap_ratio, unprotected, endpoints, modulus = classify_speeds(label, speeds)
        print(
            f"   {label:<29} {classification:<13} {float(gap_ratio):>7.6f} "
            f"{unprotected:>6} {endpoints:>9} {modulus:>7}"
        )
    print()
    print("   The nine-cover plus a 16-gate still needs five of six even breakers.")
    print("   Every simple five-breaker choice stays positive-gap; high 16-multiple")
    print("   Egyptian probes also stay positive-gap with exposed endpoints.")
    print()


def print_synthesis() -> None:
    print("5. Synthesis")
    print("   Egyptian fractions and n=16 LRC share the same hidden question:")
    print("     can reciprocal mass be split without leaving a residue?")
    print("   For k/N = 1/b + 1/c the residue is a divisor congruence")
    print("     d | N^2, d == -N mod k.")
    print("   For an LRC endpoint, the residue is an inequality")
    print("     |p*(16m+eps) - 16*a*u| < u.")
    print("   The S394 q-counterterm is the finite-part version of the same")
    print("   move: replace an infinite reciprocal object by a finite triangular")
    print("   one and keep track of the anti-reciprocal defect x-1/x.")
    print()
    print("   Creative n=16 proof angle:")
    print("     prove that any all-protected endpoint system with a 16-gate")
    print("     induces an Egyptian split ledger across dyadic layers; the")
    print("     unsplit residues must be paid by even breaker speeds, and those")
    print("     breakers create a q-like negative counterterm visible as a")
    print("     positive LRC gap or a private endpoint leaf.")


def main() -> None:
    print("LRC/Egyptian/q bridge (codex-2026-05-31 S397)")
    print("n=16 Lonely Runner threshold = 1/16")
    print()
    print_coordinate_table()
    print_harmonic_shells()
    print_breaker_probes()
    print_synthesis()


if __name__ == "__main__":
    main()
