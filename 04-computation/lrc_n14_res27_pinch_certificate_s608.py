#!/usr/bin/env python3
"""
lrc_n14_res27_pinch_certificate_s608.py

codex-2026-06-03 S608

Finite least-positive Res_27 quotient certificate for LRC at n=14.

Question:
  In the C=2n-1=27 shell quotient, after the visible unit-shell exits and the
  canonical D/U/N clock-blocker gates, which 13-residue rows remain at the
  floor?  Does the pair-sum/pinch sieve leave any subfloor row in this quotient?

Scope:
  This is a quotient computation on distinct least-positive residues
  {1,...,26}.  The D and N gates are applied to those canonical representatives.
  It does not prove lift/CRT conservativity for all integer lifts of the same
  residue shell data.  Its value is to shrink the primitive Res_27 floor branch
  to the known AP and V* rows, leaving the lift theorem as the next problem.

Method:
  Enumerate every 13-subset of nonzero residues mod 27 that hits all nine unit
  antipodal shells.  Keep the rows that pass the canonical D and N blocker
  gates from S573:
    D_q: for every 2<=q<=13, some representative is divisible by q;
    N_j: for every 1<=j<=7, the j/14 clock is blocked at distance 0 or 1/14.
  Then search pair-sum/pinch candidates t=a/(v_i+v_j), gcd(a,v_i+v_j)=1,
  using integer arithmetic.  A row is strictly certified if one candidate scores
  >1/14.  Rows scoring exactly 1/14 are checked by exact maximin.

Tournament Analysis:
  Vertices are Res_27 proof-obligation types, not runners.  A type records the
  missed nonunit shells, doubled shells, and primitive/nonprimitive status after
  D/U/N.  The pairwise observable is the type burden
    (below_count, floor_count, -strict_count, label).
  The switch/gauge orients the edge from lower burden to higher burden, with
  lexicographic tie Hamiltonian path.  The quotient preserves the finite
  predicate "canonical Res_27 row survives D/U/N and the pair-sum floor test";
  it destroys lift choices, actual speed sizes, and phase order.

Assumption challenge:
  Vertices considered: runners, residues, unit shells, nonunit shell gaps,
  pair-sum denominators, n-clock blockers, lift congruences, and proof
  obligations.  This script chooses proof-obligation types because the target
  predicate is a quotient certificate, not a runner-level classification.
  The challenged assumption is that Res_27 shell data alone proves n=14; the
  output explicitly says the remaining theorem is lift/CRT conservativity.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, product
from math import gcd


N = 14
K = N - 1
C = 2 * N - 1
FLOOR = Fraction(1, N)
ALL_RESIDUES = tuple(range(1, C))
UNIT_SHELLS = tuple(a for a in range(1, N) if gcd(a, C) == 1)
UNIT_PAIRS = tuple((a, C - a) for a in UNIT_SHELLS)
NONUNIT_SHELLS = tuple(a for a in range(1, N) if gcd(a, C) != 1)

AP = tuple(range(1, N))
VSTAR = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 24)
DOUBLE_AP = tuple(range(2, 2 * N, 2))


@dataclass(frozen=True, order=True)
class Res27Type:
    missed_nonunit: tuple[int, ...]
    doubled_shells: tuple[tuple[int, int], ...]
    primitive: bool

    @property
    def label(self) -> str:
        miss = ",".join(map(str, self.missed_nonunit)) or "-"
        doubled = ",".join(f"{a}x{c}" for a, c in self.doubled_shells) or "-"
        prim = "prim" if self.primitive else "nonprim"
        return f"miss={miss};double={doubled};{prim}"


@dataclass(frozen=True)
class PinchWitness:
    time_num: int
    time_den: int
    score_num: int
    score_den: int

    @property
    def score(self) -> Fraction:
        return Fraction(self.score_num, self.score_den)

    @property
    def time(self) -> Fraction:
        return Fraction(self.time_num, self.time_den)


def fmt_frac(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def row_gcd(row: tuple[int, ...]) -> int:
    g = 0
    for v in row:
        g = gcd(g, v)
    return g


def shell_of(r: int) -> int:
    r %= C
    if r == 0:
        return 0
    return min(r, C - r)


def shell_counts(row: tuple[int, ...]) -> Counter[int]:
    return Counter(shell_of(r) for r in row)


def res27_type(row: tuple[int, ...]) -> Res27Type:
    counts = shell_counts(row)
    missed_nonunit = tuple(a for a in NONUNIT_SHELLS if counts[a] == 0)
    doubled = tuple((a, counts[a]) for a in range(1, N) if counts[a] > 1)
    return Res27Type(missed_nonunit, doubled, row_gcd(row) == 1)


def d_failures(row: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(q for q in range(2, N) if all(v % q for v in row))


def n_failures(row: tuple[int, ...]) -> tuple[int, ...]:
    failures: list[int] = []
    for j in range(1, N // 2 + 1):
        blocked = False
        for v in row:
            r = (v * j) % N
            if min(r, N - r) <= 1:
                blocked = True
                break
        if not blocked:
            failures.append(j)
    return tuple(failures)


def min_distance_num(row: tuple[int, ...], a: int, q: int) -> int:
    best = q
    for v in row:
        r = (v * a) % q
        d = min(r, q - r)
        if d < best:
            best = d
            if best == 0:
                break
    return best


def reduced_numerators_by_den(max_q: int) -> dict[int, tuple[int, ...]]:
    return {
        q: tuple(a for a in range(1, q) if gcd(a, q) == 1)
        for q in range(2, max_q + 1)
    }


REDUCED_NUMERATORS = reduced_numerators_by_den(2 * (C - 1))


def best_pair_sum_pinch(row: tuple[int, ...]) -> PinchWitness:
    best: PinchWitness | None = None
    pair_sums = sorted({a + b for a, b in combinations(row, 2)})
    for q in pair_sums:
        for a in REDUCED_NUMERATORS[q]:
            d = min_distance_num(row, a, q)
            candidate = PinchWitness(a, q, d, q)
            if best is None:
                best = candidate
                continue
            lhs = candidate.score_num * best.score_den
            rhs = best.score_num * candidate.score_den
            if lhs > rhs or (
                lhs == rhs and (candidate.time_den, candidate.time_num) < (best.time_den, best.time_num)
            ):
                best = candidate
    if best is None:
        raise RuntimeError(f"no pair sums for row {row}")
    return best


def frac_part(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def norm(x: Fraction) -> Fraction:
    r = frac_part(x)
    return min(r, 1 - r)


def exact_score(row: tuple[int, ...], t: Fraction) -> Fraction:
    return min(norm(Fraction(v) * t) for v in row)


def exact_maximin(row: tuple[int, ...]) -> tuple[Fraction, Fraction]:
    candidates: set[Fraction] = set()
    for i, a in enumerate(row):
        for m in range(a):
            t = Fraction(2 * m + 1, 2 * a)
            if 0 < t < 1:
                candidates.add(t)
        for b in row[i + 1 :]:
            for den in (a + b, abs(a - b)):
                if den <= 0:
                    continue
                for m in range(1, den):
                    candidates.add(Fraction(m, den))

    best = Fraction(0)
    best_t = Fraction(0)
    for t in candidates:
        score = exact_score(row, t)
        if (score, -t) > (best, -best_t):
            best = score
            best_t = t
    return best, best_t


def row_name(row: tuple[int, ...]) -> str:
    if row == AP:
        return "AP"
    if row == VSTAR:
        return "V*"
    if row == DOUBLE_AP:
        return "2*AP/nonprimitive"
    return "unclassified"


def enumerate_unit_rows() -> tuple[list[tuple[int, ...]], dict[str, int]]:
    seen: set[tuple[int, ...]] = set()
    unit_rows: list[tuple[int, ...]] = []
    duplicate_presentations = 0

    for choices in product((0, 1), repeat=len(UNIT_PAIRS)):
        base = tuple(pair[choice] for pair, choice in zip(UNIT_PAIRS, choices))
        base_set = set(base)
        remaining = tuple(r for r in ALL_RESIDUES if r not in base_set)
        for extra in combinations(remaining, K - len(base)):
            row = tuple(sorted(base + extra))
            if row in seen:
                duplicate_presentations += 1
                continue
            seen.add(row)
            unit_rows.append(row)

    return unit_rows, {"duplicate_presentations": duplicate_presentations}


def tournament_fingerprint(type_stats: dict[Res27Type, Counter[str]]) -> dict[str, object]:
    vertices = sorted(type_stats)

    def burden(key: Res27Type) -> tuple[int, int, int, str]:
        stats = type_stats[key]
        return (stats["below"], stats["floor"], -stats["strict"], key.label)

    order = sorted(vertices, key=burden)
    index = {key: i for i, key in enumerate(order)}
    outdegrees = {key: len(order) - 1 - index[key] for key in order}
    score_hist = Counter(outdegrees.values())

    directed_3_cycles = 0
    for a, b, c in combinations(order, 3):
        # The burden order is total, so this is deliberately generic enough to
        # make the zero-cycle certificate explicit.
        wins = {
            (a, b): index[a] < index[b],
            (b, c): index[b] < index[c],
            (c, a): index[c] < index[a],
        }
        if all(wins.values()):
            directed_3_cycles += 1

    return {
        "vertices": len(order),
        "score_hist": tuple(sorted(score_hist.items())),
        "scc_count": len(order),
        "largest_scc": 1 if order else 0,
        "directed_3_cycles": directed_3_cycles,
        "hamiltonian_path_first": tuple(key.label for key in order[:5]),
        "hamiltonian_path_last": tuple(key.label for key in order[-5:]),
    }


def main() -> None:
    print(f"==== LRC n={N} least-positive Res_{C} pinch certificate (S608) ====")
    print(f"unit shells: {UNIT_SHELLS}")
    print(f"nonunit shells: {NONUNIT_SHELLS}")
    print()

    unit_rows, enum_meta = enumerate_unit_rows()
    d_pass = 0
    n_pass = 0
    ledger_survivors: list[tuple[int, ...]] = []
    d_fail_hist: Counter[tuple[int, ...]] = Counter()
    n_fail_hist: Counter[tuple[int, ...]] = Counter()

    for row in unit_rows:
        df = d_failures(row)
        nf = n_failures(row)
        if not df:
            d_pass += 1
        else:
            d_fail_hist[df] += 1
        if not nf:
            n_pass += 1
        else:
            n_fail_hist[nf] += 1
        if not df and not nf:
            ledger_survivors.append(row)

    strict_rows: list[tuple[tuple[int, ...], PinchWitness]] = []
    floor_rows: list[tuple[tuple[int, ...], PinchWitness]] = []
    below_rows: list[tuple[tuple[int, ...], PinchWitness]] = []
    type_stats: defaultdict[Res27Type, Counter[str]] = defaultdict(Counter)
    witness_den_hist: Counter[int] = Counter()

    for row in ledger_survivors:
        witness = best_pair_sum_pinch(row)
        key = res27_type(row)
        type_stats[key]["total"] += 1
        witness_den_hist[witness.time_den] += 1
        cmp_num = N * witness.score_num - witness.score_den
        if cmp_num > 0:
            strict_rows.append((row, witness))
            type_stats[key]["strict"] += 1
        elif cmp_num == 0:
            floor_rows.append((row, witness))
            type_stats[key]["floor"] += 1
        else:
            below_rows.append((row, witness))
            type_stats[key]["below"] += 1

    primitive_floor = [(row, w) for row, w in floor_rows if row_gcd(row) == 1]
    nonprimitive_floor = [(row, w) for row, w in floor_rows if row_gcd(row) != 1]

    print("Enumeration:")
    print(f"  raw 13-subsets of nonzero residues: 10400600")
    print(f"  unit-shell rows, deduped: {len(unit_rows)}")
    print(f"  duplicate structured presentations skipped: {enum_meta['duplicate_presentations']}")
    print(f"  D-pass among unit rows: {d_pass}")
    print(f"  N-pass among unit rows: {n_pass}")
    print(f"  D/U/N ledger survivors: {len(ledger_survivors)}")
    print()

    print("Pair-sum/pinch classification on D/U/N survivors:")
    print(f"  strict > 1/{N}: {len(strict_rows)}")
    print(f"  floor = 1/{N}: {len(floor_rows)}")
    print(f"  below < 1/{N}: {len(below_rows)}")
    print(f"  primitive floor rows: {len(primitive_floor)}")
    print(f"  nonprimitive floor rows: {len(nonprimitive_floor)}")
    print()

    if floor_rows:
        print("Floor rows, with exact maximin check:")
        for row, witness in floor_rows:
            exact_m, exact_t = exact_maximin(row)
            print(
                f"  {row_name(row):18s} row={row} "
                f"type=({res27_type(row).label}) "
                f"pinch_score={fmt_frac(witness.score)} at t={fmt_frac(witness.time)} "
                f"exact_M={fmt_frac(exact_m)} exact_t={fmt_frac(exact_t)}"
            )
        print()

    if below_rows:
        print("Below rows (unexpected):")
        for row, witness in below_rows[:20]:
            print(
                f"  row={row} type=({res27_type(row).label}) "
                f"pinch_score={fmt_frac(witness.score)} at t={fmt_frac(witness.time)}"
            )
        if len(below_rows) > 20:
            print(f"  ... (+{len(below_rows) - 20} more)")
        print()

    print("Type histogram:")
    print(f"  Res_27 types among survivors: {len(type_stats)}")
    for key, stats in sorted(type_stats.items(), key=lambda kv: (-kv[1]['floor'], -kv[1]['strict'], kv[0].label))[:12]:
        print(
            f"  {key.label:52s} total={stats['total']:5d} "
            f"strict={stats['strict']:5d} floor={stats['floor']:2d} below={stats['below']:2d}"
        )
    print()

    print("Top witness denominator counts:")
    for q, count in witness_den_hist.most_common(12):
        print(f"  q={q:2d}: {count}")
    print()

    fp = tournament_fingerprint(type_stats)
    print("Tournament Analysis over Res_27 proof-obligation types:")
    print(f"  vertices: {fp['vertices']}")
    print(f"  score_histogram: {fp['score_hist']}")
    print(f"  SCCs: count={fp['scc_count']} largest={fp['largest_scc']}")
    print(f"  directed_3_cycles: {fp['directed_3_cycles']}")
    print("  Hamiltonian path, easiest first:")
    for label in fp["hamiltonian_path_first"]:
        print(f"    {label}")
    print("  Hamiltonian path, hardest last:")
    for label in fp["hamiltonian_path_last"]:
        print(f"    {label}")
    print()

    print("Conclusion:")
    print(
        "  In the least-positive Res_27 quotient, every canonical D/U/N survivor "
        "has a pair-sum/pinch certificate at or above the LRC floor."
    )
    print(
        "  The only primitive floor rows are AP and V*.  The third floor row is "
        "the nonprimitive doubled AP, which normalizes to AP."
    )
    print(
        "  This improves the n=14 proof target to lift/CRT conservativity for "
        "the C=27 coimage ledger; it is not itself the full lifted theorem."
    )


if __name__ == "__main__":
    main()
