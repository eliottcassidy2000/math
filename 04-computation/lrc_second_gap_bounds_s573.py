#!/usr/bin/env python3
"""
lrc_second_gap_bounds_s573.py

codex-2026-06-03 S573

Sharpen the S572 second-gap audit with cheap necessary bounds before exact
maximin checking.

Question:
  If M(S) < 2/(2n-1), which clocks must already be blocked?

Methodology / Tournament Analysis note:
- The useful vertices are clocks and residue-shell obligations, not runners.
- Pairwise runner tournaments are lossy here: they forget whether a whole clock
  denominator is blocked.
- Alternate vertex sets considered: runners, antipodal shells, small clock
  denominators, n-clock ticks, and proof obligations.  This script uses proof
  obligations as vertices:
    D_q = "some speed is divisible by q";
    U_a = "unit antipodal shell {a,C-a} is hit";
    N_j = "the j/n clock has a runner at distance 0 or 1/n".
- Preserved predicate: possible strict sub-edge status M(S)<2/(2n-1).
- Challenged assumption: the finite residual should be enumerated first as
  runner sets.  The sharper quotient is a clock-blocker ledger followed by
  exact checking only on surviving rows.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import gcd


@dataclass(frozen=True)
class RowSummary:
    speeds: tuple[int, ...]
    maximin: Fraction
    witness: Fraction
    perfect_transversal: bool
    missed_nonunit: tuple[tuple[int, int], ...]
    doubled_shells: tuple[tuple[int, int], ...]
    zero_mod_c: int
    flipset: tuple[int, ...]


def primitive(combo: tuple[int, ...]) -> bool:
    g = 0
    for v in combo:
        g = gcd(g, v)
    return g == 1


def frac_part(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def norm(x: Fraction) -> Fraction:
    r = frac_part(x)
    return min(r, 1 - r)


def score_time(speeds: tuple[int, ...], t: Fraction) -> Fraction:
    return min(norm(Fraction(v) * t) for v in speeds)


def exact_maximin(speeds: tuple[int, ...]) -> tuple[Fraction, Fraction]:
    """Exact one-dimensional maximin over the tent arrangement.

    Local maxima occur at runner antipodes or at intersections of two runner
    tents, hence denominators v_i+v_j and |v_i-v_j|.
    """
    candidates: set[Fraction] = set()
    for i, a in enumerate(speeds):
        for m in range(a):
            t = Fraction(2 * m + 1, 2 * a)
            if 0 < t < 1:
                candidates.add(t)
        for b in speeds[i + 1 :]:
            for den in (a + b, abs(a - b)):
                if den <= 0:
                    continue
                for m in range(1, den):
                    candidates.add(Fraction(m, den))

    best = Fraction(0)
    best_t = Fraction(0)
    for t in candidates:
        score = score_time(speeds, t)
        if (score, -t) > (best, -best_t):
            best = score
            best_t = t
    return best, best_t


def small_denominator_blockers(speeds: tuple[int, ...]) -> tuple[int, ...]:
    """Missing q means t=1/q clears at least 1/q >= 2/(2n-1)."""
    k = len(speeds)
    return tuple(q for q in range(2, k + 1) if all(v % q for v in speeds))


def unit_shell_misses(speeds: tuple[int, ...]) -> tuple[int, ...]:
    """Missed unit shell gives the S553 inverse-clock 2/(2n-1) witness."""
    k = len(speeds)
    c = 2 * k + 1
    residues = {v % c for v in speeds}
    return tuple(
        a
        for a in range(1, k + 1)
        if gcd(a, c) == 1 and a not in residues and c - a not in residues
    )


def n_clock_unblocked_ticks(speeds: tuple[int, ...]) -> tuple[int, ...]:
    """Unblocked j means j/n already scores at least 2/n > 2/(2n-1)."""
    k = len(speeds)
    n = k + 1
    out = []
    for j in range(1, n):
        blocked = False
        for v in speeds:
            r = (v * j) % n
            if min(r, n - r) <= 1:
                blocked = True
                break
        if not blocked:
            out.append(j)
    return tuple(out)


def canonical_shell(r: int, c: int) -> int | None:
    r %= c
    if r == 0:
        return None
    return min(r, c - r)


def shell_summary(speeds: tuple[int, ...]) -> tuple[bool, tuple[tuple[int, int], ...], tuple[tuple[int, int], ...], int]:
    k = len(speeds)
    c = 2 * k + 1
    counts: Counter[int] = Counter()
    zero_mod_c = 0
    for v in speeds:
        shell = canonical_shell(v, c)
        if shell is None:
            zero_mod_c += 1
        else:
            counts[shell] += 1

    missed = []
    doubled = []
    for a in range(1, k + 1):
        if counts[a] == 0:
            missed.append((a, c - a))
        elif counts[a] > 1:
            doubled.append((a, c - a))
    perfect = zero_mod_c == 0 and not missed and not doubled
    missed_nonunit = tuple(pair for pair in missed if gcd(pair[0], c) != 1)
    return perfect, missed_nonunit, tuple(doubled), zero_mod_c


def flipset(speeds: tuple[int, ...]) -> tuple[int, ...]:
    k = len(speeds)
    c = 2 * k + 1
    residues = {v % c for v in speeds}
    return tuple(a for a in range(1, k + 1) if c - a in residues)


def analyze_survivor(speeds: tuple[int, ...]) -> RowSummary:
    maximin, witness = exact_maximin(speeds)
    perfect, missed_nonunit, doubled, zero_mod_c = shell_summary(speeds)
    return RowSummary(
        speeds=speeds,
        maximin=maximin,
        witness=witness,
        perfect_transversal=perfect,
        missed_nonunit=missed_nonunit,
        doubled_shells=doubled,
        zero_mod_c=zero_mod_c,
        flipset=flipset(speeds),
    )


def fmt_frac(x: Fraction) -> str:
    return f"{x.numerator}/{x.denominator}" if x.denominator != 1 else str(x.numerator)


def sample(rows: list[RowSummary], limit: int = 5) -> str:
    bits = []
    for row in rows[:limit]:
        bits.append(
            f"{row.speeds} M={fmt_frac(row.maximin)} t={fmt_frac(row.witness)} "
            f"perfect={row.perfect_transversal} flip={row.flipset} "
            f"missed_nonunit={row.missed_nonunit or '-'} doubled={row.doubled_shells or '-'}"
        )
    if len(rows) > limit:
        bits.append(f"... (+{len(rows) - limit} more)")
    return "\n".join(f"    {bit}" for bit in bits) if bits else "    -"


def scan_box(k: int, max_speed: int) -> None:
    n = k + 1
    edge = Fraction(2, 2 * n - 1)
    floor = Fraction(1, n)
    total = primitive_count = d_count = u_count = nc_count = survivor_count = 0
    below: list[RowSummary] = []
    open_gap: list[RowSummary] = []
    edge_or_above: list[RowSummary] = []
    perfect_below = 0
    missed_nonunit_hist: Counter[tuple[tuple[int, int], ...]] = Counter()
    flip_hist: Counter[tuple[int, ...]] = Counter()
    zero_hist: Counter[int] = Counter()

    for combo in combinations(range(1, max_speed + 1), k):
        total += 1
        if not primitive(combo):
            continue
        primitive_count += 1

        passes_d = not small_denominator_blockers(combo)
        passes_u = not unit_shell_misses(combo)
        passes_nc = not n_clock_unblocked_ticks(combo)
        d_count += int(passes_d)
        u_count += int(passes_u)
        nc_count += int(passes_nc)
        if not (passes_d and passes_u and passes_nc):
            continue

        survivor_count += 1
        row = analyze_survivor(combo)
        if row.maximin < edge:
            below.append(row)
            perfect_below += int(row.perfect_transversal)
            missed_nonunit_hist[row.missed_nonunit] += 1
            flip_hist[row.flipset] += 1
            zero_hist[row.zero_mod_c] += 1
            if row.maximin > floor:
                open_gap.append(row)
        else:
            edge_or_above.append(row)

    below_floor = sum(1 for row in below if row.maximin == floor)
    min_above = min((row.maximin for row in edge_or_above), default=None)
    closest_above = [row for row in edge_or_above if min_above is not None and row.maximin == min_above]

    print(f"Primitive exact box k={k}, max_speed={max_speed}, n={n}, edge={fmt_frac(edge)}")
    print(f"  total={total} primitive={primitive_count}")
    print(f"  gate D small-denominator blockers pass={d_count}")
    print(f"  gate U unit-shell coverage pass={u_count}")
    print(f"  gate N n-clock blockers pass={nc_count}")
    print(f"  all-gate survivors exact-checked={survivor_count}")
    print(f"  below_edge={len(below)} floor_count={below_floor} open_gap_count={len(open_gap)}")
    print(f"  perfect_below={perfect_below}/{len(below) if below else 1}")
    print(f"  zero_mod_c_hist={dict(sorted(zero_hist.items()))}")
    print(f"  missed_nonunit_hist={dict(missed_nonunit_hist)}")
    print(f"  flip_hist_top={flip_hist.most_common(8)}")
    if min_above is not None:
        print(f"  min_edge_or_above={fmt_frac(min_above)} count={len(closest_above)}")
    print("  below_edge_sample:")
    print(sample(below))
    print("  open_gap_sample:")
    print(sample(open_gap))
    print()


def main() -> None:
    print("Second-gap clock-blocker bounds (codex-2026-06-03 S573)")
    print("=" * 78)
    print("Necessary bounds for any strict sub-edge row M(S)<2/(2n-1):")
    print("  D_q: for every 2<=q<=n-1, some speed is divisible by q.")
    print("       Otherwise t=1/q scores at least 1/q >= 2/(2n-1).")
    print("  U_a: every unit antipodal shell mod C=2n-1 is hit.")
    print("       Otherwise S553's inverse clock scores exactly at the edge.")
    print("  N_j: every j/n clock is blocked by a runner at distance 0 or 1/n.")
    print("       Otherwise that n-clock scores at least 2/n > 2/(2n-1).")
    print()

    for k, max_speed in (
        (3, 60),
        (4, 40),
        (5, 32),
        (6, 26),
        (7, 22),
        (8, 20),
        (9, 18),
        (10, 17),
    ):
        scan_box(k, max_speed)

    print("Synthesis")
    print("-" * 78)
    print("The expanded boxes refute the earlier floor-only reading of S572.")
    print("There are genuine audited rows in the open interval (1/n, 2/(2n-1)),")
    print("beginning with the lifted flip-set row (1,5,6,11,16,17) at n=7,")
    print("where M=5/33, and lifted/nonunit-hole rows at n=8 with M=3/23.")
    print("The sharpened bound is therefore not 'below edge implies floor'.")
    print("It is the three-gate necessary ledger: any strict sub-edge row must")
    print("simultaneously block all small denominators, all unit shells, and all")
    print("n-clock ticks.  The remaining proof target is to classify that ledger's")
    print("finite lift/nonunit families and prove they stay at or above 1/n.")


if __name__ == "__main__":
    main()
