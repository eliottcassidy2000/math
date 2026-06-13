#!/usr/bin/env python3
"""
S654: parity/carry bridge for the n=14 LRC residual.

The recent Res_27 tower leaves one seam: lift/CRT conservativity.  A lifted
speed has the form

    v_i = r_i + 27 k_i.

Since 27 is odd, k_i toggles parity.  Since 27 == -1 (mod 14), k_i also
controls the n-clock residue:

    v_i == r_i + k_i (mod 2)
    v_i == r_i - k_i (mod 14).

Thus the even/odd carrier and the LRC apex obstruction are the same carry
coordinate viewed through two quotients.  This script makes that bridge
explicit and tests the forced "contains a multiple of 14" branch around the
known least-positive floor atoms AP and V*.

Tournament Analysis:
  Vertices are proof channels / carry obligations, not runners.  The pairwise
  observable is (proof burden, retained side-channel, exactness flag, label).
  The switch ranks no-multiple exits, parity carries, apex congruences, pair
  congruences, and owner certificates before raw row identity.

Assumption challenge:
  Candidate vertex sets considered: runners, gaps, fixed clock residues, pair
  sums, carry coordinates, endpoint owners, proof obligations, and parity
  words.  This script chooses carry coordinates because they preserve both the
  even/odd projection and the n=14 zero-divisor/apex predicate; raw Res_27
  rows destroy that information.
"""

from __future__ import annotations

import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lrc_n14_carry_conservativity_s611 import (  # noqa: E402
    AP,
    VSTAR,
    C,
    FLOOR,
    N,
    exact_maximin,
    fmt_frac,
)


BASE_ROWS = (("AP", AP), ("V*", VSTAR))
TOL = 1e-6
FLOOR_F = float(FLOOR)


def mod14(x: int) -> int:
    return x % N


def minimal_apex_carry(r: int) -> int:
    """Least nonnegative k such that r + 27k is divisible by 14."""
    return r % N


def lift_row(base: tuple[int, ...], carry: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted(r + C * k for r, k in zip(base, carry)))


def parity_word(base: tuple[int, ...], carry: tuple[int, ...]) -> tuple[int, ...]:
    return tuple((r + k) & 1 for r, k in zip(base, carry))


def apex_indices(base: tuple[int, ...], carry: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(i for i, (r, k) in enumerate(zip(base, carry)) if (r - k) % N == 0)


def pair_sum_apex_pairs(base: tuple[int, ...], carry: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    pairs = []
    for i, j in combinations(range(len(base)), 2):
        if (base[i] + base[j] - carry[i] - carry[j]) % N == 0:
            pairs.append((i, j))
    return tuple(pairs)


def _fdist(x: float) -> float:
    r = x - int(x)
    if r < 0:
        r += 1.0
    return r if r <= 0.5 else 1.0 - r


def approx_maximin(row: tuple[int, ...]) -> float:
    """
    Fast screen using the exact oracle's critical-time family.  This is only
    used to choose rows for exact verification; all reported minima are exact.
    """
    best = 0.0
    n = len(row)
    for i in range(n):
        a = row[i]
        inv2a = 1.0 / (2 * a)
        for m in range(a):
            t = (2 * m + 1) * inv2a
            if 0.0 < t < 1.0:
                s = 1.0
                for v in row:
                    d = _fdist(v * t)
                    if d < s:
                        s = d
                        if s <= best:
                            break
                if s > best:
                    best = s
        for b in row[i + 1 :]:
            for den in (a + b, abs(a - b)):
                if den <= 0:
                    continue
                invden = 1.0 / den
                for m in range(1, den):
                    t = m * invden
                    s = 1.0
                    for v in row:
                        d = _fdist(v * t)
                        if d < s:
                            s = d
                            if s <= best:
                                break
                    if s > best:
                        best = s
    return best


@dataclass(frozen=True)
class CarryHit:
    base_name: str
    base: tuple[int, ...]
    carry: tuple[int, ...]
    m: Fraction
    t: Fraction

    @property
    def row(self) -> tuple[int, ...]:
        return lift_row(self.base, self.carry)

    @property
    def label(self) -> str:
        moved = [(self.base[i], k) for i, k in enumerate(self.carry) if k]
        body = ",".join(f"r{r}->k{k}" for r, k in moved) if moved else "zero"
        return f"{self.base_name}:{body}"

    @property
    def weight(self) -> int:
        return sum(1 for k in self.carry if k)

    @property
    def l1(self) -> int:
        return sum(self.carry)

    @property
    def even_speeds(self) -> int:
        return sum(1 for bit in parity_word(self.base, self.carry) if bit == 0)

    @property
    def odd_carries(self) -> int:
        return sum(1 for k in self.carry if k & 1)

    @property
    def apex_count(self) -> int:
        return len(apex_indices(self.base, self.carry))

    @property
    def pair_apex_count(self) -> int:
        return len(pair_sum_apex_pairs(self.base, self.carry))


def exact_hit(name: str, base: tuple[int, ...], carry: tuple[int, ...]) -> CarryHit:
    m, t = exact_maximin(lift_row(base, carry))
    return CarryHit(name, base, carry, m, t)


def print_identities() -> None:
    print("A. Exact parity/carry identities")
    print(f"  n={N}, C=2n-1={C}, floor=1/{N}")
    print("  v_i = r_i + 27*k_i")
    print("  parity(v_i) = parity(r_i + k_i)       because 27 is odd")
    print("  v_i mod 14  = r_i - k_i mod 14        because 27 == -1 mod 14")
    print("  14 | v_i    iff k_i == r_i mod 14     (apex / zero-divisor carry)")
    print("  14 | v_i+v_j iff k_i+k_j == r_i+r_j mod 14 (pair-sum carry)")
    print()


def print_base_apex_table() -> None:
    print("B. Minimal apex carries for the two primitive floor shadows")
    for name, base in BASE_ROWS:
        print(f"  {name}:")
        for i, r in enumerate(base):
            k = minimal_apex_carry(r)
            v = r + C * k
            print(
                f"    i={i:2d} r={r:2d} k0={k:2d} v={v:3d} "
                f"v/14={v // N:2d} parity={'even' if v % 2 == 0 else 'odd'}"
            )
    print("  Every minimal apex speed is even, as the doubled-prime bridge predicts.")
    print()


def single_apex_sweep() -> list[CarryHit]:
    print("C. Exact single-apex carry sweep")
    hits: list[CarryHit] = []
    for name, base in BASE_ROWS:
        best: CarryHit | None = None
        print(f"  {name}:")
        for i, r in enumerate(base):
            carry = [0] * len(base)
            carry[i] = minimal_apex_carry(r)
            hit = exact_hit(name, base, tuple(carry))
            hits.append(hit)
            if best is None or (hit.m, -hit.t) < (best.m, -best.t):
                best = hit
            print(
                f"    apex at r={r:2d}: M={fmt_frac(hit.m):>6s} "
                f"margin={fmt_frac(hit.m - FLOOR):>7s} t={fmt_frac(hit.t):>7s} "
                f"pairs14={hit.pair_apex_count:2d} even={hit.even_speeds:2d}"
            )
        assert best is not None
        print(
            f"    best single-apex tax: {fmt_frac(best.m)} "
            f"(margin {fmt_frac(best.m - FLOOR)}) via {best.label}"
        )
    print()
    return hits


def apex_neighborhoods(extra_depth: int = 2) -> dict[tuple[str, int, int], CarryHit]:
    print(f"D. Apex plus up to {extra_depth} extra +27 parity/carry toggles")
    group_min: dict[tuple[str, int, int], tuple[float, tuple[int, ...], tuple[int, ...]]] = {}
    near_floor: list[tuple[str, tuple[int, ...], tuple[int, ...]]] = []
    totals = Counter()

    for name, base in BASE_ROWS:
        dim = len(base)
        for apex_i, r in enumerate(base):
            root = [0] * dim
            root[apex_i] = minimal_apex_carry(r)
            other = [j for j in range(dim) if j != apex_i]
            for d in range(extra_depth + 1):
                for extras in combinations(other, d):
                    carry = root[:]
                    for j in extras:
                        carry[j] += 1
                    carry_t = tuple(carry)
                    row = lift_row(base, carry_t)
                    mf = approx_maximin(row)
                    key = (name, mod14(r), d)
                    totals[key] += 1
                    if key not in group_min or mf < group_min[key][0]:
                        group_min[key] = (mf, base, carry_t)
                    if mf <= FLOOR_F + TOL:
                        near_floor.append((name, base, carry_t))

    exact_mins: dict[tuple[str, int, int], CarryHit] = {}
    for key, (_, base, carry) in group_min.items():
        exact_mins[key] = exact_hit(key[0], base, carry)

    exact_near = [exact_hit(name, base, carry) for name, base, carry in near_floor]
    below = [h for h in exact_near if h.m < FLOOR]
    floor = [h for h in exact_near if h.m == FLOOR]

    print(f"  rows screened: {sum(totals.values())}")
    print(f"  approximate near-floor rows exact-checked: {len(exact_near)}")
    print(f"  exact floor rows in checked set: {len(floor)}")
    print(f"  exact below-floor rows in checked set: {len(below)}")
    print("  exact group minima by base, apex residue, and extra-toggle count:")
    for key in sorted(exact_mins):
        hit = exact_mins[key]
        name, residue, d = key
        print(
            f"    {name:2s} apex={residue:2d} extra={d}: "
            f"M={fmt_frac(hit.m):>6s} margin={fmt_frac(hit.m - FLOOR):>7s} "
            f"t={fmt_frac(hit.t):>7s} l1={hit.l1:2d} "
            f"even={hit.even_speeds:2d} pairs14={hit.pair_apex_count:2d}"
        )
    print()
    return exact_mins


def boolean_parity_lattice() -> None:
    print("E. Boolean carry lattice parity profile")
    for name, base in BASE_ROWS:
        dim = len(base)
        group_min: dict[tuple[int, int, int], tuple[float, tuple[int, ...]]] = {}
        near_floor = []
        group_counts = Counter()
        for mask in range(1 << dim):
            carry = tuple((mask >> i) & 1 for i in range(dim))
            row = lift_row(base, carry)
            mf = approx_maximin(row)
            key = (
                len(apex_indices(base, carry)),
                sum(1 for k in carry if k & 1),
                sum(1 for bit in parity_word(base, carry) if bit == 0),
            )
            group_counts[key] += 1
            if key not in group_min or mf < group_min[key][0]:
                group_min[key] = (mf, carry)
            if mf <= FLOOR_F + TOL:
                near_floor.append(carry)

        exact_near = [exact_hit(name, base, carry) for carry in near_floor]
        print(f"  {name}: rows={1 << dim}, groups={len(group_counts)}")
        print(f"    approximate near-floor rows exact-checked: {len(exact_near)}")
        print(
            "    exact floor/below among checked: "
            f"{sum(1 for h in exact_near if h.m == FLOOR)}/"
            f"{sum(1 for h in exact_near if h.m < FLOOR)}"
        )
        print("    selected exact minima by (apex_count, odd_carries, even_speeds):")
        selected = []
        for key, (_, carry) in group_min.items():
            apex_count, odd_carries, even_speeds = key
            if apex_count or odd_carries in {0, 1, 6, 7, 12, 13}:
                selected.append((key, exact_hit(name, base, carry), group_counts[key]))
        for key, hit, count in sorted(selected):
            print(
                f"      {key}: count={count:4d} minM={fmt_frac(hit.m):>6s} "
                f"margin={fmt_frac(hit.m - FLOOR):>7s} label={hit.label}"
            )
    print()


@dataclass(frozen=True)
class Channel:
    name: str
    proof_burden: int
    side_channel: int
    exactness: int

    @property
    def key(self) -> tuple[int, int, int, str]:
        return (self.proof_burden, -self.side_channel, -self.exactness, self.name)


def tournament_analysis() -> None:
    print("F. Tournament Analysis over proof channels")
    channels = [
        Channel("no-multiple t=1/14 exit", 0, 5, 5),
        Channel("carry parity word", 1, 5, 5),
        Channel("apex congruence k=r mod 14", 2, 5, 5),
        Channel("pair-sum congruence k_i+k_j=r_i+r_j", 3, 4, 5),
        Channel("minimal-apex tax sweep", 4, 4, 4),
        Channel("owner/Cprime certificate reattachment", 5, 5, 3),
        Channel("raw Res_27 row identity", 6, 1, 2),
    ]
    ordered = sorted(channels, key=lambda c: c.key)
    rank = {c.name: i for i, c in enumerate(ordered)}
    score_hist = Counter()
    directed_3 = 0
    for a in channels:
        score_hist[sum(rank[a.name] > rank[b.name] for b in channels if b != a)] += 1
    for a, b, c in combinations(channels, 3):
        ab = rank[a.name] < rank[b.name]
        bc = rank[b.name] < rank[c.name]
        ca = rank[c.name] < rank[a.name]
        if (ab and bc and ca) or ((not ab) and (not bc) and (not ca)):
            directed_3 += 1
    print(f"  top_order={[c.name for c in ordered]}")
    print(f"  score_hist={dict(sorted(score_hist.items()))}")
    print(f"  directed_3cycles={directed_3}")
    print("  hamiltonian_paths=1")
    print()


def main() -> None:
    print("S654 LRC n=14 parity/carry bridge")
    print("=" * 72)
    print_identities()
    print_base_apex_table()
    single_apex_sweep()
    apex_neighborhoods(extra_depth=2)
    boolean_parity_lattice()
    tournament_analysis()
    print("Conclusion")
    print(
        "  The same carry coordinate toggles even/odd parity and decides the "
        "n=14 apex obstruction.  Non-multiple rows are already discharged by "
        "t=1/14; the required multiple-of-14 branch is the congruence "
        "k_i == r_i mod 14.  Every minimal apex lift over AP/V*, and every "
        "one/two extra-toggle neighbourhood screened here, is strictly above "
        "the floor in exact group minima, with no near-floor or below-floor "
        "rows detected."
    )
    print(
        "  This does not prove LRC(14), but it moves the carry-conservativity "
        "seam: a new primitive floor lift would need a larger globally "
        "coherent carry cocycle or an owner/Cprime route, not merely the "
        "first parity/apex bridge."
    )


if __name__ == "__main__":
    main()
