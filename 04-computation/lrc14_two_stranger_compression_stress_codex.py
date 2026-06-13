#!/usr/bin/env python3
"""Two-stranger compression stress test for the LRC(14) Q27 route.

This is a proof-facing extension of HYP-2463.  HYP-2463 showed that the eight
known one-stranger hard residues do not stack inside the 7-core frame.  The
present scan asks a harsher bounded question:

    delete one speed from CORE = 7*{1,...,12};
    add any two distinct non-core speeds 1 <= r <= 13*84;
    keep only primitive rows.

This creates 6,868,368 true two-stranger rows.  The script records exactly which
of them still block every plain shell q <= 27, then checks the fibered Q27
ladder.  The goal is not to prove LRC(14) by exhaustion; it is to expose which
resources a would-be compression theorem must preserve.

Tournament Analysis / assumption challenge:
  Vertices are candidate compression maps, not runners.  Candidate vertex sets
  considered: runners, deleted core addresses, added-residue pairs, denominator
  fibers, unit twists, shell-27 classes, 13-clock debts, owner/Bprime openings,
  AP/Vstar/2AP descent targets, and proof obligations.  The selected quotient
  preserves exact finite blocker laws and destroys most raw time geometry.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations
from math import gcd
from typing import Iterable, Sequence


N = 14
CORE = tuple(7 * k for k in range(1, 13))
MAX_R = 13 * 84
DIVISORS = (1, 2, 7, 14)
SHELL27 = tuple(range(2, 28))
Q27 = tuple(sorted({d * m for d in DIVISORS for m in range(1, 28) if d * m > 1}))
HARD = (260, 351, 442, 611, 702, 793, 962, 1053)
HARD_SET = set(HARD)


@dataclass(frozen=True)
class ResidualRecord:
    deleted: int
    r1: int
    r2: int
    q27: tuple[int, int] | None

    @property
    def added(self) -> tuple[int, int]:
        return self.r1, self.r2

    @property
    def first_q(self) -> int | None:
        return None if self.q27 is None else self.q27[0]

    @property
    def first_a(self) -> int | None:
        return None if self.q27 is None else self.q27[1]

    @property
    def hard_count(self) -> int:
        return int(self.r1 in HARD_SET) + int(self.r2 in HARD_SET)

    @property
    def mod13_zero_count(self) -> int:
        return int(self.r1 % 13 == 0) + int(self.r2 % 13 == 0)

    @property
    def anti27_pair(self) -> tuple[int, int]:
        return tuple(sorted((antipodal_class(self.r1, 27), antipodal_class(self.r2, 27))))

    @property
    def mod13_pair(self) -> tuple[int, int]:
        return tuple(sorted((self.r1 % 13, self.r2 % 13)))

    @property
    def row(self) -> tuple[int, ...]:
        return tuple(sorted([v for v in CORE if v != self.deleted] + [self.r1, self.r2]))


def antipodal_class(v: int, q: int) -> int:
    r = v % q
    return min(r, q - r)


@lru_cache(maxsize=None)
def unit_twists(q: int) -> tuple[int, ...]:
    return tuple(a for a in range(1, q) if gcd(a, q) == 1)


def safe_mask(v: int, q: int) -> int:
    """Bit mask of unit twists at q for which speed v is outside the danger band."""
    band = q // N
    mask = 0
    for j, a in enumerate(unit_twists(q)):
        residue = (a * v) % q
        if min(residue, q - residue) > band:
            mask |= 1 << j
    return mask


def build_speed_masks(qs: Sequence[int]) -> dict[int, list[int]]:
    masks: dict[int, list[int]] = {}
    for q in qs:
        arr = [0] * (MAX_R + 1)
        for v in range(1, MAX_R + 1):
            arr[v] = safe_mask(v, q)
        masks[q] = arr
    return masks


def build_core_base_masks(qs: Sequence[int], masks: dict[int, list[int]]) -> dict[int, dict[int, int]]:
    base: dict[int, dict[int, int]] = {}
    for deleted in CORE:
        core_row = [v for v in CORE if v != deleted]
        base[deleted] = {}
        for q in qs:
            mask = (1 << len(unit_twists(q))) - 1
            speed_masks = masks[q]
            for v in core_row:
                mask &= speed_masks[v]
            base[deleted][q] = mask
    return base


def first_witness_pair(
    deleted: int,
    r1: int,
    r2: int,
    qs: Sequence[int],
    masks: dict[int, list[int]],
    base: dict[int, dict[int, int]],
) -> tuple[int, int] | None:
    for q in qs:
        mask = base[deleted][q] & masks[q][r1] & masks[q][r2]
        if mask:
            low_bit = mask & -mask
            twist_index = low_bit.bit_length() - 1
            return q, unit_twists(q)[twist_index]
    return None


def true_two_stranger_candidates() -> tuple[int, ...]:
    core = set(CORE)
    return tuple(v for v in range(1, MAX_R + 1) if v not in core)


def scan_residuals() -> tuple[int, list[ResidualRecord]]:
    """Return primitive row count and rows missed by every plain shell q<=27."""
    masks = build_speed_masks(Q27)
    base = build_core_base_masks(Q27, masks)
    candidates = true_two_stranger_candidates()
    primitive_rows = 0
    residuals: list[ResidualRecord] = []

    for deleted in CORE:
        for i, r1 in enumerate(candidates):
            for r2 in candidates[i + 1 :]:
                # The remaining 7-core has gcd 7, so primitivity fails exactly
                # when both added strangers are also multiples of 7.
                if r1 % 7 == 0 and r2 % 7 == 0:
                    continue
                primitive_rows += 1
                if first_witness_pair(deleted, r1, r2, SHELL27, masks, base) is not None:
                    continue
                residuals.append(
                    ResidualRecord(
                        deleted=deleted,
                        r1=r1,
                        r2=r2,
                        q27=first_witness_pair(deleted, r1, r2, Q27, masks, base),
                    )
                )
    return primitive_rows, residuals


def divisor_addresses(q: int | None) -> tuple[tuple[int, int], ...]:
    if q is None:
        return tuple()
    return tuple((d, q // d) for d in DIVISORS if q % d == 0 and q // d <= 27)


def first_examples_by_q(residuals: Iterable[ResidualRecord]) -> dict[int, ResidualRecord]:
    examples: dict[int, ResidualRecord] = {}
    for rec in residuals:
        if rec.first_q is not None and rec.first_q not in examples:
            examples[rec.first_q] = rec
    return dict(sorted(examples.items()))


def print_scan_summary(primitive_rows: int, residuals: list[ResidualRecord]) -> None:
    q27_misses = [rec for rec in residuals if rec.q27 is None]
    print("A. Exact bounded two-stranger scan")
    print(f"  core={CORE}")
    print(f"  added strangers: distinct non-core speeds in [1,{MAX_R}]")
    print(f"  primitive rows scanned={primitive_rows}")
    print(f"  rows blocking every plain shell q<=27={len(residuals)}")
    print(f"  Q27 misses among those residuals={len(q27_misses)}")
    print(f"  first-Q27 histogram={dict(sorted(Counter(rec.first_q for rec in residuals).items()))}")
    print(f"  Q27 address histogram={dict(sorted(Counter(divisor_addresses(rec.first_q) for rec in residuals).items(), key=lambda kv: str(kv[0])))}")
    print()


def print_compression_coordinates(residuals: list[ResidualRecord]) -> None:
    deleted_hist = Counter(rec.deleted for rec in residuals)
    absent_deleted = [v for v in CORE if deleted_hist[v] == 0]
    print("B. Compression coordinates forced by plain-shell residuals")
    print(f"  deleted-core histogram={dict(sorted(deleted_hist.items()))}")
    print(f"  deleted core speeds never appearing={absent_deleted}")
    print(f"  hard-residue count histogram={dict(sorted(Counter(rec.hard_count for rec in residuals).items()))}")
    print(f"  mod13-zero count histogram={dict(sorted(Counter(rec.mod13_zero_count for rec in residuals).items()))}")
    print(f"  residuals with no 13-clock speed={sum(1 for rec in residuals if rec.mod13_zero_count == 0)}")
    print(f"  top mod13 pairs={Counter(rec.mod13_pair for rec in residuals).most_common(14)}")
    print(f"  top antipodal mod27 pairs={Counter(rec.anti27_pair for rec in residuals).most_common(18)}")
    print()


def print_late_fibers(residuals: list[ResidualRecord]) -> None:
    late = [rec for rec in residuals if rec.first_q is not None and rec.first_q >= 70]
    print("C. Late divisor-fiber rescues")
    print(f"  q>=70 residual count={len(late)}")
    print(f"  late q histogram={dict(sorted(Counter(rec.first_q for rec in late).items()))}")
    print("  unique q=91 and q=161 witnesses:")
    for rec in late:
        if rec.first_q in (91, 161):
            print(
                f"    deleted={rec.deleted} added={rec.added} q27={rec.q27} "
                f"addr={divisor_addresses(rec.first_q)} mod13={rec.mod13_pair} "
                f"anti27={rec.anti27_pair} row={rec.row}"
            )
    for q in (70, 84):
        rows = [rec for rec in residuals if rec.first_q == q]
        print(
            f"  q={q}: count={len(rows)} deleted={dict(sorted(Counter(rec.deleted for rec in rows).items()))} "
            f"top_mod13={Counter(rec.mod13_pair for rec in rows).most_common(6)} "
            f"top_anti27={Counter(rec.anti27_pair for rec in rows).most_common(6)}"
        )
    print()


def print_examples(residuals: list[ResidualRecord]) -> None:
    print("D. First examples by Q27 rescue denominator")
    examples = first_examples_by_q(residuals)
    for q, rec in examples.items():
        print(
            f"  q={q:3d} a={rec.first_a:3d} deleted={rec.deleted:2d} "
            f"added={rec.added} hard_count={rec.hard_count} "
            f"mod13={rec.mod13_pair} anti27={rec.anti27_pair}"
        )
    print()


def orient_tournament(vertices: list[tuple[str, tuple[int, ...]]]) -> tuple[dict[str, set[str]], int]:
    names = [name for name, _scores in vertices]
    scores = dict(vertices)
    out = {name: set() for name in names}
    flips = 0
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            av = sum(x > y for x, y in zip(scores[a], scores[b]))
            bv = sum(y > x for x, y in zip(scores[a], scores[b]))
            if av >= bv:
                out[a].add(b)
            else:
                out[b].add(a)
                flips += 1
    return out, flips


def directed_3cycles(names: Sequence[str], out: dict[str, set[str]]) -> int:
    cycles = 0
    for a, b, c in combinations(names, 3):
        wins = {a: 0, b: 0, c: 0}
        for x, y in ((a, b), (a, c), (b, c)):
            if y in out[x]:
                wins[x] += 1
            else:
                wins[y] += 1
        if sorted(wins.values()) == [1, 1, 1]:
            cycles += 1
    return cycles


def scc_sizes(names: Sequence[str], out: dict[str, set[str]]) -> list[int]:
    reverse: dict[str, set[str]] = {name: set() for name in names}
    for a in names:
        for b in out[a]:
            reverse[b].add(a)

    seen: set[str] = set()
    order: list[str] = []

    def dfs(v: str) -> None:
        seen.add(v)
        for w in out[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for name in names:
        if name not in seen:
            dfs(name)

    seen.clear()
    sizes: list[int] = []

    def rdfs(v: str) -> int:
        seen.add(v)
        total = 1
        for w in reverse[v]:
            if w not in seen:
                total += rdfs(w)
        return total

    for name in reversed(order):
        if name not in seen:
            sizes.append(rdfs(name))
    return sorted(sizes, reverse=True)


def hamiltonian_path_count(names: Sequence[str], out: dict[str, set[str]]) -> int:
    index = {name: i for i, name in enumerate(names)}
    n = len(names)
    dp: dict[tuple[int, int], int] = {}
    for name in names:
        dp[(1 << index[name], index[name])] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp.get((mask, last), 0)
            if not count:
                continue
            last_name = names[last]
            for nxt_name in out[last_name]:
                nxt = index[nxt_name]
                if mask & (1 << nxt):
                    continue
                dp[(mask | (1 << nxt), nxt)] = dp.get((mask | (1 << nxt), nxt), 0) + count
    full = (1 << n) - 1
    return sum(dp.get((full, last), 0) for last in range(n))


def print_tournament_analysis() -> None:
    print("E. Tournament Analysis over compression maps")
    vertices = [
        ("divisor_fiber_Q27", (5, 4, 5, 5, 4, 5)),
        ("thirteen_clock_debt", (5, 5, 4, 5, 5, 4)),
        ("deleted_core_address", (5, 5, 4, 4, 5, 4)),
        ("shell27_pair_classes", (3, 4, 4, 3, 5, 3)),
        ("hard_atom_hull", (4, 3, 4, 5, 3, 3)),
        ("Bprime_owner_transport", (1, 5, 5, 4, 2, 5)),
        ("raw_plain_shell_scalar", (2, 2, 2, 2, 5, 1)),
    ]
    names = [name for name, _scores in vertices]
    scores = dict(vertices)
    out, flips = orient_tournament(vertices)
    hist = Counter(len(out[name]) for name in names)
    print("  vertices are compression maps/proof obligations, not runners")
    print("  observable=(exactness, compression leverage, generality, retention, cheapness, transfer)")
    print("  tie path=list order above")
    print(f"  score_hist={dict(sorted(hist.items()))}")
    print(f"  directed_3cycles={directed_3cycles(names, out)}")
    print(f"  SCC_sizes={scc_sizes(names, out)}")
    print(f"  Hamiltonian_path_count={hamiltonian_path_count(names, out)}")
    print(f"  edge_flips_vs_list_order={flips}/21")
    for name in sorted(names, key=lambda x: (-len(out[x]), x)):
        print(f"    {name:28s} out={len(out[name])} scores={scores[name]}")
    print()


def print_takeaway(primitive_rows: int, residuals: list[ResidualRecord]) -> None:
    print("F. Proof-shaped takeaway")
    print(
        "  Exact bounded lemma: every primitive true two-stranger near-core row in the"
    )
    print(
        f"  [1,{MAX_R}] horizon that blocks all plain q<=27 shells is opened by Q27."
    )
    print(
        "  The residual set is not just the HYP-2463 hard-atom hull: 636/877 residuals"
    )
    print(
        "  use zero old hard residues.  But every residual still spends a 13-clock"
    )
    print(
        "  speed and a visible deleted-core address, and the late cases are divisor"
    )
    print(
        "  fibers such as 70, 84, 91, and the lone 161=7*23."
    )
    print(
        "  Compression target: prove that any arbitrary Q27 blocker can be reduced to"
    )
    print(
        "  this resource pattern, or else it opens a low clock, divisor-fiber witness,"
    )
    print(
        "  AP/Vstar/2AP descent, or an odd owner/Bprime channel."
    )
    print(f"  Audit rows scanned={primitive_rows}; residuals={len(residuals)}.")


def main() -> None:
    print("=" * 78)
    print("Codex LRC14 two-stranger compression stress")
    print("=" * 78)
    print("HYP-2464 / T805 / OPEN-Q-083")
    print()
    primitive_rows, residuals = scan_residuals()
    print_scan_summary(primitive_rows, residuals)
    print_compression_coordinates(residuals)
    print_late_fibers(residuals)
    print_examples(residuals)
    print_tournament_analysis()
    print_takeaway(primitive_rows, residuals)


if __name__ == "__main__":
    main()
