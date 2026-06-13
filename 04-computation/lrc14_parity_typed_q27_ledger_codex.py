#!/usr/bin/env python3
"""Parity-typed Q27 ledger for LRC(14).

This extends HYP-2459's parity projector protocol into the live LRC14 Q27
program.  The computation focuses on the shell-27/13-clock one-stranger
obstruction from HYP-2444, then asks whether those "hard" resource atoms can be
stacked by replacing more of the 7-core.

The key test is finite and deliberately modest:

  core = 7*{1,...,12}
  hard residues = the eight one-stranger rows that block all plain q<=27 shells

For k hard residues, delete k-1 core speeds so the row still has 13 moving
runners.  The stored run scans the whole hard replacement hull `k<=8`: all
77,520 primitive ways of stacking these hard one-stranger resource atoms inside
the 7-core frame.  If all rows in this hull have a Q27 witness, then the
hard-resource stacking attempt is forced to open other clocks.  The proof lesson
is not that this proves LRC14; it says which resource-independence lemma the
full proof should try to generalize.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import lru_cache, reduce
from itertools import combinations
from math import gcd
from pathlib import Path
import sys
from typing import Iterable, Sequence


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[0]
sys.path.append(str(HERE))

from lrc14_marked_ladder_support_codex import block_ledger, denominators
from lrc14_pisano_band_ladder_codex import (
    CORE,
    N,
    antipodal_class,
    bprime_any,
    q_lattice,
)


HARD_RS = (260, 351, 442, 611, 702, 793, 962, 1053)
Q27 = tuple(q_lattice(27))
SHELL27 = tuple(range(2, 28))
DIVISORS = (1, 2, 7, 14)


def row_gcd(row: Sequence[int]) -> int:
    return reduce(gcd, row, 0)


def row_for(add: Iterable[int], delete: Iterable[int]) -> tuple[int, ...]:
    deleted = set(delete)
    added = set(add)
    return tuple(sorted([v for v in CORE if v not in deleted] + list(added)))


def divisor_addresses(q: int | None) -> tuple[tuple[int, int], ...]:
    if q is None:
        return tuple()
    return tuple((d, q // d) for d in DIVISORS if q % d == 0 and q // d <= 27)


def residue_packet(r: int) -> str:
    c = antipodal_class(r, 27)
    if r % 27 == 0:
        shell = "zero27"
    elif c == 10:
        shell = "missing_pm10"
    else:
        shell = f"class{c}"
    clock = "clock13" if r % 13 == 0 else "free13"
    return f"{shell}+{clock}"


def bprime_target_type(row: Sequence[int], target: int | None) -> str:
    if target is None:
        return "none"
    if target in HARD_RS:
        return "hard_stranger"
    if target in CORE:
        return "core"
    if target % N == 0:
        return "multiple14"
    return "other"


@lru_cache(maxsize=None)
def unit_twists(q: int) -> tuple[int, ...]:
    return tuple(a for a in range(1, q) if gcd(a, q) == 1)


@lru_cache(maxsize=None)
def safe_twist_mask(q: int, v: int) -> int:
    """Bit mask of unit twists at q for which speed v is safely outside band."""
    band = q // N
    mask = 0
    for j, a in enumerate(unit_twists(q)):
        residue = (a * v) % q
        if min(residue, q - residue) > band:
            mask |= 1 << j
    return mask


def first_witness_masked(row: Sequence[int], qs: Sequence[int]) -> tuple[int, int] | None:
    """Fast first witness search using precomputed per-speed safe-twist masks."""
    for q in qs:
        mask = (1 << len(unit_twists(q))) - 1
        for v in row:
            mask &= safe_twist_mask(q, v)
            if mask == 0:
                break
        if mask:
            low_bit = mask & -mask
            twist_index = low_bit.bit_length() - 1
            return q, unit_twists(q)[twist_index]
    return None


@dataclass(frozen=True)
class ReplacementRecord:
    hard_count: int
    deleted: tuple[int, ...]
    added: tuple[int, ...]
    row: tuple[int, ...]
    plain27: tuple[int, int] | None
    q27: tuple[int, int] | None

    @property
    def first_q(self) -> int | None:
        return None if self.q27 is None else self.q27[0]

    @property
    def first_plain_q(self) -> int | None:
        return None if self.plain27 is None else self.plain27[0]


def replacement_records(max_hard_count: int = 4) -> list[ReplacementRecord]:
    records: list[ReplacementRecord] = []
    for hard_count in range(1, max_hard_count + 1):
        delete_count = hard_count - 1
        for deleted in combinations(CORE, delete_count):
            for added in combinations(HARD_RS, hard_count):
                row = row_for(added, deleted)
                if len(row) != 13 or len(set(row)) != 13 or row_gcd(row) != 1:
                    continue
                plain27 = first_witness_masked(row, SHELL27)
                q27 = first_witness_masked(row, Q27)
                records.append(
                    ReplacementRecord(
                        hard_count=hard_count,
                        deleted=tuple(deleted),
                        added=tuple(added),
                        row=row,
                        plain27=plain27,
                        q27=q27,
                    )
                )
    return records


def marked_pressure(row: Sequence[int]) -> dict[str, object]:
    ledgers = [block_ledger(tuple(row), q) for q in denominators(27)]
    blocked = [ledger for ledger in ledgers if ledger.blocked]
    cover_load: Counter[int] = Counter()
    for ledger in blocked:
        for v in ledger.cover_speeds:
            cover_load[v] += 1
    return {
        "blocked_q": len(blocked),
        "universal_q": sum(1 for ledger in blocked if ledger.has_universal),
        "cover_q": sum(1 for ledger in blocked if ledger.blocked and not ledger.has_universal),
        "max_min_blockers": max((ledger.min_blockers for ledger in blocked), default=0),
        "top_cover_load": tuple(cover_load.most_common(6)),
    }


def print_one_stranger_packets() -> None:
    print("A. One-stranger shell-27/13-clock packets")
    print("  hard residues are exactly the HYP-2444 plain-shell q<=27 blockers.")
    for r in HARD_RS:
        row = row_for([r], [])
        plain = first_witness_masked(row, SHELL27)
        q27 = first_witness_masked(row, Q27)
        b_ok, b_target, b_margin = bprime_any(list(row), preferred=r)
        print(
            f"  r={r:4d} packet={residue_packet(r):22s} "
            f"plain27={plain} q27={q27} addr={divisor_addresses(None if q27 is None else q27[0])} "
            f"Bprime={b_ok}:{bprime_target_type(row,b_target)}:{b_target} margin={b_margin}"
        )
    print()


def summarize_replacement_hull(records: Sequence[ReplacementRecord]) -> None:
    print("B. Hard-resource replacement hull")
    print("  Replace k-1 core runners by k hard shell-27/13-clock residues.")
    by_k: dict[int, list[ReplacementRecord]] = defaultdict(list)
    for rec in records:
        by_k[rec.hard_count].append(rec)
    for k in sorted(by_k):
        rows = by_k[k]
        q_hist = Counter(rec.first_q for rec in rows)
        plain_miss = sum(1 for rec in rows if rec.plain27 is None)
        q27_miss = sum(1 for rec in rows if rec.q27 is None)
        low_hits = sum(1 for rec in rows if rec.first_q is not None and rec.first_q <= 13)
        print(
            f"  k={k}: rows={len(rows):5d} Q27_miss={q27_miss:3d} "
            f"plain27_miss={plain_miss:5d} q<=13_hits={low_hits:5d}"
        )
        print(f"       first-Q27 histogram={dict(sorted(q_hist.items(), key=lambda kv: (9999 if kv[0] is None else kv[0])))}")
    print()


def print_plain_shell_residuals(records: Sequence[ReplacementRecord]) -> None:
    print("C. Plain shell residuals after hard-resource replacement")
    residuals = [rec for rec in records if rec.plain27 is None]
    print(
        "  These are the only rows in the complete replacement hull still missed by plain q<=27."
    )
    for rec in residuals:
        print(
            f"  k={rec.hard_count} q27={rec.q27} "
            f"deleted={rec.deleted} added={rec.added}"
        )
    print()


def print_hardest_rows(records: Sequence[ReplacementRecord], limit: int = 12) -> None:
    print("D. Rows with latest first Q27 witness")
    ranked = sorted(
        records,
        key=lambda rec: (
            -(-1 if rec.first_q is None else rec.first_q),
            rec.plain27 is not None,
            -rec.hard_count,
            rec.deleted,
            rec.added,
        ),
    )
    for rec in ranked[:limit]:
        pressure = marked_pressure(rec.row)
        b_ok, b_target, b_margin = bprime_any(list(rec.row))
        b_type = bprime_target_type(rec.row, b_target)
        print(
            f"  k={rec.hard_count} q27={rec.q27} plain27={rec.plain27} "
            f"deleted={rec.deleted} added={rec.added}"
        )
        print(
            f"     addr={divisor_addresses(rec.first_q)} "
            f"Bprime={b_ok}:{b_type}:{b_target} margin={b_margin}"
        )
        print(
            f"     pressure blocked={pressure['blocked_q']} universal={pressure['universal_q']} "
            f"cover={pressure['cover_q']} max_tau={pressure['max_min_blockers']} "
            f"top_load={pressure['top_cover_load']}"
        )
    print()


def print_deleted_core_effect(records: Sequence[ReplacementRecord]) -> None:
    print("E. Deleted-core low-clock effect")
    # Only k>=2 rows delete core runners.  Count which deleted core speeds tend to
    # create a very early witness.
    exposure: Counter[int] = Counter()
    exposure_low: Counter[int] = Counter()
    exposure_plain: Counter[int] = Counter()
    for rec in records:
        for v in rec.deleted:
            exposure[v] += 1
            if rec.first_q is not None and rec.first_q <= 13:
                exposure_low[v] += 1
            if rec.plain27 is not None:
                exposure_plain[v] += 1
    print("  deleted core speed -> low-Q27-hit fraction / plain27-hit fraction")
    for v in CORE:
        if exposure[v] == 0:
            continue
        low = f"{exposure_low[v]}/{exposure[v]}"
        plain = f"{exposure_plain[v]}/{exposure[v]}"
        print(f"    {v:3d}: low={low:>11s} plain={plain:>11s}")
    print()


def print_tournament_analysis() -> None:
    print("F. Tournament Analysis over proof obligations")
    vertices = [
        ("typed_Q27_ledger", (5, 5, 5, 5, 4, 5)),
        ("hard_resource_replacement_hull", (5, 5, 4, 5, 5, 4)),
        ("deleted_core_low_clock", (4, 5, 5, 4, 5, 4)),
        ("shell27_x_13_packet", (4, 4, 5, 5, 3, 5)),
        ("Bprime_target_transport", (4, 4, 3, 4, 4, 4)),
        ("resonance_bonus_analogy", (3, 3, 4, 3, 3, 5)),
        ("raw_plain_shell_q27", (2, 2, 3, 2, 5, 2)),
    ]
    names = [v[0] for v in vertices]
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
    hist = Counter(len(out[name]) for name in names)

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

    print("  vertices are proof obligations, not runners or denominators")
    print("  observable=(parity typing, exactness, LRC leverage, support retention, computation, transfer)")
    print(f"  score_hist={dict(sorted(hist.items()))}")
    print(f"  directed_3cycles={cycles}; edge_flips_vs_list_order={flips}/21")
    for name in sorted(names, key=lambda x: (-len(out[x]), x)):
        print(f"    {name:32s} out={len(out[name])} scores={scores[name]}")
    print()


def main() -> None:
    print("=" * 78)
    print("Codex LRC14 parity-typed Q27 ledger")
    print("=" * 78)
    print("HYP-2459 -> HYP-2444/HYP-2443 extension")
    print(f"n={N}; hard residues={HARD_RS}; Q27 size={len(Q27)}")
    print()

    print_one_stranger_packets()
    max_hard_count = len(HARD_RS)
    records = replacement_records(max_hard_count=max_hard_count)
    summarize_replacement_hull(records)
    print_plain_shell_residuals(records)
    print_hardest_rows(records)
    print_deleted_core_effect(records)
    print_tournament_analysis()

    q27_misses = [rec for rec in records if rec.q27 is None]
    print("G. Proof-shaped takeaway")
    print(f"  hard-resource replacement rows checked with k<= {max_hard_count}: {len(records)}")
    print(f"  Q27 misses inside replacement hull: {len(q27_misses)}")
    if not q27_misses:
        print(
            "  Exact finite lemma for this complete hard-replacement hull: stacking the one-stranger shell-27/13-clock"
        )
        print(
            "  obstruction never escapes Q27.  Replacing core runners opens low clocks or"
        )
        print(
            "  divisor-fiber addresses before the shell-27 resource can accumulate."
        )
    print(
        "  LRC proof target: generalize this from the hard replacement hull to arbitrary"
    )
    print(
        "  primitive multi-stranger rows by proving resource independence: a runner can spend"
    )
    print(
        "  shell-27 class, 13-clock, low-divisor coverage, or Bprime width, but not all of them"
    )
    print("  indefinitely across the Q27 ladder.")


if __name__ == "__main__":
    main()
