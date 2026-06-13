#!/usr/bin/env python3
"""
lrc_even_clean_lane_fibre_collapse_s577.py

codex-2026-06-03 S577

Merge the live even-ladder ideas from several agents:
  * S576/HYP-2092 (Codex): HYP-2091 gives a four-lane proof router.
  * S568/HYP-2093 (Opus): the even lane target is measure-zero => floor-tight.
  * S576o/HYP-2094 (Oracle): on the even lane, worry collapses from the
    converse-merged round seam to 2^((n-2)/2) self-converse round classes.
  * S574/HYP-2090 and THM-397: the remaining finite target must keep D/U/N
    private-pivot labels and pair-sum endpoint-owner labels.

This script computes the exact class-count funnel and clock-burden labels for
even LRC n.  It is a bookkeeping/proof-planning computation, not a proof.

Tournament Analysis / assumption challenge:
  Vertices are even-LRC rungs, not runners.
  Observable is the finite residual burden after successive quotient collapses:
      all tournaments -> round -> converse-merged -> self-converse fixed,
      plus nonunit shell count and D/U/N obligation count.
  Switch orients toward larger residual burden.  Tie path is increasing n.
  This preserves the size and type of the remaining proof obligation; it
  destroys individual speed realizations, which must return through per-class
  owner-labelled certificates.
"""

from __future__ import annotations

from collections import Counter, deque
from dataclasses import dataclass
from itertools import combinations
from math import factorial, gcd, log10


@dataclass(frozen=True)
class EvenRung:
    n: int
    m: int
    c: int
    factors: str
    all_tournaments: int
    round_classes: int
    converse_merged: int
    self_converse_round: int
    unit_shells: int
    nonunit_shells: int
    d_obligations: int
    u_obligations: int
    n_obligations: int
    bounded_tight_count: int | None
    bounded_nontransversal_tight: int | None

    @property
    def total_dun(self) -> int:
        return self.d_obligations + self.u_obligations + self.n_obligations

    @property
    def burden(self) -> tuple[int, int, int, int]:
        return (
            self.self_converse_round,
            self.nonunit_shells,
            self.total_dun,
            self.converse_merged,
        )


S576O_BOUNDED_TIGHT = {
    4: (1, 0),
    6: (2, 0),
    8: (3, 2),
    10: (1, 0),
    12: (1, 0),
    14: (1, 0),
}


def divisors(n: int) -> list[int]:
    return [d for d in range(1, n + 1) if n % d == 0]


def totient(n: int) -> int:
    return sum(1 for a in range(1, n + 1) if gcd(a, n) == 1)


def factorization(n: int) -> str:
    out: list[str] = []
    d = 2
    x = n
    while d * d <= x:
        if x % d:
            d += 1
            continue
        e = 0
        while x % d == 0:
            x //= d
            e += 1
        out.append(str(d) if e == 1 else f"{d}^{e}")
        d += 1
    if x > 1:
        out.append(str(x))
    return "*".join(out)


def partitions(n: int, max_part: int | None = None) -> list[tuple[int, ...]]:
    if n == 0:
        return [()]
    if max_part is None or max_part > n:
        max_part = n
    out: list[tuple[int, ...]] = []
    for first in range(min(max_part, n), 0, -1):
        for rest in partitions(n - first, first):
            out.append((first,) + rest)
    return out


def a000568_tournaments(n: int) -> int:
    total = 0
    for parts in partitions(n):
        if any(p % 2 == 0 for p in parts):
            continue
        counts = Counter(parts)
        autom = 1
        for length, multiplicity in counts.items():
            autom *= (length ** multiplicity) * factorial(multiplicity)
        exponent = sum((p - 1) // 2 for p in parts)
        for i, a in enumerate(parts):
            for b in parts[i + 1 :]:
                exponent += gcd(a, b)
        total += factorial(n) // autom * 2**exponent
    return total // factorial(n)


def round_classes(m: int) -> int:
    return sum(totient(d) * 2 ** (m // d) for d in divisors(m) if d % 2 == 1) // (2 * m)


def shell_counts(c: int) -> tuple[int, int]:
    k = (c - 1) // 2
    unit = sum(1 for a in range(1, k + 1) if gcd(a, c) == 1)
    return unit, k - unit


def rung(n: int) -> EvenRung:
    m = n - 1
    c = 2 * n - 1
    round_count = round_classes(m)
    sc_round = 2 ** ((n - 2) // 2)
    unit, nonunit = shell_counts(c)
    tight = S576O_BOUNDED_TIGHT.get(n)
    return EvenRung(
        n=n,
        m=m,
        c=c,
        factors=factorization(c),
        all_tournaments=a000568_tournaments(m),
        round_classes=round_count,
        converse_merged=(round_count + sc_round) // 2,
        self_converse_round=sc_round,
        unit_shells=unit,
        nonunit_shells=nonunit,
        d_obligations=n - 2,
        u_obligations=totient(c) // 2,
        n_obligations=n - 1,
        bounded_tight_count=None if tight is None else tight[0],
        bounded_nontransversal_tight=None if tight is None else tight[1],
    )


def fmt_count(x: int) -> str:
    if x < 1_000_000:
        return str(x)
    return f"10^{log10(x):.2f}"


def fmt_ratio(num: int, den: int) -> str:
    if num == 0:
        return "0"
    return f"10^{log10(num / den):.2f}"


def route(row: EvenRung) -> str:
    if row.nonunit_shells:
        return "fixed-SC target + nonunit gcd descent"
    return "fixed-round certificates + D/U/N exchange"


def tournament_fingerprint(rows: list[EvenRung]) -> dict[str, object]:
    n = len(rows)
    adj = [[False] * n for _ in range(n)]
    for i, a in enumerate(rows):
        for j, b in enumerate(rows):
            if i == j:
                continue
            adj[i][j] = (a.burden, a.n) > (b.burden, b.n)
    scores = [sum(row) for row in adj]
    c3 = 0
    for i, j, k in combinations(range(n), 3):
        if (adj[i][j] and adj[j][k] and adj[k][i]) or (
            adj[i][k] and adj[k][j] and adj[j][i]
        ):
            c3 += 1

    def reach(start: int) -> set[int]:
        seen = {start}
        q = deque([start])
        while q:
            v = q.popleft()
            for w, edge in enumerate(adj[v]):
                if edge and w not in seen:
                    seen.add(w)
                    q.append(w)
        return seen

    remaining = set(range(n))
    sccs: list[int] = []
    while remaining:
        start = next(iter(remaining))
        forward = reach(start)
        comp = {v for v in remaining if v in forward and start in reach(v)}
        sccs.append(len(comp))
        remaining -= comp
    return {
        "score_hist": dict(sorted(Counter(scores).items())),
        "directed_3_cycles": c3,
        "sccs": sorted(sccs, reverse=True),
        "hamiltonian_path": [f"n={r.n}" for r in sorted(rows, key=lambda x: (x.burden, x.n), reverse=True)],
    }


def main() -> None:
    rows = [rung(n) for n in range(4, 20, 2)]

    print("S577 even clean-lane fibre collapse")
    print("=" * 78)
    print("Agent merge: S576 route + HYP-2093 floor-tight target + HYP-2094")
    print("self-converse worry collapse + HYP-2090/THM-397 owner labels.")
    print()

    print("n  C       all_T        round merged fixedSC  fixed/all fixed/merged U/nonU DUN  bounded tight route")
    for row in rows:
        tight = "-" if row.bounded_tight_count is None else f"{row.bounded_tight_count}/{row.bounded_nontransversal_tight}"
        print(
            f"{row.n:2d} {row.c:<5d} "
            f"{fmt_count(row.all_tournaments):>10s} "
            f"{row.round_classes:6d} {row.converse_merged:6d} {row.self_converse_round:7d} "
            f"{fmt_ratio(row.self_converse_round, row.all_tournaments):>9s} "
            f"{row.self_converse_round / row.converse_merged:11.4f} "
            f"{row.unit_shells:2d}/{row.nonunit_shells:<2d} {row.total_dun:3d} "
            f"{tight:>6s}  {route(row)}"
        )

    n14 = next(r for r in rows if r.n == 14)
    print("\nN=14 funnel")
    print(f"  all runner tournaments A000568(13): {n14.all_tournaments}")
    print(f"  open round body A000016(13):       {n14.round_classes}")
    print(f"  converse-merged round seam:        {n14.converse_merged}")
    print(f"  self-converse fixed worry nodes:   {n14.self_converse_round}")
    print(f"  D/U/N labels per speed set:        D={n14.d_obligations}, U={n14.u_obligations}, N={n14.n_obligations}")
    print(f"  pair-sum endpoint cells:           C(13,2)=78")
    print(f"  nonunit summand shells mod 27:      {n14.nonunit_shells}")
    print(f"  fixed/all compression:             {n14.self_converse_round}/{n14.all_tournaments} = {fmt_ratio(n14.self_converse_round, n14.all_tournaments)}")
    print(f"  fixed/merged compression:          {n14.self_converse_round}/{n14.converse_merged} = {n14.self_converse_round / n14.converse_merged:.4f}")

    print("\nMerged lemma queue")
    print("  E1 containment: counterexample optimal class is self-converse round.")
    print("  E2 certificate table: each fixed n=14 class gets a pinch or n-clock witness.")
    print("  E3 realization independence: the class certificate survives all speed realizations.")
    print("  E4 nonunit descent: gcd-3/gcd-9 shell defects lift to a second clock or endpoint owner.")
    print("  E5 owner compatibility: D/U/N private pivots and THM-397 endpoint owners agree.")

    print("\nWhat changed from S576")
    print("  S576 used the 190 converse-merged round seam as the clean-lane target.")
    print("  HYP-2094 sharpens the actual even-lane worry target to the 64 fixed nodes.")
    print("  Therefore the next n=14 computation should not label all 190 equally.")
    print("  It should label the 64 fixed nodes first, then use the 126 paired nodes as")
    print("  controls for why generic round classes are already lonely.")

    print("\nTournament Analysis")
    print("  vertices: even LRC rungs")
    print("  observable: fixed SC nodes, nonunit shell count, D/U/N total, merged seam")
    print("  switch: orient toward larger finite residual burden")
    print(f"  fingerprints: {tournament_fingerprint(rows)}")


if __name__ == "__main__":
    main()
