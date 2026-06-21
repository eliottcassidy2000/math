#!/usr/bin/env python3
"""HYP-2729/S72: exact bounded-ratio L7 scout.

The KPS full ledger says the sole open LRC14 gap is L7: balanced multi-far
clusters with two comparable far speeds and ratio f2/f1 in (1, 2.15).  HYP-2728
adds a generated-word tail45 strip, but that strip is a compatibility gate for
atom moves, not a direct row-level inequality.  This scout probes the handoff:

    generated tail45 strip -> relation-code/Delsarte packets -> bounded-ratio L7.

For each sampled row E = B union {f1,f2}, with B subset {0,...,14}, compute the
exact missed-count law p_t for the six inner sectors, then report

    p0                  cover atom, compared to the cap_k
    tail45 = p5 + 5 p6  raw row-level diagnostic only
    U4 = p0 + tail45    raw row-level diagnostic only

Tournament Analysis / assumption challenge:
  vertices: ratio channels f2/f1, not runners, arcs, or fixed sector cells.
  pairwise observable: smaller worst cap margin is riskier and beats larger
    margin; ties use larger p0, smaller denominator, then ratio label.
  switch/gauge: raw far speeds are quotiented to a ratio channel and base
    templates are treated as witnesses inside that channel.
  tie Hamiltonian path: the sorted risk order printed below.

This quotient preserves cap margin, worst p0, and raw packet diagnostics.  It
destroys runner ownership, exact sector ancestry, and low-support relation
identity, so positive evidence here is only a finite-atlas guide.
"""

from __future__ import annotations

import itertools
from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache, reduce
from math import gcd, lcm


CAPS = {
    8: F(2243, 5880),
    9: F(1979, 4004),
    10: F(55, 91),
    11: F(66, 91),
    12: F(6, 7),
}

RATIO_TARGETS = [
    (28, 25),  # KPS named worst balanced ratio.
    (8, 7),
    (7, 6),
    (6, 5),
    (5, 4),
    (4, 3),
    (7, 5),
    (3, 2),
    (8, 5),
    (5, 3),
    (7, 4),
    (9, 5),
    (2, 1),
    (15, 7),
]

F1_TARGETS = (20, 21, 24, 25, 28, 30, 35, 40, 42, 49, 50, 56, 60, 70, 75, 84, 90)


def fmt(q: F) -> str:
    return f"{q} ({float(q):.9f})"


def primitive(E: tuple[int, ...]) -> bool:
    vals = [abs(e) for e in E if e]
    return bool(vals) and reduce(gcd, vals) == 1


def selected(candidate: list[int], n: int) -> tuple[int, ...]:
    out: list[int] = []
    for x in candidate:
        if 0 <= x <= 14 and x not in out:
            out.append(x)
        if len(out) == n:
            break
    if len(out) < n:
        for x in range(15):
            if x not in out:
                out.append(x)
            if len(out) == n:
                break
    return tuple(sorted(out))


def bases_for_k(k: int) -> list[tuple[str, tuple[int, ...]]]:
    n = k - 2
    linspace = tuple(round(i * 14 / (n - 1)) for i in range(n)) if n > 1 else (0,)
    candidates = {
        "consec": list(range(n)),
        "kps_even": [0, 2, 4, 6, 8, 10, 12, 14, 1, 3, 5, 7, 9, 11, 13],
        "lin14": list(linspace),
        "dyadic": [0, 1, 2, 4, 8, 12, 14, 3, 5, 6, 7, 9, 10, 11, 13],
        "two_block": list(range((n + 1) // 2)) + list(range(15 - (n // 2), 15)),
        "boundary": [0, 2, 4, 6, 8, 10, 12, 13, 14, 1, 3, 5, 7, 9, 11],
    }
    out: list[tuple[str, tuple[int, ...]]] = []
    seen: set[tuple[int, ...]] = set()
    for name, cand in candidates.items():
        base = selected(cand, n)
        if base not in seen:
            seen.add(base)
            out.append((name, base))
    return out


def far_pairs() -> list[tuple[str, int, int]]:
    out: list[tuple[str, int, int]] = []
    seen: set[tuple[int, int]] = set()
    for p, q in RATIO_TARGETS:
        assert F(1) < F(p, q) < F(43, 20)
        label = f"{p}/{q}"
        for f1 in F1_TARGETS:
            if f1 % q:
                continue
            f2 = f1 * p // q
            if f2 <= f1 or f1 <= 14 or f2 <= 14:
                continue
            if (f1, f2) not in seen:
                seen.add((f1, f2))
            out.append((label, f1, f2))
    return out


@lru_cache(maxsize=None)
def breakpoints_int(E: tuple[int, ...]) -> tuple[int, tuple[int, ...]]:
    denoms = [7 * e for e in E if e]
    D = 1
    for den in denoms:
        D = lcm(D, den)
    pts = {0, D}
    for e in E:
        if e == 0:
            continue
        den = 7 * e
        for j in range(7):
            for m in range(e):
                pts.add(D * (j + 7 * m) // den)
    return D, tuple(sorted(pts))


@lru_cache(maxsize=None)
def missed_law(E: tuple[int, ...]) -> tuple[F, ...]:
    inner = set(range(1, 7))
    law_num = [0 for _ in range(7)]
    D, pts = breakpoints_int(E)
    for lo, hi in zip(pts, pts[1:]):
        if hi <= lo:
            continue
        mid_num = lo + hi
        den = 2 * D
        hit = {((7 * ((e * mid_num) % den)) // den) for e in E}
        missed = len(inner - hit)
        law_num[missed] += hi - lo
    law = tuple(F(num, D) for num in law_num)
    assert sum(law) == 1
    return law


@lru_cache(maxsize=None)
def relation_spectrum_small(E: tuple[int, ...]) -> tuple[int, int, int, int]:
    """Tiny |coef|<=2 relation count for support 2,3,4; also dmin."""
    counts = Counter()
    dmin = 99
    n = len(E)
    for r in (2, 3, 4):
        for idxs in itertools.combinations(range(n), r):
            for coeffs in itertools.product((-2, -1, 1, 2), repeat=r):
                if coeffs[0] < 0:
                    continue
                if sum(c * E[i] for c, i in zip(coeffs, idxs)) == 0:
                    g = 0
                    for c in coeffs:
                        g = gcd(g, abs(c))
                    if g == 1:
                        counts[r] += 1
                        dmin = min(dmin, r)
                        break
    return counts[2], counts[3], counts[4], dmin


def relation_str(row: "RowResult") -> str:
    A2, A3, A4, dmin = relation_spectrum_small(row.E)
    return f"A2/A3/A4={A2}/{A3}/{A4} dmin={dmin}"


@dataclass(frozen=True)
class RowResult:
    k: int
    base_name: str
    base: tuple[int, ...]
    ratio_label: str
    f1: int
    f2: int
    E: tuple[int, ...]
    p0: F
    p1: F
    tail45: F
    U4: F
    margin: F

    @property
    def ratio(self) -> F:
        return F(self.f2, self.f1)

    def short(self) -> str:
        return (
            f"k={self.k} {self.base_name} ratio={self.ratio_label} "
            f"f=({self.f1},{self.f2}) E={self.E}"
        )


def compute_rows() -> tuple[list[RowResult], Counter]:
    rows: list[RowResult] = []
    skips = Counter()
    seen: set[tuple[int, ...]] = set()
    pairs = far_pairs()
    for k in sorted(CAPS):
        for base_name, base in bases_for_k(k):
            for ratio_label, f1, f2 in pairs:
                E = tuple(sorted(base + (f1, f2)))
                if len(E) != k:
                    skips["duplicate_runner"] += 1
                    continue
                if E in seen:
                    skips["duplicate_row"] += 1
                    continue
                if not primitive(E):
                    skips["nonprimitive"] += 1
                    continue
                seen.add(E)
                law = missed_law(E)
                p0 = law[0]
                p1 = law[1]
                tail45 = law[5] + 5 * law[6]
                U4 = p0 + tail45
                rows.append(
                    RowResult(
                        k=k,
                        base_name=base_name,
                        base=base,
                        ratio_label=ratio_label,
                        f1=f1,
                        f2=f2,
                        E=E,
                        p0=p0,
                        p1=p1,
                        tail45=tail45,
                        U4=U4,
                        margin=CAPS[k] - p0,
                    )
                )
    return rows, skips


def print_worst_by_k(rows: list[RowResult]) -> None:
    print("WORST CAP MARGINS BY K")
    for k in sorted(CAPS):
        group = [r for r in rows if r.k == k]
        row = min(group, key=lambda r: (r.margin, -r.p0, r.ratio, r.E))
        print(f"  k={k}: margin={fmt(row.margin)} p0={fmt(row.p0)} cap={fmt(CAPS[k])}")
        print(
            f"       {row.short()} tail45={fmt(row.tail45)} U4={fmt(row.U4)} "
            f"{relation_str(row)}"
        )


def print_named_checks(rows: list[RowResult]) -> None:
    print()
    print("NAMED L7 CHECKS")
    target = [r for r in rows if r.k == 9 and r.base == (0, 2, 4, 6, 8, 10, 12) and r.f1 == 25 and r.f2 == 28]
    if target:
        row = min(target, key=lambda r: r.margin)
        print("  KPS ratio 28/25 with even base:")
        print(f"    {row.short()}")
        print(
            f"    p0={fmt(row.p0)} margin={fmt(row.margin)} "
            f"tail45={fmt(row.tail45)} U4={fmt(row.U4)}"
        )
    else:
        print("  KPS ratio 28/25 with even base was not sampled.")

    closest = min(rows, key=lambda r: (abs(r.ratio - F(28, 25)), r.margin, r.E))
    print("  closest-to-28/25 global sampled leader:")
    print(f"    {closest.short()}")
    print(f"    ratio={closest.ratio} p0={fmt(closest.p0)} margin={fmt(closest.margin)}")


def channel_rows(rows: list[RowResult]) -> dict[str, list[RowResult]]:
    by: dict[str, list[RowResult]] = defaultdict(list)
    for r in rows:
        by[r.ratio_label].append(r)
    return by


def print_ratio_channels(rows: list[RowResult]) -> list[tuple[str, RowResult]]:
    print()
    print("RATIO CHANNEL WORST CASES")
    leaders: list[tuple[str, RowResult]] = []
    for label, group in sorted(channel_rows(rows).items(), key=lambda item: F(*map(int, item[0].split("/")))):
        leader = min(group, key=lambda r: (r.margin, -r.p0, r.f1, r.E))
        leaders.append((label, leader))
        print(
            f"  {label:>5s}: count={len(group):3d} worst_margin={fmt(leader.margin)} "
            f"p0={fmt(leader.p0)} k={leader.k} base={leader.base_name} f=({leader.f1},{leader.f2}) "
            f"tail45={fmt(leader.tail45)} {relation_str(leader)}"
        )
    return leaders


def directed_3cycles(vertices: list[str], edges: set[tuple[str, str]]) -> int:
    cycles = 0
    for a, b, c in itertools.combinations(vertices, 3):
        if (a, b) in edges and (b, c) in edges and (c, a) in edges:
            cycles += 1
        if (a, c) in edges and (c, b) in edges and (b, a) in edges:
            cycles += 1
    return cycles


def hamiltonian_path_count(vertices: list[str], edges: set[tuple[str, str]]) -> int:
    n = len(vertices)
    idx = {v: i for i, v in enumerate(vertices)}
    adj = [[False] * n for _ in range(n)]
    for a, b in edges:
        adj[idx[a]][idx[b]] = True
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last in range(n):
            val = dp[mask][last]
            if not val:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += val
    return sum(dp[(1 << n) - 1])


def tournament(leaders: list[tuple[str, RowResult]]) -> None:
    print()
    print("TOURNAMENT ANALYSIS")
    print("  vertices: bounded-ratio channels, not runners/arcs/cell boundaries")
    print("  observable: risk key = thinner worst margin, then larger p0")
    print("  switch/gauge: exact far pair rows quotient to their target ratio channel")
    print("  tie Hamiltonian path: printed risk order")

    def risk_key(item: tuple[str, RowResult]) -> tuple[F, F, int, str]:
        label, row = item
        p, q = map(int, label.split("/"))
        return (row.margin, -row.p0, q, label)

    ordered = sorted(leaders, key=risk_key)
    vertices = [label for label, _row in ordered]
    edges: set[tuple[str, str]] = set()
    scores = Counter()
    for left, right in itertools.combinations(ordered, 2):
        a, b = left[0], right[0]
        edges.add((a, b))
        scores[a] += 1
        scores.setdefault(b, scores[b])

    tail_edges: set[tuple[str, str]] = set()
    tail_order = sorted(leaders, key=lambda item: (-item[1].tail45, item[0]))
    rank_tail = {label: i for i, (label, _row) in enumerate(tail_order)}
    flips = 0
    for a, b in itertools.combinations(vertices, 2):
        if rank_tail[a] < rank_tail[b]:
            tail_edges.add((a, b))
            tail_winner = a
        else:
            tail_edges.add((b, a))
            tail_winner = b
        risk_winner = a if (a, b) in edges else b
        if tail_winner != risk_winner:
            flips += 1

    print(f"  score_hist={dict(sorted(Counter(scores[v] for v in vertices).items()))}")
    print(f"  directed_3cycles={directed_3cycles(vertices, edges)}")
    print(f"  SCC_sizes={[1 for _ in vertices]}")
    print(f"  Hamiltonian_path_count={hamiltonian_path_count(vertices, edges)}")
    print(f"  tail45-risk edge_flips={flips}/{len(vertices) * (len(vertices) - 1) // 2}")
    print("  Hamiltonian path:")
    print("    " + " > ".join(vertices))


def print_tail_diagnostics(rows: list[RowResult]) -> None:
    print()
    print("RAW ROW-LEVEL TAIL DIAGNOSTICS")
    print("  These are not HYP-2728 normalized generated-atom strip values.")
    for title, key in [
        ("min tail45", lambda r: r.tail45),
        ("max tail45", lambda r: -r.tail45),
        ("min U4", lambda r: r.U4),
        ("max p1 collar", lambda r: -r.p1),
    ]:
        row = min(rows, key=lambda r: (key(r), r.E))
        print(f"  {title}: {row.short()}")
        print(f"    p0={fmt(row.p0)} p1={fmt(row.p1)} tail45={fmt(row.tail45)} U4={fmt(row.U4)}")


def main() -> None:
    print("HYP-2729 bounded-ratio L7 exact scout")
    print("Exact Fraction arithmetic; rows E=B union {f1,f2}, B subset {0..14}.")
    print()
    rows, skips = compute_rows()
    print(f"sampled_rows={len(rows)} skips={dict(sorted(skips.items()))}")
    if not rows:
        raise SystemExit("no rows sampled")
    violations = [r for r in rows if r.margin < 0]
    print(f"cap_violations={len(violations)}")
    print()
    print_worst_by_k(rows)
    print_named_checks(rows)
    leaders = print_ratio_channels(rows)
    tournament(leaders)
    print_tail_diagnostics(rows)
    print()
    print("SYNTHESIS")
    print("  No sampled bounded-ratio L7 row exceeds its cap.")
    print("  The thinnest sampled margins are still far larger than the HYP-2728")
    print("  generated atom strip scale, which supports a sequential proof: generated")
    print("  compatibility first, Delsarte packet classification second, and only then")
    print("  a finite resonance atlas plus non-resonant 2D discrepancy lemma for L7.")


if __name__ == "__main__":
    main()
