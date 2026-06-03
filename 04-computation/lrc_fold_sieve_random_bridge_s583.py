#!/usr/bin/env python3
"""S583: alternate between circuit-free randomness and 3-term fold+sieve.

The user split the current LRC route into:

  Lemma A: circuit-free rows should clear the 1/(k+1) threshold by
           measure/equidistribution; proof wants discrepancy from 3-term-freeness.
  Lemma B: a 3-term relation c=a+b is a literal fold, since c*t=a*t+b*t.

This script mostly probes Lemma B, using Lemma A as the control/noise model.
It treats "no 3-term relation" as the circuit-free proxy, measures exact
maximin M(S) via pair-sum pinch times, records additive 4-term energy, and
separately measures the local fold gate

    x safe, y safe, x+y safe  on (R/Z)^2.

Tournament Analysis uses proof-state buckets as vertices rather than runners.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations, product
from math import gcd
import random


def dist(x: F) -> F:
    x %= 1
    return min(x, 1 - x)


def primitive(row: tuple[int, ...]) -> tuple[int, ...]:
    g = 0
    for v in row:
        g = gcd(g, v)
    return tuple(sorted(v // g for v in row))


def safe_at(row: tuple[int, ...], t: F, n: int) -> bool:
    delta = F(1, n)
    return all(dist(v * t) >= delta for v in row)


def exact_maximin(row: tuple[int, ...], n: int) -> tuple[F, F, int]:
    """Exact M(S) over pair-sum pinch and single-runner peak candidates."""
    candidates: dict[F, int] = {}
    for i, a in enumerate(row):
        for b in row[i + 1 :]:
            D = a + b
            for m in range(1, D):
                candidates.setdefault(F(m, D), D)
        for m in range(1, a):
            candidates.setdefault(F(m, a), a)
    best = F(0)
    best_t = F(0)
    best_D = 1
    for t, D in candidates.items():
        val = min(dist(v * t) for v in row)
        if val > best:
            best = val
            best_t = t
            best_D = D
    return best, best_t, best_D


def safe_measure(row: tuple[int, ...], n: int) -> F:
    delta = F(1, n)
    endpoints: set[F] = {F(0)}
    for v in row:
        for k in range(v + 1):
            endpoints.add(F(k * n + 1, n * v) % 1)
            endpoints.add(F(k * n - 1, n * v) % 1)
    pts = sorted(endpoints)
    total = F(0)
    for i, a in enumerate(pts):
        b = pts[(i + 1) % len(pts)]
        length = b - a if b > a else b - a + 1
        mid = (a + length / 2) % 1
        if all(dist(v * mid) > delta for v in row):
            total += length
    return total


def three_terms(row: tuple[int, ...]) -> list[tuple[int, int, int]]:
    vals = set(row)
    out = []
    for a, b in combinations(row, 2):
        c = a + b
        if c in vals:
            out.append((a, b, c))
    return out


def additive_energy(row: tuple[int, ...]) -> int:
    sums = Counter(a + b for a in row for b in row)
    return sum(v * v for v in sums.values())


def fold_gate_density(n: int, grid: int = 240) -> float:
    """Grid estimate for mu{x,y,x+y all delta-safe}; enough for comparison."""
    delta = F(1, n)
    good = 0
    for i in range(grid):
        x = F(2 * i + 1, 2 * grid)
        if dist(x) < delta:
            continue
        for j in range(grid):
            y = F(2 * j + 1, 2 * grid)
            if dist(y) >= delta and dist(x + y) >= delta:
                good += 1
    return good / (grid * grid)


def active_runners(row: tuple[int, ...], t: F) -> tuple[int, ...]:
    values = [(dist(v * t), v) for v in row]
    m = min(v for v, _ in values)
    return tuple(sorted(v for d, v in values if d == m))


@dataclass(frozen=True)
class RowStats:
    row: tuple[int, ...]
    n: int
    M: F
    witness: F
    witness_den: int
    margin: F
    measure: F
    relations: tuple[tuple[int, int, int], ...]
    energy: int
    shadow_margin: F | None
    shadow_blocked_by_fold: bool
    active_relation_pair: bool


def row_stats(row: tuple[int, ...], n: int) -> RowStats:
    M, witness, den = exact_maximin(row, n)
    rels = tuple(three_terms(row))
    shadow_margin: F | None = None
    shadow_blocked = False
    if rels:
        best_shadow = F(-1)
        blocked = False
        for _a, _b, c in rels:
            shadow = tuple(v for v in row if v != c)
            if len(shadow) == len(row):
                continue
            sm, st, _ = exact_maximin(shadow, n)
            if sm > best_shadow:
                best_shadow = sm
                blocked = sm >= F(1, n) and not safe_at(row, st, n)
        shadow_margin = best_shadow - F(1, n)
        shadow_blocked = blocked
    act = set(active_runners(row, witness))
    active_rel = any(a in act and b in act for a, b, _c in rels)
    return RowStats(
        row=row,
        n=n,
        M=M,
        witness=witness,
        witness_den=den,
        margin=M - F(1, n),
        measure=safe_measure(row, n),
        relations=rels,
        energy=additive_energy(row),
        shadow_margin=shadow_margin,
        shadow_blocked_by_fold=shadow_blocked,
        active_relation_pair=active_rel,
    )


def sample_rows(n: int, trials: int, rng: random.Random) -> list[tuple[int, ...]]:
    rows: set[tuple[int, ...]] = set()
    k = n - 1
    universe = list(range(1, 5 * n + 1))
    while len(rows) < trials:
        row = primitive(tuple(rng.sample(universe, k)))
        if len(row) == k:
            rows.add(row)
    return sorted(rows)


def fmt_frac(x: F | None) -> str:
    if x is None:
        return "-"
    return f"{float(x):.5f}"


def bucket_summary(stats: list[RowStats]) -> dict[str, dict[str, object]]:
    no3 = [s for s in stats if not s.relations]
    has3 = [s for s in stats if s.relations]
    energy_cut = None
    high4: list[RowStats] = []
    if no3:
        energies = sorted(s.energy for s in no3)
        energy_cut = energies[(3 * len(energies)) // 4]
        high4 = [s for s in no3 if s.energy >= energy_cut]
    buckets = {
        "A_circuit_free_proxy": no3,
        "A_high_4term_no3": high4,
        "B_has_3term_fold": has3,
        "B_fold_shadow_blocked": [s for s in has3 if s.shadow_blocked_by_fold],
    }
    out: dict[str, dict[str, object]] = {}
    for name, bucket in buckets.items():
        if not bucket:
            out[name] = {"count": 0}
            continue
        margins = [s.margin for s in bucket]
        measures = [s.measure for s in bucket]
        out[name] = {
            "count": len(bucket),
            "min_margin": min(margins),
            "avg_margin": sum(margins, F(0)) / len(margins),
            "min_measure": min(measures),
            "min_row": min(bucket, key=lambda s: (s.margin, max(s.row), s.row)).row,
            "energy_cut": energy_cut,
        }
    return out


def tournament_fingerprint(route_counts: Counter[str]) -> str:
    vertices = [
        ("A_circuit_free_proxy", 4, route_counts["A_circuit_free_proxy"]),
        ("A_high_4term_no3", 3, route_counts["A_high_4term_no3"]),
        ("B_has_3term_fold", 2, route_counts["B_has_3term_fold"]),
        ("B_fold_shadow_blocked", 1, route_counts["B_fold_shadow_blocked"]),
    ]

    def beats(a, b) -> bool:
        # Pair observable: (proof-cleanliness rank, sample coverage, label).
        # Switch favors the cleaner/easier proof state; fold-blocked is hardest.
        return (a[1], a[2], a[0]) > (b[1], b[2], b[0])

    scores = Counter({v[0]: 0 for v in vertices})
    cycles = 0
    for a, b in combinations(vertices, 2):
        scores[(a if beats(a, b) else b)[0]] += 1
    for a, b, c in combinations(vertices, 3):
        ab, bc, ca = beats(a, b), beats(b, c), beats(c, a)
        if (ab and bc and ca) or ((not ab) and (not bc) and (not ca)):
            cycles += 1
    path = " -> ".join(v[0] for v in sorted(vertices, key=lambda x: (x[1], x[2], x[0]), reverse=True))
    return f"score_hist={dict(sorted(Counter(scores.values()).items()))}; directed_3_cycles={cycles}; path={path}"


def targeted_rows(ns: list[int]) -> None:
    print("TARGETED 4-TERM-RICH CONTROL ROWS")
    print("AP has many 3-term folds and is tight; shifted_AP has AP-like 4-term energy but no a+b=c.")
    for n in ns:
        k = n - 1
        ap = tuple(range(1, n))
        shift = k + 1
        shifted = tuple(range(shift + 1, shift + k + 1))
        rows = [("AP_fold_rich", ap), ("shifted_AP_no3_high4", shifted)]
        if n == 14:
            rows.append(("Vstar_fold_wall", (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 24)))
        print(f"n={n:2d} k={k:2d}")
        for label, row in rows:
            stats = row_stats(row, n)
            print(
                f"  {label:22s}: M={fmt_frac(stats.M)} margin={fmt_frac(stats.margin)} "
                f"mu={fmt_frac(stats.measure)} rels={len(stats.relations):2d} "
                f"energy={stats.energy:4d} row={row}"
            )
    print()


def main() -> None:
    rng = random.Random(583)
    ns = [6, 7, 8, 9, 10, 11, 14]
    trials_by_n = {6: 260, 7: 240, 8: 220, 9: 180, 10: 150, 11: 120, 14: 70}
    print("S583 fold+sieve versus circuit-free randomness")
    print()
    print("LOCAL FOLD GATE")
    print("delta=1/n; compare independent 3-runner clearance (1-2d)^3 with fold gate x,y,x+y safe")
    for n in ns:
        d = 1 / n
        indep = (1 - 2 * d) ** 3
        fold = fold_gate_density(n)
        print(f"  n={n:2d}: independent={indep:.5f}; fold_gate~={fold:.5f}; fold_penalty={fold-indep:+.5f}")
    print()

    targeted_rows(ns)

    route_totals: Counter[str] = Counter()
    all_stats: list[RowStats] = []
    print("SAMPLE EXACT MAXIMIN BUCKETS")
    print("circuit-free proxy = no relation a+b=c; high-4term = top quartile additive energy inside no3")
    for n in ns:
        rows = sample_rows(n, trials_by_n[n], rng)
        stats = [row_stats(row, n) for row in rows]
        all_stats.extend(stats)
        summary = bucket_summary(stats)
        print(f"n={n:2d} k={n-1:2d} rows={len(stats)} delta={1/n:.5f}")
        for key in ["A_circuit_free_proxy", "A_high_4term_no3", "B_has_3term_fold", "B_fold_shadow_blocked"]:
            item = summary[key]
            route_totals[key] += int(item["count"])
            if item["count"] == 0:
                print(f"  {key:24s}: count=0")
                continue
            print(
                f"  {key:24s}: count={item['count']:3d}; "
                f"min(M-delta)={fmt_frac(item['min_margin'])}; avg={fmt_frac(item['avg_margin'])}; "
                f"min_mu={fmt_frac(item['min_measure'])}; min_row={item['min_row']}"
            )
        print()

    relation_rows = [s for s in all_stats if s.relations]
    denominator_shield = sum(1 for s in relation_rows if s.witness_den in {c for _a, _b, c in s.relations})
    active_relation = sum(1 for s in relation_rows if s.active_relation_pair)
    shadow_blocked = sum(1 for s in relation_rows if s.shadow_blocked_by_fold)
    print("FOLD/SIEVE DIAGNOSTICS")
    print(f"  3-term rows: {len(relation_rows)}")
    print(f"  exact witness denominator equals a folded c=a+b shield: {denominator_shield}")
    print(f"  exact active bottleneck uses a relation pair (a,b): {active_relation}")
    print(f"  best deletion-shadow witness blocked when c is reattached: {shadow_blocked}")
    if relation_rows:
        hard = sorted(relation_rows, key=lambda s: (s.margin, -len(s.relations), max(s.row)))[:5]
        print("  hardest fold rows:")
        for s in hard:
            print(
                f"    n={s.n} row={s.row} M={fmt_frac(s.M)} margin={fmt_frac(s.margin)} "
                f"rels={s.relations[:3]} witness={s.witness} den={s.witness_den} "
                f"shadow_margin={fmt_frac(s.shadow_margin)} blocked={s.shadow_blocked_by_fold}"
            )
    print()

    print("TOURNAMENT ANALYSIS")
    print("vertices: proof-state buckets, not runners")
    print("pair observable: (cleanliness rank, coverage); switch favors easier proof state")
    print(tournament_fingerprint(route_totals))
    print()
    print("ASSUMPTION CHALLENGE")
    print("Vertices considered: runners, 3-term relations, fold gates, additive-energy fibres, exact witnesses, proof buckets.")
    print("Chosen quotient: proof buckets plus relation diagnostics.")
    print("Preserved predicate: whether exact M(S) clears delta and which proof route plausibly explains the clearance.")
    print("Destroyed information: exact endpoint-cover components and prime-fibre certificates; those re-enter via HYP-2112/Phi and the paper sieve.")
    print("Challenged assumption: 4-term additive structure is a hard signal by itself; in this sample it is noise unless it carries a 3-term fold.")


if __name__ == "__main__":
    main()
