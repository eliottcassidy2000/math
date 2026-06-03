#!/usr/bin/env python3
"""S587: Lemma B as a denominator-divisibility sieve.

Main lens:

  Lemma B is stronger than the visible fold a+b=c.

At a pair-pinch clock t=m/(a+b), the denominator D=a+b is the local object.
If any runner v is divisible by D, then v*t is an integer for every m.  Thus
the direct fold c=a+b is only the first divisibility shield; the V* apex 24
shielding D=12 is the same mechanism one level higher.

This script alternates that structure with the Lemma-A control: when the low
denominators are not killed by observer-coupled divisibility, rows should sit
in the randomness/margin regime rather than at the delta floor.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations, permutations
from math import gcd, sqrt
import random
import statistics


def dist(x: F) -> F:
    x %= 1
    return min(x, 1 - x)


def primitive(row: tuple[int, ...]) -> tuple[int, ...]:
    g = 0
    for v in row:
        g = gcd(g, v)
    return tuple(sorted(v // g for v in row))


def exact_maximin(row: tuple[int, ...]) -> tuple[F, F]:
    candidates: set[F] = set()
    for a in row:
        candidates.add(F(1, 2 * a))
    for a, b in combinations(row, 2):
        for denom in (a + b, abs(a - b)):
            if denom:
                for m in range(1, denom):
                    candidates.add(F(m, denom))
    best = F(0)
    witness = F(0)
    for t in candidates:
        val = min(dist(v * t) for v in row)
        if val > best:
            best = val
            witness = t
    return best, witness


def safe_at(row: tuple[int, ...], t: F, n: int) -> bool:
    delta = F(1, n)
    return all(dist(v * t) >= delta for v in row)


def pair_sums(row: tuple[int, ...]) -> dict[int, tuple[tuple[int, int], ...]]:
    out: dict[int, list[tuple[int, int]]] = {}
    for a, b in combinations(sorted(row), 2):
        out.setdefault(a + b, []).append((a, b))
    return {d: tuple(pairs) for d, pairs in out.items()}


def visible_folds(row: tuple[int, ...]) -> tuple[tuple[int, int, int], ...]:
    vals = set(row)
    return tuple((a, b, a + b) for a, b in combinations(sorted(row), 2) if a + b in vals)


def unbalanced_count(row: tuple[int, ...]) -> int:
    vals = set(row)
    folds = sum(1 for a, b in combinations(sorted(row), 2) if a + b in vals)
    doubles = sum(1 for a in row if 2 * a in vals)
    return folds + doubles


def balanced_count(row: tuple[int, ...]) -> int:
    vals = set(row)
    aps = 0
    for a in row:
        for b in row:
            c = 2 * a - b
            if b < a < c and c in vals:
                aps += 1
    fibers = Counter(a + b for a, b in combinations(sorted(row), 2))
    energy = sum(m * (m - 1) // 2 for m in fibers.values())
    return aps + energy


def best_pinch_clearance(row: tuple[int, ...], denom: int) -> tuple[F, F]:
    best = F(0)
    best_t = F(0)
    for m in range(1, denom):
        t = F(m, denom)
        val = min(dist(v * t) for v in row)
        if val > best:
            best = val
            best_t = t
    return best, best_t


@dataclass(frozen=True)
class DenomCell:
    denom: int
    pairs: tuple[tuple[int, int], ...]
    shields: tuple[int, ...]
    best_clearance: F
    best_t: F

    @property
    def has_pair(self) -> bool:
        return bool(self.pairs)

    @property
    def shielded(self) -> bool:
        return bool(self.shields)


def low_denominator_cells(row: tuple[int, ...], n: int) -> tuple[DenomCell, ...]:
    sums = pair_sums(row)
    cells = []
    for denom in range(2, n + 1):
        pairs = sums.get(denom, tuple())
        shields = tuple(v for v in row if v % denom == 0)
        if pairs:
            best, best_t = best_pinch_clearance(row, denom)
        else:
            best, best_t = F(0), F(0)
        cells.append(DenomCell(denom, pairs, shields, best, best_t))
    return tuple(cells)


@dataclass(frozen=True)
class RowAudit:
    label: str
    row: tuple[int, ...]
    n: int
    M: F
    witness: F
    margin: F
    unbalanced: int
    balanced: int
    low_unshielded: tuple[int, ...]
    low_shielded: tuple[int, ...]
    n_pairs: tuple[tuple[int, int], ...]
    n_shields: tuple[int, ...]
    best_low_survivor_margin: F | None


def audit_row(label: str, row: tuple[int, ...], n: int) -> RowAudit:
    M, witness = exact_maximin(row)
    delta = F(1, n)
    cells = low_denominator_cells(row, n)
    low_unshielded = tuple(c.denom for c in cells if c.denom < n and c.has_pair and not c.shielded)
    low_shielded = tuple(c.denom for c in cells if c.denom < n and c.has_pair and c.shielded)
    n_cell = cells[-1]
    survivor_margins = [
        c.best_clearance - delta
        for c in cells
        if c.denom < n and c.has_pair and not c.shielded
    ]
    best_low = max(survivor_margins) if survivor_margins else None
    return RowAudit(
        label=label,
        row=row,
        n=n,
        M=M,
        witness=witness,
        margin=M - delta,
        unbalanced=unbalanced_count(row),
        balanced=balanced_count(row),
        low_unshielded=low_unshielded,
        low_shielded=low_shielded,
        n_pairs=n_cell.pairs,
        n_shields=n_cell.shields,
        best_low_survivor_margin=best_low,
    )


def residual_denominator(t: F, shield: int) -> int:
    """Denominator left after multiplying witness t by shield speed."""
    q = t.denominator
    p = t.numerator % q
    return q // gcd((shield * p) % q, q)


def deletion_shadow_audit(row: tuple[int, ...], n: int) -> tuple[int, int, int, Counter[int]]:
    """For visible folds a+b=c, delete c and see whether reattaching c blocks."""
    total = blocked = exact_zero = 0
    residuals: Counter[int] = Counter()
    delta = F(1, n)
    for _a, _b, c in visible_folds(row):
        shadow = tuple(v for v in row if v != c)
        if len(shadow) == len(row):
            continue
        sm, st = exact_maximin(shadow)
        if sm < delta:
            continue
        total += 1
        if not safe_at(row, st, n):
            blocked += 1
            res = residual_denominator(st, c)
            residuals[res] += 1
            if res == 1:
                exact_zero += 1
    return total, blocked, exact_zero, residuals


def target_rows() -> list[tuple[str, tuple[int, ...], int]]:
    return [
        ("AP_n14", tuple(range(1, 14)), 14),
        ("Vstar_n14", (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 24), 14),
        ("unit_shift_AP_n14", tuple(range(2, 15)), 14),
        ("far_shift_AP_n14", tuple(range(15, 28)), 14),
        ("doubled_apex_n14", tuple(range(1, 13)) + (26,), 14),
        ("hidden_S584_k9", (6, 11, 14, 15, 16, 18, 19, 23, 28), 10),
        ("AP_n10", tuple(range(1, 10)), 10),
        ("unit_shift_AP_n10", tuple(range(2, 11)), 10),
    ]


def fmt(x: F | None) -> str:
    if x is None:
        return "-"
    return f"{float(x):+.5f}"


def target_report() -> list[str]:
    lines = []
    lines.append("TARGETED DENOMINATOR-SIEVE LEDGER")
    lines.append(
        "  label                 n   M-delta   u  b   low_unshielded<D=n   low_shielded<D=n     D=n pairs/shields       best_low_survivor"
    )
    for label, row, n in target_rows():
        audit = audit_row(label, row, n)
        n_status = f"pairs={len(audit.n_pairs)},shields={list(audit.n_shields)}"
        lines.append(
            f"  {label:21s} {n:2d} {fmt(audit.margin):>9s} {audit.unbalanced:3d} {audit.balanced:3d} "
            f"{list(audit.low_unshielded)!s:20s} {list(audit.low_shielded)!s:20s} "
            f"{n_status:22s} {fmt(audit.best_low_survivor_margin):>9s}"
        )
    lines.append("")
    lines.append("  reading: AP and V* kill every pair denominator D<n that appears, but leave D=n unshielded.")
    lines.append("  V* explains the missing D=12 direct fold by divisibility: 24 is a shield for every m/12 pinch.")
    lines.append("  unit-shift AP kills D=n by the multiple 14, so the delta-clock survivor is gone and the row is loose.")
    return lines


def deletion_report() -> list[str]:
    lines = []
    lines.append("VISIBLE-FOLD DELETION SHADOWS")
    lines.append("  A deleted fold c is prime/divisibility-useful when the shadow witness has c*t in the danger band.")
    for label, row, n in target_rows():
        total, blocked, exact_zero, residuals = deletion_shadow_audit(row, n)
        common = dict(residuals.most_common(5))
        lines.append(
            f"  {label:21s}: shadow_tests={total:3d}, blocked={blocked:3d}, "
            f"exact_zero={exact_zero:3d}, residual_denoms={common}"
        )
    return lines


def sample_rows(n: int, trials: int, rng: random.Random) -> list[tuple[int, ...]]:
    rows: set[tuple[int, ...]] = set()
    k = n - 1
    universe = list(range(1, 5 * n + 1))
    attempts = 0
    while len(rows) < trials and attempts < 300 * trials:
        attempts += 1
        row = primitive(tuple(rng.sample(universe, k)))
        if len(row) == k:
            rows.add(row)
    return sorted(rows)


def pearson(xs: list[float], ys: list[float]) -> float | None:
    if len(xs) < 2:
        return None
    mx = statistics.mean(xs)
    my = statistics.mean(ys)
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    if vx == 0 or vy == 0:
        return None
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sqrt(vx * vy)


def sample_report() -> list[str]:
    rng = random.Random(587)
    lines = []
    lines.append("ALTERNATING A/B SAMPLE AUDIT")
    lines.append("  B-axis: low denominators killed by divisibility; A-axis: u=0/balanced-only margin control.")
    for n, trials in ((8, 160), (10, 140), (14, 60)):
        audits = [audit_row(f"sample_{i}", row, n) for i, row in enumerate(sample_rows(n, trials, rng))]
        a_control = [a for a in audits if a.unbalanced == 0]
        b_gate = [a for a in audits if not a.low_unshielded and a.n_pairs and not a.n_shields]
        n_killed = [a for a in audits if a.n_pairs and a.n_shields]
        fold_rich = sorted(audits, key=lambda a: a.unbalanced, reverse=True)[: max(1, len(audits) // 5)]

        def min_desc(bucket: list[RowAudit]) -> str:
            if not bucket:
                return "count=0"
            hard = min(bucket, key=lambda a: (a.margin, len(a.low_unshielded), max(a.row)))
            return (
                f"count={len(bucket):3d}; min_margin={fmt(hard.margin)}; "
                f"u={hard.unbalanced}; low_unshielded={list(hard.low_unshielded)}; "
                f"Dn_shields={list(hard.n_shields)}; row={hard.row}"
            )

        hard = [-float(a.margin) for a in audits]
        u = [float(a.unbalanced) for a in audits]
        b = [float(a.balanced) for a in audits]
        low_survivors = [float(len(a.low_unshielded)) for a in audits]
        n_shield = [1.0 if a.n_shields else 0.0 for a in audits]

        lines.append(f"  n={n:2d}, rows={len(audits)}")
        lines.append(f"    A_control_u0          {min_desc(a_control)}")
        lines.append(f"    B_sieve_complete      {min_desc(b_gate)}")
        lines.append(f"    D_eq_n_killed         {min_desc(n_killed)}")
        lines.append(f"    fold_rich_top20pct    {min_desc(fold_rich)}")
        lines.append(
            "    corr(-margin, unbalanced/balanced/low_survivors/n_shield)="
            f"{pearson(hard, u):+.3f}/{pearson(hard, b):+.3f}/"
            f"{pearson(hard, low_survivors):+.3f}/{pearson(hard, n_shield):+.3f}"
        )
    return lines


@dataclass(frozen=True)
class Lens:
    name: str
    lemma_b_delivery: int
    lemma_a_control: int
    divisibility_power: int
    observer_coupling: int
    maturity: int

    def key(self) -> tuple[int, int, int, int, int]:
        return (
            self.lemma_b_delivery,
            self.lemma_a_control,
            self.divisibility_power,
            self.observer_coupling,
            self.maturity,
        )


LENSES = [
    Lens("D_denominator_divisibility_shield", 5, 2, 5, 5, 4),
    Lens("visible_fold_a_plus_b_eq_c", 4, 2, 4, 5, 5),
    Lens("D_eq_n_unshielded_delta_clock", 4, 1, 3, 5, 4),
    Lens("Phi_endpoint_prime_residual", 4, 1, 5, 4, 5),
    Lens("augmentation_nonzero_count", 3, 3, 3, 5, 3),
    Lens("circuit_free_A_margin", 1, 5, 1, 1, 3),
    Lens("balanced_energy_background", 1, 3, 1, 0, 4),
]


def tournament_report() -> list[str]:
    names = [lens.name for lens in LENSES]
    edge: dict[tuple[str, str], bool] = {}
    for a, b in combinations(LENSES, 2):
        winner, loser = (a, b) if a.key() > b.key() else (b, a)
        edge[(winner.name, loser.name)] = True
        edge[(loser.name, winner.name)] = False

    scores = Counter()
    for u in names:
        scores[sum(1 for v in names if u != v and edge[(u, v)])] += 1

    c3 = 0
    for triple in combinations(names, 3):
        for a, b, c in permutations(triple):
            if edge[(a, b)] and edge[(b, c)] and edge[(c, a)]:
                c3 += 1
                break

    hpaths = 0
    first = None
    for path in permutations(names):
        if all(edge[(path[i], path[i + 1])] for i in range(len(path) - 1)):
            hpaths += 1
            if first is None:
                first = path

    ranking = [lens.name for lens in sorted(LENSES, key=lambda x: x.key(), reverse=True)]
    return [
        "TOURNAMENT ANALYSIS",
        "  vertices: proof lenses/denominator gates, not runners.",
        "  pair observable: (Lemma-B delivery, Lemma-A control, divisibility power, observer coupling, maturity).",
        "  switch: lexicographically larger observable gets the arc; declaration order is the tie Hamiltonian path.",
        f"  ranking: {ranking}",
        f"  score_histogram: {dict(sorted(scores.items()))}",
        f"  directed_3_cycles: {c3}",
        f"  Hamiltonian_path_count: {hpaths}",
        f"  first_Hamiltonian_path: {first}",
    ]


def assumption_challenge() -> list[str]:
    return [
        "ASSUMPTION CHALLENGE",
        "  Candidate vertices considered: runners, pair denominators D=a+b, shield speeds v with D|v, fold clauses, deleted-fold shadows, residues, endpoint components, Phi terms, balanced relations, and proof obligations.",
        "  Chosen quotient: denominator gates plus A/B proof buckets.",
        "  Preserved predicate: which low pair-pinches are killed before they can become lonely clocks, and whether D=n survives as the delta-clock.",
        "  Destroyed information: full endpoint-owner geometry and full prime-fiber languages; those re-enter at the Phi/endpoint residual.",
        "  Challenged assumption: Lemma B is only the direct relation a+b=c.  The useful structural unit is D|v at a pair denominator D=a+b.",
    ]


def main() -> None:
    print("S587 Lemma B fold-divisibility sieve, with Lemma A as control")
    print()
    print("PAIR-PINCH DIVISIBILITY FACT")
    print("  If D=a+b and t=m/D, then a*t+b*t=m.")
    print("  Any speed v with D|v has v*t in Z, so the D-pinch is shielded.")
    print("  The visible fold c=a+b is the special shield v=D; V*'s 24 shielding D=12 is the same gate.")
    print()
    for block in (
        target_report(),
        [""],
        deletion_report(),
        [""],
        sample_report(),
        [""],
        tournament_report(),
        [""],
        assumption_challenge(),
    ):
        for line in block:
            print(line)


if __name__ == "__main__":
    main()
