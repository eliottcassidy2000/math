#!/usr/bin/env python3
"""Exact audit for scale-separated inductive LRC reductions.

Notation in this workspace: LRC(n) means n-1 speeds and threshold 1/n.
If a seed B has positive safe set at threshold 1/n, then adding a sufficiently
large speed v preserves a safe point.  The only input is the finite-comb
discrepancy of the v-periodic unsafe set against a finite union of intervals.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from math import ceil


Interval = tuple[Fraction, Fraction]


def merge_intervals(arcs: list[Interval]) -> list[Interval]:
    if not arcs:
        return []
    arcs = sorted(arcs)
    merged: list[Interval] = []
    lo, hi = arcs[0]
    for a, b in arcs[1:]:
        if a <= hi:
            hi = max(hi, b)
        else:
            merged.append((lo, hi))
            lo, hi = a, b
    merged.append((lo, hi))
    return [(a, b) for a, b in merged if a < b]


def complement_intervals(arcs: list[Interval]) -> list[Interval]:
    arcs = merge_intervals(arcs)
    out: list[Interval] = []
    cursor = Fraction(0)
    for lo, hi in arcs:
        if cursor < lo:
            out.append((cursor, lo))
        cursor = max(cursor, hi)
    if cursor < 1:
        out.append((cursor, Fraction(1)))
    return out


def intersect_intervals(a: list[Interval], b: list[Interval]) -> list[Interval]:
    a = merge_intervals(a)
    b = merge_intervals(b)
    out: list[Interval] = []
    i = j = 0
    while i < len(a) and j < len(b):
        lo = max(a[i][0], b[j][0])
        hi = min(a[i][1], b[j][1])
        if lo < hi:
            out.append((lo, hi))
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return out


def measure(arcs: list[Interval]) -> Fraction:
    return sum((hi - lo for lo, hi in arcs), Fraction(0))


def unsafe_intervals_for_speed(speed: int, n: int) -> list[Interval]:
    """Return {t: ||speed*t|| < 1/n} as intervals in [0,1)."""
    half_width = Fraction(1, n * speed)
    arcs: list[Interval] = []
    for k in range(speed):
        center = Fraction(k, speed)
        lo = center - half_width
        hi = center + half_width
        if lo < 0:
            arcs.append((lo + 1, Fraction(1)))
            arcs.append((Fraction(0), hi))
        elif hi > 1:
            arcs.append((lo, Fraction(1)))
            arcs.append((Fraction(0), hi - 1))
        else:
            arcs.append((lo, hi))
    return merge_intervals(arcs)


def safe_set(speeds: tuple[int, ...], n: int) -> list[Interval]:
    arcs: list[Interval] = []
    for speed in speeds:
        arcs.extend(unsafe_intervals_for_speed(speed, n))
    return complement_intervals(arcs)


def fmt(q: Fraction) -> str:
    if q.denominator == 1:
        return str(q.numerator)
    return f"{q.numerator}/{q.denominator}"


@dataclass(frozen=True)
class InductionRow:
    name: str
    seed: tuple[int, ...]
    n: int
    added: tuple[int, ...]


def lower_bound_after_one(seed_safe: list[Interval], n: int, v: int) -> Fraction:
    """Comb lower bound for measure(seed_safe cap safe(v))."""
    mu = measure(seed_safe)
    c = len(seed_safe)
    return (Fraction(1) - Fraction(2, n)) * mu - Fraction(2 * c, v)


def sufficient_v(seed_safe: list[Interval], n: int) -> int | None:
    mu = measure(seed_safe)
    c = len(seed_safe)
    if mu <= 0:
        return None
    denom = (Fraction(1) - Fraction(2, n)) * mu
    return (2 * c * denom.denominator) // denom.numerator + 1


def audit_row(row: InductionRow) -> None:
    print(f"\n{row.name}")
    print("-" * 78)
    print(f"  LRC size parameter n={row.n}; threshold=1/{row.n}")
    print(f"  seed={row.seed}")
    seed_safe = safe_set(row.seed, row.n)
    mu = measure(seed_safe)
    print(f"  seed safe intervals={len(seed_safe)}")
    print(f"  seed safe measure={fmt(mu)} ~= {float(mu):.8f}")
    least_v = sufficient_v(seed_safe, row.n)
    print(f"  least integer v certified by comb bound={least_v or 'inf'}")
    for v in row.added:
        exact = intersect_intervals(seed_safe, safe_set((v,), row.n))
        exact_mu = measure(exact)
        predicted = (Fraction(1) - Fraction(2, row.n)) * mu
        bound = lower_bound_after_one(seed_safe, row.n, v)
        print(f"  add v={v}:")
        print(f"    seed_mu={fmt(mu)} intervals={len(seed_safe)}")
        print(f"    exact_after={fmt(exact_mu)} ~= {float(exact_mu):.8f}")
        print(f"    equidist_main={(fmt(predicted))} ~= {float(predicted):.8f}")
        print(f"    comb_lower={fmt(bound)} ~= {float(bound):.8f}")
        print(f"    positive_exact={exact_mu > 0} positive_by_comb={bound > 0}")


def print_scaling_guardrail(base: tuple[int, ...], n: int, scales: list[int]) -> None:
    print("\nDilation guardrail: pure size-only induction cannot be uniform")
    print("=" * 78)
    print(f"  base={base}, n={n}; dilation preserves measure but raises component count")
    print("  scale  safe_measure  intervals  least_certified_v")
    for q in scales:
        seed = tuple(q * b for b in base)
        ss = safe_set(seed, n)
        print(f"  {q:5d}  {fmt(measure(ss)):>12s}  {len(ss):9d}  {sufficient_v(ss, n):12d}")
    print()
    print("  Consequence: LRC(n-1) implies positive seed measure at threshold 1/n,")
    print("  but it does not give a size-only threshold for the next speed.  An")
    print("  induction needs scale normalization, a component/arc budget, or a")
    print("  bounded-core Node-2 reduction.")


def print_tournament_analysis() -> None:
    print("\nTournament Analysis over induction proof carriers")
    print("=" * 78)
    features = {
        "finite_comb_budget": {"preserves_predicate", "measure", "components", "threshold"},
        "scale_normalized_seed": {"components", "gcd", "bounded_core", "threshold"},
        "node3_large_speed": {"equidistribution", "threshold", "large_speed"},
        "node2_bounded_core": {"bounded_core", "finite_atlas", "cap"},
        "pure_size_induction": {"threshold"},
        "raw_runner_vertices": set(),
    }
    names = list(features)
    scores = Counter({name: 0 for name in names})
    adj = {name: set() for name in names}
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            if j <= i:
                continue
            key_a = (len(features[a]), "preserves_predicate" in features[a], -i)
            key_b = (len(features[b]), "preserves_predicate" in features[b], -j)
            if key_a >= key_b:
                adj[a].add(b)
                scores[a] += 1
            else:
                adj[b].add(a)
                scores[b] += 1
    cycles3 = 0
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            if j <= i:
                continue
            for k, c in enumerate(names):
                if k <= j:
                    continue
                if b in adj[a] and c in adj[b] and a in adj[c]:
                    cycles3 += 1
                if c in adj[a] and b in adj[c] and a in adj[b]:
                    cycles3 += 1
    path = sorted(names, key=lambda name: (scores[name], len(features[name]), name), reverse=True)
    print("  vertices are proof carriers, not runners")
    print("  observable=(LRC predicate retained, component budget, scale data)")
    print(f"  score_hist={dict(sorted(Counter(scores.values()).items()))}")
    print(f"  directed_3cycles={cycles3}")
    print("  Hamiltonian path=" + " > ".join(path))


def main() -> None:
    print("Scale-separated induction audit for LRC")
    print("=" * 78)
    print("Lemma audited: if A is the seed safe set at threshold 1/n,")
    print("  mu(A cap Safe_v) >= (1-2/n)mu(A) - 2*components(A)/v.")
    print("Thus adding a sufficiently large speed reduces LRC(n) to a smaller")
    print("seed problem, but the sufficient threshold depends on the seed's")
    print("component budget, not only on n.")

    rows = [
        InductionRow(
            "S46 AP-core seed with committed 30030 multiples",
            tuple(list(range(1, 12)) + [13]),
            14,
            (30030, 60060, 510510),
        ),
        InductionRow(
            "Sub-14 training: AP seed for n=10",
            tuple(range(1, 9)),
            10,
            (210, 420),
        ),
        InductionRow(
            "Small induction seed for n=8",
            tuple(range(1, 7)),
            8,
            (840,),
        ),
    ]
    for row in rows:
        audit_row(row)

    print_scaling_guardrail(tuple(list(range(1, 12)) + [13]), 14, [1, 2, 5, 10, 25, 50])
    print_tournament_analysis()

    print("\nProof-route conclusion")
    print("=" * 78)
    print("  Valid induction branch: LRC(n-1) gives nonempty seed-safe set at")
    print("  threshold 1/n; finite-comb equidistribution adds any sufficiently")
    print("  large next speed.")
    print("  Guardrail: dilation keeps the measure but increases interval count,")
    print("  so pure size-only induction is not uniform.  The missing effective")
    print("  input is a scale-normalized component bound or a bounded Node-2 atlas.")


if __name__ == "__main__":
    main()
