#!/usr/bin/env python3
"""S585: index-stratified extremal families for the LRC additive branch.

This probe extends the S578/S584 picture.

S578: translating an arithmetic progression preserves 4-term energy while
destroying visible 3-term folds.

S584: the destroyed folds are not gone; a 4-term pair-sum collision is two
3-term folds after adding the hidden sum node.

S585: make the translation index explicit.  For an interval
I(k,q)={q,...,q+k-1}, the pair-sum multiplicity profile is translation
invariant, but the original window clips it at a different index.  Visible
folds are exactly the clipped prefix; hidden virtual folds are the same tent
profile shifted into shells above the row.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations, permutations
from math import comb, gcd


def dist(x: F) -> F:
    x %= 1
    return min(x, 1 - x)


def exact_maximin(values: tuple[int, ...]) -> tuple[F, tuple[F, ...]]:
    candidates: set[F] = set()
    for v in values:
        candidates.add(F(1, 2 * v))
    for a, b in combinations(values, 2):
        for denom in (a + b, abs(a - b)):
            if denom == 0:
                continue
            for m in range(1, denom):
                candidates.add(F(m, denom))

    best = F(0)
    witnesses: list[F] = []
    for t in candidates:
        margin = min(dist(v * t) for v in values)
        if margin > best:
            best = margin
            witnesses = [t]
        elif margin == best:
            witnesses.append(t)
    return best, tuple(sorted(witnesses)[:6])


def pair_sum_profile(values: tuple[int, ...]) -> Counter[int]:
    return Counter(a + b for a, b in combinations(values, 2))


def visible_folds(values: tuple[int, ...]) -> list[tuple[int, int, int]]:
    present = set(values)
    return [(a, b, a + b) for a, b in combinations(values, 2) if a + b in present]


def four_term_count(values: tuple[int, ...]) -> int:
    profile = pair_sum_profile(values)
    return sum(comb(count, 2) for count in profile.values())


def hidden_shell_profile(values: tuple[int, ...]) -> Counter[int]:
    """Shell d means pair sum = max(values)+d.  Shell 0 aggregates visible folds."""
    max_v = max(values)
    present = set(values)
    shells: Counter[int] = Counter()
    for s, count in pair_sum_profile(values).items():
        if s in present:
            shells[0] += count
        elif s > max_v:
            shells[s - max_v] += count
        else:
            # This can occur for non-interval rows.  It is not used below, but keep
            # a signed shell so the function remains honest.
            shells[s - max_v] += count
    return shells


def top_hidden_sum(values: tuple[int, ...]) -> tuple[int, int, int]:
    """Return (sum, multiplicity, shell) for the strongest hidden virtual node."""
    max_v = max(values)
    present = set(values)
    hidden = [
        (s, count, s - max_v)
        for s, count in pair_sum_profile(values).items()
        if s not in present and s > max_v
    ]
    if not hidden:
        return (0, 0, 0)
    return max(hidden, key=lambda item: (item[1], -item[2], -item[0]))


def formula_visible_folds(k: int, q: int) -> int:
    """Count pairs 0<=i<j<k with i+j <= k-1-q."""
    limit = k - 1 - q
    if limit < 1:
        return 0
    return sum(1 for i, j in combinations(range(k), 2) if i + j <= limit)


@dataclass(frozen=True)
class IntervalRow:
    q: int
    folds: int
    formula_folds: int
    four_terms: int
    maximin: F
    witness: F
    top_hidden: int
    top_mult: int
    top_shell: int
    augmented_maximin: F


def interval_rows(k: int, q_values: range) -> list[IntervalRow]:
    rows: list[IntervalRow] = []
    for q in q_values:
        values = tuple(range(q, q + k))
        maximin, witnesses = exact_maximin(values)
        s, mult, shell = top_hidden_sum(values)
        augmented_maximin = F(0)
        if s:
            augmented = tuple(sorted(values + (s,)))
            augmented_maximin, _ = exact_maximin(augmented)
        rows.append(
            IntervalRow(
                q=q,
                folds=len(visible_folds(values)),
                formula_folds=formula_visible_folds(k, q),
                four_terms=four_term_count(values),
                maximin=maximin,
                witness=witnesses[0],
                top_hidden=s,
                top_mult=mult,
                top_shell=shell,
                augmented_maximin=augmented_maximin,
            )
        )
    return rows


def distinct_pair_completion_widths(width: int, levels: int) -> list[int]:
    """Interval support width under distinct-pair summand completion.

    If a level is a full interval of width w>=3, then pair sums of distinct
    nodes form a full interval of width 2w-3.
    """
    out = [width]
    for _ in range(levels):
        width = 2 * width - 3
        out.append(width)
    return out


def gcd_shells(C: int) -> dict[int, int]:
    counts: Counter[int] = Counter()
    for a in range(1, (C - 1) // 2 + 1):
        counts[gcd(a, C)] += 1
    return dict(sorted(counts.items()))


@dataclass(frozen=True)
class Lens:
    name: str
    certificate_rank: int
    recursion_rank: int
    preserves_labels: int
    implementation_cost: int
    maturity: int

    def key(self) -> tuple[int, int, int, int, int]:
        return (
            self.certificate_rank,
            self.recursion_rank,
            self.preserves_labels,
            -self.implementation_cost,
            self.maturity,
        )


def tournament_analysis() -> list[str]:
    lenses = [
        Lens("Phi_gap_value", 5, 3, 5, 3, 4),
        Lens("visible_fold_index", 4, 4, 4, 2, 4),
        Lens("C_gcd_shell_index", 4, 5, 5, 3, 3),
        Lens("hidden_sum_shell_index", 3, 5, 4, 2, 2),
        Lens("dyadic_debt_depth", 3, 5, 3, 3, 3),
        Lens("raw_4term_energy", 1, 2, 1, 1, 4),
    ]
    names = [lens.name for lens in lenses]
    edge: dict[tuple[str, str], bool] = {}
    for a, b in combinations(lenses, 2):
        winner, loser = (a, b) if a.key() > b.key() else (b, a)
        edge[(winner.name, loser.name)] = True
        edge[(loser.name, winner.name)] = False

    scores = Counter()
    for u in names:
        scores[sum(1 for v in names if v != u and edge[(u, v)])] += 1

    directed_3_cycles = 0
    for a, b, c in combinations(names, 3):
        cyc = False
        for x, y, z in permutations((a, b, c)):
            if edge[(x, y)] and edge[(y, z)] and edge[(z, x)]:
                cyc = True
                break
        directed_3_cycles += int(cyc)

    hpaths = 0
    first_path: tuple[str, ...] | None = None
    for path in permutations(names):
        if all(edge[(path[i], path[i + 1])] for i in range(len(path) - 1)):
            hpaths += 1
            if first_path is None:
                first_path = path

    ranked = [lens.name for lens in sorted(lenses, key=lambda lens: lens.key(), reverse=True)]
    return [
        "Tournament Analysis: index lenses, not runners",
        "  observable=(certificate_rank, recursion_rank, preserves_labels, -cost, maturity)",
        "  switch=lexicographically larger observable; tie path is declaration order",
        f"  ranking={ranked}",
        f"  score_hist={dict(sorted(scores.items()))}",
        f"  directed_3_cycles={directed_3_cycles}",
        f"  Hamiltonian_path_count={hpaths}",
        f"  first_path={first_path}",
    ]


def print_interval_table(k: int) -> None:
    delta = F(1, k + 1)
    rows = interval_rows(k, range(1, k + 1))
    print(f"Translated AP index ladder I(k,q)={{q,...,q+k-1}} with k={k}, delta=1/{k+1}")
    print("4-term count is translation-invariant; visible folds are a clipped prefix of the same pair-sum tent.")
    print(
        "  q  folds formula  4term   M       M/delta  witness    hidden(s,mult,shell)  M(+s)   drop"
    )
    for row in rows:
        drop = row.maximin - row.augmented_maximin
        print(
            f" {row.q:2d} {row.folds:6d} {row.formula_folds:7d} {row.four_terms:6d} "
            f"{float(row.maximin):7.4f} {float(row.maximin / delta):8.3f} "
            f"{str(row.witness):>8s} "
            f"({row.top_hidden:2d},{row.top_mult:2d},{row.top_shell:2d}) "
            f"{float(row.augmented_maximin):7.4f} {float(drop):7.4f}"
        )
    mismatches = [(row.q, row.folds, row.formula_folds) for row in rows if row.folds != row.formula_folds]
    four_terms = sorted({row.four_terms for row in rows})
    print(f"  formula_check_mismatches={mismatches}")
    print(f"  distinct_4term_counts_across_translation={four_terms}")


def print_shell_examples(k: int) -> None:
    print()
    print("Hidden shell profiles for selected translated APs")
    for q in (1, 2, 6, 12, 13):
        values = tuple(range(q, q + k))
        shells = hidden_shell_profile(values)
        visible = shells.get(0, 0)
        first_shells = [(d, shells[d]) for d in sorted(shells) if d > 0][:12]
        print(f"  q={q:2d}: visible_shell0={visible:2d}; first_hidden_shells={first_shells}")


def print_completion_recursion(k: int) -> None:
    print()
    print("Recursive summand completion widths")
    widths = distinct_pair_completion_widths(k, 5)
    print(f"  distinct-pair interval support widths from k={k}: {widths}")
    print("  recurrence: w_{r+1}=2*w_r-3, so w_r=2^r*(k-3)+3 for interval supports")
    print("  reading: each virtual layer is the same clipped tent mechanism at doubled scale")


def print_cross_index_table() -> None:
    print()
    print("Cross-index coordinates already present in the repo")
    print("  additive translation q: clips visible folds, preserves 4-term energy")
    print(f"  C=2n-1 shell index at n=14: C=27, antipodal shell gcd counts={gcd_shells(27)}")
    print("  HYP-2091 parity index: even n -> clean polygon; odd n -> wall/tie mesh")
    print("  HYP-2079 dyadic depth h: gap halves, boundary debt doubles, gap*debt is conserved")
    print("  S584 virtual depth: 4-term collisions become 3-folds after adjoining hidden sum nodes")


def main() -> None:
    print("=" * 88)
    print("S585: index-stratified extremal families and recursive hidden folds")
    print("=" * 88)
    print_interval_table(13)
    print_shell_examples(13)
    print_completion_recursion(13)
    print_cross_index_table()
    print()
    for line in tournament_analysis():
        print(line)
    print()
    print("Synthesis")
    print("  The extremal AP is not one point but the index-0 face of a translation family.")
    print("  Increasing q moves the same pair-sum tent from visible 3-folds into hidden shells.")
    print("  Raw 4-term energy cannot see q; hidden-shell and visible-fold indices can.")
    print("  Recursing S -> S+S turns the hidden shells into visible folds one layer higher.")
    print("  The proof object should therefore be an indexed summand sheaf over")
    print("  (parity lane, C-gcd shell, visible fold clip, hidden sum shell, dyadic debt depth).")


if __name__ == "__main__":
    main()
