#!/usr/bin/env python3
"""S584: four-term additive structure as hidden virtual 3-term folds.

The user question is the S577 fork:

* Lemma A: no 3-term additive circuit should imply LRC margin by
  randomness/equidistribution.
* Lemma B: a 3-term relation v_c=v_a+v_b is a literal phase fold.

This probe tests the in-between case: rows with no visible 3-term fold but
with 4-term pair-sum collisions a+b=c+d.  Each such collision is encoded by a
hidden virtual speed s=a+b=c+d, producing two 3-term folds (a,b,s) and
(c,d,s) after augmentation.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations, permutations
from math import gcd
import random


def dist(x: F) -> F:
    x = x % 1
    return min(x, 1 - x)


def prim(values: tuple[int, ...]) -> tuple[int, ...]:
    g = 0
    for value in values:
        g = gcd(g, value)
    return tuple(sorted(value // g for value in values))


def mexact_witnesses(values: tuple[int, ...]) -> tuple[F, tuple[F, ...]]:
    candidates: set[F] = set()
    for vi in values:
        candidates.add(F(1, 2 * vi))
    for vi, vj in combinations(values, 2):
        for denom in (vi + vj, abs(vi - vj)):
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
    return best, tuple(sorted(witnesses)[:8])


def three_terms(values: tuple[int, ...]) -> list[tuple[int, int, int]]:
    present = set(values)
    return [
        (a, b, a + b)
        for a, b in combinations(sorted(values), 2)
        if a + b in present
    ]


def pair_sum_fibers(values: tuple[int, ...]) -> dict[int, list[tuple[int, int]]]:
    fibers: dict[int, list[tuple[int, int]]] = {}
    for a, b in combinations(sorted(values), 2):
        fibers.setdefault(a + b, []).append((a, b))
    return {s: pairs for s, pairs in fibers.items() if len(pairs) >= 2}


def four_term_count(values: tuple[int, ...]) -> int:
    return sum(len(pairs) * (len(pairs) - 1) // 2 for pairs in pair_sum_fibers(values).values())


def best_hidden_sum(values: tuple[int, ...]) -> tuple[int, list[tuple[int, int]]] | None:
    hidden = [
        (s, pairs)
        for s, pairs in pair_sum_fibers(values).items()
        if s not in set(values)
    ]
    if not hidden:
        return None
    return max(hidden, key=lambda item: (len(item[1]), item[0]))


@dataclass
class RowWitness:
    values: tuple[int, ...]
    margin: F
    delta: F
    hidden_sum: int
    hidden_pairs: tuple[tuple[int, int], ...]
    virtual_margin: F
    virtual_delta: F
    virtual_at_witness: F

    @property
    def drop(self) -> F:
        return self.margin - self.virtual_margin


def virtual_witness(values: tuple[int, ...], margin: F, witnesses: tuple[F, ...]) -> RowWitness | None:
    hidden = best_hidden_sum(values)
    if hidden is None:
        return None
    s, pairs = hidden
    augmented = tuple(sorted(values + (s,)))
    virtual_margin, _ = mexact_witnesses(augmented)
    virtual_at_witness = max(dist(s * t) for t in witnesses)
    k = len(values)
    return RowWitness(
        values=values,
        margin=margin,
        delta=F(1, k + 1),
        hidden_sum=s,
        hidden_pairs=tuple(pairs),
        virtual_margin=virtual_margin,
        virtual_delta=F(1, k + 2),
        virtual_at_witness=virtual_at_witness,
    )


def fixed_sum_deformation_report() -> list[str]:
    lines = []
    lines.append("Fixed hidden-sum deformation families: choose r disjoint pairs x+(s-x).")
    lines.append("Rows with visible 3-term folds are filtered out.")
    for r in (2, 3):
        for s in (10, 12, 14, 16, 18, 20):
            pairs = [(x, s - x) for x in range(1, s // 2 + 1) if x < s - x]
            rows = []
            for chosen in combinations(pairs, r):
                values = tuple(sorted(v for pair in chosen for v in pair))
                if len(set(values)) != 2 * r or three_terms(values):
                    continue
                margin, _ = mexact_witnesses(values)
                rows.append((margin, values))
            if not rows:
                continue
            rows.sort()
            delta = F(1, 2 * r + 1)
            lines.append(
                f"  r={r}, hidden_sum={s:2d}, rows={len(rows):2d}, "
                f"M range={float(rows[0][0]):.4f}..{float(rows[-1][0]):.4f}, "
                f"min margin vs delta={float(rows[0][0] - delta):+.4f}, "
                f"min={rows[0][1]}, max={rows[-1][1]}"
            )
    return lines


def random_audit() -> tuple[list[str], list[RowWitness]]:
    rng = random.Random(584)
    lines = []
    hidden_examples: list[RowWitness] = []
    lines.append("Random audit by additive type; B=3k+4 as in S577.")
    header = (
        "  k  attempts  cf_seen  rich_seen  three_seen  min_cf_margin  "
        "min_rich_margin  min_three_margin  min_virtual_aug_margin"
    )
    lines.append(header)
    for k in range(4, 11):
        delta = F(1, k + 1)
        attempts = 0
        target_cf = 90 if k <= 8 else 45
        max_attempts = 26000
        three_eval_cap = 120
        three_evals = 0
        B = 3 * k + 4
        cf_seen = rich_seen = three_seen = 0
        min_cf: tuple[F, tuple[int, ...]] | None = None
        min_rich: tuple[F, tuple[int, ...]] | None = None
        min_three: tuple[F, tuple[int, ...]] | None = None
        min_virtual: RowWitness | None = None
        while cf_seen < target_cf and attempts < max_attempts:
            attempts += 1
            values = prim(tuple(sorted(rng.sample(range(1, B + 1), k))))
            triples = three_terms(values)
            if triples:
                three_seen += 1
                if three_evals < three_eval_cap:
                    three_evals += 1
                    margin, _ = mexact_witnesses(values)
                    if min_three is None or margin < min_three[0]:
                        min_three = (margin, values)
                continue

            margin, witnesses = mexact_witnesses(values)
            cf_seen += 1
            if min_cf is None or margin < min_cf[0]:
                min_cf = (margin, values)
            if four_term_count(values) >= 2:
                rich_seen += 1
                if min_rich is None or margin < min_rich[0]:
                    min_rich = (margin, values)
                vw = virtual_witness(values, margin, witnesses)
                if vw is not None:
                    hidden_examples.append(vw)
                    if min_virtual is None or vw.virtual_margin - vw.virtual_delta < min_virtual.virtual_margin - min_virtual.virtual_delta:
                        min_virtual = vw

        def fmt_margin(item: tuple[F, tuple[int, ...]] | None) -> str:
            if item is None:
                return "n/a"
            return f"{float(item[0] - delta):+.4f} {item[1]}"

        virtual_text = "n/a"
        if min_virtual is not None:
            virtual_text = (
                f"{float(min_virtual.virtual_margin - min_virtual.virtual_delta):+.4f} "
                f"s={min_virtual.hidden_sum} V={min_virtual.values}"
            )
        lines.append(
            f"  {k:2d} {attempts:9d} {cf_seen:8d} {rich_seen:10d} {three_seen:11d} "
            f"{fmt_margin(min_cf):>29s} {fmt_margin(min_rich):>31s} "
            f"{fmt_margin(min_three):>31s} {virtual_text}"
        )
    hidden_examples.sort(
        key=lambda w: (
            w.virtual_margin - w.virtual_delta,
            -w.drop,
            w.margin - w.delta,
        )
    )
    return lines, hidden_examples[:8]


@dataclass(frozen=True)
class Lens:
    name: str
    certificate_safety: int
    proof_payoff: int
    hidden_exposure: int
    implementation_cost: int
    maturity: int

    def key(self) -> tuple[int, int, int, int, int]:
        return (
            self.certificate_safety,
            self.proof_payoff,
            self.hidden_exposure,
            -self.implementation_cost,
            self.maturity,
        )


LENSES = [
    Lens("Phi_gap_gate", 5, 5, 3, 3, 4),
    Lens("visible_3_fold_shield", 4, 5, 5, 2, 3),
    Lens("hidden_4_virtual_sum", 3, 4, 5, 3, 2),
    Lens("circuit_free_randomness", 2, 5, 2, 4, 2),
    Lens("fixed_sum_deformation_fiber", 2, 3, 4, 2, 2),
    Lens("raw_additive_energy", 1, 2, 2, 1, 4),
]


def tournament_fingerprint() -> list[str]:
    names = [lens.name for lens in LENSES]
    matrix = {u: {v: False for v in names if v != u} for u in names}
    for a, b in combinations(LENSES, 2):
        winner, loser = (a, b) if a.key() > b.key() else (b, a)
        matrix[winner.name][loser.name] = True
        matrix[loser.name][winner.name] = False

    scores = Counter()
    for u in names:
        scores[sum(1 for v in names if v != u and matrix[u][v])] += 1

    c3 = 0
    for triple in combinations(names, 3):
        cyclic = False
        for a, b, c in permutations(triple):
            if matrix[a][b] and matrix[b][c] and matrix[c][a]:
                cyclic = True
                break
        c3 += int(cyclic)

    hpaths = 0
    example = None
    for path in permutations(names):
        if all(matrix[path[i]][path[i + 1]] for i in range(len(path) - 1)):
            hpaths += 1
            if example is None:
                example = path

    ranked = sorted(LENSES, key=lambda lens: lens.key(), reverse=True)
    lines = []
    lines.append("Tournament Analysis on explanatory lenses")
    lines.append("  vertices: proof/search lenses, not runners")
    lines.append("  observable: (certificate_safety, proof_payoff, hidden_exposure, -cost, maturity)")
    lines.append("  switch: lexicographically larger observable wins; tie path is declaration order")
    lines.append(f"  ranking: {[lens.name for lens in ranked]}")
    lines.append(f"  score histogram: {dict(sorted(scores.items()))}")
    lines.append(f"  directed 3-cycles: {c3}")
    lines.append(f"  Hamiltonian path count: {hpaths}")
    lines.append(f"  first Hamiltonian path: {example}")
    return lines


def main() -> None:
    print("S584 four-term hidden-fold probe")
    print("Question: can 4-term structure with no visible 3-term fold hide information that reappears as virtual 3-folds?")
    print()

    for line in fixed_sum_deformation_report():
        print(line)
    print()

    random_lines, examples = random_audit()
    for line in random_lines:
        print(line)
    print()

    print("Most constrained virtual encodings found in the random audit")
    for ex in examples:
        encoded = [(a, b, ex.hidden_sum) for a, b in ex.hidden_pairs]
        print(
            f"  V={ex.values}, M={float(ex.margin):.4f}, margin={float(ex.margin - ex.delta):+.4f}; "
            f"hidden s={ex.hidden_sum}, pairs={list(ex.hidden_pairs)}, encoded_3folds={encoded}; "
            f"M(V+s)={float(ex.virtual_margin):.4f}, aug_margin={float(ex.virtual_margin - ex.virtual_delta):+.4f}, "
            f"drop={float(ex.drop):.4f}, best_virtual_dist_at_original_witness={float(ex.virtual_at_witness):.4f}"
        )
    print()

    print("Interpretation")
    print("  1. A 4-term collision is two 3-term folds sharing a hidden sum node s.")
    print("  2. No visible 3-term fold means s is outside the speed row, so raw additive energy can be invisible to the exact LRC witness.")
    print("  3. Fixed hidden-sum deformation families move the individual speeds while preserving the pair-sum fiber; M varies inside that fiber, so the 4-term label alone is not a certificate.")
    print("  4. The useful certificate object is labelled: original speeds + pair edges + hidden sum nodes + virtual fold pressure.")
    print("  5. This matches HYP-2113: keep middle labels until a Phi, fold/shield, CRT, or endpoint gate consumes them.")
    print()

    for line in tournament_fingerprint():
        print(line)


if __name__ == "__main__":
    main()
