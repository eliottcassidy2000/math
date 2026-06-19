#!/usr/bin/env python3
"""HYP-2639/T887: relation-covered rows, additive energy, and GAP structure.

This scout tests the tempting proof slogan:

    every element in a small relation -> high additive energy -> Freiman GAP.

The old summand/multiplicand graph work makes that slogan too coarse.  Pair-sum
collisions are balanced and translation-invariant; visible folds are
observer-coupled and translation-sensitive.  Shifted APs have maximal additive
energy but no visible folds, so the useful invariant must remember where the
collision node lands.

Tournament Analysis vertices are proof obligations rather than runners.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, replace
from fractions import Fraction as F
from itertools import combinations, permutations
from math import gcd
import random


def primitive(row: tuple[int, ...]) -> tuple[int, ...]:
    g = 0
    for v in row:
        g = gcd(g, abs(v))
    if g <= 1:
        return tuple(sorted(row))
    return tuple(sorted(v // g for v in row))


def dist(x: F) -> F:
    x %= 1
    return min(x, 1 - x)


def exact_maximin(row: tuple[int, ...]) -> F:
    """Exact maximin over standard single and pair pinch candidates."""
    values = tuple(v for v in row if v != 0)
    candidates: set[F] = {F(0)}
    for v in values:
        for m in range(v):
            candidates.add(F(2 * m + 1, 2 * v))
    for a, b in combinations(values, 2):
        for denom in (a + b, abs(a - b)):
            if denom == 0:
                continue
            for m in range(1, denom):
                candidates.add(F(m, denom))
    best = F(0)
    for t in candidates:
        best = max(best, min(dist(v * t) for v in values))
    return best


def full_sumset(row: tuple[int, ...]) -> Counter[int]:
    return Counter(a + b for a in row for b in row)


def restricted_pair_sums(row: tuple[int, ...]) -> Counter[int]:
    return Counter(a + b for a, b in combinations(row, 2))


def visible_folds(row: tuple[int, ...]) -> list[tuple[int, int, int]]:
    present = set(row)
    return [(a, b, a + b) for a, b in combinations(row, 2) if a + b in present]


def three_ap_triples(row: tuple[int, ...]) -> list[tuple[int, int, int]]:
    present = set(row)
    triples = []
    for a, c in combinations(row, 2):
        s = a + c
        if s % 2 == 0:
            b = s // 2
            if b in present and b != a and b != c:
                triples.append((a, b, c))
    return triples


def collision_fibers(row: tuple[int, ...]) -> dict[int, tuple[tuple[int, int], ...]]:
    fibers: dict[int, list[tuple[int, int]]] = {}
    for a, b in combinations(row, 2):
        fibers.setdefault(a + b, []).append((a, b))
    return {s: tuple(ps) for s, ps in sorted(fibers.items()) if len(ps) >= 2}


def collision_relation_count(fibers: dict[int, tuple[tuple[int, int], ...]]) -> int:
    return sum(len(ps) * (len(ps) - 1) // 2 for ps in fibers.values())


def relation_covered(row: tuple[int, ...]) -> tuple[set[int], dict[str, set[int]]]:
    folds = visible_folds(row)
    aps = three_ap_triples(row)
    fibers = collision_fibers(row)
    by_kind = {
        "fold": {x for rel in folds for x in rel},
        "midpoint": {x for rel in aps for x in rel},
        "pair_collision": {x for pairs in fibers.values() for pair in pairs for x in pair},
    }
    all_covered = set().union(*by_kind.values()) if by_kind else set()
    return all_covered, by_kind


def multiplicand_profile(row: tuple[int, ...]) -> tuple[int, int, int, int]:
    """Return coarse/fine pair-sum nodes and divisor-obstructed node counts."""
    nonzero = tuple(v for v in row if v != 0)
    if not nonzero:
        return (0, 0, 0, 0)
    max_v = max(nonzero)
    nodes = sorted(restricted_pair_sums(nonzero))
    coarse = [c for c in nodes if c <= max_v]
    fine = [c for c in nodes if c > max_v]
    coarse_blocked = sum(1 for c in coarse if any(v % c == 0 for v in nonzero))
    fine_blocked = sum(1 for c in fine if any(v % c == 0 for v in nonzero))
    return (len(coarse), len(fine), coarse_blocked, fine_blocked)


@dataclass(frozen=True)
class RowReport:
    name: str
    row: tuple[int, ...]
    n: int
    maximin: F
    full_sumset_size: int
    restricted_sumset_size: int
    ordered_energy: int
    sidon_excess: int
    pair_collision_relations: int
    visible_fold_count: int
    midpoint_count: int
    hidden_collision_nodes: int
    visible_collision_nodes: int
    covered_nonzero: int
    nonzero_count: int
    all_nonzero_covered: bool
    freiman_3k4_window: bool
    coarse_nodes: int
    fine_nodes: int
    coarse_blocked: int
    fine_blocked: int


def row_report(name: str, row: tuple[int, ...], n: int) -> RowReport:
    row = primitive(tuple(sorted(row)))
    k = len(row)
    nonzero = tuple(v for v in row if v != 0)
    fs = full_sumset(row)
    rs = restricted_pair_sums(row)
    fibers = collision_fibers(row)
    present = set(row)
    hidden_nodes = sum(1 for s in fibers if s not in present)
    visible_nodes = sum(1 for s in fibers if s in present)
    covered, _by_kind = relation_covered(row)
    covered_nonzero = len(set(nonzero) & covered)
    coarse, fine, coarse_blocked, fine_blocked = multiplicand_profile(row)
    energy = sum(c * c for c in fs.values())
    return RowReport(
        name=name,
        row=row,
        n=n,
        maximin=exact_maximin(row),
        full_sumset_size=len(fs),
        restricted_sumset_size=len(rs),
        ordered_energy=energy,
        sidon_excess=energy - (2 * k * k - k),
        pair_collision_relations=collision_relation_count(fibers),
        visible_fold_count=len(visible_folds(row)),
        midpoint_count=len(three_ap_triples(row)),
        hidden_collision_nodes=hidden_nodes,
        visible_collision_nodes=visible_nodes,
        covered_nonzero=covered_nonzero,
        nonzero_count=len(nonzero),
        all_nonzero_covered=(covered_nonzero == len(nonzero)),
        freiman_3k4_window=(len(fs) <= 3 * k - 4),
        coarse_nodes=coarse,
        fine_nodes=fine,
        coarse_blocked=coarse_blocked,
        fine_blocked=fine_blocked,
    )


def sample_relation_covered_rows(trials: int = 12000) -> list[RowReport]:
    rng = random.Random(2637)
    found: list[RowReport] = []
    seen: set[tuple[int, ...]] = set()
    while len(found) < 8 and trials > 0:
        trials -= 1
        row = primitive((0,) + tuple(rng.sample(range(1, 42), 7)))
        if row in seen or len(row) != 8 or max(row) <= 20:
            continue
        seen.add(row)
        report = row_report("sample", row, n=9)
        if report.all_nonzero_covered and not report.freiman_3k4_window:
            found.append(report)
    return sorted(
        found,
        key=lambda r: (
            -r.visible_fold_count,
            -r.pair_collision_relations,
            r.full_sumset_size,
            r.row,
        ),
    )


@dataclass(frozen=True)
class Lane:
    name: str
    certificate_power: int
    decoy_resistance: int
    preserves_sign: int
    graph_bridge: int
    maturity: int

    def key(self) -> tuple[int, int, int, int, int]:
        return (
            self.certificate_power,
            self.decoy_resistance,
            self.preserves_sign,
            self.graph_bridge,
            self.maturity,
        )


def tournament_analysis() -> tuple[Counter[int], int, int, tuple[str, ...]]:
    lanes = [
        Lane("observer_coupled_visible_folds", 5, 5, 5, 4, 4),
        Lane("low_hidden_summand_shells", 4, 4, 4, 5, 3),
        Lane("multiplicand_clearance_sieve", 4, 3, 3, 5, 4),
        Lane("relation_coverage_hypergraph", 3, 4, 4, 3, 3),
        Lane("freiman_small_doubling_GAP", 3, 2, 2, 3, 5),
        Lane("balanced_pair_energy", 2, 1, 1, 3, 5),
        Lane("raw_runner_vertices", 1, 1, 1, 1, 5),
    ]
    names = [lane.name for lane in lanes]
    edge: dict[tuple[str, str], bool] = {}
    for a, b in combinations(lanes, 2):
        winner, loser = (a, b) if a.key() >= b.key() else (b, a)
        edge[(winner.name, loser.name)] = True
        edge[(loser.name, winner.name)] = False
    scores = Counter(sum(1 for v in names if v != u and edge[(u, v)]) for u in names)
    cycles = 0
    for a, b, c in combinations(names, 3):
        if any(edge[(x, y)] and edge[(y, z)] and edge[(z, x)] for x, y, z in permutations((a, b, c))):
            cycles += 1
    hpaths = 0
    first_path: tuple[str, ...] | None = None
    for path in permutations(names):
        if all(edge[(path[i], path[i + 1])] for i in range(len(path) - 1)):
            hpaths += 1
            if first_path is None:
                first_path = path
    return scores, cycles, hpaths, first_path or tuple()


def print_report(report: RowReport) -> None:
    k = len(report.row)
    gap = report.full_sumset_size - (2 * k - 1)
    print(f"{report.name:22s} k={k:2d} span={max(report.row)-min(report.row):3d} "
          f"M={float(report.maximin):.6f} M*n={float(report.maximin * report.n):.3f}")
    print(
        f"  sumset |S+S|={report.full_sumset_size:2d} gap={gap:2d} "
        f"Freiman3k-4={report.freiman_3k4_window} "
        f"restricted={report.restricted_sumset_size:2d}"
    )
    print(
        f"  energy={report.ordered_energy:4d} sidon_excess={report.sidon_excess:4d} "
        f"pair_collision_rel={report.pair_collision_relations:3d} "
        f"visible_folds={report.visible_fold_count:2d} midpoints={report.midpoint_count:2d}"
    )
    print(
        f"  collision_nodes visible={report.visible_collision_nodes:2d} "
        f"hidden={report.hidden_collision_nodes:2d}; "
        f"covered_nonzero={report.covered_nonzero}/{report.nonzero_count} "
        f"all={report.all_nonzero_covered}"
    )
    print(
        f"  multiplicand pair-sum nodes coarse/fine={report.coarse_nodes}/{report.fine_nodes}; "
        f"blocked coarse/fine={report.coarse_blocked}/{report.fine_blocked}"
    )


def main() -> None:
    rows = [
        ("AP13", tuple(range(1, 14)), 14),
        ("shiftAP13", tuple(range(13, 26)), 14),
        ("Vstar", tuple(list(range(1, 12)) + [13, 24]), 14),
        ("KPS_third_pocket_A", (0, 3, 5, 16, 28, 30, 33, 35), 9),
        ("KPS_third_pocket_B", (0, 4, 12, 15, 20, 21, 25, 31), 9),
    ]
    reports = [row_report(name, row, n) for name, row, n in rows]
    samples = sample_relation_covered_rows()

    print("=" * 88)
    print("HYP-2639/T887 relation-covered GAP scout")
    print("=" * 88)
    print("Small relation motifs tracked:")
    print("  fold:       a+b=c       coefficient sum 1  -> observer-coupled / translation-sensitive")
    print("  midpoint:   a+c=2b      coefficient sum 0  -> balanced but 2-adic midpoint")
    print("  collision:  a+b=c+d     coefficient sum 0  -> balanced pair energy")
    print("Multiplicand layer: pair-sum denominator C clears nonzero w unless C divides w.")
    print()

    print("CALIBRATION ROWS")
    for report in reports:
        print_report(report)

    print("\nSAMPLED WIDE ROWS WITH EVERY NONZERO ELEMENT IN A SMALL ADDITIVE MOTIF")
    for i, report in enumerate(samples[:6], 1):
        print_report(replace(report, name=f"sample_{i}"))

    print("\nREADING")
    print("  1. Direct Freiman small-doubling is too strong for the third pocket:")
    print("     KPS relation-covered examples have |S+S| well above 3k-4, so they are not")
    print("     automatically one-dimensional GAP rows by the elementary 3k-4 window.")
    print("  2. The useful condition is relation coverage plus shell location.  Shifted AP")
    print("     has maximal additive energy but zero visible folds; AP has the same energy")
    print("     and many observer-coupled folds.  Energy without the summand-node address")
    print("     is a decoy.")
    print("  3. Addition supplies candidate denominators C in S+S; multiplication tests")
    print("     clearance by divisibility C|w.  Positive/negative signs in a relation say")
    print("     whether the relation is balanced (even scalar, translation-invariant) or")
    print("     observer-coupled (odd marked, translation-sensitive).")
    print("  4. The proof target should be a relation-hypergraph regularity lemma:")
    print("     either a dissociated stranger peels, or the small-relation hypergraph has")
    print("     enough low visible/hidden summand-shell payload that the signed reciprocal")
    print("     tail can be channelized before absolute values.")

    scores, cycles, hpaths, path = tournament_analysis()
    print("\nTOURNAMENT ANALYSIS")
    print(f"  score_hist={dict(sorted(scores.items()))}")
    print(f"  directed_3cycles={cycles}")
    print(f"  hamiltonian_paths={hpaths}")
    print("  tie_hamiltonian_path=" + " > ".join(path))


if __name__ == "__main__":
    main()
