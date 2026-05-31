#!/usr/bin/env python3
"""
chirality_perspective_atlas_s370.py

codex-2026-05-31 S370

Apply the S369 chirality/perspective lens beyond the first 56 coincidence.

The working grammar:

  projection count  ->  inherited/core layer + mirror/chiral residue
  boundary layer    ->  outer caps + interior cells
  42                ->  doubled 21-boundary, or the 2*3*7 CRT/Fano core

This script collects exact small ledgers for:

1. Six-vertex tournament classes.
2. S367 LRC missed cells.
3. Paley/Fano T7 directed odd cycles and support excess.
4. Base-42 residue classes from the Erdos-Straus thread.

It is intentionally an atlas: each row is exact, but the proposed bridges are
marked as synthesis, not proof.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations
import importlib.util
import math
import sys


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


bridge = load_module("s369_bridge", "04-computation/lrc_tournament_56_bridge_s369.py")
s367 = load_module("s367_scalar_gauge", "04-computation/lonely_runner_k13_scalar_gauge_s367.py")


def hamiltonian_paths(bits: int, n: int) -> int:
    adj = bridge.adjacency(bits, n)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for v in range(n):
            val = dp[mask][v]
            if not val:
                continue
            for u in adj[v]:
                if not (mask >> u) & 1:
                    dp[mask | (1 << u)][u] += val
    return sum(dp[(1 << n) - 1])


def summarize_values(values: list[int]) -> str:
    return (
        f"count={len(values)} min={min(values)} max={max(values)} "
        f"mean={sum(values) / len(values):.2f} values={sorted(set(values))}"
    )


def tournament_six_atlas() -> None:
    data6 = bridge.class_data(6)
    by_mirror = {
        "self_converse": [item for item in data6 if item.self_converse],
        "chiral": [item for item in data6 if not item.self_converse],
    }
    by_strong = {
        "strong": [item for item in data6 if item.strongly_connected],
        "non_strong": [item for item in data6 if not item.strongly_connected],
    }
    h_by_rep = {item.rep: hamiltonian_paths(item.rep, 6) for item in data6}

    print("Six-vertex tournament atlas")
    print(f"  classes={len(data6)}")
    print(
        "  mirror split="
        f"{len(by_mirror['self_converse'])} self-converse + "
        f"{len(by_mirror['chiral'])} chiral"
    )
    print(
        "  strong split="
        f"{len(by_strong['strong'])} strong + "
        f"{len(by_strong['non_strong'])} non-strong"
    )
    print(f"  score_sequence_count={len({item.score_sequence for item in data6})}")
    print(f"  aut_size_hist={Counter(item.aut_size for item in data6)}")
    print(f"  aut3_count={sum(1 for item in data6 if item.aut_size == 3)}")
    for label, items in by_mirror.items():
        print(f"  H {label}: {summarize_values([h_by_rep[item.rep] for item in items])}")
    for label, items in by_strong.items():
        print(f"  H {label}: {summarize_values([h_by_rep[item.rep] for item in items])}")
    print("  structural readings:")
    print("    56 = 12 self-converse + 44 chiral")
    print("    44 = 2 * 22 chiral converse pairs / score-sequence count")
    print("    42 = 2 * 21 non-strong classes")
    print("    12 = T(5) = phi(42) = self-converse layer at n=6")
    print()


def lrc_atlas() -> None:
    system = s367.build_pattern_system(14)
    vector = (0, 0, 0, 0, 0, 7, 0, 0, 0, 0, 0, 0, 0)
    missed = s367.missed_candidates(system, vector)
    pattern_ids = sorted({pattern_idx for _, pattern_idx in missed})
    widths = Counter(system.patterns[p].hi - system.patterns[p].lo for p in pattern_ids)
    mirror_pairs = []
    used = set()
    for p in pattern_ids:
        if p in used:
            continue
        bins = system.patterns[p].bins
        mirror = tuple((system.n - 1 - value) % system.n for value in bins)
        mate = next(q for q in pattern_ids if system.patterns[q].bins == mirror)
        width = system.patterns[p].hi - system.patterns[p].lo
        mirror_pairs.append((p, mate, width))
        used.add(p)
        used.add(mate)

    outer_width = max(widths)
    outer = [p for p in pattern_ids if system.patterns[p].hi - system.patterns[p].lo == outer_width]
    inner = [p for p in pattern_ids if p not in outer]

    print("S367 LRC missed-cell atlas")
    print(f"  missed_cells={len(missed)}")
    print(f"  shifts={sorted({shift for shift, _ in missed})}")
    print(f"  stencils={len(pattern_ids)} mirror_pairs={len(mirror_pairs)}")
    print(
        "  widths="
        f"{[(s367.fmt_frac(width), count) for width, count in sorted(widths.items())]}"
    )
    print(
        "  denominator_factorization="
        f"{[(width.denominator, prime_factorization(width.denominator)) for width in sorted(widths)]}"
    )
    print(
        "  outer/interior="
        f"{len(outer)}*7={len(outer) * 7}; {len(inner)}*7={len(inner) * 7}"
    )
    print("  structural readings:")
    print("    56 = 7 odd shifts * 8 mirror-stencils")
    print("    14 = 7 shifts * 2 outer caps")
    print("    42 = 7 shifts * 6 interior stencil positions")
    print("    the three interior denominators are all multiples of 42")
    print()


def prime_factorization(n: int) -> str:
    out = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            power = 0
            while n % d == 0:
                n //= d
                power += 1
            out.append(f"{d}^{power}" if power > 1 else str(d))
        d += 1 if d == 2 else 2
    if n > 1:
        out.append(str(n))
    return "*".join(out)


def paley_tournament(p: int) -> list[list[int]]:
    residues = {pow(x, 2, p) for x in range(1, p)}
    adj = [[0] * p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in residues:
                adj[i][j] = 1
    return adj


def directed_cycles_by_length(adj: list[list[int]]) -> tuple[Counter[int], Counter[int]]:
    n = len(adj)
    cycle_counts: Counter[int] = Counter()
    support_counts: Counter[int] = Counter()
    for length in range(3, n + 1, 2):
        for subset in combinations(range(n), length):
            support_has_cycle = False
            first = min(subset)
            rest = [v for v in subset if v != first]
            for perm in permutations(rest):
                cyc = (first,) + perm
                if all(adj[cyc[i]][cyc[(i + 1) % length]] for i in range(length)):
                    cycle_counts[length] += 1
                    support_has_cycle = True
            if support_has_cycle:
                support_counts[length] += 1
    return cycle_counts, support_counts


def paley_fano_atlas() -> None:
    adj = paley_tournament(7)
    cycles, supports = directed_cycles_by_length(adj)
    total_cycles = sum(cycles.values())
    total_supports = sum(supports.values())
    print("Paley/Fano T7 atlas")
    print(f"  directed_odd_cycles={dict(cycles)} total={total_cycles}")
    print(f"  odd_cycle_supports={dict(supports)} total={total_supports}")
    print(f"  support_excess={total_cycles - total_supports}")
    print("  structural readings:")
    print("    80 = 14 triangles + 42 pentagons + 24 heptagons")
    print("    36 supports = 14 triangle supports + 21 pentagon supports + 1 full support")
    print("    44 = 80 - 36 support multiplicity excess")
    print("    42 = 2 directed pentagons on each of 21 five-subsets")
    print()


def base42_atlas() -> None:
    units = [r for r in range(42) if math.gcd(r, 42) == 1]
    hard = [r for r in units if r % 12 == 1]
    easy = [r for r in units if r not in hard]
    print("Base-42 residue atlas")
    print(f"  phi(42)={len(units)} units={units}")
    print(f"  hard_p_eq_1_mod_12={hard} mod7={[r % 7 for r in hard]}")
    print(f"  easy_units={len(easy)}")
    print("  structural readings:")
    print("    12 = phi(42) = T(5) = self-converse T(6) count")
    print("    8 easy unit classes echo the 8 LRC stencils / 5-to-6 perspective gap")
    print("    4 hard classes echo the 4 LRC mirror-stencil pairs")
    print()


def synthesis() -> None:
    print("Synthesis: applying the perspective")
    print("  A. 12 is the inherited symmetric core:")
    print("     T(5), phi(42), self-converse T(6), and order-3 aut count at T(6).")
    print("  B. 44 is the chiral/support residue:")
    print("     T(6) chiral classes and Paley T7 odd-cycle support excess.")
    print("  C. 42 is the doubled boundary/interior:")
    print("     LRC interior cells, Paley pentagons, and 2*non-strong T(6).")
    print("  D. 8 is the projection failure:")
    print("     LRC stencils, 5-to-6 perspective gap, and easy base-42 unit classes.")
    print("  Proposed proof tactic for HYP-1823:")
    print("     split the normalized vector problem into outer caps (14),")
    print("     interior 42, and mirror/chiral pairs before attempting any full SAT.")


def main() -> None:
    print("Chirality / perspective atlas (S370)")
    print()
    tournament_six_atlas()
    lrc_atlas()
    paley_fano_atlas()
    base42_atlas()
    synthesis()


if __name__ == "__main__":
    main()
