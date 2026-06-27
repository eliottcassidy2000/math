#!/usr/bin/env python3
"""Exact scout for the paper Conjecture 7.1 bridge in LRC(14).

This script separates two claims that were easy to conflate:

1. The paper's literal Conjecture 7.1 asks for an absolute denominator d
   working for every non-tight primitive tuple.
2. The project witness route needs a uniform arc in the normalized
   apex/ruler coordinate, after the large speed has been peeled.

The divisor-loaded covering rows from THM-566 refute (1) as stated, but the
same rows explain why (2) is the viable repair.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import gcd
from collections import Counter


THRESHOLD = Fraction(1, 14)


def lcm(a: int, b: int) -> int:
    return a // gcd(a, b) * b


def lcm_upto(n: int) -> int:
    out = 1
    for v in range(1, n + 1):
        out = lcm(out, v)
    return out


def frac(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def norm(x: Fraction) -> Fraction:
    r = frac(x)
    return min(r, 1 - r)


def min_distance(speeds: tuple[int, ...], t: Fraction) -> Fraction:
    return min(norm(Fraction(v) * t) for v in speeds)


def is_witness(speeds: tuple[int, ...], t: Fraction, strict: bool = False) -> bool:
    md = min_distance(speeds, t)
    return md > THRESHOLD if strict else md >= THRESHOLD


def has_grid_witness(speeds: tuple[int, ...], d: int) -> bool:
    for a in range(d):
        if is_witness(speeds, Fraction(a, d)):
            return True
    return False


def covering_row(B: int) -> tuple[int, tuple[int, ...]]:
    n = 84 * lcm_upto(B)
    return n, tuple(list(range(1, 12)) + [13, n])


def explicit_non_tight_witness(n: int) -> Fraction:
    # Since n is a multiple of 84, n/12 is integral and n*(1/(2n)) = 1/2.
    return Fraction(1, 12) + Fraction(1, 2 * n)


def direct_safe_components(speeds: tuple[int, ...]) -> list[tuple[Fraction, Fraction]]:
    points = {Fraction(0), Fraction(1)}
    for s in speeds:
        for m in range(s):
            points.add(Fraction(14 * m + 1, 14 * s))
            points.add(Fraction(14 * m + 13, 14 * s))
    xs = sorted(points)
    comps: list[tuple[Fraction, Fraction]] = []
    for a, b in zip(xs, xs[1:]):
        if a == b:
            continue
        mid = (a + b) / 2
        if is_witness(speeds, mid, strict=True):
            comps.append((a, b))
    return comps


def strict_grid_bound(length: Fraction) -> int | None:
    if length <= 0:
        return None
    return length.denominator // length.numerator + 1


def fmt_frac(x: Fraction) -> str:
    return f"{x.numerator}/{x.denominator}"


def fmt_float(x: Fraction) -> str:
    return f"{float(x):.12g}"


@dataclass(frozen=True)
class TournamentVertex:
    name: str
    preserves_lrc: int
    survives_divisor_loading: int
    reduces_paper_bottleneck: int
    proof_scope: int
    denominator_assumption_safe: int

    def vector(self) -> tuple[int, int, int, int, int, str]:
        return (
            self.preserves_lrc,
            self.survives_divisor_loading,
            self.reduces_paper_bottleneck,
            self.proof_scope,
            self.denominator_assumption_safe,
            self.name,
        )


def tournament() -> None:
    vertices = [
        TournamentVertex("normalized_ruler_arc_floor", 1, 1, 1, 4, 1),
        TournamentVertex("level7_c7_lift_sieve_THM573", 1, 1, 1, 3, 1),
        TournamentVertex("two_tier_CRT_I_k_p_1_sieve", 1, 1, 1, 2, 1),
        TournamentVertex("THM566_divisor_loaded_obstruction", 1, 1, 0, 2, 1),
        TournamentVertex("direct_t_largest_arc_floor", 1, 0, 0, 1, 0),
        TournamentVertex("raw_Conjecture_7_1_denominator_floor", 0, 0, 0, 1, 0),
        TournamentVertex("prime_field_Prop4_1_shortcut", 0, 0, 1, 1, 1),
    ]
    scores = Counter({v.name: 0 for v in vertices})
    flips = []
    for i, a in enumerate(vertices):
        for b in vertices[i + 1 :]:
            winner = a if a.vector() > b.vector() else b
            loser = b if winner is a else a
            scores[winner.name] += 1
            if winner is b:
                flips.append((b.name, a.name))

    ordered = sorted(vertices, key=lambda v: (scores[v.name], v.vector()), reverse=True)
    hist = Counter(scores.values())
    print("TOURNAMENT_ANALYSIS")
    print("  vertices_are: proof obligations / coordinate carriers, not runners")
    print("  pairwise_observable: lexicographic retention of LRC predicate, survival under divisor loading, removal of paper I(k,p,1)/c=14 bottleneck, proof scope, denominator-assumption safety")
    print("  score_histogram:", dict(sorted(hist.items())))
    print("  hamiltonian_path:", " > ".join(v.name for v in ordered))
    print("  edge_flips_against_input_order:", len(flips))
    if flips:
        print("  first_flips:", "; ".join(f"{a}>{b}" for a, b in flips[:5]))


def component_scout() -> None:
    rows = [
        ("tight_AP_1_to_13", tuple(range(1, 14))),
        ("base_divisor_row_N84", tuple(list(range(1, 12)) + [13, 84])),
        ("first_large_THM566_row_B6", covering_row(6)[1]),
    ]
    print("DIRECT_COMPONENT_SCOUT")
    for name, speeds in rows:
        comps = direct_safe_components(speeds)
        measure = sum((b - a) for a, b in comps)
        largest = max((b - a for a, b in comps), default=Fraction(0))
        print(
            f"  {name}: components={len(comps)} measure={fmt_frac(measure)} "
            f"largest={fmt_frac(largest)} largest_float={fmt_float(largest)} "
            f"strict_grid_D={strict_grid_bound(largest)} mult7_count={sum(v % 7 == 0 for v in speeds)}"
        )
    print("  readout: direct t-components shrink with the loaded apex; this is the failure that points to the normalized ruler coordinate.")


def refutation_scout() -> None:
    print("CONJECTURE_7_1_LITERAL_REFUTATION_SCOUT")
    print("  family: S_B = {1,2,...,11,13,84*lcm(1..B)}")
    print("  obstruction: every d<=B divides the last speed, so that runner is at 0 for every a/d")
    print("  non_tight_witness: t = 1/12 + 1/(2N), valid once N>462; for B>=6, N>=5040")
    for B in (6, 10, 14, 26, 41, 67):
        n, speeds = covering_row(B)
        t = explicit_non_tight_witness(n)
        md = min_distance(speeds, t)
        margin = md - THRESHOLD
        grid_hits = [d for d in range(1, B + 1) if has_grid_witness(speeds, d)]
        apex_upper = Fraction(6, 7 * n)
        print(
            f"  B={B}: N_digits={len(str(n))} gcd={gcd_all(speeds)} "
            f"mult7_count={sum(v % 7 == 0 for v in speeds)} "
            f"explicit_margin={fmt_frac(margin)} no_grid_hits_d<=B={not grid_hits} "
            f"direct_component_upper_by_apex={fmt_frac(apex_upper)} "
            f"direct_strict_D_at_least={strict_grid_bound(apex_upper)}"
        )
    print("  conclusion: for any proposed universal D, choose B>=max(D,6) and d=max(D,6).")
    print("  Then S_B is primitive and non-tight, but has no witness in (1/d)Z.")


def gcd_all(values: tuple[int, ...]) -> int:
    out = 0
    for v in values:
        out = gcd(out, v)
    return out


def assumption_challenge() -> None:
    print("ASSUMPTION_CHALLENGE")
    print("  considered_vertices: runners, residues mod 14, residues mod 7, denominators, direct t-components, normalized ruler arcs, danger endpoints, lift factors, proof obligations")
    print("  selected_vertices: proof obligations and coordinate carriers")
    print("  preserved_predicate: existence of level-1/14 witness after legal lifting/normalization")
    print("  destroyed_by_raw_Conj71: apex scale; absolute denominator d is sensitive to divisor loading")
    print("  repaired_obligation: prove a uniform normalized arc/component-count floor in the THM-573 residual, then translate through the finite-ruler/equidistribution gate")


def main() -> None:
    refutation_scout()
    print()
    component_scout()
    print()
    tournament()
    print()
    assumption_challenge()


if __name__ == "__main__":
    main()
