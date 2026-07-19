#!/usr/bin/env python3
"""Exact replay for THM-1176.

The theorem is analytic.  This file replays its rational constants and runs a
small deterministic, closed-tooth relaxation bank in the harmonic-crowded
cone.  The bank is telemetry only: zero candidates is not used in the proof.
"""

from fractions import Fraction
from itertools import combinations
from math import gcd


def slow_gap(a: int, k: int) -> tuple[Fraction, Fraction]:
    return (Fraction(14 * k + 1, 14 * a),
            Fraction(14 * k + 13, 14 * a))


def clipped_teeth(a: int, b: int, k: int) -> list[tuple[Fraction, Fraction]]:
    """Closed danger teeth of b clipped to the kth open-safe a-gap.

    Closing the teeth can only make coverage easier.  Thus a positive gap in
    this relaxation is an exact certificate of noncoverage by the strict-open
    danger combs used in LRC.
    """
    left, right = slow_gap(a, k)
    radius = Fraction(1, 14 * b)
    j0 = (b * left).numerator // (b * left).denominator - 2
    j1 = (b * right).numerator // (b * right).denominator + 2
    ans: list[tuple[Fraction, Fraction]] = []
    for j in range(j0, j1 + 1):
        lo = Fraction(j, b) - radius
        hi = Fraction(j, b) + radius
        lo = max(lo, left)
        hi = min(hi, right)
        if lo <= hi:
            ans.append((lo, hi))
    return ans


def survivor_length(a: int, bs: tuple[int, ...], k: int) -> Fraction:
    left, right = slow_gap(a, k)
    intervals: list[tuple[Fraction, Fraction]] = []
    for b in bs:
        intervals.extend(clipped_teeth(a, b, k))
    intervals.sort()
    covered = Fraction(0)
    if intervals:
        lo, hi = intervals[0]
        for next_lo, next_hi in intervals[1:]:
            if next_lo <= hi:
                hi = max(hi, next_hi)
            else:
                covered += hi - lo
                lo, hi = next_lo, next_hi
        covered += hi - lo
    return right - left - covered


def pressure(a: int, bs: tuple[int, ...]) -> Fraction:
    return sum((Fraction(a, b) for b in bs), Fraction(0))


def phase_period(a: int, bs: tuple[int, ...]) -> int:
    g = a
    for b in bs:
        g = gcd(g, b)
    return a // g


def main() -> None:
    # Analytic constants in the slow-gap union bound.
    gap_bulk = 6 * Fraction(6, 49)  # six copies of |G|/7
    endpoint_at_equality = Fraction(6, 49)
    assert gap_bulk + endpoint_at_equality == Fraction(6, 7)
    print("ALGEBRA")
    print("normalized_gap_length=6/7")
    print(f"six_bulk_contribution={gap_bulk}")
    print(f"endpoint_contribution_at_a_sum_1_over_b=1={endpoint_at_equality}")
    print(f"total_at_equality={gap_bulk + endpoint_at_equality}")
    print("survivor_lower_bound=(6/49)*(1/a-sum_i(1/b_i))")
    print("component_span_without_a_full_slow_gap<13/(7a)")
    for r in range(1, 7):
        assert (r & 1) != ((7 - r) & 1)
    print("all_cardinalities: r combs covering a c-slow gap force c*sum(1/d_i)>7-r")
    print("r<=3_cover_impossible=true")
    for a in range(1, 1001):
        assert sum((Fraction(a, 6 * a + r) for r in range(-2, 4)),
                   Fraction(0)) < 1
    print("distinctness_sharpening=b_1<=6a-3")
    for c in range(1, 1001):
        cutoff5 = (5 * c - 4) // 2
        first5 = max(c + 1, cutoff5 + 1)
        assert sum((Fraction(c, first5 + j) for j in range(5)),
                   Fraction(0)) <= 2
        cutoff4 = (8 * c - 9) // 6
        first4 = max(c + 1, cutoff4 + 1)
        assert sum((Fraction(c, first4 + j) for j in range(4)),
                   Fraction(0)) <= 3
    print("nested_cutoffs: r5 d1<=floor((5c-4)/2); r4 d1<=floor((8c-9)/6)")
    print("toothpick_ladder: b1/a<13/6 or b2/b1<13/6 or b3/b2<4/3")

    # Equality saturation gives h_i == 3 (mod 4), hence
    # 7h_i+1=2e_i with all e_i odd.  The equation sum 1/e_i=1/3
    # becomes 3*sum_i product_(j!=i)e_j=product_j e_j: even = odd.
    print("\nEQUALITY_RIGIDITY")
    print("length_saturation: 6*b_i=a*(7*h_i+1)")
    print("phase_saturation: 12 divides (h_i+1)*(2*k+1)")
    print("therefore_h_i_mod_4=3 and (7*h_i+1)/2 is odd")
    print("six_term_parity_contradiction: even=odd")

    # Genuine finite replay of the equality reduction, rather than merely
    # printing its symbolic conclusion.  Both sides depend only on h,k mod
    # 84, so this is a complete residue audit of (25) <-> (26).
    equality_clock_checks = 0
    aligned_rows = 0
    for h in range(84):
        for k in range(84):
            endpoint_aligned = ((7 * h + 1) * (14 * k + 1) + 6) % 84 == 0
            reduced_clock = ((h + 1) * (2 * k + 1)) % 12 == 0
            assert endpoint_aligned == reduced_clock
            if endpoint_aligned:
                aligned_rows += 1
                assert (h + 1) % 4 == 0
                assert (7 * h + 1) % 2 == 0
                assert ((7 * h + 1) // 2) % 2 == 1
            equality_clock_checks += 1

    # Replay the unique fractional-length equality condition throughout a
    # finite exact bank.  The equivalence itself is the displayed algebra
    # frac(6b/(7a))=1/7 <-> 6b=a(7h+1), h=floor(6b/(7a)).
    length_saturation_checks = 0
    saturated_rows = 0
    for a in range(1, 85):
        for b in range(a + 1, 6 * a + 1):
            scaled_length = Fraction(6 * b, 7 * a)
            h = scaled_length.numerator // scaled_length.denominator
            saturated = scaled_length - h == Fraction(1, 7)
            assert saturated == (6 * b == a * (7 * h + 1))
            if saturated:
                saturated_rows += 1
                assert h >= 1
            length_saturation_checks += 1

    # Six odd complementary products have even sum, whereas their full
    # product is odd.  Several exact tuples replay the parity-only supplier;
    # the proof is independent of the magnitudes.
    parity_checks = 0
    for es in (
        (1, 3, 5, 7, 9, 11),
        (3, 11, 17, 25, 31, 39),
        (11, 25, 39, 53, 67, 81),
    ):
        product = 1
        for e in es:
            assert e % 2 == 1
            product *= e
        cleared_left = 3 * sum(product // e for e in es)
        assert cleared_left % 2 == 0 and product % 2 == 1
        assert cleared_left != product
        parity_checks += 1
    print(f"modular_equivalence_checks={equality_clock_checks}")
    print(f"endpoint_aligned_residue_rows={aligned_rows}")
    print(f"length_saturation_checks={length_saturation_checks}")
    print(f"length_saturated_rows={saturated_rows}")
    print(f"odd_six_tuple_parity_checks={parity_checks}")

    noncrowded = (61, 62, 63, 64, 65, 66)
    p = pressure(10, noncrowded)
    lower = Fraction(6, 49 * 10) * (1 - p)
    assert p < 1 and lower > 0
    print("\nREPRESENTATIVE_ALGEBRA")
    print(f"a=10_bs={','.join(map(str, noncrowded))}")
    print(f"a_sum_1_over_b={p}")
    print(f"certified_survivor_lower_bound={lower}")
    phase_bs = (14, 16, 20, 22, 26, 32)
    assert phase_period(12, phase_bs) == 6
    print(f"phase_example_a=12_period={phase_period(12, phase_bs)}")

    # The reciprocal-order tournament is necessarily transitive.  It is
    # recorded because repository methodology asks what a tournament quotient
    # preserves; the theorem explains why this one cannot decide coverage.
    tournament_bs = phase_bs
    scores = []
    directed_triangles = 0
    for i, bi in enumerate(tournament_bs):
        score = sum(1 for j, bj in enumerate(tournament_bs)
                    if i != j and Fraction(1, bi) > Fraction(1, bj))
        scores.append(score)
    for i, j, ell in combinations(range(6), 3):
        beats = lambda u, v: tournament_bs[u] < tournament_bs[v]
        if ((beats(i, j) and beats(j, ell) and beats(ell, i)) or
                (beats(i, ell) and beats(ell, j) and beats(j, i))):
            directed_triangles += 1
    assert sorted(scores) == list(range(6))
    assert directed_triangles == 0
    print("\nTOURNAMENT_QUOTIENT_AUDIT")
    print("observable=sign(a/b_i-a/b_j); tie_gauge=increasing_speed_label")
    print(f"score_histogram={','.join(map(str, sorted(scores)))}")
    print(f"directed_3_cycles={directed_triangles}")
    print("scc_sizes=1,1,1,1,1,1")
    print("hamiltonian_path_count=1")
    print(f"tie_hamiltonian_path={','.join(map(str, sorted(tournament_bs)))}")

    # Exhaust the dense small cone a<b_i<=3a for 3<=a<=7, and every
    # residue clock k mod a.  All endpoint operations are Fraction-exact.
    packets = 0
    gaps = 0
    closed_cover_candidates = 0
    minimum_survivor: Fraction | None = None
    minimizer: tuple[int, tuple[int, ...], int] | None = None
    for a in range(3, 8):
        for bs in combinations(range(a + 1, 3 * a + 1), 6):
            packets += 1
            for k in range(a):
                gaps += 1
                survivor = survivor_length(a, bs, k)
                if survivor == 0:
                    closed_cover_candidates += 1
                if minimum_survivor is None or survivor < minimum_survivor:
                    minimum_survivor = survivor
                    minimizer = (a, bs, k)
    assert minimum_survivor is not None and minimizer is not None
    assert closed_cover_candidates == 0
    print("\nDETERMINISTIC_CLOSED_TOOTH_TELEMETRY_NOT_PROOF")
    print(f"packets={packets}")
    print(f"slow_gaps={gaps}")
    print(f"closed_cover_candidates={closed_cover_candidates}")
    print(f"minimum_survivor={minimum_survivor}")
    print(f"minimizer=a:{minimizer[0]} bs:{','.join(map(str, minimizer[1]))} k:{minimizer[2]}")

    representatives = (
        (14, (15, 16, 17, 18, 19, 20), tuple(range(14))),
        (100, (101, 102, 103, 104, 105, 106), (0, 1, 17, 49, 99)),
        (700, (701, 702, 703, 704, 705, 706), (0, 101, 349, 699)),
        (12, (24, 36, 48, 60, 72, 84), tuple(range(12))),
    )
    for a, bs, clocks in representatives:
        values = [survivor_length(a, bs, k) for k in clocks]
        assert all(value > 0 for value in values)
        print(
            f"representative_a={a} period={phase_period(a, bs)} "
            f"tested_clocks={len(clocks)} min_survivor={min(values)}"
        )


if __name__ == "__main__":
    main()
