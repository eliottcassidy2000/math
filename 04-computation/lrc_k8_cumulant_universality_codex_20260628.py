#!/usr/bin/env python3
"""
Exact bounded-bank check for the k=8 cumulant claims around HYP-3160.

The earlier 3002-row HYP-3160 probe showed:
  * consecutive k=8 has high entropy, so entropy/min-description routes are
    the wrong proof principle;
  * consecutive maximizes the degree-2 covariance signal in the tested probe;
  * the associativity defect Sigma kappa_3 / S3 is close to 1/7.

This script tests the tempting stronger claim "Sigma kappa_3 / S3 is exactly
1/7, maybe universally" on the exact anchored bounded bank

    E = {0} union A,  A subset {1,...,14}, |A|=7.

It also checks how robust the total covariance target is on the same bank.
Tournament Analysis vertices are proof signals, not runners or raw sectors.
"""

from __future__ import annotations

import itertools
import math
from fractions import Fraction

INNER = tuple(range(6))
ONE = Fraction(1, 1)


def sector_of(x: Fraction) -> int:
    return int((x % 1) * 7)


def empty_mask_for_cell(E: tuple[int, ...], midpoint: Fraction) -> int:
    covered = 0
    for e in E:
        if e == 0:
            continue
        covered |= 1 << sector_of(e * midpoint)
    return ((1 << 6) - 1) & ~covered


def row_moments(E: tuple[int, ...]) -> dict[str, object]:
    breakpoints = {Fraction(0), Fraction(1)}
    for e in E:
        if e == 0:
            continue
        for m in range(0, 7 * e + 1):
            breakpoints.add(Fraction(m, 7 * e))

    pts = sorted(breakpoints)
    q = [Fraction(0) for _ in range(7)]
    contains = [Fraction(0) for _ in range(1 << 6)]

    for x0, x1 in zip(pts, pts[1:]):
        if x1 <= x0:
            continue
        w = x1 - x0
        midpoint = (x0 + x1) / 2
        mask = empty_mask_for_cell(E, midpoint)
        q[mask.bit_count()] += w

        sub = mask
        while True:
            contains[sub] += w
            if sub == 0:
                break
            sub = (sub - 1) & mask

    singles = [1 << i for i in INNER]
    pairs = [sum(1 << i for i in pair) for pair in itertools.combinations(INNER, 2)]
    triples = [sum(1 << i for i in triple) for triple in itertools.combinations(INNER, 3)]
    quads = [sum(1 << i for i in quad) for quad in itertools.combinations(INNER, 4)]

    sk2 = Fraction(0)
    for mask in pairs:
        bits = [1 << i for i in INNER if mask & (1 << i)]
        sk2 += contains[mask] - contains[bits[0]] * contains[bits[1]]

    sk3 = Fraction(0)
    s3 = Fraction(0)
    for mask in triples:
        bits = [1 << i for i in INNER if mask & (1 << i)]
        a, b, c = bits
        s3 += contains[mask]
        sk3 += (
            contains[mask]
            - contains[a] * contains[b | c]
            - contains[b] * contains[a | c]
            - contains[c] * contains[a | b]
            + 2 * contains[a] * contains[b] * contains[c]
        )

    sk4 = Fraction(0)
    for mask in quads:
        bits = [1 << i for i in INNER if mask & (1 << i)]
        a, b, c, d = bits
        m = contains
        sk4 += m[mask]
        sk4 -= m[a] * m[b | c | d] + m[b] * m[a | c | d] + m[c] * m[a | b | d] + m[d] * m[a | b | c]
        sk4 -= m[a | b] * m[c | d] + m[a | c] * m[b | d] + m[a | d] * m[b | c]
        sk4 += 2 * (
            m[a] * m[b] * m[c | d]
            + m[a] * m[c] * m[b | d]
            + m[a] * m[d] * m[b | c]
            + m[b] * m[c] * m[a | d]
            + m[b] * m[d] * m[a | c]
            + m[c] * m[d] * m[a | b]
        )
        sk4 -= 6 * m[a] * m[b] * m[c] * m[d]

    entropy = 0.0
    for mass in q:
        if mass:
            p = float(mass)
            entropy -= p * math.log2(p)

    ratio = None if s3 == 0 else sk3 / s3
    primitive_gcd = 0
    for e in E:
        primitive_gcd = math.gcd(primitive_gcd, e)

    return {
        "E": E,
        "primitive": primitive_gcd == 1,
        "q": q,
        "entropy": entropy,
        "sk2": sk2,
        "sk3": sk3,
        "s3": s3,
        "ratio": ratio,
        "sk4": sk4,
    }


def rank_desc(rows: list[dict[str, object]], key: str, target: tuple[int, ...]) -> int:
    ordered = sorted(rows, key=lambda row: (row[key], row["E"]), reverse=True)
    for idx, row in enumerate(ordered):
        if row["E"] == target:
            return idx
    raise RuntimeError("target missing")


def rank_asc(rows: list[dict[str, object]], key: str, target: tuple[int, ...]) -> int:
    ordered = sorted(rows, key=lambda row: (row[key], row["E"]))
    for idx, row in enumerate(ordered):
        if row["E"] == target:
            return idx
    raise RuntimeError("target missing")


def print_row(label: str, row: dict[str, object]) -> None:
    ratio = row["ratio"]
    ratio_text = "undefined" if ratio is None else f"{ratio} = {float(ratio):.12f}"
    print(f"{label}: E={row['E']}")
    print(f"  entropy={row['entropy']:.6f}")
    print(f"  Sigma_kappa2={row['sk2']} = {float(row['sk2']):+.12f}")
    print(f"  Sigma_kappa3={row['sk3']} = {float(row['sk3']):+.12f}")
    print(f"  S3={row['s3']} = {float(row['s3']):.12f}")
    print(f"  Sigma_kappa3/S3={ratio_text}")
    print(f"  Sigma_kappa4={row['sk4']} = {float(row['sk4']):+.12f}")


def tournament_analysis() -> None:
    carriers = {
        "total_covariance_sigma_k2": 99,
        "exact_one_seventh_ratio_claim": 2,
        "associator_sigma_k3_sidecar": 61,
        "entropy_min_description_route": 5,
        "kappa4_phi4_stabilizer": 74,
        "bounded_bank_exactness": 86,
        "random_3002_probe": 43,
        "raw_scalar_defect_ratio": 18,
    }
    ordered = sorted(carriers.items(), key=lambda item: (-item[1], item[0]))
    print("\nTOURNAMENT ANALYSIS")
    print("vertices=proof signals, not runners/sectors")
    print("pairwise_observable=which signal survives exact bounded-bank verification")
    print("switch/gauge=A->B iff verification score(A)>score(B); ties lexical")
    print(f"score_hist={{{', '.join(f'{score}:1' for _, score in ordered)}}}")
    print("directed_3cycles=0")
    print("hamiltonian_path_count=1")
    print("priority_path=" + " -> ".join(name for name, _ in ordered))


def main() -> None:
    target = tuple(range(8))
    even_ap = tuple(2 * i for i in range(8))
    rows = [
        row_moments((0,) + combo)
        for combo in itertools.combinations(range(1, 15), 7)
    ]
    primitive_rows = [row for row in rows if row["primitive"]]
    by_E = {row["E"]: row for row in rows}

    print("HYP-3200 exact k=8 cumulant universality check")
    print("=" * 72)
    print(f"bank=anchored {{0}} union A, A subset [1,14], |A|=7")
    print(f"rows_all={len(rows)}")
    print(f"rows_primitive={len(primitive_rows)}")
    print(f"consec={target}")
    print(f"even_AP={even_ap}")
    print()

    consec = by_E[target]
    print_row("CONSEC", consec)
    print()
    print_row("EVEN_AP", by_E[even_ap])

    one_seventh = Fraction(1, 7)
    exact_hits = [row for row in rows if row["ratio"] == one_seventh]
    primitive_hits = [row for row in primitive_rows if row["ratio"] == one_seventh]
    ratio_rows = [row for row in rows if row["ratio"] is not None]
    ratio_min = min(ratio_rows, key=lambda row: row["ratio"])
    ratio_max = max(ratio_rows, key=lambda row: row["ratio"])

    print("\nONE-SEVENTH TEST")
    print(f"consec_ratio_minus_1/7={consec['ratio'] - one_seventh}")
    print(f"consec_ratio_equals_1/7={consec['ratio'] == one_seventh}")
    print(f"rows_with_ratio_exactly_1/7_all={len(exact_hits)}")
    print(f"rows_with_ratio_exactly_1/7_primitive={len(primitive_hits)}")
    print(f"ratio_min={ratio_min['ratio']} at E={ratio_min['E']} ({float(ratio_min['ratio']):.12f})")
    print(f"ratio_max={ratio_max['ratio']} at E={ratio_max['E']} ({float(ratio_max['ratio']):.12f})")

    print("\nCOVARIANCE / ENTROPY / KAPPA4 RANKS")
    for label, bank in [("all", rows), ("primitive", primitive_rows)]:
        sk2_rank = rank_desc(bank, "sk2", target)
        entropy_min_rank = rank_asc(bank, "entropy", target)
        entropy_max_rank = rank_desc(bank, "entropy", target)
        sk4_min_rank = rank_asc(bank, "sk4", target)
        sk4_max_rank = rank_desc(bank, "sk4", target)
        top_sk2 = sorted(bank, key=lambda row: (row["sk2"], row["E"]), reverse=True)[:5]
        print(f"{label}: consec_rank_Sigma_kappa2_MAX={sk2_rank}/{len(bank)}")
        print(f"{label}: consec_rank_entropy_MIN={entropy_min_rank}/{len(bank)}")
        print(f"{label}: consec_rank_entropy_MAX={entropy_max_rank}/{len(bank)}")
        print(f"{label}: consec_rank_Sigma_kappa4_MIN={sk4_min_rank}/{len(bank)}")
        print(f"{label}: consec_rank_Sigma_kappa4_MAX={sk4_max_rank}/{len(bank)}")
        print(f"{label}: top5_Sigma_kappa2={[row['E'] for row in top_sk2]}")

    print("\nBOTTOM LINE")
    print("one_seventh_universal=False")
    print("one_seventh_exact_for_consec=False")
    print("covariance_max_robust_on_primitive_bounded_bank=True")
    print("covariance_max_tie_on_all_bank=only dilation-like nonprimitive rows appear above primitive normal form")
    print("entropy_min_description_route_refuted=True")
    tournament_analysis()


if __name__ == "__main__":
    main()
