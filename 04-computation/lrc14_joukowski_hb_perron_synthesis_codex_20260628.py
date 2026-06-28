#!/usr/bin/env python3
"""Joukowski/Hermite-Biehler and Perron-Frobenius synthesis scout.

This is a small proof-carrier scout for the LRC14 k=8 node.  It does not try
to re-enumerate the bounded bank.  Instead it turns the incoming Perron,
Toeplitz, Joukowski, and Hermite-Biehler ideas into exact local certificates
and quotient guardrails that can be attached to HYP-3200/HYP-3202.

Tournament Analysis vertices are proof obligations and certificates, not
runners or raw arcs.
"""

from __future__ import annotations

import math
from collections import Counter
from fractions import Fraction


D1 = Fraction(308509, 1080450)
D2 = Fraction(547577, 2160900)
D3 = Fraction(225577, 1234800)


def fmt(frac: Fraction) -> str:
    return f"{frac} = {float(frac):+.12f}"


def poly_eval(coeffs: tuple[int, ...], x: float) -> float:
    out = 0.0
    for c in coeffs:
        out = out * x + c
    return out


def hb_leg_certificate() -> None:
    """Print the exact k=8 toy Hermite-Biehler interlacing certificate."""
    print("HERMITE-BIEHLER LEG CERTIFICATE")
    print("=" * 72)
    print("even raw fold: v^2 - 5v + 4 has roots 1,4")
    print("negative-axis even leg: E(x)=x^2+5x+4 has roots -4,-1")
    print("odd Eulerian/Worpitzky leg: O(x)=A_3(x)=x^2+4x+1")
    odd_roots = (-2 - math.sqrt(3), -2 + math.sqrt(3))
    print(f"odd roots=(-2-sqrt(3), -2+sqrt(3)) ~= ({odd_roots[0]:+.12f}, {odd_roots[1]:+.12f})")
    print("strict_interlacing: -4 < -2-sqrt(3) < -1 < -2+sqrt(3)")
    print("Wronskian orientation:")
    print("  W=E*O'-E'*O=x^2+6x+11=(x+3)^2+2 > 0 on R")
    print("readout=the minimal k=8 even biquadratic leg and odd A_3 leg form a strict HB pair")
    print("guardrail=this is a leg-level certificate; full miss-PGF still owes the self-inversive defect")
    print()


def joukowski_certificate() -> None:
    print("JOUKOWSKI / DE MOIVRE CERTIFICATE")
    print("=" * 72)
    print("cyclotomic real-axis polynomial: c(x)=x^3+x^2-2x-1")
    roots = [2 * math.cos(2 * math.pi * j / 7) for j in (1, 2, 3)]
    roots = sorted(roots)
    residuals = [abs(poly_eval((1, 1, -2, -1), r)) for r in roots]
    print("de_moivre_roots_2cos(2pi*j/7)=" + ", ".join(f"{r:+.12f}" for r in roots))
    print(f"max_float_residual_in_cubic={max(residuals):.3e}")
    print("readout=the Joukowski image of the 7th-root ideal is the real cubic axis")
    print("guardrail=off-circle Im(w) and q_t R^t - q_{6-t} R^{6-t} remain named sidecars")
    print()


def circulant_pf_eigenvalues(w1: Fraction, w2: Fraction, w3: Fraction) -> list[Fraction]:
    """Eigenvalues of first row [0,w1,w2,w3,w2,w1] on C6."""
    return [
        2 * w1 + 2 * w2 + w3,
        w1 - w2 - w3,
        -w1 - w2 + w3,
        -2 * w1 + 2 * w2 - w3,
        -w1 - w2 + w3,
        w1 - w2 - w3,
    ]


def antiferro_eigenvalues(w1: Fraction, w2: Fraction, w3: Fraction) -> list[Fraction]:
    """Eigenvalues of signed first row [0,-w1,w2,-w3,w2,-w1]."""
    return [
        -2 * w1 + 2 * w2 - w3,
        -w1 - w2 + w3,
        w1 - w2 - w3,
        2 * w1 + 2 * w2 + w3,
        w1 - w2 - w3,
        -w1 - w2 + w3,
    ]


def pf_alignment_certificate() -> None:
    print("PERRON-FROBENIUS DISTANCE-QUOTIENT CERTIFICATE")
    print("=" * 72)
    print("input=HYP-3202 consecutive distance-layer sums D1,D2,D3")
    print(f"D1={fmt(D1)}")
    print(f"D2={fmt(D2)}")
    print(f"D3={fmt(D3)}")
    print(f"Sigma_kappa2={fmt(D1 + D2 + D3)}")
    print()

    # On an ideal C6 distance quotient there are 6 d1 edges, 6 d2 edges,
    # and 3 antipodal d3 edges.
    w1 = D1 / 6
    w2 = D2 / 6
    w3 = D3 / 3
    eig = circulant_pf_eigenvalues(w1, w2, w3)
    rayleigh_ones = (D1 + D2 + D3) / 3
    top = max(eig)
    print("ideal_C6_layer_weights:")
    print(f"  w1=D1/6={fmt(w1)}")
    print(f"  w2=D2/6={fmt(w2)}")
    print(f"  w3=D3/3={fmt(w3)}")
    print("ideal_nonnegative_circulant_eigenvalues:")
    for i, val in enumerate(eig):
        print(f"  lambda_{i}={fmt(val)}")
    print(f"ones_Rayleigh=(1^T C 1)/6={fmt(rayleigh_ones)}")
    print(f"pf_alignment_defect=lambda_max-ones_Rayleigh={fmt(top - rayleigh_ones)}")
    print("row_sum_variance=0 on the quotient")
    print("readout=after distance-layer quotienting, all-ones is exactly the Perron vector")
    print("guardrail=the actual empty-sector covariance matrix has boundary/apex deviation; this quotient cannot replace that sidecar")
    print()

    afm = antiferro_eigenvalues(w1, w2, w3)
    afm_top = max(range(len(afm)), key=lambda i: afm[i])
    print("AFM SIGN-CONTRAST TOY")
    print("signed_first_row=[0,-w1,+w2,-w3,+w2,-w1]")
    for i, val in enumerate(afm):
        print(f"  signed_lambda_{i}={fmt(val)}")
    print(f"signed_top_mode=k{afm_top}")
    print(f"signed_ones_mode_lambda_0={fmt(afm[0])}")
    print("readout=coherence/sign pattern, not raw layer magnitudes alone, decides whether ones is Perron")
    print()


def quotient_guardrails() -> None:
    print("COMPRESSION / INFORMATION-THEORY GUARDRAILS")
    print("=" * 72)
    rows = [
        (
            "radius_only",
            "q0=q6*R^6",
            "root angles, Im(w), self-inversive defect",
            "illegal unless the Lee-Yang circle defect is zero or retained",
        ),
        (
            "total_covariance_only",
            "Sigma_kappa2",
            "D1/D2/D3 layer profile and Perron alignment",
            "too lossy for the HYP-3202 proof route",
        ),
        (
            "commutative_pair_mass",
            "Pascal/binomial cap",
            "odd Worpitzky associator and bracket/order sidecar",
            "captures cap but not dip",
        ),
        (
            "positive_association",
            "nonnegative pair covariances",
            "location of covariance mass and finite trap identity",
            "FKG alone misses 19 primitive nonnegative decoys",
        ),
        (
            "hb_leg_interlacing",
            "E/O real-rooted interlacing",
            "full miss-PGF self-inversive defect",
            "promising only with Joukowski defect sidecar",
        ),
    ]
    for name, preserves, destroys, verdict in rows:
        print(f"{name}: preserves={preserves}; destroys={destroys}; verdict={verdict}")
    print()


def directed_3cycles(adj: list[list[bool]]) -> int:
    n = len(adj)
    total = 0
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                if adj[i][j] and adj[j][k] and adj[k][i]:
                    total += 1
                if adj[i][k] and adj[k][j] and adj[j][i]:
                    total += 1
    return total


def hamiltonian_path_count(adj: list[list[bool]]) -> int:
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last in range(n):
            value = dp[mask][last]
            if not value:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += value
    return sum(dp[-1])


def tournament_analysis() -> None:
    print("TOURNAMENT ANALYSIS")
    print("=" * 72)
    vertices = [
        ("hermite_biehler_interlacing_certificate", 92),
        ("perron_alignment_certificate", 90),
        ("toeplitz_caratheodory_margin", 84),
        ("distance_layer_pf_quotient", 79),
        ("self_inversive_defect_sidecar", 73),
        ("random_current_coupling_order", 65),
        ("law_defect_entropy_meter", 59),
        ("raw_total_covariance_scalar", 41),
        ("plain_positive_association", 29),
        ("row_entropy_scalar", 3),
    ]
    n = len(vertices)
    adj = [[False] * n for _ in range(n)]
    for i, (_, si) in enumerate(vertices):
        for j, (_, sj) in enumerate(vertices):
            if i == j:
                continue
            adj[i][j] = (si, vertices[j][0]) > (sj, vertices[i][0])
    scores = [sum(1 for j in range(n) if adj[i][j]) for i in range(n)]
    path = [name for name, _ in sorted(vertices, key=lambda item: (-item[1], item[0]))]
    print("vertices=proof certificates/obligations, not runners, raw arcs, or sectors")
    print("pairwise_observable=which carrier preserves more LRC14 proof payload toward the k=8 node")
    print("switch/gauge=A beats B iff route_score(A)>route_score(B); ties lexical")
    print(f"score_hist={dict(sorted(Counter(scores).items()))}")
    print(f"directed_3cycles={directed_3cycles(adj)}")
    print("scc_sizes=[1,1,1,1,1,1,1,1,1,1]")
    print(f"hamiltonian_path_count={hamiltonian_path_count(adj)}")
    print("priority_path=" + " -> ".join(path))
    print("challenged_assumption=the best vertices are not runners; they are certificates plus sidecars")
    print()


def main() -> None:
    print("HYP-3222 Joukowski/Hermite-Biehler + Perron-Frobenius synthesis")
    print("=" * 72)
    print("status=exact local certificates and quotient guardrails; not an LRC14 proof")
    print()
    hb_leg_certificate()
    joukowski_certificate()
    pf_alignment_certificate()
    quotient_guardrails()
    tournament_analysis()
    print("NEXT PROOF OBLIGATIONS")
    print("=" * 72)
    print("1. Lift the exact E/O interlacing certificate through the Joukowski map while measuring self-inversive defect.")
    print("2. Replace the ideal C6 PF quotient by a boundary-aware covariance matrix inequality.")
    print("3. Join Toeplitz lambda_min, Perron alignment, distance layers, and random-current order into one spectral packet.")


if __name__ == "__main__":
    main()
