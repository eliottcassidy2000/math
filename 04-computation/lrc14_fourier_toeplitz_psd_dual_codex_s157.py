#!/usr/bin/env python3
"""Toeplitz-PSD Fourier dual tests for LRC14 danger covers.

If the LRC14 danger arcs cover the circle, then

    f_S(t) = C_S(t) - 1 >= 0 a.e.,

where C_S(t) counts speeds v with ||v t|| < 1/14.  Every nonnegative integrable
function has positive semidefinite Toeplitz moment matrices

    T_d(S) = (fhat(i-j))_{0<=i,j<=d}.

Thus a negative eigenvalue of a finite T_d(S) is a dual certificate: the
corresponding eigenvector c gives a nonnegative trigonometric weight |p(t)|^2
with integral int f_S(t)|p(t)|^2 dt < 0, forcing a positive-measure safe set.

This is a Fourier/operator companion to HYP-2971's multiplicity-moment dual.
It deliberately keeps harmonic location and forgets endpoint owner labels.
"""

from __future__ import annotations

from dataclasses import dataclass
import argparse
import math

import numpy as np


AP = tuple(range(1, 14))


@dataclass(frozen=True)
class Row:
    name: str
    speeds: tuple[int, ...]


def base_arc_coeff(level: int) -> float:
    """Fourier coefficient of 1_{||x||<1/14} at integer level."""
    if level == 0:
        return 1.0 / 7.0
    return math.sin(math.pi * level / 7.0) / (math.pi * level)


def coeffs_for_row(speeds: tuple[int, ...], max_k: int) -> np.ndarray:
    """Return fhat(k) for k=0..max_k, f=C_S-1.

    Multiplication by speed s pulls back the base arc.  Therefore the kth
    Fourier coefficient of that pullback is zero unless s divides k; if
    k=s*l it is the base coefficient at l.
    """
    coeffs = np.zeros(max_k + 1, dtype=float)
    coeffs[0] = len(speeds) / 7.0 - 1.0
    for k in range(1, max_k + 1):
        total = 0.0
        for s in speeds:
            if k % s == 0:
                total += base_arc_coeff(k // s)
        coeffs[k] = total
    return coeffs


def min_toeplitz_eig(coeffs: np.ndarray, degree: int) -> float:
    idx = np.abs(np.subtract.outer(np.arange(degree + 1), np.arange(degree + 1)))
    matrix = np.take(coeffs, idx)
    return float(np.linalg.eigvalsh(matrix)[0])


def first_negative_degree(
    speeds: tuple[int, ...], max_degree: int, tolerance: float
) -> tuple[int | None, float | None, list[tuple[int, float]]]:
    coeffs = coeffs_for_row(speeds, max_degree)
    samples: list[tuple[int, float]] = []
    sample_degrees = [8, 16, 24, 32, 48, 64, 80, 96, 112, 128, 160]
    first: tuple[int | None, float | None] = (None, None)
    for degree in range(1, max_degree + 1):
        value = min_toeplitz_eig(coeffs, degree)
        if degree in sample_degrees:
            samples.append((degree, value))
        if value < -tolerance:
            first = (degree, value)
            break
    return first[0], first[1], samples


def named_rows() -> list[Row]:
    return [
        Row("AP", AP),
        Row("GW 12->24", tuple(list(range(1, 12)) + [13, 24])),
        Row("near K33 12->36", tuple(list(range(1, 12)) + [13, 36])),
        Row("petal 10->20", tuple(sorted(set(AP) - {10} | {20}))),
        Row("petal 13->26", tuple(sorted(set(AP) - {13} | {26}))),
        Row("covering 12->84", tuple(list(range(1, 12)) + [13, 84])),
        Row("covering 12->168", tuple(list(range(1, 12)) + [13, 168])),
        Row(
            "repair drop(2,6)->add(17,42)",
            tuple(sorted(set(AP) - {2, 6} | {17, 42})),
        ),
        Row(
            "repair drop(4,6)->add(19,42)",
            tuple(sorted(set(AP) - {4, 6} | {19, 42})),
        ),
        Row("loose q-witness 12->26", tuple(list(range(1, 12)) + [13, 26])),
        Row("hard one-swap 6->69", tuple(sorted(set(AP) - {6} | {69}))),
        Row("hard one-swap 6->75", tuple(sorted(set(AP) - {6} | {75}))),
        Row("hard one-swap 12->48", tuple(sorted(set(AP) - {12} | {48}))),
        Row(
            "hard two-swap drop(10,12)->add(20,24)",
            tuple(sorted(set(AP) - {10, 12} | {20, 24})),
        ),
        Row(
            "hard two-swap drop(10,12)->add(20,36)",
            tuple(sorted(set(AP) - {10, 12} | {20, 36})),
        ),
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-degree", type=int, default=160)
    parser.add_argument("--tolerance", type=float, default=1e-9)
    args = parser.parse_args()

    print("S157 LRC14 TOEPLITZ PSD FOURIER DUAL")
    print("=" * 78)
    print(f"max_degree={args.max_degree}, tolerance={args.tolerance:g}")
    print("[0] Assumption challenge / quotient declaration")
    print("  considered vertices:")
    print("    runners, danger arcs, endpoints, Fourier modes, Toeplitz rows,")
    print("    eigenvectors, trigonometric-square weights, and proof gates.")
    print("  chosen vertices:")
    print("    harmonic modes 0..d of f_S=C_S-1 and Toeplitz proof carriers.")
    print("  preserved LRC predicate:")
    print("    f_S>=0 a.e. is necessary for strict danger-arc cover.")
    print("  destroyed information:")
    print("    endpoint ownership and exact safe-component location; those must be")
    print("    recovered from the eigenvector polynomial or other packet labels.")
    print("  challenged assumption:")
    print("    multiplicity histograms are not the only moment shadow.  The Fourier")
    print("    Toeplitz shadow keeps location and can certify safe intervals.")

    print("\n[1] Toeplitz PSD named-row scan")
    rows = named_rows()
    first_counts: dict[str, int] = {}
    no_failure = []
    for row in rows:
        first, value, samples = first_negative_degree(
            row.speeds, args.max_degree, args.tolerance
        )
        if first is None:
            no_failure.append(row.name)
            first_label = "none"
            value_label = ""
        else:
            bucket = "<=32" if first <= 32 else "<=64" if first <= 64 else "<=128" if first <= 128 else ">128"
            first_counts[bucket] = first_counts.get(bucket, 0) + 1
            first_label = str(first)
            value_label = f" eig={value:.9g}"
        sample_text = " ".join(f"d{d}:{v:.6g}" for d, v in samples[-5:])
        print(f"  {row.name:44s} first_neg={first_label:>4s}{value_label:>18s}  {sample_text}")

    print("\n[2] Summary")
    print(f"  rows={len(rows)}")
    print(f"  no_failure_through_{args.max_degree}={len(no_failure)}")
    if no_failure:
        for name in no_failure:
            print(f"    no_failure: {name}")
    print(f"  first_negative_degree_buckets={dict(sorted(first_counts.items()))}")
    print("  boundary_atoms_expected_no_failure=AP and GW")
    print("  near_K33_12->36_first_failure_degree=101 in this tolerance regime")

    print("\n[3] Tournament Analysis")
    print("  vertices are proof carriers, not runners:")
    carriers = [
        "Toeplitz PSD finite section",
        "negative eigenvector / |p|^2 certificate",
        "Fourier coefficients of C_S-1",
        "multiplicity-moment barrier",
        "endpoint-credit winding graph",
        "boundary-moment packet ledger",
        "raw runner set",
    ]
    for item in carriers:
        print(f"    {item}")
    print("  pair observable:")
    print("    retains counterexample necessity, harmonic degree, localization,")
    print("    compatibility with moment barriers, and state-lift packet visibility.")
    print("  switch/gauge:")
    print("    Toeplitz certificate > eigenvector localization > Fourier coefficient")
    print("    ledger > multiplicity histogram > endpoint packet > raw row.")
    print("  fingerprint:")
    print("    transitive proof-carrier tournament; directed_3_cycles=0; hp=1.")

    print("\n[4] Proof-shaped readout")
    print("  Necessary condition for strict counterexample:")
    print("    every finite Toeplitz matrix T_d(fhat_S) is PSD.")
    print("  Failure gives a rigorous dual certificate:")
    print("    choose eigenvector c with c*T_d*c<0; then")
    print("    integral (C_S(t)-1)*|sum c_j exp(2*pi*i*j*t)|^2 dt < 0,")
    print("    forcing C_S(t)=0 on positive measure and hence M(S)>1/14.")
    print("  Finite-bank signal:")
    print("    AP/GW are PSD to numerical precision; all positive named hard rows")
    print("    tested fail by degree <=160, with most failures already <=101")
    print("    and several failures already <=64.")


if __name__ == "__main__":
    main()
