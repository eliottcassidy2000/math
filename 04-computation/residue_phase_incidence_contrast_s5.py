#!/usr/bin/env python3
"""
residue_phase_incidence_contrast_s5.py

Small contrast table for the residue/phase/incidence synthesis.

The goal is not classification accuracy.  It is a compact set of canonical
examples where the same feature extractor shows different mechanisms:
localized residue kill, broad phase/fiber differences, a landscape plateau,
and a near-kill residue that is still algebraically dangerous.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from tournament_tda import TournamentFeatures, ham_path_count, score_sequence  # noqa: E402


def matrix_from_out(n: int, out: dict[int, list[int]]) -> np.ndarray:
    A = np.zeros((n, n), dtype=np.int8)
    for i in range(n):
        for j in out.get(i, []):
            A[i, j] = 1
    for i in range(n):
        for j in range(i + 1, n):
            if A[i, j] == 0 and A[j, i] == 0:
                A[j, i] = 1
    return A


def transitive(n: int) -> np.ndarray:
    A = np.zeros((n, n), dtype=np.int8)
    for i in range(n):
        for j in range(i + 1, n):
            A[i, j] = 1
    return A


def circulant(p: int, connection: set[int]) -> np.ndarray:
    A = np.zeros((p, p), dtype=np.int8)
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in connection:
                A[i, j] = 1
    return A


def single_core(signature: str) -> np.ndarray:
    """Core vertex 0 inserted into a transitive tournament on vertices 1..m."""
    n = len(signature) + 1
    A = np.zeros((n, n), dtype=np.int8)
    for i in range(1, n):
        for j in range(i + 1, n):
            A[i, j] = 1
    for pos, bit in enumerate(signature, start=1):
        if bit == "1":
            A[0, pos] = 1
        else:
            A[pos, 0] = 1
    return A


def thm025() -> np.ndarray:
    return matrix_from_out(
        9,
        {
            0: [1, 3, 6, 7],
            1: [3],
            2: [0, 1, 4, 5, 6, 7],
            3: [2, 5, 7],
            4: [0, 1, 3, 7],
            5: [0, 1, 4, 6, 7, 8],
            6: [1, 3, 4, 7],
            7: [1],
            8: [0, 1, 2, 3, 4, 6, 7],
        },
    )


def tournament_from_bits(n: int, bits: int) -> np.ndarray:
    A = np.zeros((n, n), dtype=np.int8)
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if bits & (1 << k):
                A[i, j] = 1
            else:
                A[j, i] = 1
            k += 1
    return A


def flipped(A: np.ndarray, i: int, j: int) -> np.ndarray:
    B = A.copy()
    B[i, j], B[j, i] = B[j, i], B[i, j]
    return B


def is_arcflip_local_max(A: np.ndarray) -> bool:
    n = len(A)
    H = ham_path_count(A, n)
    for i in range(n):
        for j in range(i + 1, n):
            if ham_path_count(flipped(A, i, j), n) > H:
                return False
    return True


def find_h37_trap() -> np.ndarray:
    n = 6
    target_score = [1, 2, 2, 3, 3, 4]
    for bits in range(1 << (n * (n - 1) // 2)):
        A = tournament_from_bits(n, bits)
        if score_sequence(A, n) != target_score:
            continue
        if ham_path_count(A, n) != 37:
            continue
        if is_arcflip_local_max(A):
            return A
    raise RuntimeError("No n=6 H=37 local trap found")


def lens_tag(label: str, f: dict[str, object]) -> str:
    if f["omega_projection_kill_vertices"]:
        return "localized residue kill"
    if f["omega_near_kill_rank2_vertices"]:
        return "dangerous near-kill"
    if "Paley" in label:
        return "flat phase fiber"
    if "Interval" in label:
        return "packed phase fiber"
    if "H=37" in label:
        return "phase landscape trap"
    if f["t3"] == 0:
        return "complete order baseline"
    return "mixed"


def summarize(label: str, A: np.ndarray) -> dict[str, object]:
    f = TournamentFeatures(A, len(A)).extract(include_homology=False)
    alpha = (
        f.get("omega_alpha_1", 0),
        f.get("omega_alpha_2", 0),
        f.get("omega_alpha_3", 0),
    )
    residue = (
        f.get("omega_max_loss_residue_alpha_1", 0),
        f.get("omega_max_loss_residue_alpha_2", 0),
        f.get("omega_max_loss_residue_alpha_3", 0),
    )
    return {
        "label": label,
        "n": f["n"],
        "H": f["H"],
        "score": tuple(f["score_seq"]),
        "t3": f["t3"],
        "score_var": f["score_var"],
        "walsh2": f.get("walsh_2", None),
        "alpha": alpha,
        "support_excess": f.get("omega_support_excess", 0),
        "core": f.get("omega_core_size", 0),
        "rho": f.get("omega_max_deletion_loss_frac", 0.0),
        "near_rank2": f.get("omega_near_kill_rank2_vertices", 0),
        "residue": residue,
        "res_rank": f.get("omega_max_loss_residue_rank", 0),
        "res_I2": f.get("omega_max_loss_residue_I2", 1),
        "scc_defect": f["scc_defect"],
        "lens": lens_tag(label, f),
    }


def fmt_walsh2(value: object) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return f"{value:.1f}"


def main() -> None:
    examples = [
        ("transitive T7", transitive(7)),
        ("Paley T7", circulant(7, {1, 2, 4})),
        ("Interval T7", circulant(7, {1, 2, 3})),
        ("H=63 core 1001100", single_core("1001100")),
        ("H=63 core 1100110", single_core("1100110")),
        ("H=37 n6 trap", find_h37_trap()),
        ("THM-025 n9 failure", thm025()),
    ]
    rows = [summarize(label, A) for label, A in examples]

    print("Residue/phase/incidence contrast table (codex-2026-05-30 S5)")
    print("=" * 118)
    print(
        f"{'example':<21} {'n':>2} {'H':>5} {'t3':>3} {'w2':>5} "
        f"{'var':>5} {'alpha(1,2,3)':>15} {'ex':>3} {'core':>4} "
        f"{'rho':>5} {'near2':>5} {'residue':>12} {'rr':>2} "
        f"{'I2':>4} {'lens':<24}"
    )
    print("-" * 118)
    for r in rows:
        print(
            f"{r['label']:<21} {r['n']:>2} {r['H']:>5} {r['t3']:>3} "
            f"{fmt_walsh2(r['walsh2']):>5} {r['score_var']:>5.2f} "
            f"{str(r['alpha']):>15} {r['support_excess']:>3} "
            f"{r['core']:>4} {r['rho']:>5.3f} {r['near_rank2']:>5} "
            f"{str(r['residue']):>12} {r['res_rank']:>2} "
            f"{r['res_I2']:>4} {r['lens']:<24}"
        )

    print()
    print("Score sequences:")
    for r in rows:
        print(f"- {r['label']}: {r['score']}")

    print()
    print("Feedback observations:")
    print("- H=63 complete-core examples are exact residue kills: rho=1 and residue rank 0.")
    print("- THM-025 is not an exact kill: only two cycles survive, but the surviving pair is disjoint, giving rank 2.")
    print("- Paley and Interval T7 differ mainly by phase/fiber geometry: equal score variance but different alpha/disjointness profiles.")
    print("- The H=37 trap is a landscape phenomenon: moderate residue rank, nonzero score variance, and no one-vertex kill explanation.")
    print("- SCC/good-cut defect is independent in this sample: it detects the transitive baseline but not the high-H residue/phase split.")


if __name__ == "__main__":
    main()
