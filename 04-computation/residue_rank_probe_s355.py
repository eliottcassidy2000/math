#!/usr/bin/env python3
"""
residue_rank_probe_s355.py

codex-2026-05-30

Small feedback-loop probe for HYP-1780/HYP-1785.  It compares familiar
examples by the new `omega_*` deletion-residue rank features in
`tournament_tda.py`.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from tournament_tda import TournamentFeatures  # noqa: E402


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


def summarize(label: str, A: np.ndarray) -> dict[str, object]:
    f = TournamentFeatures(A, len(A)).extract()
    return {
        "label": label,
        "n": f["n"],
        "H": f["H"],
        "scc": f["scc_count"],
        "scc_defect": f["scc_defect"],
        "a1": f["omega_alpha_1"],
        "a2": f["omega_alpha_2"],
        "core": f["omega_core_size"],
        "excess": f["omega_support_excess"],
        "rho": f["omega_max_deletion_loss_frac"],
        "min_keep": f["omega_min_deletion_residue_cycles"],
        "near": f["omega_near_kill_vertices"],
        "near2": f["omega_near_kill_rank2_vertices"],
        "rv": f["omega_max_loss_vertex"],
        "ra1": f["omega_max_loss_residue_alpha_1"],
        "ra2": f["omega_max_loss_residue_alpha_2"],
        "rrank": f["omega_max_loss_residue_rank"],
        "rI2": f["omega_max_loss_residue_I2"],
    }


def main() -> None:
    examples = [
        ("transitive T7", transitive(7)),
        ("Paley T7", circulant(7, {1, 2, 4})),
        ("Interval T7", circulant(7, {1, 2, 3})),
        ("H=63 core 1001100", single_core("1001100")),
        ("H=63 core 1100110", single_core("1100110")),
        ("THM-025 n9 failure", thm025()),
    ]
    rows = [summarize(label, A) for label, A in examples]

    print("Residue-rank feature probe (codex-2026-05-30 S355)")
    print("=" * 88)
    print(
        f"{'example':<22} {'n':>2} {'H':>5} {'SCC':>3} {'def':>3} "
        f"{'a1':>4} {'a2':>4} {'core':>4} {'ex':>4} {'rho':>5} "
        f"{'minK':>5} {'near':>4} {'near2':>5} {'rv':>3} "
        f"{'res(a1,a2)':>12} {'rank':>4} {'I2':>4}"
    )
    print("-" * 88)
    for r in rows:
        print(
            f"{r['label']:<22} {r['n']:>2} {r['H']:>5} {r['scc']:>3} "
            f"{r['scc_defect']:>3} {r['a1']:>4} {r['a2']:>4} "
            f"{r['core']:>4} {r['excess']:>4} {r['rho']:>5.3f} "
            f"{r['min_keep']:>5} {r['near']:>4} {r['near2']:>5} "
            f"{r['rv']:>3} {str((r['ra1'], r['ra2'])):>12} "
            f"{r['rrank']:>4} {r['rI2']:>4}"
        )

    print()
    print("Feedback observations:")
    print("- H=63 examples are exact kills: max-loss deletion residue has rank 0.")
    print("- THM-025 is a near-kill: only two cycles survive, but they are disjoint, so rank 2 and I2=9.")
    print("- Paley/Interval have broad residues instead of near-kills; their contrast is fiber/disjointness, not small-residue rank.")
    print("- SCC defect is orthogonal here: all nontransitive examples are strongly connected, while transitive T7 has defect 0.")


if __name__ == "__main__":
    main()
