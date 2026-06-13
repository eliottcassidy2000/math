#!/usr/bin/env python3
"""
erdos_moser_support_gate_codex_p5.py

Codex 2026-06-11 P5: a small support-gate atlas for the Erdos-Moser
Mersenne tower after THM-483 (the zigzag law).

The point is to separate two things that were easy to conflate:

  * the scalar recurrence a_{k+2} = 2 a_k + 1 suggested by the exact
    tower data, and
  * the support realization needed to produce the extra +1.

Pure maximum chains lose the +1 under two skew doublings:

    trans(D(D(TT_t))) = 2t.

But small full tower cores keep it:

    trans(D(D(T_m))) = 2 trans(T_m) + 1,  m = 3,7,15,31.

For the open T511 frontier, the recent t511 witness-lift computation found
only the chain-only value 30 from a maximum TT15 chain in T127, so the next
proof target is a marked support packet, not a bigger brute-force search.
"""

from __future__ import annotations

import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

from erdos_moser_trans_tower_kps2 import (  # noqa: E402
    TransSolver,
    build_tower,
    masks_of,
    trans_of,
    transitive_T,
)
from skew_doubling_core_kps1 import D_skew, H_count  # noqa: E402

OUT = ROOT / "05-knowledge" / "results" / "erdos_moser_support_gate_codex_p5.out"


def dd(A: np.ndarray) -> np.ndarray:
    """Two skew doublings."""
    return D_skew(D_skew(A)[0])[0]


def induced(A: np.ndarray, verts: list[int]) -> np.ndarray:
    return A[np.ix_(verts, verts)]


def max_chain(A: np.ndarray) -> tuple[int, list[int]]:
    """Return trans(A) and one source-ordered maximum transitive chain."""
    masks, n = masks_of(A)
    solver = TransSolver(masks)
    S = (1 << n) - 1
    target = solver.trans(S)
    need = target
    chain: list[int] = []
    while need:
        items = []
        T = S
        while T:
            b = T & -T
            u = b.bit_length() - 1
            T ^= b
            items.append(((masks[u] & S).bit_count(), u))
        items.sort(reverse=True)
        for _, u in items:
            nxt = masks[u] & S
            if 1 + solver.trans(nxt) == need:
                chain.append(u)
                S = nxt
                need -= 1
                break
        else:
            raise RuntimeError("failed to reconstruct maximum chain")
    return target, chain


def route_tournament() -> tuple[np.ndarray, list[tuple[str, tuple[int, int, int, int]]]]:
    """Tournament Analysis over proof routes, not over tower vertices."""
    routes = [
        ("marked_support_gate", (5, 5, 4, 4)),
        ("two_step_recurrence", (3, 5, 3, 5)),
        ("chain_only_lift", (1, 2, 5, 5)),
        ("border_twin_accounting", (3, 3, 4, 4)),
        ("direct_T511_sat", (1, 1, 2, 1)),
        ("random_baseline", (0, 1, 2, 3)),
    ]
    n = len(routes)
    A = np.zeros((n, n), dtype=np.int64)
    for i, (ri, vi) in enumerate(routes):
        for j, (rj, vj) in enumerate(routes):
            if i == j:
                continue
            if vi > vj or (vi == vj and ri < rj):
                A[i, j] = 1
    return A, routes


def ham_path_vertices(A: np.ndarray) -> list[int]:
    """Insertion Hamiltonian path, used only as a tie-path fingerprint."""
    path = [0]
    for v in range(1, A.shape[0]):
        k = 0
        while k < len(path) and A[path[k], v]:
            k += 1
        path.insert(k, v)
    return path


def recurrence_table() -> list[tuple[int, int, str, int, str]]:
    exact = {
        2: ("2", "anchor T3"),
        3: ("3", "anchor T7"),
        4: ("5", "stored exact T15"),
        5: ("7", "stored exact T31"),
        6: ("11", "stored exact T63"),
        7: ("15", "tower_step_structure_kps2.out"),
        8: ("23", "tower_t255_probe_kps2.out"),
        9: ("[30,47]", "t511_witness/direct capped; prediction open"),
    }
    pred: dict[int, int] = {2: 2, 3: 3}
    for k in range(4, 10):
        pred[k] = 2 * pred[k - 2] + 1
    rows = []
    for k in range(2, 10):
        n = (1 << k) - 1
        known, source = exact[k]
        rows.append((k, n, known, pred[k], source))
    return rows


def main() -> None:
    lines: list[str] = []
    lines.append("=== erdos_moser_support_gate_codex_p5 ===")
    lines.append("Question: where does the two-step +1 in the Mersenne tower live?")
    lines.append("")

    lines.append("--- 1. Exact tower table after recent agents ---")
    lines.append("k  n=2^k-1  known trans(T_n)  two-step prediction  source")
    for k, n, known, pred, source in recurrence_table():
        marker = "OK" if known == str(pred) else "OPEN"
        lines.append(f"{k:1d}  {n:7d}  {known:>15}  {pred:>19}  {marker}: {source}")
    lines.append("")
    lines.append(
        "Progress: HYP-2413's old next datum T127 is now known. "
        "The live frontier is T511: prove or refute the predicted value 31."
    )
    lines.append("")

    lines.append("--- 2. Support gate: full core vs chain-only two-step lift ---")
    lines.append(
        "m   trans(T_m)  trans(D(D(T_m)))  full_bonus  "
        "max_chain_size  trans(D(D(chain)))  chain_bonus  seconds"
    )
    cores = build_tower()
    for m in [3, 7, 15, 31]:
        A = cores[m]
        t, chain = max_chain(A)
        tic = time.perf_counter()
        full_dd = dd(A)
        t_full_dd, sv_full = trans_of(full_dd)
        chain_dd = dd(induced(A, chain))
        t_chain_dd, sv_chain = trans_of(chain_dd)
        sec = time.perf_counter() - tic
        lines.append(
            f"{m:2d}  {t:10d}  {t_full_dd:17d}  {t_full_dd - 2*t:10d}  "
            f"{len(chain):14d}  {t_chain_dd:18d}  {t_chain_dd - 2*t:11d}  "
            f"{sec:7.2f}  nodes(full={sv_full.nodes}, chain={sv_chain.nodes})"
        )
    lines.append("")

    lines.append("--- 3. Pure-chain control ---")
    lines.append("t  trans(D(D(TT_t)))  defect")
    for t in [2, 3, 5, 7, 11, 15]:
        A = transitive_T(t)
        t_dd, _ = trans_of(dd(A))
        lines.append(f"{t:2d}  {t_dd:18d}  {t_dd - 2*t:6d}")
    lines.append("")

    lines.append("--- 4. T511 reading ---")
    lines.append(
        "The stored T511 witness lift used a maximum TT15 chain C in T127 and found "
        "trans(D(D(C))) = 30. Twelve one-extra packets and eight two-extra packets "
        "with trans(X)=15 also stayed at 30."
    )
    lines.append(
        "Therefore the predicted T511=31 cannot be proved by a naked maximum chain. "
        "One must find a marked T127 support packet X with trans(X)=15 and "
        "trans(D(D(X)))>=31, or prove no such packet exists and look for a nonlocal witness."
    )
    lines.append("")

    lines.append("--- 5. Tournament Analysis of proof routes ---")
    lines.append(
        "Vertices: proof routes/support packets, not runners. Observable vector = "
        "(support locality, proof leverage, computational tractability, reuse). "
        "Switch: lexicographic domination. Tie Hamiltonian path: insertion path."
    )
    A, routes = route_tournament()
    scores = [int(x) for x in A.sum(axis=1)]
    tri = 0
    for i in range(A.shape[0]):
        for j in range(A.shape[0]):
            for k in range(A.shape[0]):
                if i < j < k:
                    if A[i, j] and A[j, k] and A[k, i]:
                        tri += 1
                    if A[i, k] and A[k, j] and A[j, i]:
                        tri += 1
    path = ham_path_vertices(A)
    lines.append(f"score histogram: {dict(sorted(Counter(scores).items()))}")
    lines.append(f"directed triangles: {tri}")
    lines.append(f"Hamiltonian path count H: {H_count(A)}")
    lines.append("tie path: " + " -> ".join(routes[i][0] for i in path))
    lines.append("")

    lines.append("--- 6. Challenged assumption ---")
    lines.append(
        "Do not assume tournament vertices are the original vertices of T_m. "
        "For this proof, the useful vertices are support packets: chains, marked chains, "
        "border-twin packets, orbit packets, and two-step-lift obligations."
    )
    lines.append(
        "Preserved predicate: existence of a transitive witness after two skew doublings. "
        "Destroyed data: the full tower arc algebra outside the chosen support packet. "
        "The T511 failure says that destroyed data may contain exactly the missing +1."
    )
    lines.append("")
    lines.append("Conclusion: replace the next target by a support-realization problem:")
    lines.append("  classify packets X subset T127 with trans(X)=15 by q(X)=trans(D(D(X)))-30.")
    lines.append("  A packet with q(X)>=1 proves trans(T511)>=31; q(X)=0 for every local packet")
    lines.append("  would refute the naive recurrence mechanism and force a global witness search.")

    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
