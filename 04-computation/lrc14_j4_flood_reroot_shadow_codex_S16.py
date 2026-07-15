#!/usr/bin/env python3
r"""Exact re-root descent and five-small shadow for the THM-741 flood tail.

Write H={8,...,14} and E_e=H union e for an edge e of K_7.  A completed
13-speed family is not attached to a unique root edge: if its small-label set
contains a triple {a,b,c}, the same family is a four-speed extension of each
of E_ab, E_ac, and E_bc.  This script verifies the resulting transition
graphs and the exact finite tail needed below.

The three already-certified bodies E_56, E_57, and E_67 therefore close every
final family containing one of those anchors.  For a final family containing
five labels from {1,...,7}, 18 of the 21 possible five-sets contain an anchor.
For each of the three residual five-sets K, put P=H union K.  P has twelve speeds, so
only one speed w>=15 remains.  If G(P) has r components and measure m, the
THM-732 discrepancy estimate gives

    |G(P) \ D_w| >= 6m/7 - sqrt(2) r/(7w)
                  > 6m/7 - (99/70) r/(7w).

Thus only 15 <= w <= floor((99/70)r/(6m)) need exact interval sweeps.  The
three finite ranges contain 260 values in total and all are positive.  Every
six-set contains an anchor.  This proves the strict shadow:
every 13-speed family containing H and at least five small labels is lonely.
It does not close any whole unresolved flood body, whose extensions may have
only two, three, or four small labels.

Tournament Analysis deliberately changes vertices from the 21 root edges to
the three residual five-label bases.  Its pair observable is the difference
of exact one-speed tail caps C(K)=(99/70)r(K)/(6m(K)); the proof-work gauge
points from lower to higher cap, with lexicographic order as the tie gauge.
This preserves the finite search horizon but destroys the interval comb and
the exact clearance.  Root edges are rejected as vertices for transport:
the non-Fano re-root triangles connect all 21 of them, so any edge-only scalar
that descends to completed-family data is constant.
"""

from __future__ import annotations

import hashlib
import importlib.util
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CORE_PATH = ROOT / "04-computation/lrc14_thm741_2002_body_j4_tree_kps_S128c5.py"
CORE_SHA256 = "5aa81d9d78273c8f9e3e7a6574091a3bc3f64ab6086c7024c15f9420c99dac96"
ANCHOR_OUTPUTS = {
    (5, 6): (
        ROOT / "05-knowledge/results/lrc14_j4_flood_56_exact_codex_S16.out",
        "f75cce66766f68e4c44c3a5d68a17136f135cdf775f396591517ac55e793a233",
    ),
    (5, 7): (
        ROOT / "05-knowledge/results/lrc14_j4_flood_57_exact_codex_S14.out",
        "582097f06cb0f66f5deedd7814e127e22f48b22051bed612305ecd7ea3c062b0",
    ),
    (6, 7): (
        ROOT / "05-knowledge/results/lrc14_j4_flood_67_exact_codex_S15.out",
        "5f8e5dc2108c15ba4f7d80ce7d77805c58eb5c292b326f704ab3a28279d11694",
    ),
}
H = frozenset(range(8, 15))
SMALL = frozenset(range(1, 8))
ANCHORS = tuple(frozenset(edge) for edge in ANCHOR_OUTPUTS)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_core():
    digest = sha256(CORE_PATH)
    require(digest == CORE_SHA256, f"THM-741 dependency changed: {digest}")
    spec = importlib.util.spec_from_file_location("thm741_reroot_dependency", CORE_PATH)
    require(spec is not None and spec.loader is not None, "cannot load THM-741 core")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def transition_components(
    edges: tuple[tuple[int, int], ...], triples: tuple[tuple[int, int, int], ...]
) -> tuple[tuple[int, ...], tuple[int, ...], int]:
    """Component sizes, degree set, and undirected edge count of chart transitions."""
    index = {edge: i for i, edge in enumerate(edges)}
    adjacency = [set() for _ in edges]
    for triple in triples:
        charts = tuple(combinations(triple, 2))
        for left, right in combinations(charts, 2):
            i, j = index[left], index[right]
            adjacency[i].add(j)
            adjacency[j].add(i)
    unseen = set(range(len(edges)))
    components = []
    while unseen:
        seed = min(unseen)
        component = {seed}
        frontier = [seed]
        unseen.remove(seed)
        while frontier:
            vertex = frontier.pop()
            for neighbour in adjacency[vertex]:
                if neighbour in unseen:
                    unseen.remove(neighbour)
                    component.add(neighbour)
                    frontier.append(neighbour)
        components.append(len(component))
    return (
        tuple(sorted(components)),
        tuple(sorted({len(neighbours) for neighbours in adjacency})),
        sum(len(neighbours) for neighbours in adjacency) // 2,
    )


def main() -> None:
    core = load_core()
    require(core.S2 == F(99, 70) and core.S2 * core.S2 > 2, "bad sqrt(2) majorant")

    for edge, (path, digest) in ANCHOR_OUTPUTS.items():
        require(sha256(path) == digest, f"anchor output {edge} changed")
        payload = path.read_text()
        require(
            "EXACTLY CLOSED REQUESTED FLOOD BODIES=1" in payload
            and "covering failures=0" in payload,
            f"anchor output {edge} lacks its closure verdict",
        )

    # Re-root descent.  Every triple gives three presentations of one actual
    # family H union triple union {15,16,17}.
    edges = tuple(combinations(range(1, 8), 2))
    triples = tuple(combinations(range(1, 8), 3))
    fano = tuple(sorted({tuple(sorted((a, b, a ^ b))) for a, b in edges}))
    nonfano = tuple(triple for triple in triples if triple not in set(fano))
    require((len(edges), len(triples), len(fano), len(nonfano)) == (21, 35, 7, 28), "bad atlas")
    fano_components = transition_components(edges, fano)
    nonfano_components = transition_components(edges, nonfano)
    all_components = transition_components(edges, triples)
    require(fano_components == ((3, 3, 3, 3, 3, 3, 3), (2,), 21), "bad Fano transitions")
    require(nonfano_components == ((21,), (8,), 84), "bad non-Fano transitions")
    require(all_components == ((21,), (10,), 105), "bad full transitions")

    five_sets = tuple(frozenset(K) for K in combinations(range(1, 8), 5))
    anchored_five = tuple(K for K in five_sets if any(anchor <= K for anchor in ANCHORS))
    residual_five = tuple(K for K in five_sets if K not in anchored_five)
    require((len(five_sets), len(anchored_five), len(residual_five)) == (21, 18, 3), "bad five-set split")

    rows = []
    total_sweeps = 0
    global_minimum: tuple[F, tuple[int, ...], int] | None = None
    for K in residual_five:
        P = tuple(sorted(H | K))
        require(len(P) == 12, "residual base is not a twelve-speed family")
        G, r, m = core.good_norm(P)
        require(m > 0, f"empty residual base K={tuple(sorted(K))}")
        cap = core.S2 * r / (6 * m)
        bmax = cap.numerator // cap.denominator
        minimum: tuple[F, int] | None = None
        for w in range(15, bmax + 1):
            value = core.subtract_sparse(G, w)
            require(value > 0, f"zero exact sweep K={tuple(sorted(K))}, w={w}")
            candidate = (value, w)
            if minimum is None or candidate < minimum:
                minimum = candidate
            total_sweeps += 1
            global_candidate = (value, tuple(sorted(K)), w)
            if global_minimum is None or global_candidate < global_minimum:
                global_minimum = global_candidate
        first_tail = max(15, bmax + 1)
        tail_lower = F(6, 7) * m - core.S2 * r / (7 * first_tail)
        require(first_tail > cap and tail_lower > 0, f"tail did not close K={tuple(sorted(K))}")
        rows.append(
            {
                "K": tuple(sorted(K)),
                "r": r,
                "m": m,
                "cap": cap,
                "bmax": bmax,
                "sweeps": max(0, bmax - 14),
                "minimum": minimum,
                "first_tail": first_tail,
                "tail_lower": tail_lower,
            }
        )
    require(total_sweeps == 260, f"unexpected sweep total {total_sweeps}")
    require(global_minimum is not None, "no exact sweeps")

    six_sets = tuple(frozenset(K) for K in combinations(range(1, 8), 6))
    anchored_six = tuple(K for K in six_sets if any(anchor <= K for anchor in ANCHORS))
    residual_six = tuple(K for K in six_sets if K not in anchored_six)
    require((len(six_sets), len(anchored_six), len(residual_six)) == (7, 7, 0), "bad six-set split")

    # Proof-work tournament on residual bases.  All caps are distinct, hence
    # the gauge is transitive and has a unique Hamiltonian path.
    lex_rows = sorted(rows, key=lambda row: row["K"])
    require(len({row["cap"] for row in rows}) == len(rows), "tail-cap tie")
    cap_path = tuple(row["K"] for row in sorted(rows, key=lambda row: row["cap"]))
    flips = sum(
        1
        for i, left in enumerate(lex_rows)
        for right in lex_rows[i + 1 :]
        if left["cap"] > right["cap"]
    )

    print("THM-741 FLOOD RE-ROOT DESCENT AND FIVE-SMALL SHADOW")
    print("=" * 88)
    print(f"dependency_sha256={CORE_SHA256}")
    print(
        "anchor_output_sha256="
        + ",".join(f"{edge}:{digest}" for edge, (_, digest) in ANCHOR_OUTPUTS.items())
    )
    print("root charts=21 K7 edges; transition families=35 small-label triples")
    print(
        "Fano transitions: triples=7 components=7x3 degree=2 edges=21 "
        "equality_rank=14 quotient_dim=7"
    )
    print(
        "non-Fano transitions: triples=28 components=1x21 degree=8 edges=84 "
        "equality_rank=20 quotient_dim=1"
    )
    print(
        "all transitions: triples=35 components=1x21 degree=10 edges=105 "
        "equality_rank=20 quotient_dim=1"
    )
    print("descent no-go: every root-edge-only completed-family invariant is constant")
    print("five-small sets total/anchor-shadow/residual=21/18/3")
    print("residual rows K ; r ; m ; cap ; exact_w ; minimum ; first_tail ; tail_lower")
    for row in rows:
        minimum = row["minimum"]
        require(minimum is not None, "missing row minimum")
        print(
            f"  {row['K']} ; {row['r']} ; {row['m']} ; {row['cap']} ; "
            f"15..{row['bmax']} ({row['sweeps']}) ; {minimum[0]}@{minimum[1]} ; "
            f"{row['first_tail']} ; {row['tail_lower']}"
        )
    print(
        f"exact residual sweeps={total_sweeps} zeros=0 ; "
        f"global_exact_minimum={global_minimum[0]} at K={global_minimum[1]},w={global_minimum[2]}"
    )
    print(
        "six-small sets total/anchor-shadow/residual=7/7/0 ; "
        "the former residual K=(1,2,3,4,5,6) is the (5,6)-anchor minimum"
    )
    print("Tournament Analysis vertices: three residual five-label bases (not runners or root edges)")
    print("pair observable: C(B)-C(A), C=(99/70)r/(6m); gauge: lower C -> higher C; tie: lex")
    print(
        f"fingerprint: score_hist={{0..2:1}}, directed_3cycles=0, SCC_sizes=3x1, "
        f"edge_flips_vs_lex={flips}, Hamiltonian_paths=1"
    )
    print(f"tie Hamiltonian path={cap_path}")
    print("kept: exact one-speed horizon; destroyed: interval comb and exact clearance")
    print("PROVED SHADOW: H plus at least five small labels is strictly lonely")
    print("SCOPE: unresolved flood families have at most four small labels; THM-741 remains CLAIMED")
    print(f"source_sha256={sha256(Path(__file__).resolve())}")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
