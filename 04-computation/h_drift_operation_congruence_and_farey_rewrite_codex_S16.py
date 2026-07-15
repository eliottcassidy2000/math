#!/usr/bin/env python3
"""Exact operation-congruence audit for H-drift and the Farey ladder.

This file makes two distinctions which are easy to blur under the word
"self-similar".

1.  The Hamiltonian-path forward polynomial is a finite Koopman state for
    averaged H-observables.  It is *not* a state-transition quotient.  We
    enumerate converse-merged tournament classes through n=5 and locate the
    first strong-lumpability failure.

2.  The exact recursive carrier of the interval-core Farey ladder is the
    ordered denominator pair of a gap, with the Stern--Brocot mediant rewrite.
    The coarser pair of divisor-chain lengths already fails at order six.

No theorem or hypothesis identifier is claimed here.  The self-line section
hash-checks the three-way-certified n=8 obstruction already frozen in the
repository; it does not present that census as a new computation.
"""

from __future__ import annotations

import hashlib
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import permutations
from math import comb
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(str(message))


# ---------------------------------------------------------------------------
# Tournament / forward-polynomial audit


def edges(n: int) -> tuple[tuple[int, int], ...]:
    return tuple((i, j) for i in range(n) for j in range(i + 1, n))


def adjacency(mask: int, n: int, left: int, right: int) -> int:
    """Return 1 iff left -> right in the lexicographic full-arc mask."""
    if left == right:
        return 0
    if left > right:
        return 1 - adjacency(mask, n, right, left)
    index = edges(n).index((left, right))
    return (mask >> index) & 1


def encode_relabel(mask: int, n: int, image: tuple[int, ...]) -> int:
    """Relabel old vertex v as image[v]."""
    answer = 0
    index = {edge: k for k, edge in enumerate(edges(n))}
    for old_left, old_right in edges(n):
        winner, loser = (
            (old_left, old_right)
            if adjacency(mask, n, old_left, old_right)
            else (old_right, old_left)
        )
        new_winner, new_loser = image[winner], image[loser]
        low, high = sorted((new_winner, new_loser))
        if new_winner == low:
            answer |= 1 << index[(low, high)]
    return answer


def merged_code(mask: int, n: int) -> int:
    """Canonical code modulo relabelling and full converse."""
    full = (1 << comb(n, 2)) - 1
    return min(
        candidate
        for image in permutations(range(n))
        for candidate in (
            encode_relabel(mask, n, image),
            encode_relabel(mask ^ full, n, image),
        )
    )


def forward_counts(mask: int, n: int) -> tuple[int, ...]:
    """a_j = number of vertex orders with j forward consecutive arcs."""
    answer = [0] * n
    for order in permutations(range(n)):
        count = sum(
            adjacency(mask, n, order[k], order[k + 1])
            for k in range(n - 1)
        )
        answer[count] += 1
    return tuple(answer)


def krawtchouk(degree: int, position: int, dimension: int) -> int:
    return sum(
        (-1) ** chosen
        * comb(position, chosen)
        * comb(dimension - position, degree - chosen)
        for chosen in range(max(0, degree - dimension + position), min(degree, position) + 1)
    )


def walsh_stalk(mask: int, n: int) -> tuple[F, ...]:
    """The nonconstant even Krawtchouk/Walsh layers of H."""
    values = forward_counts(mask, n)
    dimension = n - 1
    return tuple(
        F(
            sum(
                krawtchouk(degree, position, dimension) * values[position]
                for position in range(n)
            ),
            1 << dimension,
        )
        for degree in range(2, dimension + 1, 2)
    )


def cyclic_triangles(mask: int, n: int) -> int:
    answer = 0
    for a in range(n):
        for b in range(a + 1, n):
            for c in range(b + 1, n):
                answer += int(
                    adjacency(mask, n, a, b)
                    and adjacency(mask, n, b, c)
                    and adjacency(mask, n, c, a)
                )
                answer += int(
                    adjacency(mask, n, a, c)
                    and adjacency(mask, n, c, b)
                    and adjacency(mask, n, b, a)
                )
    return answer


def score_sequence(mask: int, n: int) -> tuple[int, ...]:
    return tuple(
        sorted(sum(adjacency(mask, n, vertex, other) for other in range(n)) for vertex in range(n))
    )


def target_stalk_histogram(mask: int, n: int) -> tuple[tuple[tuple[F, ...], int], ...]:
    histogram = Counter(walsh_stalk(mask ^ (1 << edge), n) for edge in range(comb(n, 2)))
    return tuple(sorted(histogram.items()))


def h_gradient(mask: int, n: int) -> tuple[int, ...]:
    source_h = forward_counts(mask, n)[-1]
    return tuple(
        sorted(
            forward_counts(mask ^ (1 << edge), n)[-1] - source_h
            for edge in range(comb(n, 2))
        )
    )


def merged_rows(n: int) -> list[dict[str, object]]:
    representatives: dict[int, int] = {}
    for mask in range(1 << comb(n, 2)):
        code = merged_code(mask, n)
        representatives.setdefault(code, mask)
    rows = []
    for code, mask in sorted(representatives.items()):
        polynomial = forward_counts(mask, n)
        rows.append(
            {
                "code": code,
                "mask": mask,
                "A": polynomial,
                "H": polynomial[-1],
                "c3": cyclic_triangles(mask, n),
                "score": score_sequence(mask, n),
                "W": walsh_stalk(mask, n),
                "targetW": target_stalk_histogram(mask, n),
            }
        )
    return rows


def h_drift_audit() -> tuple[list[dict[str, object]], dict[str, object]]:
    summary = []
    rows_by_n = {}
    for n in range(2, 6):
        rows = merged_rows(n)
        rows_by_n[n] = rows
        fibres: dict[tuple[F, ...], list[dict[str, object]]] = defaultdict(list)
        for row in rows:
            fibres[row["W"]].append(row)
        split = sum(
            len({row["targetW"] for row in fibre}) > 1
            for fibre in fibres.values()
        )
        summary.append(
            {
                "n": n,
                "nodes": len(rows),
                "stalk_cells": len(fibres),
                "transition_split_cells": split,
            }
        )

    require([row["transition_split_cells"] for row in summary] == [0, 0, 0, 2], "split census changed")
    require([row["nodes"] for row in summary] == [1, 2, 3, 10], "merged-node census changed")

    # These are the full-arc masks corresponding to fixed-path tiling masks
    # n5-a01=1 and n5-a05=5 in the THM-848 atlas.
    left_mask, right_mask = 8, 10
    require(merged_code(left_mask, 5) != merged_code(right_mask, 5), "witness classes merged")
    left_A, right_A = forward_counts(left_mask, 5), forward_counts(right_mask, 5)
    left_W, right_W = walsh_stalk(left_mask, 5), walsh_stalk(right_mask, 5)
    left_targets = target_stalk_histogram(left_mask, 5)
    right_targets = target_stalk_histogram(right_mask, 5)
    require(left_A == right_A == (9, 30, 42, 30, 9), "forward-polynomial witness changed")
    require(left_W == right_W == (F(3, 2), F(0)), "Walsh witness changed")
    require(cyclic_triangles(left_mask, 5) == cyclic_triangles(right_mask, 5) == 3, "c3 witness changed")
    require(score_sequence(left_mask, 5) == score_sequence(right_mask, 5) == (1, 1, 2, 3, 3), "score witness changed")
    require(left_targets != right_targets, "target-W histograms no longer split")
    left_labelled_target = walsh_stalk(left_mask ^ 1, 5)
    right_labelled_target = walsh_stalk(right_mask ^ 1, 5)
    require(left_labelled_target == (F(-3, 2), F(-1)), "left labelled target changed")
    require(right_labelled_target == (F(-9, 2), F(0)), "right labelled target changed")

    left_gradient = h_gradient(left_mask, 5)
    right_gradient = h_gradient(right_mask, 5)
    require(sum(left_gradient) == sum(right_gradient) == -6, "mean drift changed")
    require(sum(x * x for x in left_gradient) == 140, "left diffusion changed")
    require(sum(x * x for x in right_gradient) == 156, "right diffusion changed")

    witness = {
        "left_mask": left_mask,
        "right_mask": right_mask,
        "A": left_A,
        "W": left_W,
        "H": 9,
        "c3": 3,
        "score": (1, 1, 2, 3, 3),
        "left_gradient": left_gradient,
        "right_gradient": right_gradient,
        "left_diffusion": F(14),
        "right_diffusion": F(78, 5),
        "left_targets": left_targets,
        "right_targets": right_targets,
        "left_labelled_target": left_labelled_target,
        "right_labelled_target": right_labelled_target,
    }
    return summary, {"rows": rows_by_n[5], "witness": witness}


# ---------------------------------------------------------------------------
# Farey / toothpick operation audit


def farey(order: int) -> tuple[F, ...]:
    return tuple(sorted({F(a, d) for d in range(1, order + 1) for a in range(d)}))


def farey_gaps(order: int) -> tuple[tuple[F, F, tuple[int, int]], ...]:
    fractions = farey(order)
    answer = []
    for left, right in zip(fractions, fractions[1:] + (F(1),)):
        pair = (left.denominator, 1 if right == 1 else right.denominator)
        answer.append((left, right, pair))
    return tuple(answer)


def pair_rewrite(pair: tuple[int, int], next_order: int) -> tuple[tuple[int, int], ...]:
    left, right = pair
    if left + right == next_order:
        return ((left, left + right), (left + right, right))
    if left + right > next_order:
        return (pair,)
    raise AssertionError("a Farey-neighbour pair cannot lag the current order")


def farey_rewrite_audit(max_order: int = 64) -> dict[str, object]:
    gap_steps = 0
    split_steps = 0
    duplicate_pair_fibres = 0
    for order in range(1, max_order):
        gaps = farey_gaps(order)
        expected_next = set(farey(order))
        pair_fibres: dict[tuple[int, int], set[tuple[tuple[int, int], ...]]] = defaultdict(set)
        for left, right, pair in gaps:
            i, j = pair
            require(right - left == F(1, i * j), "Farey gap length changed")
            rewrite = pair_rewrite(pair, order + 1)
            pair_fibres[pair].add(rewrite)
            gap_steps += 1
            if len(rewrite) == 2:
                split_steps += 1
                expected_next.add(F(left.numerator + right.numerator, i + j))
        require(tuple(sorted(expected_next)) == farey(order + 1), "Farey rewrite failed")
        require(all(len(signatures) == 1 for signatures in pair_fibres.values()), "pair rewrite became ambiguous")
        duplicate_pair_fibres += sum(
            sum(1 for _, _, seen_pair in gaps if seen_pair == pair) > 1
            for pair in pair_fibres
        )

    # The divisor-chain length code already fails to commute with refinement.
    order = 6
    examples = {}
    for left, right, pair in farey_gaps(order):
        chain_code = (order // pair[0], order // pair[1])
        if chain_code == (1, 2) and pair in ((4, 3), (5, 3)):
            examples[pair] = (left, right, pair_rewrite(pair, order + 1))
    require(examples == {
        (4, 3): (F(1, 4), F(1, 3), ((4, 7), (7, 3))),
        (5, 3): (F(3, 5), F(2, 3), ((5, 3),)),
    }, "chain-code counterexample changed")
    return {
        "orders": (1, max_order),
        "gap_steps": gap_steps,
        "split_steps": split_steps,
        "duplicate_pair_fibres": duplicate_pair_fibres,
        "chain_counterexample": examples,
    }


# ---------------------------------------------------------------------------
# Tournament Analysis on competing proof-facing quotients


def partition_stats(rows: list[dict[str, object]], key) -> tuple[int, int]:
    fibres = Counter(key(row) for row in rows)
    total = comb(len(rows), 2)
    separated = total - sum(comb(size, 2) for size in fibres.values())
    return len(fibres), separated


def tournament_fingerprint(values: list[F], tie_order: list[int]) -> dict[str, object]:
    n = len(values)
    rank = {vertex: position for position, vertex in enumerate(tie_order)}
    arcs = set()
    for left in range(n):
        for right in range(left + 1, n):
            if values[left] > values[right] or (
                values[left] == values[right] and rank[left] < rank[right]
            ):
                arcs.add((left, right))
            else:
                arcs.add((right, left))
    scores = [sum((v, w) in arcs for w in range(n) if w != v) for v in range(n)]
    score_hist = dict(sorted(Counter(scores).items()))
    triangles = sum(
        ((a, b) in arcs and (b, c) in arcs and (c, a) in arcs)
        or ((a, c) in arcs and (c, b) in arcs and (b, a) in arcs)
        for a in range(n)
        for b in range(a + 1, n)
        for c in range(b + 1, n)
    )

    # Tarjan is overkill here; direct mutual reachability is transparent.
    reach = [[False] * n for _ in range(n)]
    for v in range(n):
        reach[v][v] = True
    for v, w in arcs:
        reach[v][w] = True
    for middle in range(n):
        for left in range(n):
            for right in range(n):
                reach[left][right] |= reach[left][middle] and reach[middle][right]
    unseen = set(range(n))
    scc = []
    while unseen:
        seed = min(unseen)
        component = {v for v in unseen if reach[seed][v] and reach[v][seed]}
        unseen -= component
        scc.append(len(component))

    paths = 0
    for order in permutations(range(n)):
        paths += all((order[k], order[k + 1]) in arcs for k in range(n - 1))
    return {
        "arcs": arcs,
        "score": score_hist,
        "triangles": triangles,
        "scc": sorted(scc, reverse=True),
        "hamiltonian_paths": paths,
    }


def observer_tournament(rows: list[dict[str, object]]) -> dict[str, object]:
    observers = [
        ("H", lambda row: row["H"]),
        ("(H,c3)", lambda row: (row["H"], row["c3"])),
        ("forward polynomial A", lambda row: row["A"]),
        ("Walsh stalk W", lambda row: row["W"]),
        ("(W,target-W histogram)", lambda row: (row["W"], row["targetW"])),
        ("merged tournament class", lambda row: row["code"]),
    ]
    table = []
    for name, key in observers:
        cells, separated = partition_stats(rows, key)
        # (cells-1).bit_length() is ceil(log2(cells)), without floating point.
        table.append(
            (name, cells, separated, F(separated, max(1, (cells - 1).bit_length())))
        )
    raw = tournament_fingerprint([F(row[2]) for row in table], list(range(len(table))))
    economy = tournament_fingerprint([row[3] for row in table], list(range(len(table))))
    flips = len(raw["arcs"].symmetric_difference(economy["arcs"])) // 2
    return {"table": table, "raw": raw, "economy": economy, "flips": flips}


# ---------------------------------------------------------------------------
# Frozen self-line obstruction audit


SELF_LINE_FILES = {
    "04-computation/n8_black_self_line_obstruction_codex_S15.py":
        "8a316ef13e40d65215c2c524353d2f1e423829d7662aa5b303bbec676092e737",
    "05-knowledge/results/n8_black_self_line_obstruction_codex_S15.out":
        "6fda29d74551cb7caca0f95f7b693bf7ebfeddcff0fe9575ca963872b5a49910",
    "04-computation/selfk_sc_n8_check_opus_S312.py":
        "d9a80dff243f78de303b9a101c2ef59eb9b1bfe2ffa42d2e87a959891b14511d",
    "05-knowledge/results/selfk_sc_n8_check_opus_S312.out":
        "7bb62c54944b1e278ed66403bb03839668efeaff803402a4353c19e8851ae85d",
    "04-computation/thm852_selfline_sc_bijection_kps_S128c13.py":
        "27800b8d858ab48b90e750808510456ddd87b64eefa8c63c2a6eae5e67f19fd5",
    "05-knowledge/results/thm852_selfline_sc_bijection_kps_S128c13.out":
        "8564717315c72beaf4ffc084ec8f2a87d6fa117a1d13b17cd7df7ac7bece42e1",
}


def self_line_hash_audit() -> None:
    for relative, expected in SELF_LINE_FILES.items():
        actual = hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()
        require(actual == expected, (relative, expected, actual))
    primary = (ROOT / "05-knowledge/results/n8_black_self_line_obstruction_codex_S15.out").read_text()
    for needle in (
        "Qblack  selfB  selfK",
        " 5    6       64       12    8          16   8      0       8      0      4",
        " 6   10     1024       56   12         180  16      4      12      2      6",
        " 7   15    32768      456   88        4540  88      0      88      0     44",
        " 8   21  2097152     6880  176      195004 412      8     404      4    202",
        "2*selfK(8) = Q_black(8) = 404",
    ):
        require(needle in primary, f"missing frozen self-line ledger: {needle}")


def main() -> None:
    summary, h_data = h_drift_audit()
    farey_data = farey_rewrite_audit()
    observer = observer_tournament(h_data["rows"])
    self_line_hash_audit()
    witness = h_data["witness"]

    print("H-DRIFT: KOOPMAN CLOSURE IS NOT OPERATION CONGRUENCE")
    print("=" * 78)
    print("q(T)=W_n(T); strong lumpability requires # {e:q(T^e)=w} to be q-fibre-constant")
    for row in summary:
        print(
            "n={n}: merged_nodes={nodes}, W_cells={stalk_cells}, "
            "transition-split W_cells={transition_split_cells}".format(**row)
        )
    print("first failure: n=5 (two split W-fibres; none through n=4)")
    print(
        "witness: atlas n5-a01/n5-a05; full-arc masks "
        f"{witness['left_mask']}/{witness['right_mask']}"
    )
    print(f"common A_T coefficients={witness['A']}")
    print(
        f"common W={witness['W']}, H={witness['H']}, c3={witness['c3']}, "
        f"score={witness['score']}"
    )
    print(
        "same labelled flip {1,2}: target W="
        f"{witness['left_labelled_target']} vs {witness['right_labelled_target']}"
    )
    print(f"gradient n5-a01={witness['left_gradient']}")
    print(f"gradient n5-a05={witness['right_gradient']}")
    print(
        f"common drift=-3/5; diffusion={witness['left_diffusion']} vs "
        f"{witness['right_diffusion']}"
    )
    print(f"target-W histogram n5-a01={witness['left_targets']}")
    print(f"target-W histogram n5-a05={witness['right_targets']}")
    print()
    print("exact functional-form correction:")
    print("  S B(z)=-2 z B'(z); exp(tS)B(z)=B(exp(-2t)z) (continuous Poissonization)")
    print("  P B(z)=B(z)-(2/M)zB'(z); P^t mode=(1-4r/M)^t (discrete chain)")
    print("  P is not B(z)->B(cz) once H2,H4 both occur:")
    print("  c^4 would need both (1-4/M)^2 and 1-8/M; their gap is 16/M^2")
    print("verdict: Mobius radiality is an observable/Koopman conjugacy, not a state quotient")

    print("\nFAREY LADDER: THE TRUE REWRITE CARRIER")
    print("=" * 78)
    print("ordered gap pair (i,j) at F_k:")
    print("  if i+j=k+1: (i,j) -> (i,i+j),(i+j,j)")
    print("  if i+j>k+1: (i,j) -> (i,j)")
    print(
        "exact Farey replay orders %d..%d: gap steps=%d, splits=%d, duplicate pair fibres=%d"
        % (
            farey_data["orders"][0],
            farey_data["orders"][1],
            farey_data["gap_steps"],
            farey_data["split_steps"],
            farey_data["duplicate_pair_fibres"],
        )
    )
    print("THM-841's local interval-core profile also factors through (k,i,j,lambda)")
    print("because gap length=1/(ij) and both multiplier chains are determined by i,j")
    print("coarse divisor-chain code failure at k=6:")
    for pair, (left, right, rewrite) in farey_data["chain_counterexample"].items():
        print(f"  gap ({left},{right}), pair={pair}, chain code=(1,2), next={rewrite}")
    print("verdict: Stern-Brocot ordered pairs are operation-congruent; tooth counts are not")
    print("the ordered pair is an exact address (no duplicate pair fibres through order 64),")
    print("not evidence for a smaller count-only compression")

    print("\nTOURNAMENT ANALYSIS: OBSERVER QUOTIENTS AS VERTICES")
    print("=" * 78)
    print("pairwise observable=unordered n=5 merged-node pairs separated by observer")
    print("switch/gauge=raw separation -> separation/ceil(log2(number of cells))")
    print("tie Hamiltonian path=the displayed observer order")
    for index, (name, cells, separated, economy) in enumerate(observer["table"]):
        print(f"v{index} {name:28s}: cells={cells:2d}, separated={separated:2d}, economy={economy}")
    for name in ("raw", "economy"):
        fp = observer[name]
        print(
            f"{name}: score={fp['score']}, C3={fp['triangles']}, SCC={fp['scc']}, "
            f"HP={fp['hamiltonian_paths']}"
        )
    print(f"gauge edge flips={observer['flips']}")

    print("\nALTERNATE-VERTEX / INFORMATION AUDIT")
    print("=" * 78)
    print("original tournament vertices: preserve local arcs/scores; destroy needle-boundary totals after quotient")
    print("arc-flip events/orbits: preserve one-step target incidence; destroy multi-step wall chronology")
    print("Hamiltonian needles: preserve H and the one-defect boundary U1; counts forget the failed-arc owner")
    print("path gaps/forward-count bins: preserve A_T and every averaged H mode; destroy target-W incidence")
    print("Walsh/Fourier modes: preserve all future means of H; destroy transition law already at n=5")
    print("Farey gaps (ordered denominator pairs): preserve the interval-core W-profile and mediant rewrite")
    print("fixed sections/boundaries/wall events/cover arcs: need offset+owner+side to preserve local coverage")
    print("residues: preserve modular incidence; destroy metric clearance and carry")
    print("matroid circuits: preserve dependence support; destroy coefficients, scale, and strict inequalities")
    print("proof obligations: preserve the comparison being optimized; they are observers, not runners/arcs")
    print("LRC predicate audit: none of the H/tournament quotients preserve loneliness")
    print("the Farey-pair quotient preserves only the special {1,...,k} interval-core violation predicate")
    print("and its measures; it does not preserve arbitrary-speed LRC loneliness")

    print("\nSELF-LINE ALL-n CLAIM: FROZEN EXACT OBSTRUCTION")
    print("=" * 78)
    print("black quasi-fixed endpoints 2*selfK at n=5,6,7,8: 8,12,88,404")
    print("self-converse classes SC at n=5,6,7,8:           8,12,88,176")
    print("n=6 total Q=16 includes four blue endpoints, but removing them only repairs n=6")
    print("n=8 has Qblue=8 and Qblack=404, so the 228 black excess survives the restriction")
    print("three independent source/output pairs match their frozen SHA-256 hashes")
    for relative, digest in SELF_LINE_FILES.items():
        print(f"  {digest}  {relative}")
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
