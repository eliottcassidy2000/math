#!/usr/bin/env python3
"""Exact referee and second-Krylov consequence of THM-848.

THM-848 identifies the Hamiltonian-path count H with a finite sum of even
Walsh layers and the odd-cycle-forest statistic K with (n-D)H.  This audit
pushes that dictionary one generator step farther:

    E[Delta K | T] = (a_(n-3)+a_(n-2)-C(n,2)H)/C(n,2)
                   = -4/C(n,2) sum_r r(n-2r) H_(2r).

It also checks the pointwise edge-current formula

    K(T^e)-K(T) = 4(sum_gained K(T-V(C))-sum_lost K(T-V(C))),

where the sums run over directed odd cycles using the flipped unordered
pair, and K(empty)=1/2.  Every assertion is exact.  The finite census uses
one representative of every converse-merged node for n=4,5,6,7.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import itertools
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from math import ceil, comb, log2
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
THM848_AUDIT = ROOT / "04-computation/h_drift_walsh_needle_stalk_codex_S15.py"
THM848_AUDIT_SHA256 = "dbab32bfcfb091ef8fa7f7558733127d1eeea3e8dbca63f4f622b9eff1862e6c"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_thm848_audit():
    digest = sha256(THM848_AUDIT)
    if digest != THM848_AUDIT_SHA256:
        raise RuntimeError(
            "THM-848 audit changed: expected "
            f"{THM848_AUDIT_SHA256}, got {digest}"
        )
    spec = importlib.util.spec_from_file_location("thm848_audit", THM848_AUDIT)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load THM-848 audit")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def forward_top(adjacency: list[list[int]], defects: int = 2) -> tuple[int, ...]:
    """Counts vertex orders with 0,1,...,defects backward adjacencies."""
    n = len(adjacency)
    if n == 0:
        return (1,) + (0,) * defects
    dp = [[[0] * (defects + 1) for _ in range(n)] for _ in range(1 << n)]
    for vertex in range(n):
        dp[1 << vertex][vertex][0] = 1
    for support in range(1 << n):
        for last in range(n):
            for used_defects, count in enumerate(dp[support][last]):
                if not count:
                    continue
                for nxt in range(n):
                    if (support >> nxt) & 1:
                        continue
                    new_defects = used_defects + (1 - adjacency[last][nxt])
                    if new_defects <= defects:
                        dp[support | (1 << nxt)][nxt][new_defects] += count
    return tuple(sum(dp[-1][last][r] for last in range(n)) for r in range(defects + 1))


def K_from_forward(adjacency: list[list[int]]) -> Fraction:
    """K=(a_(n-2)+(n+1)H)/2, with the natural empty convention."""
    n = len(adjacency)
    if n == 0:
        return Fraction(1, 2)
    if n == 1:
        return Fraction(1)
    h_value, one_defect, _ = forward_top(adjacency)
    return Fraction(one_defect + (n + 1) * h_value, 2)


@lru_cache(maxsize=None)
def odd_cycle_candidates(n: int) -> tuple[tuple[int, ...], ...]:
    candidates: list[tuple[int, ...]] = []
    for length in range(3, n + 1, 2):
        for vertices in itertools.combinations(range(n), length):
            anchor = vertices[0]
            for tail in itertools.permutations(vertices[1:]):
                candidates.append((anchor,) + tail)
    return tuple(candidates)


def directed_odd_cycles(adjacency: list[list[int]]) -> tuple[tuple[int, ...], ...]:
    return tuple(
        cycle
        for cycle in odd_cycle_candidates(len(adjacency))
        if all(
            adjacency[cycle[index]][cycle[(index + 1) % len(cycle)]]
            for index in range(len(cycle))
        )
    )


def cycle_uses_pair(cycle: tuple[int, ...], left: int, right: int) -> bool:
    length = len(cycle)
    return any(
        {cycle[index], cycle[(index + 1) % length]} == {left, right}
        for index in range(length)
    )


def induced(adjacency: list[list[int]], vertices: tuple[int, ...]) -> list[list[int]]:
    return [[adjacency[left][right] for right in vertices] for left in vertices]


def alpha_levels(
    adjacency: list[list[int]], cycles: tuple[tuple[int, ...], ...]
) -> tuple[tuple[int, int], ...]:
    """Return (f,alpha_f), alpha_f=sum_{I:f_I=f}2^|I|, for n<=7."""
    n = len(adjacency)
    q = n - 1
    by_length = Counter(map(len, cycles))
    triangle_masks = [sum(1 << vertex for vertex in cycle) for cycle in cycles if len(cycle) == 3]
    disjoint_triangle_pairs = sum(
        left & right == 0
        for left, right in itertools.combinations(triangle_masks, 2)
    )
    alpha = {q: 1}
    if n >= 3:
        alpha[n - 3] = 2 * by_length[3]
    if n >= 5:
        alpha[n - 5] = 2 * by_length[5] + 4 * disjoint_triangle_pairs
    if n >= 7:
        alpha[n - 7] = 2 * by_length[7]
    return tuple(sorted(alpha.items(), reverse=True))


def fibre_stats(rows: list[dict[str, object]], key, target) -> tuple[int, int, int, int]:
    fibres: dict[object, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        fibres[key(row)].append(row)
    split = sum(len({target(row) for row in fibre}) > 1 for fibre in fibres.values())
    separated = comb(len(rows), 2) - sum(comb(len(fibre), 2) for fibre in fibres.values())
    return len(fibres), max(map(len, fibres.values())), split, separated


def render_fraction(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def alpha_levels_n8(
    cycles: tuple[tuple[int, ...], ...]
) -> tuple[tuple[int, int], ...]:
    by_length = Counter(map(len, cycles))
    triangle_masks = [sum(1 << vertex for vertex in cycle) for cycle in cycles if len(cycle) == 3]
    five_masks = [sum(1 << vertex for vertex in cycle) for cycle in cycles if len(cycle) == 5]
    disjoint_triangle_pairs = sum(
        left & right == 0
        for left, right in itertools.combinations(triangle_masks, 2)
    )
    disjoint_three_five_pairs = sum(
        triangle & five == 0 for triangle in triangle_masks for five in five_masks
    )
    return (
        (7, 1),
        (5, 2 * by_length[3]),
        (3, 2 * by_length[5] + 4 * disjoint_triangle_pairs),
        (1, 2 * by_length[7] + 4 * disjoint_three_five_pairs),
    )


def n8_fibre_audit(needle, atlas_path: Path) -> None:
    representatives, digest = needle.certified_n8_representatives(atlas_path)
    fibres: dict[tuple[int, int], list[dict[str, int]]] = defaultdict(list)
    for rank, mask in enumerate(representatives):
        adjacency = needle.adjacency_from_legacy_tiling(mask, 8)
        h_value, one_defect, two_defects = forward_top(adjacency)
        k_value = (one_defect + 9 * h_value) // 2
        assert 2 * k_value == one_defect + 9 * h_value
        drift_sum = two_defects + one_defect - 28 * h_value
        fibres[(h_value, k_value)].append(
            {
                "rank": rank,
                "mask": mask,
                "b1": one_defect,
                "b2": two_defects,
                "SK": drift_sum,
            }
        )
    split = [
        (key, rows)
        for key, rows in fibres.items()
        if len({row["SK"] for row in rows}) > 1
    ]
    assert len(representatives) == 3528
    assert len(fibres) == 1727
    assert max(map(len, fibres.values())) == 13
    assert len(split) == 1 and split[0][0] == (345, 2820)
    key, witnesses = split[0]
    assert len(witnesses) == 2
    witnesses.sort(key=lambda row: row["SK"])

    print("THM-848 SECOND-KRYLOV OPTIONAL HASH-GUARDED N=8 FIBRE AUDIT")
    print("=" * 82)
    print("command: python3 04-computation/h_drift_K_second_krylov_referee_codex_S14.py \\")
    print("         --n8-atlas /tmp/n8_merged_node_rank_u16.bin")
    print("atlas producer: 04-computation/mobius_cech_n8_node_atlas_export_codex_S13.py")
    print("atlas format: 2^21 little-endian uint16 converse-merged node ranks")
    print(f"atlas SHA256: {digest}")
    print("merged nodes=3528, (H,K) cells=1727, max fibre=13")
    print("(H,K) fibres split by E[Delta K]=1 (the unique first actual split)")
    print(f"unique split key: (H,K)={key}")
    alpha_vectors = []
    for row in witnesses:
        adjacency = needle.adjacency_from_legacy_tiling(row["mask"], 8)
        levels = alpha_levels_n8(directed_odd_cycles(adjacency))
        alpha = tuple(dict(levels)[f_value] for f_value in (5, 3, 1))
        assert 1 + sum(alpha) == key[0]
        assert 128 + 32 * alpha[0] + 8 * alpha[1] + 2 * alpha[2] == key[1]
        predicted_b2 = 4293 + 189 * alpha[0] - 27 * alpha[1] + 9 * alpha[2]
        assert predicted_b2 == row["b2"]
        alpha_vectors.append(alpha)
        print(
            f"rank={row['rank']}, mask={row['mask']}, a6={row['b1']}, a5={row['b2']}, "
            f"S_K={row['SK']}, E[Delta K]={render_fraction(Fraction(row['SK'], 28))}, "
            f"alpha_(5,3,1)={alpha}"
        )
    alpha_delta = tuple(right - left for left, right in zip(*alpha_vectors))
    assert alpha_delta == (2, -10, 8)
    assert witnesses[1]["b2"] - witnesses[0]["b2"] == 720
    assert witnesses[1]["SK"] - witnesses[0]["SK"] == 720
    print("alpha difference=(2,-10,8)=2*(1,-5,4): the n=7 ambient kernel is realized")
    print("Delta a5=Delta S_K=720; drift gap=180/7")
    print("scope: finite-exact conditional on the hash-verified scratch atlas; no all-n claim")
    print("ALL ASSERTIONS PASSED")


def main() -> None:
    if not __debug__:
        raise RuntimeError("exact audit requires assertions; do not run with python -O")
    parser = argparse.ArgumentParser()
    parser.add_argument("--n8-atlas", type=Path)
    args = parser.parse_args()
    needle = load_thm848_audit()
    if args.n8_atlas is not None:
        n8_fibre_audit(needle, args.n8_atlas)
        return
    small = __import__("json").loads(needle.DEFAULT_SMALL_ATLAS.read_text())
    n7 = __import__("json").loads(needle.DEFAULT_N7_ATLAS.read_text())
    atlas_by_n = {entry["n"]: entry for entry in small["sizes"]}
    atlas_by_n[7] = n7
    support_bank = {
        n: {
            level: needle.even_path_forests(n, level)
            for level in range(1, (n - 1) // 2 + 1)
        }
        for n in range(4, 8)
    }

    print("THM-848 LIVE-PULL REFEREE: K IS THE SECOND KRYLOV NEEDLE")
    print("=" * 82)
    print(f"guarded THM-848 audit SHA256={THM848_AUDIT_SHA256}")
    print("generator: S f=sum_e(f(T^e)-f(T)), M=C(n,2)")
    print("S a_j=(j+1)a_(j+1)+(n-j)a_(j-1)-(n-1)a_j")
    print("K=(a_(n-2)+(n+1)H)/2")
    print("S K=a_(n-3)+a_(n-2)-M H")
    print("Walsh form: E[Delta K|T]=-(4/M) sum_r r(n-2r)H_(2r)")

    all_rows: dict[int, list[dict[str, object]]] = {}
    edge_formula_checks = 0
    for n in range(4, 8):
        arc_count = comb(n, 2)
        rows: list[dict[str, object]] = []
        for node in atlas_by_n[n]["nodes"]:
            adjacency = needle.adjacency_from_tiling(node["tiling_masks"][0], n)
            top = forward_top(adjacency)
            h_value, one_defect, two_defects = top
            k_value = Fraction(one_defect + (n + 1) * h_value, 2)
            assert k_value.denominator == 1
            layers = tuple(
                needle.needle_layer(adjacency, support_bank[n][level])
                for level in range(1, (n - 1) // 2 + 1)
            )
            k_walsh = n * h_value - 2 * sum(
                level * layer for level, layer in enumerate(layers, start=1)
            )
            assert k_value == k_walsh
            predicted_h_sum = one_defect - (n - 1) * h_value
            predicted_k_sum = two_defects + one_defect - arc_count * h_value
            walsh_k_sum = -4 * sum(
                level * (n - 2 * level) * layer
                for level, layer in enumerate(layers, start=1)
            )
            assert predicted_k_sum == walsh_k_sum

            old_cycles = directed_odd_cycles(adjacency)
            levels = alpha_levels(adjacency, old_cycles)
            assert h_value == sum(alpha for _, alpha in levels)
            assert k_value == sum((1 << f_value) * alpha for f_value, alpha in levels)
            complement_k: dict[int, Fraction] = {}

            def deleted_k(cycle: tuple[int, ...]) -> Fraction:
                retained_mask = ((1 << n) - 1) ^ sum(1 << vertex for vertex in cycle)
                if retained_mask not in complement_k:
                    vertices = tuple(v for v in range(n) if (retained_mask >> v) & 1)
                    complement_k[retained_mask] = K_from_forward(induced(adjacency, vertices))
                return complement_k[retained_mask]

            direct_k_sum = Fraction(0)
            k_gradient: list[int] = []
            for left in range(n):
                for right in range(left + 1, n):
                    lost = tuple(
                        cycle for cycle in old_cycles if cycle_uses_pair(cycle, left, right)
                    )
                    adjacency[left][right], adjacency[right][left] = (
                        adjacency[right][left],
                        adjacency[left][right],
                    )
                    gained_cycles = directed_odd_cycles(adjacency)
                    gained_levels = alpha_levels(adjacency, gained_cycles)
                    new_k = Fraction(
                        sum((1 << f_value) * alpha for f_value, alpha in gained_levels)
                    )
                    gained = tuple(
                        cycle
                        for cycle in gained_cycles
                        if cycle_uses_pair(cycle, left, right)
                    )
                    current = new_k - k_value
                    cycle_current = 4 * (
                        sum((deleted_k(cycle) for cycle in gained), Fraction(0))
                        - sum((deleted_k(cycle) for cycle in lost), Fraction(0))
                    )
                    assert current == cycle_current
                    assert current.denominator == 1
                    direct_k_sum += current
                    k_gradient.append(current.numerator)
                    edge_formula_checks += 1
                    adjacency[left][right], adjacency[right][left] = (
                        adjacency[right][left],
                        adjacency[left][right],
                    )
            assert direct_k_sum == predicted_k_sum
            rows.append(
                {
                    "id": node["id"],
                    "mask": node["tiling_masks"][0],
                    "H": h_value,
                    "K": k_value.numerator,
                    "b1": one_defect,
                    "b2": two_defects,
                    "bH": Fraction(predicted_h_sum, arc_count),
                    "bK": Fraction(predicted_k_sum, arc_count),
                    "layers": layers,
                    "alpha": levels,
                    "Kgrad": tuple(sorted(k_gradient)),
                }
            )
        all_rows[n] = rows

        h_stats = fibre_stats(rows, lambda row: row["H"], lambda row: row["bH"])
        h_k_stats = fibre_stats(rows, lambda row: row["H"], lambda row: row["bK"])
        k_stats = fibre_stats(rows, lambda row: row["K"], lambda row: row["bK"])
        hk_stats = fibre_stats(rows, lambda row: (row["H"], row["K"]), lambda row: row["bK"])
        hkb2_stats = fibre_stats(
            rows, lambda row: (row["H"], row["K"], row["b2"]), lambda row: row["bK"]
        )
        print(
            f"n={n}: merged nodes={len(rows)}, "
            f"H cells/max/split(bH,bK)={h_stats[0]}/{h_stats[1]}/{h_stats[2]},{h_k_stats[2]}, "
            f"K cells/max/split(bK)={k_stats[0]}/{k_stats[1]}/{k_stats[2]}, "
            f"(H,K) cells/max/split(bK)={hk_stats[0]}/{hk_stats[1]}/{hk_stats[2]}, "
            f"(H,K,a_(n-3)) cells/max/split(bK)="
            f"{hkb2_stats[0]}/{hkb2_stats[1]}/{hkb2_stats[2]}"
        )

    assert all(
        fibre_stats(all_rows[n], lambda row: (row["H"], row["K"]), lambda row: row["bK"])[2]
        == 0
        for n in range(4, 8)
    )

    print("\nTHE n=7 AMBIENT KERNEL IS NOT REALIZED")
    print("-" * 82)
    n7_rows = all_rows[7]
    hk_to_alpha: dict[tuple[int, int], set[tuple[tuple[int, int], ...]]] = defaultdict(set)
    for row in n7_rows:
        hk_to_alpha[(row["H"], row["K"])].add(row["alpha"])
    assert all(len(values) == 1 for values in hk_to_alpha.values())
    hk_cells = len(hk_to_alpha)
    alpha_cells = len({row["alpha"] for row in n7_rows})
    walsh_cells = len({row["layers"] for row in n7_rows})
    assert hk_cells == alpha_cells == walsh_cells == 156
    print("alpha_6=1 and the nonempty level vector is (alpha_4,alpha_2,alpha_0)")
    print("same (H,K) would force alpha difference lambda*(1,-5,4)")
    print(
        "along that ambient kernel, a4 and S K change by 120*lambda; "
        "det[(H,K,a4) on alpha_(4,2,0)]=-360"
    )
    print(
        "exact merged-node audit: (H,K) cells = OCF-alpha cells = Walsh-stalk cells = "
        f"{hk_cells}; no nonzero kernel step is realized"
    )
    print("therefore (H,K) closes E[Delta K] through n=7, but only by finite n=7 rigidity")

    print("\nPOINTWISE ODD-CYCLE EDGE CURRENT")
    print("-" * 82)
    print(f"exact edge identities checked={edge_formula_checks}")
    print("Delta_e K=4(sum gained K(T-V(C))-sum lost K(T-V(C)))")
    print("affected cycles all share the flipped pair, so at most one occurs per OCF family")
    print("K(empty)=1/2 handles a spanning odd cycle without an exception")

    print("\nTOURNAMENT ANALYSIS / ASSUMPTION CHALLENGE")
    print("-" * 82)
    rows = all_rows[7]
    observers = [
        ("H", lambda row: row["H"]),
        ("(H,K)", lambda row: (row["H"], row["K"])),
        ("(H,K,a4)", lambda row: (row["H"], row["K"], row["b2"])),
        ("OCF alpha levels", lambda row: row["alpha"]),
        ("K-gradient multiset", lambda row: row["Kgrad"]),
    ]
    raw_values: list[Fraction] = []
    economy_values: list[Fraction] = []
    for label, key in observers:
        cells, max_fibre, _, separated = fibre_stats(rows, key, lambda row: row["bK"])
        economy = Fraction(separated, max(1, ceil(log2(cells))))
        raw_values.append(Fraction(separated))
        economy_values.append(economy)
        print(
            f"{label:24s}: cells={cells:3d}, max fibre={max_fibre:2d}, "
            f"separated={separated:5d}, economy={render_fraction(economy)}"
        )
    tie_order = list(range(len(observers)))
    raw = needle.tournament_fingerprint(raw_values, tie_order)
    economy = needle.tournament_fingerprint(economy_values, tie_order)
    edge_flips = len(raw["arcs"] ^ economy["arcs"]) // 2
    print(
        f"raw observer tournament: score={raw['score_histogram']}, C3={raw['directed_triangles']}, "
        f"SCC={raw['scc_sizes']}, HP={raw['hamiltonian_paths']}"
    )
    print(
        f"economy switch: score={economy['score_histogram']}, C3={economy['directed_triangles']}, "
        f"SCC={economy['scc_sizes']}, HP={economy['hamiltonian_paths']}, flips={edge_flips}"
    )
    print("vertices are proof-observer quotients, not runners, tournament vertices, or raw arcs")
    print("(H,K) ambiently forgets a third moment, but the n=7 realized locus retains it rigidly")
    print("edge-current carrier preserves endpoint-pair response but not labels after multiset quotient")

    print("\nFINITE-HIERARCHY BOUNDARY")
    print("-" * 82)
    print("A_T(x)=sum_j a_j x^j obeys S A=(1-x)[(1+x)A'-(n-1)A]")
    print("top coefficients form a finite tridiagonal Krylov ladder, not an all-n two-scalar law")
    print("n<=6: fixed alpha_(n-1)=1 leaves at most two nonempty f-levels; (H,K) closes algebraically")
    print("n=7: three nonempty levels occur; a third moment is ambiently independent but unrealized")
    print("the all-n hierarchy grows; the n=7 two-scalar closure is finite rigidity, not a theorem beyond n=7")
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
