#!/usr/bin/env python3
"""S140: try lifting C27 marked transfers into q=3 unital 4-blocks.

The question is whether the HYP-2937/HYP-2940 marked C=27 transfer packets can
be interpreted as blocks in the q=3 unital design 2-(28,4,1), after the
AP/Goddyn-Wong labels are attached.

This script tests three deliberately explicit models:

1. residue-pair blocks:
   H[a] -> D[d] is the four C27 residue points {a,27-a,d,27-d}.
2. global AP/GW labelled shell blocks:
   H[a] -> D[d] is {AP,GW,H_a,D_d}.
3. branch-local unital charts:
   a branch may embed a compatible family of residue-pair blocks into the
   Hermitian q=3 unital if its desired blocks do not repeat a pair.

The result is a useful negative/positive split.  One q=3 unital cannot contain
both the GW block H[12]->D[3] and the K33 near-miss block H[12]->D[9] under the
residue-pair model because they share the C27 hole pair {12,15}.  But each
branch separately embeds, and the S138 two-hole rows lift as two-block splices.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from pathlib import Path
import sys


REPO = Path(__file__).resolve().parents[1]
C = 27


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


s105 = load_module(
    "s105_design_hodge_for_s140",
    REPO / "04-computation" / "lrc14_design_hodge_carriers_codex_s105.py",
)


@dataclass(frozen=True)
class Transfer:
    name: str
    hole_shell: int
    double_shell: int
    branch: str

    @property
    def residue_block(self) -> tuple[int, int, int, int]:
        pts = {self.hole_shell, C - self.hole_shell, self.double_shell, C - self.double_shell}
        return tuple(sorted(pts))

    @property
    def global_label_block(self) -> tuple[str, str, str, str]:
        return tuple(sorted(("AP", "GW", f"H{self.hole_shell}", f"D{self.double_shell}")))


TRANSFERS = {
    "GW": Transfer("GW H12->D3", 12, 3, "tight-floor"),
    "K33": Transfer("K33 H12->D9", 12, 9, "K33-unit-excess"),
    "P10": Transfer("petal H10->D7", 10, 7, "C27-petal-two-block"),
    "P13": Transfer("petal H13->D1", 13, 1, "C27-petal-two-block"),
}

BRANCHES = {
    "tight chart: GW + both petals": ("GW", "P10", "P13"),
    "K33 chart: near-miss + both petals": ("K33", "P10", "P13"),
    "global low frontier: GW + K33 + petals": ("GW", "K33", "P10", "P13"),
    "S138 splice drop(10,12)->add(20,24)": ("P10", "GW"),
    "S138 splice drop(10,12)->add(20,36)": ("P10", "K33"),
    "12-branch superposition only": ("GW", "K33"),
}


def pair_multiplicities(blocks: list[tuple[object, ...]]) -> Counter[tuple[object, object]]:
    counts: Counter[tuple[object, object]] = Counter()
    for block in blocks:
        for a, b in combinations(block, 2):
            counts[tuple(sorted((a, b), key=str))] += 1
    return counts


def repeated_pairs(blocks: list[tuple[object, ...]]) -> dict[tuple[object, object], int]:
    return {pair: n for pair, n in sorted(pair_multiplicities(blocks).items(), key=lambda kv: (str(kv[0]), kv[1])) if n > 1}


def unital_design():
    points, blocks, sizes = s105.hermitian_unital_q3()
    block_sets = [frozenset(block) for block in blocks]
    return points, block_sets, sizes


def block_intersection_hist(blocks: list[frozenset[int]]) -> Counter[int]:
    hist: Counter[int] = Counter()
    for a, b in combinations(blocks, 2):
        hist[len(a & b)] += 1
    return hist


def find_unital_block_chart(desired_blocks: list[tuple[int, ...]], unital_blocks: list[frozenset[int]]):
    """Find unital blocks with the same intersection matrix as desired_blocks.

    For the branch-local rows used here the desired blocks are disjoint, but the
    implementation keeps the actual intersection profile visible.
    """

    desired_intersections = {}
    for i, j in combinations(range(len(desired_blocks)), 2):
        desired_intersections[(i, j)] = len(set(desired_blocks[i]) & set(desired_blocks[j]))

    chosen: list[int] = []

    def rec(k: int) -> list[int] | None:
        if k == len(desired_blocks):
            return chosen.copy()
        for idx, block in enumerate(unital_blocks):
            if idx in chosen:
                continue
            ok = True
            for prev_pos, prev_idx in enumerate(chosen):
                need = desired_intersections[(prev_pos, k)]
                got = len(block & unital_blocks[prev_idx])
                if got != need:
                    ok = False
                    break
            if not ok:
                continue
            chosen.append(idx)
            out = rec(k + 1)
            if out is not None:
                return out
            chosen.pop()
        return None

    return rec(0)


def print_assumption_challenge() -> None:
    print("[0] Assumption challenge / quotient declaration")
    print("  considered vertices:")
    print("    C27 residues, antipodal shells, hole/double statuses, AP/GW labels,")
    print("    q=3 unital points, unital 4-blocks, branch-local charts, and")
    print("    S138 two-hole proof obligations.")
    print("  chosen test predicate:")
    print("    whether a marked transfer packet can be a q=3 unital 4-block")
    print("    without violating lambda=1 pair uniqueness.")
    print("  quotient preserves:")
    print("    C27 hole/double shell pair, branch label, and unital pair incidence.")
    print("  quotient destroys:")
    print("    exact time geometry, within-shell orientation beyond the antipodal")
    print("    pair, and any global ordering of unital points.")
    print("  challenged assumption:")
    print("    all marked transfers can live in one unital chart.  The test shows")
    print("    branch-local charts work, but the GW/K33 superposition repeats a pair.")


def print_unital_sanity(blocks: list[frozenset[int]], sizes: list[int]) -> None:
    print()
    print("[1] q=3 unital design sanity")
    print(f"  points=28 blocks={len(blocks)} block_size_hist={dict(sorted(Counter(map(len, blocks)).items()))}")
    print(f"  projective-line intersection sizes={dict(sorted(Counter(sizes).items()))}")
    print(f"  block intersection histogram={dict(sorted(block_intersection_hist(blocks).items()))}")
    print("  readout: two distinct unital blocks never share two points.")


def print_transfer_blocks() -> None:
    print()
    print("[2] C27 transfer packets as residue-pair 4-blocks")
    for key, transfer in TRANSFERS.items():
        print(
            f"  {key:3s} {transfer.name:18s} branch={transfer.branch:20s} "
            f"block={transfer.residue_block}"
        )
    all_blocks = [TRANSFERS[k].residue_block for k in ("GW", "K33", "P10", "P13")]
    repeats = repeated_pairs(all_blocks)
    print()
    print(f"  simultaneous low-frontier repeated pairs={repeats}")
    print("  obstruction:")
    print("    GW and K33 both contain the hole pair (12,15).  Since a q=3")
    print("    unital has lambda=1, those two packets cannot both be blocks in")
    print("    one unital chart under the raw residue-pair lift.")


def print_global_label_model() -> None:
    print()
    print("[3] Global AP/GW-label block model")
    blocks = [TRANSFERS[k].global_label_block for k in ("GW", "K33", "P10", "P13")]
    for key, block in zip(("GW", "K33", "P10", "P13"), blocks):
        print(f"  {key:3s} block={block}")
    repeats = repeated_pairs(blocks)
    print(f"  repeated pairs={repeats}")
    print("  readout:")
    print("    treating AP and GW as two global points in every transfer block is")
    print("    worse: the pair (AP,GW) repeats immediately.  AP/GW labels must be")
    print("    chart labels or branch colors, not universal block vertices.")


def print_branch_charts(unital_blocks: list[frozenset[int]]) -> None:
    print()
    print("[4] Branch-local q=3 unital chart test")
    for name, keys in BRANCHES.items():
        desired = [TRANSFERS[k].residue_block for k in keys]
        repeats = repeated_pairs(desired)
        chart = None if repeats else find_unital_block_chart(desired, unital_blocks)
        status = "FAIL(pair-repeat)" if repeats else ("EMBEDS" if chart is not None else "FAIL(no chart)")
        print(f"  {name:45s} status={status}")
        print(f"    desired_blocks={desired}")
        if repeats:
            print(f"    repeated_pairs={repeats}")
        if chart is not None:
            print(f"    example_unital_block_indices={chart}")
            print(f"    example_unital_blocks={[tuple(sorted(unital_blocks[i])) for i in chart]}")
    print()
    print("  splice readout:")
    print("    The two S138 genuine two-hole rows are not new single unital blocks.")
    print("    Each is a two-block path: petal H10->D7 plus either the GW block")
    print("    H12->D3 or the K33 near-miss block H12->D9.")


def print_tournament_analysis() -> None:
    print()
    print("[5] Tournament Analysis on lift obligations")
    carriers = [
        ("exact M/Farey branch", (6, 6, 5, 6, 6)),
        ("C27 marked transfer", (5, 6, 6, 5, 6)),
        ("unital pair-repeat obstruction", (5, 5, 6, 6, 6)),
        ("branch-local q3 unital chart", (4, 5, 5, 5, 5)),
        ("S138 two-block splice path", (4, 5, 5, 4, 5)),
        ("global AP/GW vertex model", (1, 2, 2, 2, 5)),
        ("raw unital analogy", (0, 1, 1, 1, 1)),
    ]
    wins = Counter()
    c3 = 0
    for i, j in combinations(range(len(carriers)), 2):
        if carriers[i][1] >= carriers[j][1]:
            wins[i] += 1
        else:
            wins[j] += 1
    for i, j, k in combinations(range(len(carriers)), 3):
        def beats(a: int, b: int) -> bool:
            return carriers[a][1] >= carriers[b][1]

        if (beats(i, j) and beats(j, k) and beats(k, i)) or (
            beats(i, k) and beats(k, j) and beats(j, i)
        ):
            c3 += 1
    hist = Counter(wins[i] for i in range(len(carriers)))
    order = sorted(range(len(carriers)), key=lambda i: carriers[i][1], reverse=True)
    print("  vertices: lift obligations/proof carriers, not runners.")
    print("  pair observable, ordered lexicographically:")
    print("    theorem-scale retention, C27 predicate retention, lambda=1 incidence,")
    print("    finite checkability, anti-scalar guard.")
    print(f"  fingerprint: score_hist={dict(sorted(hist.items()))} c3={c3} hp=1")
    print("  Hamiltonian path:")
    print("    " + " > ".join(carriers[i][0] for i in order))


def print_verdict() -> None:
    print()
    print("[6] Verdict")
    print("  Negative global result:")
    print("    the raw C27 residue-pair lift cannot put GW H12->D3 and K33")
    print("    H12->D9 into the same q=3 unital, because the hole pair (12,15)")
    print("    would occur in two different blocks.")
    print("  Positive local result:")
    print("    each branch separately embeds into the q=3 unital incidence design,")
    print("    and the S138 genuine two-hole rows lift as two-block splices.")
    print("  Proof-use refinement:")
    print("    use the unital as a branch-local pair-unique chart or splice grammar,")
    print("    not as a universal atlas for all C27 marked transfers at once.")


def main() -> None:
    points, unital_blocks, sizes = unital_design()
    _ = points
    print("S140 LRC14 C27 -> q=3 UNITAL BLOCK-LIFT AUDIT")
    print("=" * 78)
    print_assumption_challenge()
    print_unital_sanity(unital_blocks, sizes)
    print_transfer_blocks()
    print_global_label_model()
    print_branch_charts(unital_blocks)
    print_tournament_analysis()
    print_verdict()


if __name__ == "__main__":
    main()
