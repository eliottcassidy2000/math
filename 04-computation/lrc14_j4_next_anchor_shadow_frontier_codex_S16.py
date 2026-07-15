#!/usr/bin/env python3
r"""Low-CPU set-system and frontier audit for the next THM-741 flood anchor.

The certified flood edges are 56,57,67.  They form a triangle A on the small
labels, and a completed family is already in their containment shadow exactly
when its small-label set contains an edge of A.  This script computes that
shadow by small-set cardinality, the marginal gain of each unresolved edge,
and exact root/E1/E2 workload data without entering the expensive v3 tree.

It also changes vertices from root edges to the thirteen currently unshadowed
four-small final bases.  Each such base P=H union K needs two external speeds
15<=a<b.  The j=2 Bonferroni threshold closes a,b above a finite V(P); below
it, exact G(P+a) data and THM-732 close b beyond a finite cap.  The remaining
finite candidate pairs are swept exactly.  Five flat-index quantiles in each
K-bank are independently replayed by the full interval-subtraction path and
serialized into a cross-check manifest.

Tournament Analysis is a proof-scheduling tournament, not an LRC quotient.
Its vertices are unresolved root edges; the pair observable is the
lexicographic vector (difference in three-small marginal gain, difference in
four-small gain, negative difference in exact E2-node count).  The gauge
prefers greater shadow gain, then smaller E2 count, with lexicographic edge
order only for a true tie.
This retains first-layer work and containment gain but destroys deeper
interval geometry and clearance margins.  For direct four-small closure the
correct alternate vertices are final small-label sets K, not runners, Fano
flags, or presentation-dependent root edges.
"""

from __future__ import annotations

import hashlib
import importlib.util
from fractions import Fraction as F
from itertools import combinations
from math import comb
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
O = frozenset(range(1, 5))
T = frozenset(range(5, 8))
ANCHORS = tuple(frozenset(edge) for edge in ANCHOR_OUTPUTS)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_core():
    require(sha256(CORE_PATH) == CORE_SHA256, "THM-741 core hash changed")
    spec = importlib.util.spec_from_file_location("thm741_next_anchor_dependency", CORE_PATH)
    require(spec is not None and spec.loader is not None, "cannot load THM-741 core")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def safe_comb(n: int, k: int) -> int:
    return comb(n, k) if 0 <= k <= n else 0


def body_root_work(core, edge: tuple[int, int]) -> dict[str, object]:
    body = tuple(sorted(H | frozenset(edge)))
    body_set = set(body)
    good, r, m = core.good_norm(body)
    threshold = 3 * m / (core.S2 * r)
    V1 = core.minV(4, threshold.numerator, threshold.denominator)
    E1 = E2 = 0
    max_V2 = 0
    for v1 in range(1, V1):
        if v1 in body_set:
            continue
        r1, m1, _ = core.subtract(good, v1)
        require(m1 > 0, f"empty E1 node edge={edge}, v1={v1}")
        E1 += 1
        threshold1 = 4 * m1 / (core.S2 * r1)
        V2 = core.minV(3, threshold1.numerator, threshold1.denominator)
        max_V2 = max(max_V2, V2)
        E2 += max(0, V2 - v1 - 1) - sum(v1 < speed < V2 for speed in body_set)
    require(E1 == V1 - 10, f"unexpected E1 count edge={edge}")
    return {"edge": edge, "r": r, "m": m, "V1": V1, "E1": E1, "E2": E2, "max_V2": max_V2}


def two_speed_frontier(core, K: tuple[int, ...]) -> dict[str, object]:
    base = tuple(sorted(H | frozenset(K)))
    require(len(base) == 11, f"K={K} is not a four-small base")
    good, r, m = core.good_norm(base)
    require(m > 0, f"empty four-small base K={K}")
    threshold = 5 * m / (core.S2 * r)
    V = core.minV(2, threshold.numerator, threshold.denominator)
    a_nodes = candidate_pairs = 0
    max_b = 0
    for a in range(15, V):
        r1, m1, _ = core.subtract(good, a)
        require(m1 > 0, f"empty one-external node K={K}, a={a}")
        a_nodes += 1
        cap = core.S2 * r1 / (6 * m1)
        bmax = cap.numerator // cap.denominator
        max_b = max(max_b, bmax)
        candidate_pairs += max(0, bmax - a)
    require(candidate_pairs > 0, f"empty candidate bank K={K}")

    sample_indices = {
        0,
        (candidate_pairs - 1) // 4,
        (candidate_pairs - 1) // 2,
        3 * (candidate_pairs - 1) // 4,
        candidate_pairs - 1,
    }
    require(len(sample_indices) == 5, f"collapsed sample indices K={K}")
    positive = zero = flat_index = 0
    minimum: tuple[F, int, int] | None = None
    failures = []
    sample_records = []
    for a in range(15, V):
        _, m1, good1 = core.subtract(good, a)
        require(m1 > 0, f"empty replay node K={K}, a={a}")
        r1 = len(good1)
        cap = core.S2 * r1 / (6 * m1)
        bmax = cap.numerator // cap.denominator
        for b in range(a + 1, bmax + 1):
            value = core.subtract_sparse(good1, b)
            if value > 0:
                positive += 1
                candidate = (value, a, b)
                if minimum is None or candidate < minimum:
                    minimum = candidate
            else:
                zero += 1
                failures.append(tuple(sorted(H | frozenset(K) | {a, b})))
            if flat_index in sample_indices:
                full_r, full_value, _ = core.subtract(good1, b)
                require(full_value == value, f"sparse/full mismatch K={K}, a={a}, b={b}")
                sample_records.append(
                    f"K={','.join(map(str, K))};index={flat_index}/{candidate_pairs};"
                    f"a={a};b={b};value={value};full_r={full_r}"
                )
            flat_index += 1
    require(flat_index == candidate_pairs, f"candidate replay count changed K={K}")
    require(len(sample_records) == 5, f"bad cross-check sample K={K}")
    require(positive + zero == candidate_pairs, f"bad sign ledger K={K}")
    return {
        "K": K,
        "r": r,
        "m": m,
        "V": V,
        "a_nodes": a_nodes,
        "candidate_pairs": candidate_pairs,
        "max_b": max_b,
        "positive": positive,
        "zero": zero,
        "minimum": minimum,
        "failures": tuple(failures),
        "sample_records": tuple(sample_records),
    }


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

    edges = tuple(combinations(range(1, 8), 2))
    unresolved = tuple(edge for edge in edges if frozenset(edge) not in ANCHORS)
    require(len(unresolved) == 18, "bad unresolved edge count")

    coverage = {}
    small_sets = {}
    for k in range(2, 7):
        sets = tuple(frozenset(K) for K in combinations(range(1, 8), k))
        shadowed = tuple(K for K in sets if any(anchor <= K for anchor in ANCHORS))
        unshadowed = tuple(K for K in sets if K not in shadowed)
        expected_unshadowed = safe_comb(4, k) + 3 * safe_comb(4, k - 1)
        require(len(unshadowed) == expected_unshadowed, f"bad unshadowed formula k={k}")
        coverage[k] = (len(sets), len(shadowed), len(unshadowed))
        small_sets[k] = (sets, shadowed, unshadowed)
    require(coverage == {2: (21, 3, 18), 3: (35, 13, 22), 4: (35, 22, 13), 5: (21, 18, 3), 6: (7, 7, 0)}, "bad coverage table")

    root_rows = []
    for edge in unresolved:
        edge_set = frozenset(edge)
        edge_type = "OO" if edge_set <= O else "OT"
        require(edge_type == "OO" or len(edge_set & O) == len(edge_set & T) == 1, "bad edge type")
        gains = {
            k: sum(edge_set <= K for K in small_sets[k][2])
            for k in range(2, 7)
        }
        expected = (
            {k: safe_comb(2, k - 2) + 3 * safe_comb(2, k - 3) for k in range(2, 7)}
            if edge_type == "OO"
            else {k: safe_comb(3, k - 2) for k in range(2, 7)}
        )
        require(gains == expected, f"bad marginal formula edge={edge}")
        row = body_root_work(core, edge)
        row.update({"type": edge_type, "gains": gains})
        root_rows.append(row)

    four_unshadowed = tuple(tuple(sorted(K)) for K in small_sets[4][2])
    tail_rows = [two_speed_frontier(core, K) for K in four_unshadowed]
    tail_by_K = {row["K"]: row for row in tail_rows}
    total_pairs = sum(int(row["candidate_pairs"]) for row in tail_rows)
    require(total_pairs == 29183, f"unexpected two-speed frontier {total_pairs}")
    total_positive = sum(int(row["positive"]) for row in tail_rows)
    total_zero = sum(int(row["zero"]) for row in tail_rows)
    failures = tuple(family for row in tail_rows for family in row["failures"])
    sample_records = tuple(record for row in tail_rows for record in row["sample_records"])
    require(len(sample_records) == 65, f"bad cross-check manifest size {len(sample_records)}")
    manifest_payload = "\n".join(sample_records).encode()
    manifest_sha256 = hashlib.sha256(manifest_payload).hexdigest()
    global_minimum = min(
        (row["minimum"][0], row["K"], row["minimum"][1], row["minimum"][2])
        for row in tail_rows
        if row["minimum"] is not None
    )

    for row in root_rows:
        edge_set = frozenset(row["edge"])
        newly_shadowed = tuple(K for K in four_unshadowed if edge_set <= frozenset(K))
        row["remaining_K4"] = len(four_unshadowed) - len(newly_shadowed)
        row["remaining_pairs"] = total_pairs - sum(
            int(tail_by_K[K]["candidate_pairs"]) for K in newly_shadowed
        )

    priority = sorted(
        root_rows,
        key=lambda row: (
            -int(row["gains"][3]),
            -int(row["gains"][4]),
            int(row["E2"]),
            row["edge"],
        ),
    )
    priority_path = tuple(row["edge"] for row in priority)
    priority_rank = {edge: index for index, edge in enumerate(priority_path)}
    lex_edges = tuple(sorted(unresolved))
    flips = sum(
        priority_rank[left] > priority_rank[right]
        for i, left in enumerate(lex_edges)
        for right in lex_edges[i + 1 :]
    )
    pareto = tuple(
        row["edge"]
        for row in root_rows
        if not any(
            other is not row
            and int(other["gains"][3]) >= int(row["gains"][3])
            and int(other["E2"]) <= int(row["E2"])
            and (
                int(other["gains"][3]) > int(row["gains"][3])
                or int(other["E2"]) < int(row["E2"])
            )
            for other in root_rows
        )
    )
    require(set(pareto) == {(3, 4), (4, 7)}, f"bad Pareto frontier {pareto}")

    matching_rows = []
    for left, right in (((1, 2), (3, 4)), ((1, 3), (2, 4)), ((1, 4), (2, 3))):
        left_row = next(row for row in root_rows if row["edge"] == left)
        right_row = next(row for row in root_rows if row["edge"] == right)
        matching_rows.append((left, right, int(left_row["E2"]) + int(right_row["E2"])))
    require(min(matching_rows, key=lambda row: row[2])[:2] == ((1, 3), (2, 4)), "bad matching optimum")

    best = priority[0]
    require(best["edge"] == (3, 4), "unexpected coverage-first recommendation")

    print("THM-741 NEXT-ANCHOR SHADOW AND TWO-SPEED FRONTIER AUDIT")
    print("=" * 92)
    print(f"dependency_sha256={CORE_SHA256}")
    print(
        "anchor_output_sha256="
        + ",".join(f"{edge}:{digest}" for edge, (_, digest) in ANCHOR_OUTPUTS.items())
    )
    print("certified anchor graph: triangle edges=((5,6),(5,7),(6,7)); vertices={5,6,7}")
    print("shadow criterion: completed small set K is closed by containment iff K contains an anchor edge")
    print("coverage k ; total ; shadowed ; unshadowed")
    for k in range(2, 7):
        print(f"  {k} ; {coverage[k][0]} ; {coverage[k][1]} ; {coverage[k][2]}")
    print("marginal formula OO: C(2,k-2)+3C(2,k-3) -> (1,5,7,3,0)")
    print("marginal formula OT: C(3,k-2) -> (1,3,3,1,0)")
    print("unresolved edge ; type ; gain(k2..k6) ; r ; m ; V1 ; E1 ; E2 ; maxV2 ; post-edge K4/pairs")
    for row in priority:
        gains = tuple(row["gains"][k] for k in range(2, 7))
        print(
            f"  {row['edge']} ; {row['type']} ; {gains} ; {row['r']} ; {row['m']} ; "
            f"{row['V1']} ; {row['E1']} ; {row['E2']} ; {row['max_V2']} ; "
            f"{row['remaining_K4']}/{row['remaining_pairs']}"
        )
    print(f"single-edge Pareto frontier (gain3 high,E2 low)={pareto}")
    print("perfect-matching OO pairs ; combined exact E2 nodes")
    for left, right, work in matching_rows:
        print(f"  {left}+{right} ; {work}")
    print("current unshadowed four-small bases ; r ; m ; V(j2) ; a-nodes ; pairs ; positive/zero ; minimum ; max-b")
    for row in tail_rows:
        minimum = row["minimum"]
        require(minimum is not None, f"no positive minimum K={row['K']}")
        print(
            f"  {row['K']} ; {row['r']} ; {row['m']} ; {row['V']} ; "
            f"{row['a_nodes']} ; {row['candidate_pairs']} ; {row['positive']}/{row['zero']} ; "
            f"{minimum[0]}@({minimum[1]},{minimum[2]}) ; {row['max_b']}"
        )
    print(
        f"two-speed reduction: bases={len(tail_rows)} a_nodes={sum(int(row['a_nodes']) for row in tail_rows)} "
        f"finite_exact_pair_candidates={total_pairs} empty_intermediate_nodes=0"
    )
    print(
        f"exact pair evaluation: positive={total_positive} zero={total_zero} "
        f"covering_failures={len(failures)} global_minimum={global_minimum[0]} "
        f"at K={global_minimum[1]},a={global_minimum[2]},b={global_minimum[3]}"
    )
    print("sparse/full cross-check rule: flat indices 0,1/4,1/2,3/4,last in every K-bank")
    print("sparse/full cross-check manifest records:")
    for record in sample_records:
        print(f"  {record}")
    print(f"crosscheck_samples={len(sample_records)} mismatches=0 manifest_sha256={manifest_sha256}")
    require(total_positive == total_pairs and total_zero == 0 and not failures, "four-small exact closure failed")
    print("PROVED FOUR-SMALL SHADOW: every exact candidate pair is positive")
    print("scope: three whole flood bodies are closed; eighteen whole bodies still remain")
    print("Tournament Analysis vertices: 18 unresolved root edges (proof scheduler only)")
    print("pair observable: lex(delta gain3, delta gain4, -delta E2); gauge: higher gain then lower E2; tie: lex edge")
    print(
        f"fingerprint: score_hist={{0..17:1}}, directed_3cycles=0, SCC_sizes=18x1, "
        f"edge_flips_vs_lex={flips}, Hamiltonian_paths=1"
    )
    print(f"tie Hamiltonian path={priority_path}")
    print("kept: containment gain and exact E2 count; destroyed: v3 geometry, bottom margins, final-family identity")
    print("challenged vertices: use 22 unshadowed three-small bases next; runners/Fano flags/root charts lose that quotient")
    print(
        "RECOMMENDATION: if the next job must close a whole edge, choose (3,4): "
        "maximum three-small gain 5 and minimum E2=62282 among all gain-5 edges"
    )
    print("STRATEGY: the direct four-small bank is closed; reserve a full (3,4) run for the <=3-small frontier")
    print(f"source_sha256={sha256(Path(__file__).resolve())}")
    print("ALL EXACT FRONTIER CHECKS PASSED")


if __name__ == "__main__":
    main()
