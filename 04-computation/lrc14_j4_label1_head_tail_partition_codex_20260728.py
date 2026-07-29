#!/usr/bin/env python3
r"""Exact head/tail closure for the six label-1 THM-741 flood roots.

Let

    H = {8,9,10,11,12,13,14},
    E_b = {1,b} union H,                  2 <= b <= 7,
    G_b = G(E_b),
    c_b(w) = |G_b intersect D_w|.

This companion treats the pure tail in which every added speed is at least
15.  If one of the four added speeds is at most 14, the resulting family has
at least three labels in `{1,...,7}` and is already in THM-741's proved
``Three-small exact closure`` addendum.  For the pure tail, this verifier
fixes W=2500, evaluates
`c_b(w)` for every `15<=w<=W`, and partitions the four added speeds by the
number above W.

For each finite coverage, four exact paths agree:

1. sparse subtraction from the root carrier;
2. full subtraction from the root carrier;
3. a direct ten-comb union;
4. a two-pointer incidence sum over the teeth of D_w.

The tail cap inherited from THM-735(ii), THM-731, and THM-732 is

    c_b(w) < u_b(w)
           := m_b/7 + (99/70) r_b/(7w).

If at least two speeds exceed W, the top two finite single coverages plus
`2u_b(W+1)` leave a positive margin.  The verifier also checks
`q_2>u_b(W+1)`.  Hence, with three tails,
`q_1+3u_b(W+1)<q_1+q_2+2u_b(W+1)`, and with four tails,
`4u_b(W+1)<q_1+q_2+2u_b(W+1)`.

If exactly one speed exceeds W, the root union bound closes every finite head
triple except those with

    c_b(a)+c_b(c)+c_b(d)+u_b(W+1) >= m_b.

For every such dangerous triple S, the exact residual carrier `G_b\D_S` has
measure m_S and r_S components and satisfies

    6m_S/7 - (99/70) r_S/(7(W+1)) > 0.

Thus its final tail speed is harmless.

If no speed exceeds W, the root union bound closes every finite quadruple
whose four individual coverages sum to less than m_b.  Every remaining
quadruple is checked by both nested subtraction and direct thirteen-comb
union, and has positive survivor.

Candidate sets are independently recovered in two ways: specialized
descending-ranked loops and a separate arity-generic upper-envelope recursion.
The expected candidate counts are load-bearing hostile controls:
each row has both dangerous triples and all-head exceptions, so the script
cannot silently replace the repaired partition by the already-failed top-four
union envelope.

This proves the pure tails of the six literal flood bodies 12,13,14,15,16,17.
Together with THM-741's earlier three-small completed-family addendum for
their non-pure branches, it closes those six whole flood bodies.  Composed
with the independently certified other fifteen literal flood bodies, it gives
21/21 flood bodies: equivalently, the THM-741 stratum of thirteen-speed
families containing H and at least two labels from `{1,...,7}`.  It does not
cover families containing H with only zero or one such label, does not promote
global THM-741, and does not prove LRC(14).
"""

from __future__ import annotations

import hashlib
import importlib.util
import sys
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CORE_PATH = ROOT / "04-computation/lrc14_thm741_2002_body_j4_tree_kps_S128c5.py"
CORE_SHA256 = "5aa81d9d78273c8f9e3e7a6574091a3bc3f64ab6086c7024c15f9420c99dac96"
H = tuple(range(8, 15))
FIRST_EXTERNAL = 15
W = 2500
TAIL_FIRST = W + 1
EXPECTED = {
    (1, 2): {
        "r": 32,
        "m": F(319927, 2522520),
        "top2": ((23, F(453587, 11603592)), (46, F(212279, 5801796))),
        "tail_cap": F(182859895, 8832351528),
        "multi_tail_margin": F(4947963661, 507860212860),
        "top4_margin": -F(3653017, 275585310),
        "dangerous_speeds": ((19, 23, 46),),
        "triple_minimum": (
            F(71621894273, 1608224007390),
            (19, 23, 46),
            F(29780131, 551170620),
            22,
        ),
        "quad_count": 2589,
        "quad_minimum": (
            F(17069726431, 693372639960),
            (17, 19, 23, 37),
            16,
        ),
        "coverage_manifest": "5011d9d975e9d6858e2f54aa4f2b3736fc0d72b330021f5f748ee9d539dff0f5",
        "triple_manifest": "b407e9ed61b369129cc81bdf2c0fd44e7336d6a72efd91366be5f4ffce3008dd",
        "quad_manifest": "2e2574f8f0eb18958feb97f1da1ebdaadbc695bc51f262ab1b695494dd2706bc",
    },
    (1, 3): {
        "r": 30,
        "m": F(3319, 25740),
        "top2": ((23, F(453587, 11603592)), (17, F(67493, 1786785))),
        "tail_cap": F(65750513, 3154411260),
        "multi_tail_margin": F(25633591579, 2466749605320),
        "top4_margin": -F(25012943, 1249320072),
        "dangerous_speeds": (
            (17, 19, 23),
            (17, 19, 46),
            (17, 21, 23),
            (17, 23, 46),
            (19, 21, 23),
            (19, 23, 46),
        ),
        "triple_minimum": (
            F(34858373223, 959291513180),
            (17, 19, 23),
            F(958105, 21917896),
            14,
        ),
        "quad_count": 7609,
        "quad_minimum": (F(48080163, 2027405380), (17, 19, 23, 37), 16),
        "coverage_manifest": "b75802caecacca8b43a4b01040e80ba79ff18c9178ea1ed60cee6228dba283dd",
        "triple_manifest": "c9864972896330c6b7218bdda6ba1d5ba400d25a4de5b2b7be6b1486025bf7fc",
        "quad_manifest": "352972a15489feee86199beba567a1d73cfc0d33823efcea11434d047c71dd22",
    },
    (1, 4): {
        "r": 30,
        "m": F(335047, 2522520),
        "top2": ((23, F(240307, 5801796)), (17, F(135481, 3573570))),
        "tail_cap": F(944979467, 44161757640),
        "multi_tail_margin": F(61555969441, 5755749079080),
        "top4_margin": -F(50908597, 3747960216),
        "dangerous_speeds": ((17, 19, 23),),
        "triple_minimum": (
            F(1865084760919, 54679616251260),
            (17, 19, 23),
            F(154798927, 3747960216),
            16,
        ),
        "quad_count": 2581,
        "quad_minimum": (F(162578519, 7415750160), (19, 23, 32, 37), 16),
        "coverage_manifest": "365498ca72dec323a18c511d18cdaa845f41e5846c97b3e6530a335d0342cc36",
        "triple_manifest": "ceccc28c4d61439997cf3ba87cf1d0234e716f47c54ccea4ff8b9d6ca3fa35c6",
        "quad_manifest": "9e02f456eda2f955a3a3cb0c462011cf6ae498b77d319bf8b9dbd1b4ab93e23d",
    },
    (1, 5): {
        "r": 26,
        "m": F(6716, 45045),
        "top2": ((23, F(41075, 892584)), (19, F(36649, 798798))),
        "tail_cap": F(527231, 22531509),
        "multi_tail_margin": F(9555282629, 918985147080),
        "top4_margin": -F(68144053, 3747960216),
        "dangerous_speeds": (
            (17, 19, 23),
            (19, 23, 25),
            (19, 23, 38),
            (19, 23, 46),
        ),
        "triple_minimum": (
            F(15348916166, 365505456225),
            (19, 23, 25),
            F(6349621, 125266050),
            18,
        ),
        "quad_count": 7290,
        "quad_minimum": (F(544970477, 18539375400), (19, 23, 25, 37), 16),
        "coverage_manifest": "7794967988d1961970f1c796c3d955e85ab4ef0cb2cc712621fe32d69afc1d43",
        "triple_manifest": "fb061c2d737aaedfb28ba26e643f7632efba1df27b82d8b15acecb844d28c168",
        "quad_manifest": "7f34243d912fb4bf39a28a3c511377b0ccdaccb445778bed893bfb38a05a7037",
    },
    (1, 6): {
        "r": 30,
        "m": F(365567, 2522520),
        "top2": ((19, F(377149, 7987980)), (23, F(47539, 1054872))),
        "tail_cap": F(1021309987, 44161757640),
        "multi_tail_margin": F(15408825262, 2412336011085),
        "top4_margin": -F(6867281, 312330018),
        "dangerous_speeds": (
            (17, 19, 23),
            (17, 19, 46),
            (19, 21, 23),
            (19, 23, 25),
            (19, 23, 38),
            (19, 23, 46),
        ),
        "triple_minimum": (
            F(156275963542, 4556634687605),
            (17, 19, 23),
            F(6484031, 156165009),
            16,
        ),
        "quad_count": 7522,
        "quad_minimum": (F(1707671939, 77041404440), (17, 19, 23, 37), 16),
        "coverage_manifest": "a67a83cd28f9f06d07d52f70944506cea9694cfdba9b70cc52a7c515036998ce",
        "triple_manifest": "ec6715015f524723afb2c7af2efe1de315b59ab3ad40c1648b8d1041f2c5155b",
        "quad_manifest": "4928a466aa690a73ed1bb8f1a2f6242e7264d917b8bc505fb46fc48e8fcf15f3",
    },
    (1, 7): {
        "r": 22,
        "m": F(384011, 2522520),
        "top2": ((19, F(128503, 2662660)), (23, F(25331, 552552))),
        "tail_cap": F(1038897919, 44161757640),
        "multi_tail_margin": F(106901554667, 9649344044340),
        "top4_margin": -F(228625, 11639628),
        "dangerous_speeds": (
            (17, 19, 23),
            (19, 23, 38),
            (19, 23, 46),
        ),
        "triple_minimum": (
            F(12865251479, 300437451930),
            (17, 19, 23),
            F(1052099, 20593188),
            12,
        ),
        "quad_count": 5087,
        "quad_minimum": (F(13420153, 455992020), (17, 19, 23, 31), 16),
        "coverage_manifest": "d7683b3db66169e0da59febc8c949994b31be3d72c04f3ba9458fa2f54ea88e9",
        "triple_manifest": "0e150e74a2b1ceab8b5d1cf82e6a58a82f3e65dd95b659d1cf12cc66b977e75c",
        "quad_manifest": "3045929547d2a7cc78a1226e588916a30f5e9662cae19a672942265658ecd61a",
    },
}
EXPECTED_COMBINED_MANIFEST = "556e997f0fb7b270943611ab560ce72a32f077c97ed3e6354da8d6176d2c76e4"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_core():
    require(sha256(CORE_PATH) == CORE_SHA256, "THM-741 core hash changed")
    spec = importlib.util.spec_from_file_location("thm741_label1_partition_dependency", CORE_PATH)
    require(spec is not None and spec.loader is not None, "cannot load THM-741 core")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def fraction_text(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def direct_tooth_coverage(good: list[tuple[F, F]], w: int) -> F:
    """Independent intersection sum over the disjoint teeth of D_w."""
    require(w >= 1, "speed must be positive")
    radius = F(1, 14 * w)
    total = F(0)
    first_component = 0
    for k in range(w + 1):
        center = F(k, w)
        left = max(F(0), center - radius)
        right = min(F(1), center + radius)
        while first_component < len(good) and good[first_component][1] <= left:
            first_component += 1
        index = first_component
        while index < len(good) and good[index][0] < right:
            overlap_left = max(left, good[index][0])
            overlap_right = min(right, good[index][1])
            if overlap_left < overlap_right:
                total += overlap_right - overlap_left
            if good[index][1] > right:
                break
            index += 1
    return total


def exact_coverage_bank(core, edge: tuple[int, int]) -> dict[str, object]:
    expected = EXPECTED[edge]
    body = tuple(sorted((*edge, *H)))
    require(len(body) == 9 and len(set(body)) == 9, f"bad root body {edge}")
    good, root_r, root_m = core.good_norm(body)
    require(root_m > 0 and root_r == len(good), f"bad root carrier {edge}")
    require((root_r, root_m) == (expected["r"], expected["m"]), f"root summary changed {edge}")

    rows: list[tuple[F, int, int]] = []
    coverage_by_speed: dict[int, F] = {}
    manifest_rows: list[str] = []
    for w in range(FIRST_EXTERNAL, W + 1):
        sparse_survivor = core.subtract_sparse(good, w)
        full_r, full_survivor, full_good = core.subtract(good, w)
        direct_good, direct_r, direct_survivor = core.good_norm((*body, w))
        coverage = root_m - sparse_survivor
        tooth_coverage = direct_tooth_coverage(good, w)
        require(sparse_survivor == full_survivor, f"sparse/full mismatch {edge},w={w}")
        require(full_survivor == direct_survivor, f"subtract/union mismatch {edge},w={w}")
        require(full_good == direct_good and full_r == direct_r, f"carrier mismatch {edge},w={w}")
        require(coverage == tooth_coverage, f"subtract/tooth mismatch {edge},w={w}")
        require(F(0) <= coverage < root_m, f"bad coverage {edge},w={w}")
        rows.append((coverage, w, full_r))
        coverage_by_speed[w] = coverage
        manifest_rows.append(
            f"{edge[0]}{edge[1]}:{w}:{fraction_text(coverage)}:{fraction_text(full_survivor)}:{full_r}"
        )

    require(len(rows) == W - FIRST_EXTERNAL + 1, f"bad finite coverage census {edge}")
    ranked = sorted(rows, key=lambda row: (-row[0], row[1]))

    def tail_cap(w: int) -> F:
        return root_m / 7 + core.S2 * root_r / (7 * w)

    tail_first_cap = tail_cap(TAIL_FIRST)
    q1, q2 = ranked[0][0], ranked[1][0]
    require(q2 > tail_first_cap, f"second finite rank does not dominate tail cap {edge}")
    multi_tail_margin = root_m - q1 - q2 - 2 * tail_first_cap
    require(multi_tail_margin > 0, f"two-tail envelope failed {edge}")
    top4_margin = root_m - sum((row[0] for row in ranked[:4]), F(0))
    require(top4_margin <= 0, f"hostile top-four control unexpectedly closed {edge}")
    observed_top2 = tuple((row[1], row[0]) for row in ranked[:2])
    coverage_manifest = hashlib.sha256("\n".join(manifest_rows).encode()).hexdigest()
    require(observed_top2 == expected["top2"], f"top two changed {edge}")
    require(tail_first_cap == expected["tail_cap"], f"tail cap changed {edge}")
    require(
        multi_tail_margin == expected["multi_tail_margin"],
        f"multi-tail margin changed {edge}",
    )
    require(top4_margin == expected["top4_margin"], f"top-four hostile changed {edge}")
    require(
        coverage_manifest == expected["coverage_manifest"],
        f"coverage manifest changed {edge}",
    )

    return {
        "edge": edge,
        "body": body,
        "good": good,
        "r": root_r,
        "m": root_m,
        "rows": ranked,
        "coverage_by_speed": coverage_by_speed,
        "tail_first_cap": tail_first_cap,
        "multi_tail_margin": multi_tail_margin,
        "top4_margin": top4_margin,
        "coverage_manifest": coverage_manifest,
    }


def ranked_triples(
    ranked: list[tuple[F, int, int]],
    threshold: F,
) -> set[tuple[int, int, int]]:
    """Descending branch-and-bound enumeration of sums at least threshold."""
    result: set[tuple[int, int, int]] = set()
    n = len(ranked)
    for i in range(n - 2):
        qi = ranked[i][0]
        if qi + ranked[i + 1][0] + ranked[i + 2][0] < threshold:
            break
        for j in range(i + 1, n - 1):
            qij = qi + ranked[j][0]
            if qij + ranked[j + 1][0] < threshold:
                break
            for k in range(j + 1, n):
                if qij + ranked[k][0] < threshold:
                    break
                result.add(tuple(sorted((ranked[i][1], ranked[j][1], ranked[k][1]))))
    return result


def ranked_quadruples(
    ranked: list[tuple[F, int, int]],
    threshold: F,
) -> set[tuple[int, int, int, int]]:
    """Descending branch-and-bound enumeration of sums at least threshold."""
    result: set[tuple[int, int, int, int]] = set()
    n = len(ranked)
    for i in range(n - 3):
        qi = ranked[i][0]
        if qi + ranked[i + 1][0] + ranked[i + 2][0] + ranked[i + 3][0] < threshold:
            break
        for j in range(i + 1, n - 2):
            qij = qi + ranked[j][0]
            if qij + ranked[j + 1][0] + ranked[j + 2][0] < threshold:
                break
            for k in range(j + 1, n - 1):
                qijk = qij + ranked[k][0]
                if qijk + ranked[k + 1][0] < threshold:
                    break
                for ell in range(k + 1, n):
                    if qijk + ranked[ell][0] < threshold:
                        break
                    result.add(
                        tuple(
                            sorted(
                                (
                                    ranked[i][1],
                                    ranked[j][1],
                                    ranked[k][1],
                                    ranked[ell][1],
                                )
                            )
                        )
                    )
    return result


def recursive_upper_envelope(
    ranked: list[tuple[F, int, int]],
    arity: int,
    threshold: F,
) -> tuple[set[tuple[int, ...]], int]:
    """Independent generic recursion with exact suffix upper envelopes."""
    require(arity in (3, 4), "only triple and quadruple controls are used")
    result: set[tuple[int, ...]] = set()
    prefix_sums = [F(0)]
    for row in ranked:
        prefix_sums.append(prefix_sums[-1] + row[0])
    visited = 0

    def best_from(start: int, count: int) -> F:
        if count == 0:
            return F(0)
        if start + count > len(ranked):
            return F(-1)
        return prefix_sums[start + count] - prefix_sums[start]

    def visit(start: int, need: int, total: F, speeds: tuple[int, ...]) -> None:
        nonlocal visited
        visited += 1
        if need == 0:
            require(total >= threshold, "recursive enumeration admitted a subthreshold tuple")
            result.add(tuple(sorted(speeds)))
            return
        if total + best_from(start, need) < threshold:
            return
        final_index = len(ranked) - need
        for index in range(start, final_index + 1):
            optimistic = total + ranked[index][0] + best_from(index + 1, need - 1)
            if optimistic < threshold:
                break
            visit(
                index + 1,
                need - 1,
                total + ranked[index][0],
                (*speeds, ranked[index][1]),
            )

    visit(0, arity, F(0), ())
    return result, visited


def subtract_sequence(core, good: list[tuple[F, F]], speeds: tuple[int, ...]):
    carrier = good
    for speed in speeds:
        _, _, carrier = core.subtract(carrier, speed)
    return carrier, len(carrier), sum((right - left for left, right in carrier), F(0))


def close_dangerous_triples(core, bank: dict[str, object]) -> dict[str, object]:
    edge = bank["edge"]
    expected = EXPECTED[edge]
    root_m = bank["m"]
    threshold = root_m - bank["tail_first_cap"]
    ranked = bank["rows"]
    candidates = ranked_triples(ranked, threshold)
    independent, recursive_nodes = recursive_upper_envelope(ranked, 3, threshold)
    require(candidates == independent, f"triple enumeration paths disagree {edge}")
    require(
        tuple(sorted(candidates)) == expected["dangerous_speeds"],
        f"dangerous triple count changed {edge}: {len(candidates)}",
    )

    minimum: tuple[F, tuple[int, int, int], F, int] | None = None
    records: list[str] = []
    for speeds in sorted(candidates):
        nested_good, residual_r, residual_m = subtract_sequence(core, bank["good"], speeds)
        direct_good, direct_r, direct_m = core.good_norm((*bank["body"], *speeds))
        reverse_good, reverse_r, reverse_m = subtract_sequence(
            core, bank["good"], tuple(reversed(speeds))
        )
        require(
            nested_good == direct_good == reverse_good,
            f"dangerous residual carrier mismatch {edge},speeds={speeds}",
        )
        require(
            (residual_r, residual_m) == (direct_r, direct_m) == (reverse_r, reverse_m),
            f"dangerous residual summary mismatch {edge},speeds={speeds}",
        )
        residual_tail_margin = (
            F(6, 7) * residual_m
            - core.S2 * residual_r / (7 * TAIL_FIRST)
        )
        require(residual_tail_margin > 0, f"dangerous tail margin failed {edge},speeds={speeds}")
        candidate = (residual_tail_margin, speeds, residual_m, residual_r)
        if minimum is None or candidate < minimum:
            minimum = candidate
        records.append(
            ":".join(
                (
                    f"{edge[0]}{edge[1]}",
                    ",".join(map(str, speeds)),
                    fraction_text(residual_m),
                    str(residual_r),
                    fraction_text(residual_tail_margin),
                )
            )
        )

    require(minimum == expected["triple_minimum"], f"dangerous minimum changed {edge}")
    triple_manifest = hashlib.sha256("\n".join(records).encode()).hexdigest()
    require(
        triple_manifest == expected["triple_manifest"],
        f"dangerous manifest changed {edge}",
    )
    return {
        "count": len(candidates),
        "recursive_nodes": recursive_nodes,
        "minimum": minimum,
        "manifest": triple_manifest,
        "speeds": tuple(sorted(candidates)),
    }


def close_head_quadruples(core, bank: dict[str, object]) -> dict[str, object]:
    edge = bank["edge"]
    expected = EXPECTED[edge]
    root_m = bank["m"]
    ranked = bank["rows"]
    candidates = ranked_quadruples(ranked, root_m)
    independent, recursive_nodes = recursive_upper_envelope(ranked, 4, root_m)
    require(candidates == independent, f"quadruple enumeration paths disagree {edge}")
    require(
        len(candidates) == expected["quad_count"],
        f"head quadruple count changed {edge}: {len(candidates)}",
    )

    cache1: dict[int, list[tuple[F, F]]] = {}
    cache2: dict[tuple[int, int], list[tuple[F, F]]] = {}
    cache3: dict[tuple[int, int, int], list[tuple[F, F]]] = {}
    minimum: tuple[F, tuple[int, int, int, int], int] | None = None
    records: list[str] = []
    for speeds in sorted(candidates):
        a, b, c, d = speeds
        if a not in cache1:
            _, _, cache1[a] = core.subtract(bank["good"], a)
        key2 = (a, b)
        if key2 not in cache2:
            _, _, cache2[key2] = core.subtract(cache1[a], b)
        key3 = (a, b, c)
        if key3 not in cache3:
            _, _, cache3[key3] = core.subtract(cache2[key2], c)
        nested_r, nested_m, nested_good = core.subtract(cache3[key3], d)
        sparse_m = core.subtract_sparse(cache3[key3], d)
        direct_good, direct_r, direct_m = core.good_norm((*bank["body"], *speeds))
        require(sparse_m == nested_m, f"terminal sparse/full mismatch {edge},speeds={speeds}")
        require(nested_good == direct_good, f"nested/direct carrier mismatch {edge},speeds={speeds}")
        require(
            (nested_r, nested_m) == (direct_r, direct_m),
            f"nested/direct summary mismatch {edge},speeds={speeds}",
        )
        require(nested_m > 0, f"nonpositive all-head survivor {edge},speeds={speeds}")
        candidate = (nested_m, speeds, nested_r)
        if minimum is None or candidate < minimum:
            minimum = candidate
        records.append(
            ":".join(
                (
                    f"{edge[0]}{edge[1]}",
                    ",".join(map(str, speeds)),
                    fraction_text(nested_m),
                    str(nested_r),
                )
            )
        )

    require(minimum == expected["quad_minimum"], f"all-head minimum changed {edge}")
    quad_manifest = hashlib.sha256("\n".join(records).encode()).hexdigest()
    require(quad_manifest == expected["quad_manifest"], f"all-head manifest changed {edge}")
    return {
        "count": len(candidates),
        "recursive_nodes": recursive_nodes,
        "minimum": minimum,
        "manifest": quad_manifest,
    }


def audit_edge(core, edge: tuple[int, int]) -> dict[str, object]:
    bank = exact_coverage_bank(core, edge)
    triples = close_dangerous_triples(core, bank)
    quadruples = close_head_quadruples(core, bank)
    return {"bank": bank, "triples": triples, "quadruples": quadruples}


def main() -> None:
    core = load_core()
    print("THM741_LABEL1_HEAD_TAIL_PARTITION")
    print("optimization_modes=normal_and_O_share_truth_gates")
    print(f"python={sys.version.split()[0]}")
    print(f"core_sha256={CORE_SHA256}")
    print(f"finite_universe={FIRST_EXTERNAL}..{W}; tail={TAIL_FIRST}..infinity")
    print("root_universe=12,13,14,15,16,17")
    print("tail_cap=m/7+(99/70)r/(7w)")
    print("coverage_paths=sparse_subtract,full_subtract,direct_union,tooth_incidence")
    print("candidate_paths=specialized_ranked_loops,generic_upper_envelope_recursion")

    combined_records: list[str] = []
    for edge in EXPECTED:
        result = audit_edge(core, edge)
        bank = result["bank"]
        triples = result["triples"]
        quadruples = result["quadruples"]
        ranked = bank["rows"]
        triple_min_margin, triple_min_speeds, triple_min_measure, triple_min_r = triples["minimum"]
        quad_min_measure, quad_min_speeds, quad_min_r = quadruples["minimum"]
        dangerous_text = ";".join(",".join(map(str, speeds)) for speeds in triples["speeds"])
        print(
            " ".join(
                (
                    f"root={edge[0]}{edge[1]}",
                    f"r={bank['r']}",
                    f"m={fraction_text(bank['m'])}",
                    f"top2={ranked[0][1]}:{fraction_text(ranked[0][0])},"
                    f"{ranked[1][1]}:{fraction_text(ranked[1][0])}",
                    f"tail_cap_2501={fraction_text(bank['tail_first_cap'])}",
                    f"two_tail_margin={fraction_text(bank['multi_tail_margin'])}",
                    f"top4_hostile={fraction_text(bank['top4_margin'])}",
                )
            )
        )
        print(
            " ".join(
                (
                    f"root={edge[0]}{edge[1]}",
                    f"dangerous_triples={triples['count']}",
                    f"recursive_nodes={triples['recursive_nodes']}",
                    f"speeds={dangerous_text}",
                    f"min_tail_margin={fraction_text(triple_min_margin)}",
                    f"min_tail_at={','.join(map(str, triple_min_speeds))}",
                    f"residual_m={fraction_text(triple_min_measure)}",
                    f"residual_r={triple_min_r}",
                )
            )
        )
        print(
            " ".join(
                (
                    f"root={edge[0]}{edge[1]}",
                    f"head_quadruples={quadruples['count']}",
                    f"recursive_nodes={quadruples['recursive_nodes']}",
                    f"min_survivor={fraction_text(quad_min_measure)}",
                    f"min_at={','.join(map(str, quad_min_speeds))}",
                    f"components={quad_min_r}",
                )
            )
        )
        print(
            " ".join(
                (
                    f"root={edge[0]}{edge[1]}",
                    f"coverage_manifest={bank['coverage_manifest']}",
                    f"triple_manifest={triples['manifest']}",
                    f"quad_manifest={quadruples['manifest']}",
                )
            )
        )
        combined_records.append(
            ":".join(
                (
                    f"{edge[0]}{edge[1]}",
                    bank["coverage_manifest"],
                    triples["manifest"],
                    quadruples["manifest"],
                )
            )
        )

    combined_manifest = hashlib.sha256("\n".join(combined_records).encode()).hexdigest()
    require(combined_manifest == EXPECTED_COMBINED_MANIFEST, "combined manifest changed")
    print(f"combined_manifest={combined_manifest}")
    print("partition_complete=6/6")
    print("combined_flood_census=21/21")
    print("global_THM741=OPEN")
    print("LRC14=OPEN")
    print("THM741_LABEL1_HEAD_TAIL_PARTITION_OK")


if __name__ == "__main__":
    main()
