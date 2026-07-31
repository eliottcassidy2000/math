#!/usr/bin/env python3
r"""Ranked-exception closure for the pure flood roots 23,25,26,27.

For a root edge `e`, let `G_e` have measure `m_e` and define the individual
root coverage

    c_e(w) = |G_e intersect D_w|.

Suppose an exact finite bank, sealed by the THM-735(ii) tail cap, has ranked
global coverages `q_1>=q_2>=...`.  If

    m_e - q_1-q_2-q_3-q_{K+1} > 0,                      (1)

then every four-speed set not wholly contained in the first `K` ranks is
closed by the ordinary union bound.  It remains only to check the finite head
`C(K,4)` exactly.

This verifier proves (1) for roots 23,25,26,27 with `K=36,22,10,23`.
The exact coverage banks end at the floor of the rational crossing where the
THM-735(ii), THM-731/732 cap falls below `q_{K+1}`.  Every coverage is compared
by sparse subtraction, full subtraction, direct ten-comb union, and an
independent two-pointer tooth-incidence sum.

The complete head banks are evaluated by cached exact nested subtraction:
full interval carriers are cached through the first three speeds and the
fourth is subtracted sparsely.  Five deterministic quantiles per row are
replayed with full terminal subtraction.  Each row's global minimum is also
rebuilt both as a direct thirteen-comb union and as four full nested
subtractions.  All `75,285` head quadruples are strictly positive.

The script hash-pins the prior ten-body top-four suite.  Its four new bodies
give a local subtotal of 14/21; composed in THM-741 with the independent
root-37 top-five exception certificate, the exact census is 15/21.  The six
remaining roots are precisely 12,13,14,15,16,17.  No Fano/chi_7 transport and
no claim of global THM-741 or LRC(14) is made.
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
PRIOR_TEN_OUTPUT = (
    ROOT / "05-knowledge/results/lrc14_j4_six_more_top_four_coverage_codex_20260728.out"
)
PRIOR_TEN_SHA256 = "8654ecafa11749504d33ecfd8ff030eae842b5365ae1462001a3a3c88bb7af55"
H = tuple(range(8, 15))
FIRST_EXTERNAL = 15

EXPECTED = {
    (2, 3): {
        "r": 30,
        "m": F(358427, 2522520),
        "K": 36,
        "finite_max": 1796,
        "head": (
            46, 23, 17, 19, 21, 110, 125, 16, 69, 22, 156, 42, 168, 50,
            100, 53, 92, 34, 154, 130, 189, 79, 212, 115, 15, 182, 161,
            25, 81, 52, 138, 75, 18, 137, 96, 243,
        ),
        "outside": (350, F(29, 1225)),
        "threshold": F(535134600, 297953),
        "gate": F(30353, 410960550),
        "previous_gate": -F(667243, 3328780455),
        "top4_margin": -F(2106724, 156165009),
        "exception_count": 79,
        "minimum": (F(1068173, 39939900), (16, 19, 22, 25), 20),
    },
    (2, 5): {
        "r": 26,
        "m": F(409261, 2522520),
        "K": 22,
        "finite_max": 1175,
        "head": (
            23, 19, 46, 17, 25, 38, 22, 50, 34, 75, 69, 16, 80, 92, 168,
            81, 154, 137, 182, 175, 111, 108,
        ),
        "outside": (117, F(317, 11466)),
        "threshold": F(92756664, 78919),
        "gate": F(37841, 551170620),
        "previous_gate": -F(2069, 91861770),
        "top4_margin": -F(5458325, 468495027),
        "exception_count": 31,
        "minimum": (F(11739671, 340723656), (17, 19, 23, 25), 14),
    },
    (2, 6): {
        "r": 30,
        "m": F(413747, 2522520),
        "K": 10,
        "finite_max": 840,
        "head": (19, 23, 46, 22, 17, 16, 110, 21, 69, 125),
        "outside": (25, F(3303, 107800)),
        "threshold": F(267567300, 318211),
        "gate": F(69547, 612411800),
        "previous_gate": -F(81687, 306205900),
        "top4_margin": -F(1226821, 122482360),
        "exception_count": 21,
        "minimum": (F(287374099, 9369900540), (17, 19, 22, 46), 20),
    },
    (2, 7): {
        "r": 22,
        "m": F(52147, 315315),
        "K": 23,
        "finite_max": 924,
        "head": (
            19, 23, 46, 17, 38, 34, 69, 25, 22, 63, 21, 16, 125, 110, 50,
            100, 79, 105, 137, 54, 26, 130, 53,
        ),
        "outside": (92, F(206, 7245)),
        "threshold": F(225648423, 244061),
        "gate": F(15721, 367447080),
        "previous_gate": -F(1269581, 58424085720),
        "top4_margin": -F(1002853, 76488984),
        "exception_count": 36,
        "minimum": (F(39526373, 1338557220), (17, 19, 26, 46), 18),
    },
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_core():
    require(sha256(CORE_PATH) == CORE_SHA256, "THM-741 core hash changed")
    require(sha256(PRIOR_TEN_OUTPUT) == PRIOR_TEN_SHA256, "prior ten-body suite changed")
    spec = importlib.util.spec_from_file_location("thm741_ranked_exception_dependency", CORE_PATH)
    require(spec is not None and spec.loader is not None, "cannot load THM-741 core")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def direct_tooth_coverage(good: list[tuple[F, F]], w: int) -> F:
    """Independent two-pointer sum of intersections with the `D_w` teeth."""
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


def fraction_text(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def exact_coverage_bank(core, edge: tuple[int, int], expected: dict[str, object]) -> dict[str, object]:
    body = tuple(sorted((*H, *edge)))
    require(len(body) == 9 and len(set(body)) == 9 and max(body) == 14, f"bad body {edge}")
    good, root_r, root_m = core.good_norm(body)
    require(root_r == expected["r"] and root_m == expected["m"], f"root changed {edge}")
    finite_max = int(expected["finite_max"])

    rows: list[tuple[F, int, int]] = []
    coverage_by_speed: dict[int, F] = {}
    manifest_rows: list[str] = []
    for w in range(FIRST_EXTERNAL, finite_max + 1):
        sparse_survivor = core.subtract_sparse(good, w)
        full_r, full_survivor, _ = core.subtract(good, w)
        _, _, union_survivor = core.good_norm((*body, w))
        coverage = root_m - sparse_survivor
        tooth_coverage = direct_tooth_coverage(good, w)
        require(sparse_survivor == full_survivor, f"sparse/full mismatch {edge},w={w}")
        require(sparse_survivor == union_survivor, f"subtract/union mismatch {edge},w={w}")
        require(coverage == tooth_coverage, f"subtract/tooth mismatch {edge},w={w}")
        require(F(0) <= coverage < root_m, f"bad coverage {edge},w={w}")
        rows.append((coverage, w, full_r))
        coverage_by_speed[w] = coverage
        manifest_rows.append(f"{edge[0]}{edge[1]}:{w}:{fraction_text(coverage)}:{full_r}")
    require(len(rows) == finite_max - FIRST_EXTERNAL + 1, f"bad coverage census {edge}")

    ranked = sorted(rows, key=lambda row: (-row[0], row[1]))
    K = int(expected["K"])
    head = tuple(w for _, w, _ in ranked[:K])
    require(head == expected["head"], f"ranked head changed {edge}")
    outside = (ranked[K][1], ranked[K][0])
    require(outside == expected["outside"], f"outside rank changed {edge}")
    require(ranked[K - 1][0] > ranked[K][0], f"head boundary is tied {edge}")

    fourth = ranked[3][0]
    top3_sum = sum((ranked[index][0] for index in range(3)), F(0))
    top4_margin = root_m - top3_sum - fourth
    previous_gate = root_m - top3_sum - ranked[K - 1][0]
    gate = root_m - top3_sum - ranked[K][0]
    require(top4_margin == expected["top4_margin"] < 0, f"bad top-four hostility {edge}")
    require(previous_gate == expected["previous_gate"] <= 0, f"bad minimal-rank hostility {edge}")
    require(gate == expected["gate"] > 0, f"bad outside gate {edge}")

    outside_coverage = ranked[K][0]

    def tail_cap(w: int) -> F:
        return root_m / 7 + core.S2 * root_r / (7 * w)

    require(outside_coverage > root_m / 7, f"outside coverage below limiting cap {edge}")
    threshold = core.S2 * root_r / (7 * (outside_coverage - root_m / 7))
    require(threshold == expected["threshold"], f"tail threshold changed {edge}")
    require(
        threshold.numerator // threshold.denominator == finite_max,
        f"finite maximum is not exact threshold floor {edge}",
    )
    require(
        tail_cap(finite_max) > outside_coverage > tail_cap(finite_max + 1),
        f"bad tail crossing {edge}",
    )

    coverage_manifest = hashlib.sha256("\n".join(manifest_rows).encode()).hexdigest()
    return {
        "edge": edge,
        "body": body,
        "good": good,
        "r": root_r,
        "m": root_m,
        "K": K,
        "finite_max": finite_max,
        "head": head,
        "outside": outside,
        "top4_margin": top4_margin,
        "previous_gate": previous_gate,
        "gate": gate,
        "threshold": threshold,
        "lower_crossing_gap": tail_cap(finite_max) - outside_coverage,
        "upper_crossing_gap": outside_coverage - tail_cap(finite_max + 1),
        "coverage_by_speed": coverage_by_speed,
        "coverage_manifest": coverage_manifest,
    }


def exact_head_bank(core, coverage: dict[str, object], expected: dict[str, object]) -> dict[str, object]:
    edge = coverage["edge"]
    body = coverage["body"]
    good = coverage["good"]
    root_m = coverage["m"]
    head_speeds = tuple(sorted(coverage["head"]))
    K = coverage["K"]
    terminal_count = comb(K, 4)
    sample_indices = {
        0,
        (terminal_count - 1) // 4,
        (terminal_count - 1) // 2,
        3 * (terminal_count - 1) // 4,
        terminal_count - 1,
    }
    require(len(sample_indices) == 5, f"collapsed sample indices {edge}")

    cache1: dict[int, list[tuple[F, F]]] = {}
    cache2: dict[tuple[int, int], list[tuple[F, F]]] = {}
    cache3: dict[tuple[int, int, int], list[tuple[F, F]]] = {}
    minimum: tuple[F, tuple[int, int, int, int]] | None = None
    exception_minimum: tuple[F, tuple[int, int, int, int]] | None = None
    positive = zero = exception_count = 0
    terminal_records: list[str] = []
    sample_records: list[str] = []

    for index, speeds in enumerate(combinations(head_speeds, 4)):
        a, b, c, d = speeds
        if a not in cache1:
            _, _, cache1[a] = core.subtract(good, a)
        key2 = (a, b)
        if key2 not in cache2:
            _, _, cache2[key2] = core.subtract(cache1[a], b)
        key3 = (a, b, c)
        if key3 not in cache3:
            _, _, cache3[key3] = core.subtract(cache2[key2], c)
        value = core.subtract_sparse(cache3[key3], d)
        if value > 0:
            positive += 1
        else:
            zero += 1
        candidate = (value, speeds)
        if minimum is None or candidate < minimum:
            minimum = candidate
        union_margin = root_m - sum(
            (coverage["coverage_by_speed"][w] for w in speeds),
            F(0),
        )
        if union_margin <= 0:
            exception_count += 1
            if exception_minimum is None or candidate < exception_minimum:
                exception_minimum = candidate
        if index in sample_indices:
            full_r, full_value, _ = core.subtract(cache3[key3], d)
            require(full_value == value, f"terminal sparse/full mismatch {edge},speeds={speeds}")
            sample_records.append(
                f"{edge}:{index}/{terminal_count}:{speeds}:{value}:{full_r}"
            )
        terminal_records.append(f"{edge}:{index}:{speeds}:{value}:{union_margin}")

    require(positive == terminal_count and zero == 0, f"nonpositive head terminal {edge}")
    require(exception_count == expected["exception_count"], f"exception count changed {edge}")
    require(minimum is not None and exception_minimum is not None, f"missing minimum {edge}")
    expected_value, expected_speeds, expected_r = expected["minimum"]
    require(minimum == (expected_value, expected_speeds), f"global head minimum changed {edge}")

    require(len(cache1) == K - 3, f"bad first-carrier cache {edge}")
    require(len(cache2) == comb(K - 2, 2), f"bad second-carrier cache {edge}")
    require(len(cache3) == comb(K - 1, 3), f"bad third-carrier cache {edge}")
    require(len(sample_records) == 5, f"bad terminal sample count {edge}")

    _, direct_r, direct_m = core.good_norm((*body, *expected_speeds))
    nested = good
    nested_r = 0
    for w in expected_speeds:
        nested_r, _, nested = core.subtract(nested, w)
    nested_m = sum((right - left for left, right in nested), F(0))
    require(direct_m == nested_m == expected_value, f"minimum direct/nested mismatch {edge}")
    require(direct_r == nested_r == expected_r, f"minimum component count changed {edge}")

    terminal_manifest = hashlib.sha256("\n".join(terminal_records).encode()).hexdigest()
    sample_manifest = hashlib.sha256("\n".join(sample_records).encode()).hexdigest()
    minimum_union_margin = root_m - sum(
        (coverage["coverage_by_speed"][w] for w in expected_speeds),
        F(0),
    )
    return {
        **coverage,
        "terminal_count": terminal_count,
        "cache_counts": (len(cache1), len(cache2), len(cache3)),
        "exception_count": exception_count,
        "minimum": minimum,
        "minimum_r": expected_r,
        "minimum_union_margin": minimum_union_margin,
        "exception_minimum": exception_minimum,
        "terminal_manifest": terminal_manifest,
        "sample_manifest": sample_manifest,
    }


def print_row(row: dict[str, object]) -> None:
    edge = "".join(map(str, row["edge"]))
    outside_w, outside_coverage = row["outside"]
    print(
        f"edge={edge} root_r={row['r']} root_m={row['m']}; "
        f"K={row['K']} finite_w=15..{row['finite_max']}"
    )
    print(
        f"  q_(K+1)=c({outside_w})={outside_coverage}; "
        f"gate={row['gate']}; previous_gate={row['previous_gate']}; "
        f"top4_margin={row['top4_margin']}"
    )
    print(
        f"  threshold={row['threshold']}; "
        f"u(W)-q_(K+1)={row['lower_crossing_gap']}; "
        f"q_(K+1)-u(W+1)={row['upper_crossing_gap']}"
    )
    minimum_value, minimum_speeds = row["minimum"]
    exception_value, exception_speeds = row["exception_minimum"]
    print(
        f"  head C(K,4)={row['terminal_count']}; caches={row['cache_counts']}; "
        f"union-bound exceptions={row['exception_count']}; all positive"
    )
    print(
        f"  global head minimum={minimum_value} at {minimum_speeds} r={row['minimum_r']}; "
        f"its union margin={row['minimum_union_margin']}; "
        f"exception minimum={exception_value} at {exception_speeds}"
    )
    print(
        f"  coverage_manifest={row['coverage_manifest']}; "
        f"terminal_manifest={row['terminal_manifest']}; "
        f"sample_manifest={row['sample_manifest']}"
    )


def main() -> None:
    core = load_core()
    require(core.S2 == F(99, 70) and core.S2 * core.S2 > 2, "bad sqrt(2) majorant")
    rows = []
    for edge, expected in EXPECTED.items():
        coverage = exact_coverage_bank(core, edge, expected)
        rows.append(exact_head_bank(core, coverage, expected))

    require(sum(row["terminal_count"] for row in rows) == 75285, "bad terminal total")
    require(sum(row["exception_count"] for row in rows) == 167, "bad exception total")
    suite_payload = "\n".join(
        f"{''.join(map(str, row['edge']))}:{row['coverage_manifest']}:"
        f"{row['terminal_manifest']}:{row['sample_manifest']}"
        for row in rows
    ).encode()
    suite_manifest = hashlib.sha256(suite_payload).hexdigest()

    print("THM-741 FLOOD TAIL: FOUR EXACT RANKED-EXCEPTION COVERAGE CERTIFICATES")
    print("=" * 92)
    print(f"dependency_sha256={CORE_SHA256}")
    print(f"prior-ten-output_sha256={PRIOR_TEN_SHA256}")
    print(
        "coverage path: sparse/full/direct-union/direct-tooth on each finite bank; "
        "tail sealed by THM-735(ii) via THM-731/732"
    )
    print(
        "head path: cached full nested carriers through three speeds, sparse fourth; "
        "five full-terminal samples and direct/full minimum replay per root"
    )
    for row in rows:
        print_row(row)
    print(f"suite-manifest SHA256={suite_manifest}")
    print("new uniformly closed roots=((2,3),(2,5),(2,6),(2,7)); count=4")
    print(
        "composition: hash-pinned prior ten roots + these four =14/21; "
        "with independent canonical root37 exception certificate =15/21"
    )
    print("remaining literal flood roots=((1,2),(1,3),(1,4),(1,5),(1,6),(1,7))")
    print(
        "scope: four literal roots, no Fano/chi_7 transport; "
        "global THM-741 and LRC(14) remain open"
    )
    print(
        "VERDICT: roots 23,25,26,27 are uniformly lonely over all "
        "15<=a<b<c<d; exact combined flood census is 15/21"
    )
    print(f"source_sha256={sha256(Path(__file__))}")
    print("ALL FOUR RANKED-EXCEPTION COVERAGE CHECKS PASSED")


if __name__ == "__main__":
    main()
