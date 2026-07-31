#!/usr/bin/env python3
r"""Six-root exact top-four coverage certificate for the THM-741 flood tail.

For a small-label edge `e`, put

    E_e = {8,9,10,11,12,13,14} union e,
    G_e = G(E_e),
    c_e(w) = |G_e intersect D_w|.

If `15<=a<b<c<d`, the final good measure is at least

    |G_e|-c_e(a)-c_e(b)-c_e(c)-c_e(d).                    (1)

For each of the six literal roots 24,35,36,45,46,47, this verifier evaluates
`c_e(w)` exactly for all 468 integers `15<=w<=482`.  THM-735(ii), through the
THM-731/732 covariance-discrepancy chain, gives

    c_e(w) <= |G_e|/7 + sqrt(2)r_e/(7w)
           < |G_e|/7 + (99/70)r_e/(7w) = u_e(w).          (2)

The fourth-largest exact finite coverage is above `u_e(483)` in every row.
Thus no omitted speed can enter the global top four.  In all six rows the sum
of those four values is below `|G_e|`, so (1) is uniformly positive.

Every finite coverage is compared through four exact paths: pinned sparse
subtraction, pinned full subtraction, a direct ten-comb union, and a local
two-pointer incidence sum over the disjoint teeth of `D_w`.  Roots 13 and 16
are run through the same complete audit as negative controls: their exact
top-four margins are negative, so this certificate correctly leaves them open
even though their literal extremal families are lonely.  The closer exceptional
root 37 is deliberately excluded for a separate fifth-rank repair.

Together with the hash-pinned prior whole-body certificates for roots
34,56,57,67, this raises the exact flood count to 10/21.  It is not a
Fano/chi_7 quotient and does not prove global THM-741 or LRC(14).
"""

from __future__ import annotations

import hashlib
import importlib.util
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CORE_PATH = ROOT / "04-computation/lrc14_thm741_2002_body_j4_tree_kps_S128c5.py"
CORE_SHA256 = "5aa81d9d78273c8f9e3e7a6574091a3bc3f64ab6086c7024c15f9420c99dac96"
H = tuple(range(8, 15))
FIRST_EXTERNAL = 15
FINITE_MAX = 482
TAIL_FIRST = FINITE_MAX + 1

PRIOR_OUTPUTS = {
    (3, 4): (
        ROOT / "05-knowledge/results/lrc14_j4_34_exact_top_four_coverage_codex_20260728.out",
        "5ac115e222c9c51205b7f32a8a06569684f0aee083eac94b32ef2dcdae14880d",
    ),
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

POSITIVE_EXPECTED = {
    (2, 4): {
        "r": 30,
        "m": F(31789, 194040),
        "top": (
            (23, F(240307, 5801796)),
            (22, F(879, 21560)),
            (17, F(10777, 274890)),
            (16, F(37445, 1009008)),
        ),
        "margin": F(3499547, 657536880),
        "threshold": F(214053840, 484061),
    },
    (3, 5): {
        "r": 24,
        "m": F(652, 3465),
        "top": (
            (17, F(345047, 7147140)),
            (19, F(38497, 798798)),
            (23, F(41075, 892584)),
            (46, F(6085, 148764)),
        ),
        "margin": F(1988211, 416440024),
        "threshold": F(10820304, 31291),
    },
    (3, 6): {
        "r": 28,
        "m": F(504167, 2522520),
        "top": (
            (19, F(30433, 614460)),
            (17, F(10649, 216580)),
            (23, F(48631, 1054872)),
            (21, F(38609, 840840)),
        ),
        "margin": F(57162379, 6246600360),
        "threshold": F(49945896, 153311),
    },
    (4, 5): {
        "r": 24,
        "m": F(514471, 2522520),
        "top": (
            (17, F(346037, 7147140)),
            (23, F(21577, 446292)),
            (19, F(3431, 72618)),
            (32, F(86851, 2018016)),
        ),
        "margin": F(6497515, 384406176),
        "threshold": F(342486144, 981901),
    },
    (4, 6): {
        "r": 28,
        "m": F(518957, 2522520),
        "top": (
            (17, F(10679, 216580)),
            (19, F(35279, 726180)),
            (23, F(12499, 263718)),
            (32, F(51859, 1261260)),
        ),
        "margin": F(2532941, 131047560),
        "threshold": F(33297264, 69023),
    },
    (4, 7): {
        "r": 20,
        "m": F(261193, 1261260),
        "top": (
            (17, F(36191, 714714)),
            (19, F(12013, 242060)),
            (23, F(13309, 276276)),
            (32, F(3467, 80080)),
        ),
        "margin": F(575561731, 37479602160),
        "threshold": F(28540512, 96835),
    },
}

NEGATIVE_EXPECTED = {
    (1, 3): {
        "r": 30,
        "m": F(3319, 25740),
        "top": (
            (23, F(453587, 11603592)),
            (17, F(67493, 1786785)),
            (19, F(148727, 3993990)),
            (46, F(67423, 1933932)),
        ),
        "margin": -F(25012943, 1249320072),
        "threshold": F(87914970, 238493),
    },
    (1, 6): {
        "r": 30,
        "m": F(365567, 2522520),
        "top": (
            (19, F(377149, 7987980)),
            (23, F(47539, 1054872)),
            (17, F(8439, 216580)),
            (46, F(1045, 29302)),
        ),
        "margin": -F(6867281, 312330018),
        "threshold": F(2461619160, 6075659),
    },
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_core():
    require(sha256(CORE_PATH) == CORE_SHA256, "THM-741 core hash changed")
    spec = importlib.util.spec_from_file_location("thm741_six_coverage_dependency", CORE_PATH)
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


def audit_edge(core, edge: tuple[int, int], expected: dict[str, object]) -> dict[str, object]:
    body = tuple(sorted((*H, *edge)))
    require(len(body) == 9 and len(set(body)) == 9 and max(body) == 14, f"bad body {edge}")
    good, root_r, root_m = core.good_norm(body)
    require(root_r == expected["r"] and root_m == expected["m"], f"root changed {edge}")

    rows: list[tuple[F, int, int]] = []
    manifest_rows: list[str] = []
    for w in range(FIRST_EXTERNAL, FINITE_MAX + 1):
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
        manifest_rows.append(f"{edge[0]}{edge[1]}:{w}:{fraction_text(coverage)}:{full_r}")
    require(len(rows) == 468, f"bad finite census {edge}")

    ranked = sorted(rows, key=lambda row: (-row[0], row[1]))
    observed_top = tuple((w, coverage) for coverage, w, _ in ranked[:4])
    require(observed_top == expected["top"], f"top four changed {edge}")
    require(ranked[3][0] > ranked[4][0], f"fourth coverage is not isolated {edge}")
    fourth = ranked[3][0]

    def tail_cap(w: int) -> F:
        return root_m / 7 + core.S2 * root_r / (7 * w)

    require(fourth > root_m / 7, f"fourth coverage does not clear limiting cap {edge}")
    threshold = core.S2 * root_r / (7 * (fourth - root_m / 7))
    require(threshold == expected["threshold"], f"threshold changed {edge}")
    threshold_floor = threshold.numerator // threshold.denominator
    require(
        tail_cap(threshold_floor) > fourth > tail_cap(threshold_floor + 1),
        f"bad exact cap crossing {edge}",
    )
    require(threshold < TAIL_FIRST and tail_cap(TAIL_FIRST) < fourth, f"tail not sealed {edge}")

    top_sum = sum((coverage for _, coverage in observed_top), F(0))
    margin = root_m - top_sum
    require(margin == expected["margin"], f"top-four margin changed {edge}")

    extremal_speeds = tuple(sorted(w for w, _ in observed_top))
    extremal_family = tuple(sorted((*body, *extremal_speeds)))
    require(len(extremal_family) == 13 and len(set(extremal_family)) == 13, f"bad extremal {edge}")
    _, extremal_r, extremal_m = core.good_norm(extremal_family)
    nested = good
    for w in extremal_speeds:
        _, _, nested = core.subtract(nested, w)
    nested_m = sum((right - left for left, right in nested), F(0))
    require(extremal_r == len(nested) and extremal_m == nested_m, f"extremal mismatch {edge}")
    require(extremal_m > 0, f"extremal family covers {edge}")
    if margin > 0:
        require(extremal_m >= margin, f"union bound exceeds exact survivor {edge}")

    manifest = hashlib.sha256("\n".join(manifest_rows).encode()).hexdigest()
    return {
        "edge": edge,
        "r": root_r,
        "m": root_m,
        "top": observed_top,
        "top_sum": top_sum,
        "margin": margin,
        "threshold": threshold,
        "threshold_floor": threshold_floor,
        "lower_crossing_gap": tail_cap(threshold_floor) - fourth,
        "upper_crossing_gap": fourth - tail_cap(threshold_floor + 1),
        "extremal_family": extremal_family,
        "extremal_r": extremal_r,
        "extremal_m": extremal_m,
        "manifest": manifest,
    }


def print_row(row: dict[str, object], verdict: str) -> None:
    edge = "".join(map(str, row["edge"]))
    top = ",".join(f"{w}:{coverage}" for w, coverage in row["top"])
    print(f"edge={edge} root_r={row['r']} root_m={row['m']}")
    print(f"  top4={top}; sum={row['top_sum']}; margin={row['margin']}")
    print(
        f"  threshold={row['threshold']} floor={row['threshold_floor']}; "
        f"u(floor)-q4={row['lower_crossing_gap']}; "
        f"q4-u(floor+1)={row['upper_crossing_gap']}"
    )
    print(
        f"  extremal={row['extremal_family']} exact_r={row['extremal_r']} "
        f"exact_m={row['extremal_m']}; manifest={row['manifest']}; verdict={verdict}"
    )


def main() -> None:
    core = load_core()
    require(core.S2 == F(99, 70) and core.S2 * core.S2 > 2, "bad sqrt(2) majorant")
    for edge, (path, digest) in PRIOR_OUTPUTS.items():
        require(sha256(path) == digest, f"prior whole-body output changed {edge}")

    positive_rows = [
        audit_edge(core, edge, expected)
        for edge, expected in POSITIVE_EXPECTED.items()
    ]
    negative_rows = [
        audit_edge(core, edge, expected)
        for edge, expected in NEGATIVE_EXPECTED.items()
    ]
    require(all(row["margin"] > 0 for row in positive_rows), "positive suite has a nonpositive row")
    require(all(row["margin"] < 0 for row in negative_rows), "negative control unexpectedly closes")
    require(
        max((row["threshold_floor"], row["edge"]) for row in positive_rows)
        == (FINITE_MAX, (4, 6)),
        "universal finite cutoff is not pinned by root 46",
    )

    new_edges = tuple(row["edge"] for row in positive_rows)
    prior_edges = tuple(PRIOR_OUTPUTS)
    combined = tuple(sorted((*prior_edges, *new_edges)))
    require(len(combined) == 10 and len(set(combined)) == 10, "bad whole-body count")
    suite_payload = "\n".join(
        f"{''.join(map(str, row['edge']))}:{row['manifest']}"
        for row in (*positive_rows, *negative_rows)
    ).encode()
    suite_manifest = hashlib.sha256(suite_payload).hexdigest()

    print("THM-741 FLOOD TAIL: SIX MORE EXACT TOP-FOUR COVERAGE ENVELOPES")
    print("=" * 88)
    print(f"dependency_sha256={CORE_SHA256}")
    print(
        "finite universe per audited root: all w=15..482 (468 speeds); "
        "sparse/full/direct-union/direct-tooth paths agree"
    )
    print("POSITIVE ROOTS")
    for row in positive_rows:
        print_row(row, "CLOSED UNIFORMLY")
    print("NEGATIVE CONTROLS")
    for row in negative_rows:
        print_row(row, "INCONCLUSIVE BY TOP-FOUR ENVELOPE")
    print(f"suite-manifest SHA256={suite_manifest}")
    print(f"hash-pinned prior whole roots={prior_edges}")
    print(f"new whole roots={new_edges}")
    print(f"combined exact whole roots={combined}; count=10/21")
    print(
        "proof partition: exact c_e(w) for 15<=w<=482; THM-735(ii) "
        "(via THM-731/732) cap for every w>=483; four-speed union bound"
    )
    print(
        "scope: six new literal root edges, no Fano/chi_7 transport; "
        "all eleven still-unclosed roots (including the two controls and "
        "root 37, which is reserved for a separate fifth-rank audit), "
        "global THM-741, and LRC(14) are not "
        "closed by this certificate"
    )
    print(
        "VERDICT: roots 24,35,36,45,46,47 are uniformly lonely over all "
        "15<=a<b<c<d; together with 34,56,57,67 this proves 10/21 flood bodies"
    )
    print(f"source_sha256={sha256(Path(__file__))}")
    print("ALL SIX-ROOT TOP-FOUR COVERAGE CHECKS PASSED")


if __name__ == "__main__":
    main()
