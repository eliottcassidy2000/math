#!/usr/bin/env python3
r"""Exact top-five-plus-one-exception certificate for the `(3,7)` flood tail.

Fix

    E = {3,7,8,9,10,11,12,13,14}

and let `G=G(E)`.  For every external speed `w`, put

    c(w)=|G intersect D_w|.

The four largest individual coverages have sum slightly above `|G|`, so the
plain top-four envelope does not close this root.  The fifth-largest coverage
repairs every quadruple except the literal top-four set: if a distinct
four-speed set is not that set, its coverage sum is at most

    c_1+c_2+c_3+c_5 < |G|.

The sole exceptional set `{17,19,21,23}` is then reconstructed exactly and
has positive survivor.

This verifier computes every `c(w)` for `15<=w<=293`.  THM-735(ii), through
the THM-731/732 covariance-discrepancy chain, gives

    c(w) <= |G|/7+sqrt(2)r(G)/(7w)
         < |G|/7+(99/70)r(G)/(7w).

The rational tail cap at `w=294` lies below the exact fifth coverage, so the
finite top five are global.  Sparse subtraction, full subtraction, direct
ten-comb union, and an independent tooth-incidence sum agree throughout.
The exceptional four-comb survivor is checked both by nested subtraction and
direct thirteen-comb union.
"""

from __future__ import annotations

import hashlib
import importlib.util
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CORE_PATH = ROOT / "04-computation/lrc14_thm741_2002_body_j4_tree_kps_S128c5.py"
CORE_SHA256 = "5aa81d9d78273c8f9e3e7a6574091a3bc3f64ab6086c7024c15f9420c99dac96"
E = (3, 7, 8, 9, 10, 11, 12, 13, 14)
FIRST_EXTERNAL = 15
FINITE_MAX = 293
TAIL_FIRST = FINITE_MAX + 1
EXPECTED_TOP_FIVE = (
    (19, F(134663, 2662660)),
    (17, F(2578, 51051)),
    (21, F(38609, 840840)),
    (23, F(25331, 552552)),
    (46, F(79435, 1933932)),
)
EXCEPTIONAL_SPEEDS = (17, 19, 21, 23)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_core():
    require(sha256(CORE_PATH) == CORE_SHA256, "THM-741 core hash changed")
    spec = importlib.util.spec_from_file_location("thm741_37_top_five_dependency", CORE_PATH)
    require(spec is not None and spec.loader is not None, "cannot load THM-741 core")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def direct_tooth_coverage(good: list[tuple[F, F]], w: int) -> F:
    """Measure `|G intersect D_w|` from literal disjoint tooth incidences."""
    require(w >= 1, "speed must be positive")
    radius = F(1, 14 * w)
    total = F(0)
    for left, right in good:
        require(F(0) <= left < right <= F(1), "bad good-set component")
        for k in range(w + 1):
            center = F(k, w)
            lo = max(left, center - radius)
            hi = min(right, center + radius)
            if lo < hi:
                total += hi - lo
    return total


def fraction_text(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def main() -> None:
    core = load_core()
    require(core.S2 == F(99, 70), "unexpected rational sqrt(2) majorant")
    require(core.S2 * core.S2 > 2, "S2 is not above sqrt(2)")

    good, root_r, root_m = core.good_norm(E)
    require(root_r == 20, "root component count changed")
    require(root_m == F(53619, 280280), "root measure changed")
    require(len(E) == 9 and len(set(E)) == 9 and max(E) == 14, "bad root body")

    rows: list[tuple[F, int, int]] = []
    manifest_rows: list[str] = []
    for w in range(FIRST_EXTERNAL, FINITE_MAX + 1):
        sparse_survivor = core.subtract_sparse(good, w)
        full_r, full_survivor, _ = core.subtract(good, w)
        _, _, direct_union_survivor = core.good_norm((*E, w))
        coverage = root_m - sparse_survivor
        direct_coverage = direct_tooth_coverage(good, w)
        require(sparse_survivor == full_survivor, f"sparse/full mismatch w={w}")
        require(sparse_survivor == direct_union_survivor, f"subtraction/union mismatch w={w}")
        require(coverage == direct_coverage, f"subtraction/tooth mismatch w={w}")
        require(F(0) <= coverage < root_m, f"bad coverage w={w}")
        rows.append((coverage, w, full_r))
        manifest_rows.append(f"{w}:{fraction_text(coverage)}:{full_r}")

    require(len(rows) == FINITE_MAX - FIRST_EXTERNAL + 1 == 279, "bad finite census")
    ranked = sorted(rows, key=lambda row: (-row[0], row[1]))
    observed_top = tuple((w, coverage) for coverage, w, _ in ranked[:5])
    require(observed_top == EXPECTED_TOP_FIVE, "top-five finite coverage list changed")
    require(ranked[4][0] > ranked[5][0], "fifth coverage is not isolated")

    fifth = EXPECTED_TOP_FIVE[4][1]

    def tail_cap(w: int) -> F:
        return root_m / 7 + core.S2 * root_r / (7 * w)

    require(tail_cap(FINITE_MAX) > fifth, "finite cutoff is not sharp below")
    require(tail_cap(TAIL_FIRST) < fifth, "tail cap does not clear fifth coverage")
    tail_threshold = core.S2 * root_r / (7 * (fifth - root_m / 7))
    require(
        F(FINITE_MAX) < tail_threshold < F(TAIL_FIRST),
        "tail threshold is not between 293 and 294",
    )

    top_three_plus_fifth = sum(
        (coverage for _, coverage in EXPECTED_TOP_FIVE[:3]),
        fifth,
    )
    generic_margin = root_m - top_three_plus_fifth
    require(
        generic_margin == F(843411, 260275015),
        "generic nonexceptional margin changed",
    )
    require(generic_margin > 0, "nonexceptional top-five gate is not positive")

    top_four_sum = sum(
        (coverage for _, coverage in EXPECTED_TOP_FIVE[:4]),
        F(0),
    )
    failed_plain_margin = root_m - top_four_sum
    require(
        failed_plain_margin == -F(3183347, 2082200120),
        "plain top-four hostile margin changed",
    )

    nested = good
    for w in EXCEPTIONAL_SPEEDS:
        _, _, nested = core.subtract(nested, w)
    nested_m = sum((right - left for left, right in nested), F(0))
    exceptional_family = tuple(sorted((*E, *EXCEPTIONAL_SPEEDS)))
    direct_good, direct_r, direct_m = core.good_norm(exceptional_family)
    require(len(direct_good) == direct_r, "direct exceptional component count mismatch")
    require(nested == direct_good, "nested/direct exceptional interval mismatch")
    require(direct_r == 8, "exceptional component count changed")
    require(direct_m == nested_m == F(137897, 2299080), "exceptional survivor changed")
    require(direct_m > 0, "exceptional top-four family is not lonely")

    manifest_sha256 = hashlib.sha256("\n".join(manifest_rows).encode()).hexdigest()
    source_sha256 = sha256(Path(__file__))
    print("THM-741 PURE (3,7) FLOOD TAIL: TOP FIVE PLUS ONE EXACT EXCEPTION")
    print("=" * 88)
    print(f"dependency_sha256={CORE_SHA256}")
    print(f"body E={E}; external universe: all integers 15<=a<b<c<d")
    print(f"root r={root_r} m={root_m}")
    print(
        "finite exact bank: w=15..293 (279 speeds); "
        "sparse/full/direct-union/direct-tooth paths agree"
    )
    for rank, (w, coverage) in enumerate(EXPECTED_TOP_FIVE, start=1):
        print(f"coverage rank {rank}: w={w} c(w)={coverage}")
    print(f"finite-coverage manifest SHA256={manifest_sha256}")
    print(
        "tail law: c(w)<u(w)=m/7+(99/70)r/(7w) for w>=294; "
        f"threshold={tail_threshold}"
    )
    print(
        f"cutoff hostility: u(293)-c(46)={tail_cap(293)-fifth}; "
        f"c(46)-u(294)={fifth-tail_cap(294)}"
    )
    print(f"plain top-four union-bound margin={failed_plain_margin}<0")
    print(
        "every non-top-four speed set has individual-cover sum at most "
        f"c1+c2+c3+c5; survivor>={generic_margin}>0"
    )
    print(
        f"sole exceptional speed set={EXCEPTIONAL_SPEEDS}; "
        f"exact survivor r={direct_r} m={direct_m}"
    )
    print(
        "proof partition: exact c(w) for 15<=w<=293; THM-735(ii) "
        "monotone cap for every w>=294; top-five ranking; one exact exception"
    )
    print(
        "scope: one literal root edge (3,7), no Fano/chi_7 transport; "
        "global THM-741 and LRC(14) remain open"
    )
    print(
        "VERDICT: E_(3,7) union {a,b,c,d} is lonely for every "
        "15<=a<b<c<d; this closes a fifth whole flood body"
    )
    print(f"source_sha256={source_sha256}")
    print("ALL EXACT TOP-FIVE/EXCEPTION CHECKS PASSED")


if __name__ == "__main__":
    main()
