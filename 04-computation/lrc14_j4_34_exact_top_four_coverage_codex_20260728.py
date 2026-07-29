#!/usr/bin/env python3
r"""Exact top-four coverage certificate for the pure `(3,4)` j=4 flood tail.

Fix

    E = {3,4,8,9,10,11,12,13,14}

and let `G=G(E)` be its exact good set.  For an external speed `w`, write

    c(w) = |G intersect D_w|.

The final good set after four distinct external speeds satisfies the elementary
union lower bound

    |G minus (D_a union D_b union D_c union D_d)|
      >= |G| - c(a)-c(b)-c(c)-c(d).                       (1)

This verifier computes `c(w)` exactly for every integer `15<=w<=472`.
THM-735(ii), through the THM-731/732 covariance-discrepancy chain, and the
rational majorant `S2=99/70>sqrt(2)` give, for every larger integer,

    c(w) <= |G|/7 + sqrt(2) r(G)/(7w)
          < |G|/7 + S2 r(G)/(7w) = u(w).                  (2)

The exact fourth-largest coverage in the finite bank is `c(32)`.  The script
checks `u(473)<c(32)`; monotonicity of `u` then proves that no omitted speed
can enter the global top four.  The sum of the four exact largest coverages is
strictly below `|G|`, so (1) closes every ordered integer tail
`15<=a<b<c<d`.

Three exact evaluation paths are compared on all 458 finite speeds:

* the pinned THM-741 sparse subtraction kernel;
* its full interval-subtraction kernel and direct ten-comb union; and
* a local endpoint-incidence sum over the disjoint `D_w` teeth.

The local path is intentionally simple and does not call either subtraction
routine.  A direct thirteen-comb reconstruction also checks the family formed
by the four extremal individual coverages.  This is a proof for one literal
flood body, not a Fano/chi_7 quotient and not global THM-741.
"""

from __future__ import annotations

import hashlib
import importlib.util
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CORE_PATH = ROOT / "04-computation/lrc14_thm741_2002_body_j4_tree_kps_S128c5.py"
CORE_SHA256 = "5aa81d9d78273c8f9e3e7a6574091a3bc3f64ab6086c7024c15f9420c99dac96"
E = (3, 4, 8, 9, 10, 11, 12, 13, 14)
FIRST_EXTERNAL = 15
FINITE_MAX = 472
TAIL_FIRST = FINITE_MAX + 1
EXPECTED_TOP = (
    (17, F(85973, 1786785)),
    (23, F(240307, 5801796)),
    (19, F(14017, 363090)),
    (32, F(11521, 315315)),
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_core():
    require(sha256(CORE_PATH) == CORE_SHA256, "THM-741 core hash changed")
    spec = importlib.util.spec_from_file_location("thm741_34_top_four_dependency", CORE_PATH)
    require(spec is not None and spec.loader is not None, "cannot load THM-741 core")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def direct_tooth_coverage(good: list[tuple[F, F]], w: int) -> F:
    """Sum `G intersect D_w` directly from the disjoint comb teeth.

    On `[0,1]`, the circular tooth at zero is represented by the two half
    teeth centered at `0` and `1`.  The `w+1` closed endpoint representatives
    have pairwise disjoint interiors, so summing their intersections with the
    disjoint components of `G` counts measure exactly once.
    """
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
    require(root_r == 28, "root component count changed")
    require(root_m == F(433607, 2522520), "root measure changed")
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

    require(len(rows) == FINITE_MAX - FIRST_EXTERNAL + 1 == 458, "bad finite census")
    ranked = sorted(rows, key=lambda row: (-row[0], row[1]))
    observed_top = tuple((w, coverage) for coverage, w, _ in ranked[:4])
    require(observed_top == EXPECTED_TOP, "top-four finite coverage list changed")
    require(ranked[3][0] > ranked[4][0], "fourth coverage is not isolated")

    fourth = EXPECTED_TOP[3][1]

    def tail_cap(w: int) -> F:
        return root_m / 7 + core.S2 * root_r / (7 * w)

    require(tail_cap(FINITE_MAX) > fourth, "finite cutoff is not sharp below")
    require(tail_cap(TAIL_FIRST) < fourth, "tail cap does not clear fourth coverage")
    tail_threshold = core.S2 * root_r / (7 * (fourth - root_m / 7))
    require(
        F(FINITE_MAX) < tail_threshold < F(TAIL_FIRST),
        "tail threshold is not between 472 and 473",
    )

    top_sum = sum((coverage for _, coverage in EXPECTED_TOP), F(0))
    uniform_margin = root_m - top_sum
    require(top_sum == F(308603791, 1873980108), "top-four sum changed")
    require(uniform_margin == F(135228493, 18739801080), "uniform margin changed")
    require(uniform_margin > 0, "top-four union lower bound is not positive")

    extremal_speeds = tuple(sorted(w for w, _ in EXPECTED_TOP))
    extremal_family = tuple(sorted((*E, *extremal_speeds)))
    require(len(extremal_family) == 13 and len(set(extremal_family)) == 13, "bad extremal family")
    extremal_good, extremal_r, extremal_m = core.good_norm(extremal_family)
    nested = good
    for w in extremal_speeds:
        _, _, nested = core.subtract(nested, w)
    nested_m = sum((right - left for left, right in nested), F(0))
    require(extremal_m == nested_m, "extremal direct/nested family mismatch")
    require(extremal_r == len(nested), "extremal component mismatch")
    require(extremal_r == 16 and extremal_m == F(805109, 16702140), "extremal control changed")
    require(extremal_m >= uniform_margin, "union lower bound exceeds exact extremal survivor")

    reciprocal_packet = sum(
        (F(1, w) for w in range(FIRST_EXTERNAL, FIRST_EXTERNAL + 4)),
        F(0),
    )
    old_root_bound = (
        root_m
        - F(4, 7) * root_m
        - core.S2 * root_r * reciprocal_packet / 7
    )
    require(old_root_bound < 0, "hostile control no longer defeats the old root envelope")

    manifest_sha256 = hashlib.sha256("\n".join(manifest_rows).encode()).hexdigest()
    source_sha256 = sha256(Path(__file__))
    print("THM-741 PURE (3,4) FLOOD TAIL: EXACT TOP-FOUR COVERAGE ENVELOPE")
    print("=" * 88)
    print(f"dependency_sha256={CORE_SHA256}")
    print(f"body E={E}; external universe: all integers 15<=a<b<c<d")
    print(f"root r={root_r} m={root_m}")
    print(
        "finite exact bank: w=15..472 (458 speeds); "
        "sparse/full/direct-union/direct-tooth paths agree"
    )
    for rank, (w, coverage) in enumerate(EXPECTED_TOP, start=1):
        print(f"coverage rank {rank}: w={w} c(w)={coverage}")
    print(f"finite-coverage manifest SHA256={manifest_sha256}")
    print(
        "tail law: c(w)<u(w)=m/7+(99/70)r/(7w) for w>=473; "
        f"threshold={tail_threshold}"
    )
    print(
        f"cutoff hostility: u(472)-c(32)={tail_cap(472)-fourth}; "
        f"c(32)-u(473)={fourth-tail_cap(473)}"
    )
    print(f"global top-four cover sum<={top_sum}")
    print(f"uniform union-bound survivor>={uniform_margin}>0")
    print(
        f"extremal individual-cover family={extremal_family}; "
        f"exact survivor r={extremal_r} m={extremal_m}"
    )
    print(
        "hostile old-envelope control: B11 at distinct starts 15,16,17,18 "
        f"gives {old_root_bound}<0"
    )
    print(
        "proof partition: exact c(w) for 15<=w<=472; THM-735(ii) "
        "(via THM-731/732) monotone cap for every w>=473; four-speed union bound"
    )
    print(
        "scope: one literal root edge (3,4), no Fano/chi_7 transport; "
        "global THM-741 and LRC(14) remain open"
    )
    print(
        "VERDICT: E_(3,4) union {a,b,c,d} is lonely for every "
        "15<=a<b<c<d; this closes a fourth whole flood body"
    )
    print(f"source_sha256={source_sha256}")
    print("ALL EXACT TOP-FOUR COVERAGE CHECKS PASSED")


if __name__ == "__main__":
    main()
