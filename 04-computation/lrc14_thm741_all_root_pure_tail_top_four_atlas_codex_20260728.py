#!/usr/bin/env python3
r"""Exact top-four pure-tail atlas for all 2,002 THM-741 root bodies.

Let ``E`` be any nine-element subset of ``{1,...,14}``, let ``G(E)`` have
measure ``m_E`` and ``r_E`` components, and put

    c_E(w) = |G(E) intersect D_w|.

This companion treats the *pure tail* in which the four speeds added to
``E`` are distinct integers at least 15.  The THM-735(ii) discrepancy bound
and ``99/70 > sqrt(2)`` give

    c_E(w) < u_E(w) := m_E/7 + (99/70) r_E/(7w).

For every one of the ``binom(14,9)=2,002`` roots, the script evaluates all
single-comb coverages for ``15<=w<=600``.  If ``q_4`` is the fourth-largest
coverage in this finite bank, it verifies

    u_E(601) < q_4.

Since ``u_E`` decreases, the finite first four ranks are therefore the
global first four over every integer ``w>=15``.  Distinctness of the added
speeds and the union bound then prove the pure tail whenever

    m_E - q_1 - q_2 - q_3 - q_4 > 0.

Exactly 584 roots pass this strict test.  This is a finite exact stratum of
THM-741, not a proof of any root with nonpositive top-four margin, not a
uniform four-extension theorem for those 584 bodies when an added speed is
at most 14, and not LRC(14).

The full 1,173,172-entry coverage atlas is bound by a canonical SHA-256
manifest.  Every root's four decisive coverages are independently rebuilt
by sparse subtraction, full subtraction, direct ten-comb union, and a
two-pointer tooth-incidence sum.  For all 584 positive rows, the literal
thirteen-comb extremal family is also rebuilt by direct union and full nested
subtraction.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
import os
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CORE_PATH = ROOT / "04-computation/lrc14_thm741_2002_body_j4_tree_kps_S128c5.py"
CORE_SHA256 = "5aa81d9d78273c8f9e3e7a6574091a3bc3f64ab6086c7024c15f9420c99dac96"
FIRST_EXTERNAL = 15
FINITE_MAX = 600
TAIL_FIRST = FINITE_MAX + 1
BODIES = tuple(combinations(range(1, 15), 9))

EXPECTED_CLOSED_COUNT = 584
EXPECTED_COVERAGE_MANIFEST = (
    "63fd4a08965c9e5f2665dde27ae8db792cd446bc4578290ab4aa14b01fc469f7"
)
EXPECTED_CLOSED_BODY_DIGEST = (
    "93ed30b15748a90ea78aa1f392e2335cb284fa683f98e8d068b8a4f0f6af7a54"
)
EXPECTED_MINIMUM_POSITIVE = (
    (1, 2, 4, 5, 6, 8, 11, 12, 14),
    F(47, 1513512),
)
EXPECTED_WORST_NEGATIVE = (
    (1, 3, 4, 5, 7, 8, 10, 11, 13),
    -F(9158777, 174053880),
)
EXPECTED_MAX_THRESHOLD = (
    (1, 2, 4, 8, 10, 11, 12, 13, 14),
    F(39783744, 67829),
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_core():
    require(file_sha256(CORE_PATH) == CORE_SHA256, "THM-741 core hash changed")
    spec = importlib.util.spec_from_file_location(
        "thm741_all_root_top_four_dependency",
        CORE_PATH,
    )
    require(spec is not None and spec.loader is not None, "cannot load THM-741 core")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


CORE = load_core()


def fraction_text(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def direct_tooth_coverage(good: list[tuple[F, F]], w: int) -> F:
    """Independent intersection sum over the disjoint teeth of ``D_w``."""

    require(w >= 1, "speed must be positive")
    radius = F(1, 14 * w)
    total = F(0)
    first_component = 0
    for k in range(w + 1):
        center = F(k, w)
        left = max(F(0), center - radius)
        right = min(F(1), center + radius)
        while (
            first_component < len(good)
            and good[first_component][1] <= left
        ):
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


def scan_body(body: tuple[int, ...]) -> dict[str, object]:
    """Compute one exact finite row and its tail certificate."""

    good, root_r, root_m = CORE.good_norm(body)
    require(root_r == len(good) and root_m > 0, f"bad root carrier {body}")

    rows: list[tuple[F, int]] = []
    row_manifest = hashlib.sha256()
    row_manifest.update(("E=" + ",".join(map(str, body)) + "\n").encode())
    for w in range(FIRST_EXTERNAL, FINITE_MAX + 1):
        coverage = root_m - CORE.subtract_sparse(good, w)
        require(F(0) <= coverage < root_m, f"bad coverage {body},w={w}")
        rows.append((coverage, w))
        row_manifest.update(f"{w}:{fraction_text(coverage)}\n".encode())

    require(
        len(rows) == FINITE_MAX - FIRST_EXTERNAL + 1,
        f"finite row length changed {body}",
    )
    ranked = sorted(rows, key=lambda item: (-item[0], item[1]))
    fourth = ranked[3][0]
    require(fourth > root_m / 7, f"fourth rank misses limiting cap {body}")
    threshold = CORE.S2 * root_r / (7 * (fourth - root_m / 7))

    def tail_cap(w: int) -> F:
        return root_m / 7 + CORE.S2 * root_r / (7 * w)

    require(
        threshold < TAIL_FIRST and tail_cap(TAIL_FIRST) < fourth,
        f"finite horizon does not seal tail {body}",
    )
    margin = root_m - sum((row[0] for row in ranked[:4]), F(0))
    return {
        "body": body,
        "r": root_r,
        "m": root_m,
        "margin": margin,
        "threshold": threshold,
        "top": tuple(ranked[:4]),
        "row_manifest": row_manifest.hexdigest(),
    }


def independent_decisive_checks(row: dict[str, object]) -> tuple[int, int]:
    """Cross-check the decisive atoms and, for positive rows, their family."""

    body = row["body"]
    require(isinstance(body, tuple), "body type changed")
    good, root_r, root_m = CORE.good_norm(body)
    require(
        root_r == row["r"] and root_m == row["m"],
        f"root replay changed {body}",
    )
    top = row["top"]
    require(isinstance(top, tuple), "top row type changed")

    atom_checks = 0
    for coverage, w in top:
        sparse_survivor = CORE.subtract_sparse(good, w)
        full_r, full_survivor, full_good = CORE.subtract(good, w)
        direct_good, direct_r, direct_survivor = CORE.good_norm((*body, w))
        tooth_coverage = direct_tooth_coverage(good, w)
        require(
            root_m - sparse_survivor == coverage,
            f"ranked sparse replay changed {body},w={w}",
        )
        require(
            sparse_survivor == full_survivor == direct_survivor,
            f"ranked survivor paths disagree {body},w={w}",
        )
        require(
            full_r == direct_r and full_good == direct_good,
            f"ranked carrier paths disagree {body},w={w}",
        )
        require(
            coverage == tooth_coverage,
            f"ranked tooth incidence disagrees {body},w={w}",
        )
        atom_checks += 1

    family_checks = 0
    margin = row["margin"]
    require(isinstance(margin, F), "margin type changed")
    if margin > 0:
        speeds = tuple(sorted(w for _, w in top))
        family = tuple(sorted((*body, *speeds)))
        require(
            len(family) == 13 and len(set(family)) == 13,
            f"bad extremal family {body}",
        )
        direct_good, direct_r, direct_m = CORE.good_norm(family)
        nested = good
        for w in speeds:
            _, _, nested = CORE.subtract(nested, w)
        nested_m = sum((right - left for left, right in nested), F(0))
        require(
            direct_good == nested
            and direct_r == len(nested)
            and direct_m == nested_m,
            f"extremal family paths disagree {body}",
        )
        require(
            direct_m >= margin > 0,
            f"extremal family violates union margin {body}",
        )
        family_checks = 1
    return atom_checks, family_checks


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--workers",
        type=int,
        default=min(10, os.cpu_count() or 1),
    )
    args = parser.parse_args()
    require(args.workers >= 1, "workers must be positive")
    require(len(BODIES) == 2002 and len(set(BODIES)) == 2002, "root census")

    if args.workers == 1:
        rows = list(map(scan_body, BODIES))
    else:
        # A file-backed main module makes spawn portable; ``map`` preserves the
        # canonical lexicographic body order independently of worker timing.
        with mp.get_context("spawn").Pool(args.workers) as pool:
            rows = pool.map(scan_body, BODIES, chunksize=2)

    require(
        tuple(row["body"] for row in rows) == BODIES,
        "parallel scan changed canonical body order",
    )
    positive = [row for row in rows if row["margin"] > 0]
    zero = [row for row in rows if row["margin"] == 0]
    negative = [row for row in rows if row["margin"] < 0]
    require(
        len(positive) == EXPECTED_CLOSED_COUNT
        and not zero
        and len(negative) == len(BODIES) - EXPECTED_CLOSED_COUNT,
        "top-four sign census changed",
    )

    coverage_manifest = hashlib.sha256()
    coverage_manifest.update(b"THM741/all-2002-pure-tail-top4-atlas/W600/v1\n")
    for row in rows:
        body = row["body"]
        coverage_manifest.update(
            (
                ",".join(map(str, body))
                + ":"
                + str(row["row_manifest"])
                + "\n"
            ).encode()
        )
    require(
        coverage_manifest.hexdigest() == EXPECTED_COVERAGE_MANIFEST,
        "coverage atlas manifest changed",
    )

    closed_body_text = "\n".join(
        ",".join(map(str, row["body"])) for row in positive
    )
    closed_body_digest = hashlib.sha256(closed_body_text.encode()).hexdigest()
    require(
        closed_body_digest == EXPECTED_CLOSED_BODY_DIGEST,
        "positive-root membership changed",
    )

    minimum_positive = min(positive, key=lambda row: row["margin"])
    worst_negative = min(negative, key=lambda row: row["margin"])
    maximum_threshold = max(rows, key=lambda row: row["threshold"])
    require(
        (
            minimum_positive["body"],
            minimum_positive["margin"],
        )
        == EXPECTED_MINIMUM_POSITIVE,
        "minimum positive row changed",
    )
    require(
        (worst_negative["body"], worst_negative["margin"])
        == EXPECTED_WORST_NEGATIVE,
        "worst negative control changed",
    )
    require(
        (maximum_threshold["body"], maximum_threshold["threshold"])
        == EXPECTED_MAX_THRESHOLD,
        "maximum cap crossing changed",
    )

    atom_checks = 0
    family_checks = 0
    for row in rows:
        atoms, families = independent_decisive_checks(row)
        atom_checks += atoms
        family_checks += families
    require(atom_checks == 4 * len(BODIES), "decisive atom check count")
    require(family_checks == len(positive), "extremal family check count")

    print("THM-741 ALL-ROOT PURE-TAIL TOP-FOUR ATLAS")
    print("status=FINITE-EXACT+GLOBAL-TAIL-SEALED")
    print("root_universe=C(14,9)=2002; added_speeds=15<=a<b<c<d")
    print("finite_horizon=15..600; tail_first=601; all_2002_tail_caps_below_q4")
    print("coverage_entries=1173172")
    print(
        "top4_signs="
        f"positive:{len(positive)},zero:{len(zero)},nonpositive:{len(negative)}"
    )
    print(
        "minimum_positive="
        + ",".join(map(str, minimum_positive["body"]))
        + ":"
        + fraction_text(minimum_positive["margin"])
    )
    print(
        "worst_negative_control="
        + ",".join(map(str, worst_negative["body"]))
        + ":"
        + fraction_text(worst_negative["margin"])
    )
    print(
        "maximum_tail_threshold="
        + ",".join(map(str, maximum_threshold["body"]))
        + ":"
        + fraction_text(maximum_threshold["threshold"])
        + "<601"
    )
    print(f"decisive_four_path_atom_checks={atom_checks}")
    print(f"positive_extremal_family_crosschecks={family_checks}")
    print(f"coverage_manifest_sha256={coverage_manifest.hexdigest()}")
    print(f"closed_body_digest_sha256={closed_body_digest}")
    print("scope=584 pure-tail roots only; global THM-741 and LRC14 remain open")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
