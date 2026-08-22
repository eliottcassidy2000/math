#!/usr/bin/env python3
"""Exact de Rham refinement of the THM-3627 non-even hostile boundary.

For the fixed hostile polynomial Q_h, stack the arbitrary-target-two-form
pullback equations with the exact closedness equations

    P_w - K_epsilon + R_c = 0

for Omega=P dc^depsilon+K dc^dw+R depsilon^dw.  The combined rational systems
remain compatible through total source degree five.  THM-3627 independently
proves that the larger arbitrary-form system fails at degree six.

Every mathematical gate is optimization-safe; there are no Python assert
statements.  This companion proves local finite-jet survival only, not a
global polynomial Darboux pair or a Jacobian counterexample.
"""

from __future__ import annotations

import ast
from contextlib import redirect_stdout
from hashlib import sha256
import importlib
import io
import json
from math import comb
from pathlib import Path
import sys

import sympy as sp
from sympy.polys.matrices import DomainMatrix


ROOT = Path(__file__).resolve().parents[1]
CANON = ROOT / "01-canon/theorems"
COMPUTATION = ROOT / "04-computation"
RESULTS = ROOT / "05-knowledge/results"

PARENT_FILES = (
    (
        "thm3624_theorem",
        CANON / "THM-3624-russell-cylinder-noneven-fold-weighted-cokernel-boundary.md",
        "d805a9318f6782f68d9e958e93a57430d3dce2159c98afa74635b4405da620e5",
    ),
    (
        "thm3624_script",
        COMPUTATION / "jc2_russell_cylinder_noneven_fold_weighted_cokernel_thm3624.py",
        "18b7bdc062909bd45d271c6e1744ef4c1ccc360385e769db419e36418dac3468",
    ),
    (
        "thm3624_output",
        RESULTS / "jc2_russell_cylinder_noneven_fold_weighted_cokernel_thm3624.out",
        "6a38a8b220bbc40b77062f92e75dc4f3272361f2fb8e39c77d00786b46baa462",
    ),
    (
        "thm3627_theorem",
        CANON / "THM-3627-russell-cylinder-noneven-hostile-degree-six-closure.md",
        "e6d9156cc6b0af1cc109749e5c31067ffb97fa01b4c7c912be28d74b9abd98de",
    ),
    (
        "thm3627_script",
        COMPUTATION / "jc2_russell_cylinder_noneven_hostile_degree_six_closure_thm3627.py",
        "fcafd6fe3f3ce42432958a4ee074ec2150b943d3dc30552213a8b9e8cd3a0689",
    ),
    (
        "thm3627_output",
        RESULTS / "jc2_russell_cylinder_noneven_hostile_degree_six_closure_thm3627.out",
        "c46ee58088554811aac00ce013e55d3d0b67f3d958e25ee55b0f24af125de0a4",
    ),
)

EXPECTED_PULLBACK_RANKS = (2, 7, 15, 26, 40, 57)
EXPECTED_CLOSED_SYSTEM_RANKS = (2, 8, 19, 36, 60, 92)
EXPECTED_AFFINE_SOLUTION_DIMENSIONS = (1, 4, 11, 24, 45, 76)
EXPECTED_SEMANTIC_SHA256 = (
    "f8ce8a774ffca01430a00519e6672d30c7458e51829037e1daf731163c138796"
)


CHECKS = 0


def require(condition, payload):
    global CHECKS
    if condition is not True:
        raise RuntimeError(payload)
    CHECKS += 1


def lf_sha256(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_json(value):
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")
    return sha256(encoded).hexdigest()


for label, path, expected in PARENT_FILES:
    observed = lf_sha256(path)
    require(observed == expected, ("parent drift", label, observed, expected))

sys.path.insert(0, str(COMPUTATION))
with redirect_stdout(io.StringIO()):
    M = importlib.import_module(
        "jc2_russell_cylinder_noneven_fold_weighted_cokernel_thm3624"
    )


def closedness_matrix(maximum_degree):
    """Exact coefficient matrix of dOmega=0 for target degree <= N."""
    exponents = M.exponent_triples(maximum_degree)
    index = {monomial: position for position, monomial in enumerate(exponents)}
    row_exponents = (
        M.exponent_triples(maximum_degree - 1) if maximum_degree else ()
    )
    width = 3 * len(exponents)
    rows = []
    for a_degree, b_degree, w_degree in row_exponents:
        row = [0] * width

        # P_w in d(P dc wedge depsilon).
        source = (a_degree, b_degree, w_degree + 1)
        if source in index:
            row[index[source]] += w_degree + 1

        # -K_epsilon in d(K dc wedge dw).
        source = (a_degree, b_degree + 1, w_degree)
        if source in index:
            row[len(exponents) + index[source]] -= b_degree + 1

        # R_c in d(R depsilon wedge dw).
        source = (a_degree + 1, b_degree, w_degree)
        if source in index:
            row[2 * len(exponents) + index[source]] += a_degree + 1

        rows.append(row)
    return sp.Matrix(rows) if rows else sp.zeros(0, width)


def main():
    local_series = tuple(
        sp.expand(M.Q_h.subs(M.x, point + M.xi)) for point in (-1, 0, 1)
    )
    records = []
    for degree in range(6):
        pullback, target = M.pullback_matrix_from_local_series(
            local_series, degree
        )
        leak = closedness_matrix(degree)
        system = pullback.col_join(leak)
        rhs = target.col_join(sp.zeros(leak.rows, 1))

        pullback_rank = DomainMatrix.from_Matrix(pullback).rank()
        leak_rank = DomainMatrix.from_Matrix(leak).rank() if leak.rows else 0
        rank = DomainMatrix.from_Matrix(system).rank()
        augmented_rank = DomainMatrix.from_Matrix(system.row_join(rhs)).rank()

        target_monomials = comb(degree + 3, 3)
        source_monomials = comb(degree + 2, 2)
        leak_rows = comb(degree + 2, 3)
        expected_pullback_shape = (3 * source_monomials, 3 * target_monomials)
        expected_leak_shape = (leak_rows, 3 * target_monomials)

        require(
            pullback.shape == expected_pullback_shape,
            ("pullback shape", degree, pullback.shape, expected_pullback_shape),
        )
        require(
            leak.shape == expected_leak_shape,
            ("leak shape", degree, leak.shape, expected_leak_shape),
        )
        require(
            pullback_rank == EXPECTED_PULLBACK_RANKS[degree],
            ("pullback rank", degree, pullback_rank),
        )
        require(leak_rank == leak_rows, ("leak full row rank", degree, leak_rank))
        require(
            rank == EXPECTED_CLOSED_SYSTEM_RANKS[degree],
            ("closed system rank", degree, rank),
        )
        require(
            rank == pullback_rank + leak_rank,
            ("closedness transverse to pullback rows", degree),
        )
        require(
            augmented_rank == rank,
            ("closed-form compatibility", degree, rank, augmented_rank),
        )

        affine_dimension = system.cols - rank
        require(
            affine_dimension == EXPECTED_AFFINE_SOLUTION_DIMENSIONS[degree],
            ("affine solution dimension", degree, affine_dimension),
        )
        records.append(
            {
                "degree": degree,
                "pullback_shape": pullback.shape,
                "leak_shape": leak.shape,
                "pullback_rank": pullback_rank,
                "leak_rank": leak_rank,
                "closed_rank": rank,
                "augmented_rank": augmented_rank,
                "solution_dimension": affine_dimension,
            }
        )

    # The N=5 target has three nonzero constant source rows, so no solution can
    # have zero target two-form at the common target point.  This is the exact
    # nonvanishing hypothesis used by local presymplectic Darboux.
    require(
        tuple(target[index * comb(7, 2)] for index in range(3)) == (12, 12, 12),
        "three normalized constant rows at N=5",
    )
    require(
        EXPECTED_CLOSED_SYSTEM_RANKS[5] < 3 * comb(8, 3),
        "positive-dimensional rational solution space at N=5",
    )

    semantic_record = {
        "parents": tuple((label, expected) for label, _path, expected in PARENT_FILES),
        "equation": "P_w-K_epsilon+R_c=0",
        "records": records,
        "N5_nonzero_constant_rows": (12, 12, 12),
        "N6_parent_boundary": (77, 78, "arbitrary forms incompatible"),
        "scope": "closed polynomial two-form/local formal Darboux through N=5; no global polynomial pair/JC2",
    }
    semantic = digest_json(semantic_record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(
            semantic == EXPECTED_SEMANTIC_SHA256,
            ("semantic", semantic, EXPECTED_SEMANTIC_SHA256),
        )

    source_path = Path(__file__).resolve()
    source_bytes = source_path.read_bytes()
    require(b"\r\n" not in source_bytes, "source must be raw LF")
    tree = ast.parse(source_bytes.decode("utf-8"))
    require(
        not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
        "Python assert node present",
    )

    print("== THM-3631 non-even closed-form order-five survival ==")
    print("parent_sha256_lf", tuple((label, expected) for label, _path, expected in PARENT_FILES))
    for record in records:
        print(
            "degree",
            record["degree"],
            "pullback",
            record["pullback_shape"],
            "leak",
            record["leak_shape"],
            "ranks",
            (record["pullback_rank"], record["leak_rank"], record["closed_rank"], record["augmented_rank"]),
            "solution_dimension",
            record["solution_dimension"],
        )
    print("closedness_equation P_w-K_epsilon+R_c=0")
    print("N5_constant_rows (12, 12, 12); target_form_at_origin=nonzero")
    print("N6_parent_boundary arbitrary_rank=77 augmented=78; all forms fail")
    print("semantic_sha256", semantic)
    print("script_sha256_lf", lf_sha256(source_path))
    print("CHECKS", CHECKS)
    print("SCOPE local closed-form/formal target-pair jet through N=5; no global polynomial pair/JC2")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
