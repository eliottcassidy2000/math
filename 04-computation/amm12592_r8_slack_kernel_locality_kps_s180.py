#!/usr/bin/env python3
"""R=8 AMM slack-kernel locality scout.

For the denominator-103 vertex F and integer entry point I, augment every
displayed inequality A y <= b by its slack.  The scaled displacement

    w=(103(I-F), 103(s_I-s_F))

is an integer kernel vector of [A I].  This script asks whether w has a
conformal decomposition into kernel moves whose y-support lies in at most L
consecutive causal rows.  The discovery LP is small; any positive candidate
is rationalized and checked with Fraction arithmetic before being reported.

This is a locality test for the actual AMM carrier.  It is not a claim about
the full Graver basis, R=512, or an AMM bound.
"""

from __future__ import annotations

import importlib.util
from fractions import Fraction
from functools import reduce
from hashlib import sha256
from math import gcd, lcm
from pathlib import Path

import highspy
import numpy as np
from scipy.optimize import Bounds, LinearConstraint, linprog, milp
from scipy.sparse import csc_matrix, csr_matrix, vstack


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "amm12592_r8_integer_entry_section_kps_s178.py"
SPEC = importlib.util.spec_from_file_location("entry", SOURCE)
ENTRY = importlib.util.module_from_spec(SPEC)
if SPEC.loader is None:
    raise RuntimeError("cannot load THM-3371 companion")
SPEC.loader.exec_module(ENTRY)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def augmented_system():
    profile, meta, inequalities, bounds = ENTRY.build_polytope()
    rows: list[tuple[int, ...]] = []
    rhs: list[int] = []
    names: list[tuple] = []
    for row, value, name in inequalities:
        rows.append(row)
        rhs.append(value)
        names.append(name)
    for j, ((lower, upper), address) in enumerate(zip(bounds, meta)):
        if lower is not None:
            row = tuple(-1 if k == j else 0 for k in range(len(meta)))
            rows.append(row)
            rhs.append(-lower)
            names.append(("bound-lower",) + address)
        if upper is not None:
            row = tuple(1 if k == j else 0 for k in range(len(meta)))
            rows.append(row)
            rhs.append(upper)
            names.append(("bound-upper",) + address)
    return profile, meta, rows, rhs, names


def points():
    integer_rows = (
        (6, -3, 2, 1),
        (-1, 0, 0, 0, 0),
        (-1, 0, 0, 0, 0),
        (0, 0, 0, 0, 0, 0),
        (0, 0, 0, 0, 0, 0, 0),
        (0, 0, 0, 0, 0, 0, 0),
        (0, 0, 0, 0, 0, 0, 0, 0),
    )
    integer = tuple(Fraction(x) for row in integer_rows for x in row)
    raw = (
        "7 -7 1 1 0 -490/103 0 0 0 "
        "-1 -169/103 128/103 -2 1 "
        "-1 140/103 -828/103 911/103 -3 1 "
        "-1 -169/103 894/103 -2150/103 751/103 7 0 "
        "0 -375/103 2270/103 -2801/103 146/103 -73/103 1 "
        "1 346/103 799/103 4969/103 419/103 -419/103 8 0"
    )
    fractional = tuple(Fraction(x) for x in raw.split())
    return integer, fractional


def matvec(rows: list[tuple[int, ...]], vector: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sum(a * x for a, x in zip(row, vector)) for row in rows)


def exact_farkas(a_eq: csr_matrix, b_eq: np.ndarray, a_ub: csr_matrix):
    """Discover with HiGHS and then verify an integer Farkas certificate.

    The system is E z=d, U z<=0, z>=0.  The returned integers satisfy
    lambda*E+mu*U>=0, mu>=0 and lambda*d<0.
    """
    matrix = csc_matrix(vstack([a_eq, a_ub], format="csr"))
    n_eq = a_eq.shape[0]
    model = highspy.HighsLp()
    model.num_col_ = matrix.shape[1]
    model.num_row_ = matrix.shape[0]
    model.col_cost_ = np.zeros(matrix.shape[1])
    model.col_lower_ = np.zeros(matrix.shape[1])
    model.col_upper_ = np.full(matrix.shape[1], highspy.kHighsInf)
    model.row_lower_ = np.concatenate(
        [b_eq, np.full(a_ub.shape[0], -highspy.kHighsInf)]
    )
    model.row_upper_ = np.concatenate([b_eq, np.zeros(a_ub.shape[0])])
    model.a_matrix_.format_ = highspy.MatrixFormat.kColwise
    model.a_matrix_.start_ = matrix.indptr.astype(np.int32)
    model.a_matrix_.index_ = matrix.indices.astype(np.int32)
    model.a_matrix_.value_ = matrix.data
    highs = highspy.Highs()
    highs.setOptionValue("output_flag", False)
    highs.setOptionValue("presolve", "off")
    require(highs.passModel(model) == highspy.HighsStatus.kOk, "passModel")
    require(highs.run() == highspy.HighsStatus.kOk, "HiGHS run")
    require(
        highs.getModelStatus() == highspy.HighsModelStatus.kInfeasible,
        "Farkas model not infeasible",
    )
    status, exists, ray = highs.getDualRay()
    require(status == highspy.HighsStatus.kOk and exists, "missing dual ray")

    fractions = [Fraction(float(-x)).limit_denominator(10**8) for x in ray]
    denominator = reduce(lcm, (x.denominator for x in fractions), 1)
    integers = [int(x * denominator) for x in fractions]
    divisor = reduce(gcd, (abs(x) for x in integers if x), 0)
    integers = [x // divisor for x in integers]
    lambdas = integers[:n_eq]
    mus = integers[n_eq:]
    require(all(x >= 0 for x in mus), "Farkas upper-row sign")

    e = np.asarray(a_eq.toarray(), dtype=np.int64)
    u = np.asarray(a_ub.toarray(), dtype=np.int64)
    coefficients = [
        sum(lambdas[i] * int(e[i, j]) for i in range(len(lambdas)))
        + sum(mus[i] * int(u[i, j]) for i in range(len(mus)))
        for j in range(matrix.shape[1])
    ]
    rhs_value = sum(x * int(y) for x, y in zip(lambdas, b_eq))
    require(all(x >= 0 for x in coefficients), "Farkas column sign")
    require(rhs_value < 0, "Farkas rhs sign")
    digest = sha256(
        ",".join(map(str, integers)).encode("ascii")
    ).hexdigest()
    return {
        "eq_nonzero": sum(x != 0 for x in lambdas),
        "ub_nonzero": sum(x != 0 for x in mus),
        "max_abs": max(abs(x) for x in integers),
        "rhs": rhs_value,
        "sha256": digest,
    }


def solve_width(
    width: int,
    meta: tuple[tuple[int, int], ...],
    rows: list[tuple[int, ...]],
    dy: tuple[int, ...],
    ds: tuple[int, ...],
):
    intervals = [
        (a, b)
        for a in range(7)
        for b in range(a, min(7, a + width))
    ]
    variables: list[tuple[int, int]] = []
    for q, (a, b) in enumerate(intervals):
        for j, (row_index, _) in enumerate(meta):
            if a <= row_index <= b and dy[j]:
                variables.append((q, j))
    index = {address: k for k, address in enumerate(variables)}
    signs = tuple((x > 0) - (x < 0) for x in dy)

    eq_rows: list[list[float]] = []
    eq_rhs: list[float] = []
    for j, target in enumerate(dy):
        if not target:
            continue
        row = [0.0] * len(variables)
        for q in range(len(intervals)):
            k = index.get((q, j))
            if k is not None:
                row[k] = 1.0
        eq_rows.append(row)
        eq_rhs.append(float(abs(target)))

    ub_rows: list[list[float]] = []
    ub_rhs: list[float] = []
    for q in range(len(intervals)):
        for i, target in enumerate(ds):
            row = [0.0] * len(variables)
            for j, coefficient in enumerate(rows[i]):
                k = index.get((q, j))
                if k is not None:
                    # slack displacement is -A*gy, gy=sign(dy)*z.
                    row[k] = -coefficient * signs[j]
            if target > 0:
                ub_rows.append([-x for x in row])
                ub_rhs.append(0.0)
            elif target < 0:
                ub_rows.append(row)
                ub_rhs.append(0.0)
            else:
                eq_rows.append(row)
                eq_rhs.append(0.0)

    a_eq = csr_matrix(np.asarray(eq_rows))
    a_ub = csr_matrix(np.asarray(ub_rows))
    result = linprog(
        np.zeros(len(variables)),
        A_ub=a_ub,
        b_ub=np.asarray(ub_rhs),
        A_eq=a_eq,
        b_eq=np.asarray(eq_rhs),
        bounds=(0, None),
        method="highs",
    )
    if not result.success:
        certificate = exact_farkas(a_eq, np.asarray(eq_rhs), a_ub)
        return intervals, variables, result, None, certificate, None

    rational = tuple(Fraction(float(x)).limit_denominator(10**6) for x in result.x)
    atoms = []
    for q, interval in enumerate(intervals):
        gy = [Fraction(0) for _ in dy]
        for j in range(len(dy)):
            k = index.get((q, j))
            if k is not None:
                gy[j] = signs[j] * rational[k]
        gs = tuple(
            -sum(Fraction(a) * x for a, x in zip(row, gy)) for row in rows
        )
        atoms.append((interval, tuple(gy), gs))

    # Exact reconstruction and conformality audit.
    require(
        all(sum(atom[1][j] for atom in atoms) == dy[j] for j in range(len(dy))),
        "rationalized y reconstruction failed",
    )
    require(
        all(sum(atom[2][i] for atom in atoms) == ds[i] for i in range(len(ds))),
        "rationalized slack reconstruction failed",
    )
    for _, gy, gs in atoms:
        for value, target in zip(gy, dy):
            require(value == 0 or value * target > 0, "y nonconformal")
        for value, target in zip(gs, ds):
            require(value == 0 or value * target > 0, "slack nonconformal")

    integer_result = milp(
        np.zeros(len(variables)),
        integrality=np.ones(len(variables)),
        bounds=Bounds(np.zeros(len(variables)), np.full(len(variables), np.inf)),
        constraints=(
            LinearConstraint(a_eq, np.asarray(eq_rhs), np.asarray(eq_rhs)),
            LinearConstraint(a_ub, np.full(len(ub_rhs), -np.inf), np.zeros(len(ub_rhs))),
        ),
        options={"time_limit": 60.0, "presolve": True},
    )
    integer_atoms = None
    if integer_result.success:
        integer_values = tuple(int(round(x)) for x in integer_result.x)
        require(
            all(abs(x - y) < 1e-6 for x, y in zip(integer_result.x, integer_values)),
            "MILP nonintegral return",
        )
        integer_atoms = []
        for q, interval in enumerate(intervals):
            gy = [0 for _ in dy]
            for j in range(len(dy)):
                k = index.get((q, j))
                if k is not None:
                    gy[j] = signs[j] * integer_values[k]
            gs = tuple(-sum(a * x for a, x in zip(row, gy)) for row in rows)
            integer_atoms.append((interval, tuple(gy), gs))
        require(
            all(
                sum(atom[1][j] for atom in integer_atoms) == dy[j]
                for j in range(len(dy))
            ),
            "integer y reconstruction",
        )
        require(
            all(
                sum(atom[2][i] for atom in integer_atoms) == ds[i]
                for i in range(len(ds))
            ),
            "integer slack reconstruction",
        )
        for _, gy, gs in integer_atoms:
            require(
                all(x == 0 or x * y > 0 for x, y in zip(gy, dy)),
                "integer y nonconformal",
            )
            require(
                all(x == 0 or x * y > 0 for x, y in zip(gs, ds)),
                "integer slack nonconformal",
            )
    return intervals, variables, result, atoms, None, integer_atoms


def main() -> None:
    profile, meta, rows, rhs, names = augmented_system()
    integer, fractional = points()
    require(len(integer) == len(fractional) == len(meta) == 42, "point size")
    require(len(rows) == len(rhs) == len(names) == 114, "augmented row count")

    dy_fraction = tuple(103 * (x - y) for x, y in zip(integer, fractional))
    require(all(x.denominator == 1 for x in dy_fraction), "scaled y lattice")
    dy = tuple(int(x) for x in dy_fraction)
    a_dy = matvec(rows, dy)
    ds = tuple(-x for x in a_dy)

    slack_integer = tuple(
        Fraction(b) - sum(Fraction(a) * x for a, x in zip(row, integer))
        for row, b in zip(rows, rhs)
    )
    slack_fractional = tuple(
        Fraction(b) - sum(Fraction(a) * x for a, x in zip(row, fractional))
        for row, b in zip(rows, rhs)
    )
    require(all(x >= 0 for x in slack_integer + slack_fractional), "feasibility")
    require(
        all(
            Fraction(ds[i]) == 103 * (slack_integer[i] - slack_fractional[i])
            for i in range(len(ds))
        ),
        "slack kernel identity",
    )

    print("AMM12592 R8 SLACK-KERNEL LOCALITY SCOUT")
    print(
        f"carrier=profile={profile};y_columns={len(meta)};"
        f"slack_columns={len(rows)};augmented_shape={len(rows)}x{len(meta)+len(rows)}"
    )
    print(
        f"displacement=y_nonzero={sum(x!=0 for x in dy)};"
        f"slack_nonzero={sum(x!=0 for x in ds)};"
        f"kernel_check={all(x+y==0 for x,y in zip(a_dy,ds))}"
    )
    first = None
    farkas_hashes = []
    integer_atom_hash = None
    for width in range(1, 8):
        intervals, variables, result, atoms, certificate, integer_atoms = solve_width(
            width, meta, rows, dy, ds
        )
        exact = atoms is not None
        nonzero_atoms = (
            sum(any(x for x in gy) for _, gy, _ in atoms) if atoms else None
        )
        print(
            f"width={width};windows={len(intervals)};allocation_variables={len(variables)};"
            f"lp_status={result.status};lp_success={result.success};"
            f"exact_rational_candidate={exact};nonzero_atoms={nonzero_atoms};"
            f"exact_farkas={certificate is not None};"
            f"exact_integer_candidate={integer_atoms is not None}"
        )
        if certificate is not None:
            farkas_hashes.append(certificate["sha256"])
            print(
                f"width_{width}_farkas=eq_nonzero={certificate['eq_nonzero']};"
                f"ub_nonzero={certificate['ub_nonzero']};"
                f"max_abs={certificate['max_abs']};rhs={certificate['rhs']};"
                f"sha256={certificate['sha256']}"
            )
        if exact and first is None:
            atom_denominator = reduce(
                lcm,
                (
                    value.denominator
                    for _, gy, _ in atoms
                    for value in gy
                ),
                1,
            )
            summaries = []
            for interval, gy, gs in atoms:
                if not any(gy):
                    continue
                denominator = reduce(
                    lcm, (value.denominator for value in gy), 1
                )
                summaries.append(
                    f"{interval[0]}-{interval[1]}:y{sum(x!=0 for x in gy)}:"
                    f"s{sum(x!=0 for x in gs)}:q{denominator}"
                )
            print(
                f"width_{width}_candidate=common_denominator={atom_denominator};"
                f"atoms={'|'.join(summaries)}"
            )
        if integer_atoms is not None and first is None:
            summaries = []
            for interval, gy, gs in integer_atoms:
                if any(gy):
                    summaries.append(
                        f"{interval[0]}-{interval[1]}:y{sum(x!=0 for x in gy)}:"
                        f"s{sum(x!=0 for x in gs)}"
                    )
            print(
                f"width_{width}_integer_candidate=atoms={'|'.join(summaries)}"
            )
            atom_text = ";".join(
                f"{a}-{b}:" + ",".join(map(str, gy)) + ":" + ",".join(map(str, gs))
                for (a, b), gy, gs in integer_atoms
            )
            integer_atom_hash = sha256(atom_text.encode("ascii")).hexdigest()
            print(f"width_{width}_integer_atom_sha256={integer_atom_hash}")
        if exact and first is None:
            first = width
    print(f"first_exact_conformal_window_width={first}")
    require(first == 5 and integer_atom_hash is not None, "locality threshold")
    semantic_text = (
        ",".join(map(str, dy))
        + ";"
        + ",".join(map(str, ds))
        + ";"
        + ";".join(farkas_hashes)
        + ";"
        + integer_atom_hash
    )
    source = Path(__file__).read_bytes().replace(b"\r\n", b"\n")
    print(f"semantic_sha256={sha256(semantic_text.encode('ascii')).hexdigest()}")
    print(f"source_lf_sha256={sha256(source).hexdigest()}")
    print(
        "scope=CONFORMAL_WINDOW_ALLOCATION_FOR_ONE_R8_DISPLACEMENT_"
        "NOT_FULL_GRAVER_BASIS_NOT_R512_NOT_AMM_CLOSURE"
    )


if __name__ == "__main__":
    main()
