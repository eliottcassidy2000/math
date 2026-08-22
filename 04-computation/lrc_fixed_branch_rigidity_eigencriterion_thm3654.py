#!/usr/bin/env python3
"""Exact companion for THM-3654's literal fixed-branch eigencriterion.

Inside THM-3636's four-space R, literal restriction to branch r=6 acts with
the alpha eigenvalue exactly on K and with the beta eigenvalue on the middle
quotient.  Its defect is an invertible reparameterization of LambdaSharp.
This is a static pinned finite-field statement, not chronology or LRC(14).
"""

from __future__ import annotations

import ast
from hashlib import sha256
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COMPUTATION = ROOT / "04-computation"
sys.path.insert(0, str(COMPUTATION))

PARENT_FILES = (
    (
        "THM3647_theorem",
        ROOT / "01-canon/theorems/THM-3647-lrc-single-reversal-paired-branch-spectral-projector.md",
        "40bd2d192cf5c9d22b41c34b2bd0ad91c29b37658bedc23923b31682371c47a5",
    ),
    (
        "THM3647_script",
        COMPUTATION / "lrc_single_reversal_paired_branch_projector_thm3647.py",
        "2e42c26d43e5d94b5a5c33c6b425dc64f72b7cd38e0759d8dba61aab3fb4c11b",
    ),
    (
        "THM3647_output",
        ROOT / "05-knowledge/results/lrc_single_reversal_paired_branch_projector_thm3647.out",
        "300a43fcc6cbc679d0ec91f571d19115c292555d1ae1fa92db64e65678002950",
    ),
)
EXPECTED_SEMANTIC_SHA256 = "4cab3b1ef03eebf75afcd7f968bf2e1de7729711de3c519261b645fc60000369"
CHECKS = 0


def require(condition: bool, payload: object) -> None:
    global CHECKS
    if condition is not True:
        raise RuntimeError(payload)
    CHECKS += 1


def raw_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def digest(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()


for label, path, expected_hash in PARENT_FILES:
    require(raw_sha256(path) == expected_hash, (label, "parent drift"))

import lrc_single_reversal_paired_branch_projector_thm3647 as T


D, A, M = T.D, T.A, T.M
MOD = A.MOD
WIDTH = M.WIDTH


def combine(coefficients, rows):
    return tuple(
        sum(coefficients[index] * rows[index][column]
            for index in range(len(rows))) % MOD
        for column in range(len(rows[0]))
    )


def multiply_two(row, matrix):
    return tuple(
        sum(row[index] * matrix[index][column] for index in range(2)) % MOD
        for column in range(2)
    )


def main() -> None:
    zeta, generators, addresses = A.build_addresses()
    identity = A.identity_matrix()
    a6 = addresses[6]
    alpha, beta = T.EXPECTED_RECORDS[6][-2:]
    pi_w = D.diagonal((1, 1, 0, 0, 1, 1))
    pi_m = A.matsub(identity, pi_w)
    require(a6 == A.matadd(A.scalar(pi_w, alpha), A.scalar(pi_m, beta)),
            "A6 spectrum")

    k2, r4, _ = M.C.reconstruct_relation_flag()
    centered_generators = M.centre_relation_sharp(generators)
    qsharp = M.C.canonical_row_basis(centered_generators, WIDTH)
    ksharp = M.lift_basis(k2, qsharp, M.mu_rows(qsharp))
    rsharp = M.lift_basis(r4, qsharp, M.mu_rows(qsharp))
    coords = lambda rows: tuple(
        tuple(M.C.coordinates_in_row_basis(centered_generators, row))
        for row in rows
    )
    k_rows = D.canonical(coords(ksharp))
    r_rows_ordered = coords(rsharp)
    r_rows = D.canonical(r_rows_ordered)
    expected_k = D.canonical(((1, 1, 0, 0, 0, 0),
                              (0, 0, 0, 0, 1, 1)))
    expected_r = D.canonical(expected_k + ((0, 0, 1, 0, 0, 0),
                                            (0, 0, 0, 1, 0, 0)))
    expected_middle = D.canonical(((0, 0, 1, 0, 0, 0),
                                   (0, 0, 0, 1, 0, 0)))
    require(k_rows == expected_k and r_rows == expected_r, "flag coordinates")

    defect_operator = A.matsub(a6, A.scalar(identity, alpha))
    beta_operator = A.matsub(a6, A.scalar(identity, beta))
    require(defect_operator == A.scalar(pi_m, (beta - alpha) % MOD),
            "spectral defect identity")
    require(beta_operator == A.scalar(pi_w, (alpha - beta) % MOD),
            "complementary spectral identity")

    defect_rows = tuple(A.row_times_matrix(row, defect_operator)
                        for row in r_rows_ordered)
    beta_rows = tuple(A.row_times_matrix(row, beta_operator)
                      for row in r_rows_ordered)
    require(M.linear_map_kernel(r_rows_ordered, defect_rows, 6) == k_rows,
            "alpha eigenspace inside R")
    require(M.linear_map_kernel(r_rows_ordered, beta_rows, 6) == expected_middle,
            "beta eigenspace inside R")
    require(D.canonical(tuple(A.row_times_matrix(row, pi_w)
                              for row in r_rows_ordered)) == k_rows,
            "PiW image on R")
    require(D.canonical(tuple(A.row_times_matrix(row, pi_m)
                              for row in r_rows_ordered)) == expected_middle,
            "PiM image on R")

    # Literal branch restriction in the parent four-way tensor.  Centering in
    # relation commutes with restriction, so A6 acts on the centered frame too.
    (_ga, _gs, _gamma, point_cores, _support, _work, _states, _branches) = (
        M.F.build_banks()
    )
    point_branch = M.F.invert_point_cores(point_cores, zeta)
    branch6_raw = tuple(
        tuple(point_branch[point][6][difference][relation]
              for difference in range(13) for relation in range(13))
        for point in range(6)
    )
    branch6_centered = M.centre_relation_sharp(branch6_raw)
    require(all(branch6_centered[point] == combine(a6[point], centered_generators)
                for point in range(6)), "literal branch-six intertwiner")

    literal_defects = tuple(
        tuple((combine(row, branch6_centered)[column]
               - alpha * combine(row, centered_generators)[column]) % MOD
              for column in range(WIDTH))
        for row in r_rows_ordered
    )
    coordinate_defects = tuple(
        combine(A.row_times_matrix(row, defect_operator), centered_generators)
        for row in r_rows_ordered
    )
    require(literal_defects == coordinate_defects, "literal spectral defects")
    require(M.linear_map_kernel(r_rows_ordered, literal_defects, WIDTH) == k_rows,
            "literal branch-six eigencriterion")

    # Rebuild LambdaSharp and prove that the literal branch defect is the same
    # quotient map up to the invertible scale beta-alpha and THM-3636's T_mid.
    psharp = M.C.canonical_row_basis(generators, WIDTH)
    centered_p = M.centre_relation_sharp(psharp)
    augmentation_images = tuple(
        tuple((raw - centered) % MOD
              for raw, centered in zip(raw_row, centered_row, strict=True))
        for raw_row, centered_row in zip(psharp, centered_p, strict=True)
    )
    augmentation_basis = M.C.canonical_row_basis(augmentation_images, WIDTH)
    lambda_on_r = tuple(
        M.C.combine_rows(
            M.C.coordinates_in_row_basis(centered_p, row), augmentation_images
        )
        for row in rsharp
    )
    lambda_coordinates = tuple(
        tuple(M.C.coordinates_in_row_basis(augmentation_basis, row))
        for row in lambda_on_r
    )
    defect_middle = tuple((row[2], row[3]) for row in defect_rows)
    scale_inverse = pow((beta - alpha) % MOD, -1, MOD)
    predicted_lambda = tuple(
        multiply_two(tuple(scale_inverse * value % MOD for value in row),
                     D.EXPECTED_T_MID)
        for row in defect_middle
    )
    require(predicted_lambda == lambda_coordinates,
            "branch defect factors LambdaSharp")

    semantic_record = {
        "parents": tuple((label, expected_hash)
                         for label, _path, expected_hash in PARENT_FILES),
        "field": MOD,
        "branch": 6,
        "spectra": (alpha, beta),
        "dimensions": {"R": 4, "K": 2, "middle": 2},
        "digests": {
            "K": D.digest(k_rows),
            "R": D.digest(r_rows),
            "defect_operator": D.digest(defect_operator),
            "literal_defects": D.digest(literal_defects),
            "lambda_coordinates": D.digest(lambda_coordinates),
        },
        "scope": "static literal branch restriction over pinned Fp; no chronology/current/char0/LRC14",
    }
    semantic = digest(semantic_record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("semantic digest", semantic))

    source = Path(__file__).resolve()
    source_bytes = source.read_bytes()
    require(b"\r\n" not in source_bytes, "source raw LF")
    require(not any(isinstance(node, ast.Assert)
                    for node in ast.walk(ast.parse(source_bytes.decode("utf-8")))),
            "Python assert node present")

    print("== THM-3654 fixed-branch rigidity eigencriterion ==")
    print("parent_sha256_raw=" + repr(semantic_record["parents"]))
    print(f"field=p:{MOD};branch=6;spectra=(alpha:{alpha},beta:{beta})")
    print("flag=R:4=K:2+middle:2")
    print("operator=A6-alpha*I=(beta-alpha)*Pi_M")
    print("criterion=c_in_K iff literal_branch6(c)=alpha*branch_sum(c)")
    print("complement=beta_eigenspace_inside_R=middle_plane")
    print("quotient=LambdaSharp=defect_middle/(beta-alpha)*T_mid;T_mid_invertible")
    print("digests=" + repr(semantic_record["digests"]))
    print("semantic_sha256=" + semantic)
    print(f"CHECKS={CHECKS};imported_address_checks={A.CHECKS}")
    print("status=FINITE-EXACT STATIC LITERAL-BRANCH EIGENCRITERION")
    print("scope=no chronology/current/characteristic-zero/row-exclusion/LRC14")


if __name__ == "__main__":
    main()
