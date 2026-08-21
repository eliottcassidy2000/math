#!/usr/bin/env python3
"""Exact eleven-cell exceptional flux gate for the LRC carry response."""

from __future__ import annotations

import ast
from collections import Counter
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
PARENT_SCRIPT = ROOT / "04-computation/lrc_exceptional_leakage_boundary_thm3660.py"
PARENT_OUTPUT = ROOT / "05-knowledge/results/lrc_exceptional_leakage_boundary_thm3660.out"
EXPECTED_PARENT_HASHES = (
    "37a7c775096d80c4d8acd9d39fc0812e32545a963ee572ea773157d0ef6804c6",
    "d24cda381808b5ec6009a47fc967f1293486e3ffb95d9153f4e4cd0bd3ed0730",
)
EXPECTED_SEMANTIC_SHA256 = "5606124428d10919903c0d7cff433e2550a5e14010c971e4cdcf261a6ff28342"

P = 13
MOD = 755373809845391722745761
EXPECTED_C = 28939885940104996879767
EXPECTED_HOSTILE_VECTOR = (
    289369708799585224299099,
    644718915560530562860929,
)
EXPECTED_TRANSVERSE = (
    (0, 10), (0, 11), (0, 12),
    (3, 8), (3, 9),
    (6, 4), (6, 5), (6, 6), (6, 7),
    (9, 2), (9, 3),
)
EXPECTED_ZERO = (
    (1, 12), (2, 12), (3, 2), (3, 12), (4, 4), (4, 12),
    (5, 4), (5, 12), (6, 12), (7, 7), (7, 12), (8, 7),
    (8, 12), (9, 9), (9, 12), (10, 12), (11, 12),
)


def require(condition: bool, payload: object) -> None:
    if condition is not True:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_json(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()


def load_parent():
    hashes = (lf_sha256(PARENT_SCRIPT), lf_sha256(PARENT_OUTPUT))
    require(hashes == EXPECTED_PARENT_HASHES, ("THM-3660 parent hashes", hashes))
    spec = importlib.util.spec_from_file_location("thm3660_parent", PARENT_SCRIPT)
    require(spec is not None and spec.loader is not None, "THM-3660 loader")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module, hashes


def split_add(left, right):
    return ((left[0] + right[0]) % P, (left[1] + right[1]) % P)


def carry(left, right):
    return int(left[0] + right[0] >= P)


def main() -> None:
    L, parent_hashes = load_parent()
    T, grandparent_hashes = L.load_parent()
    M = T.M
    tensor, reconstruction = M.reconstruct_two_current()
    require(reconstruction == M.EXPECTED_PARENT_RECONSTRUCTION_SHA256[:2],
            "two-current reconstruction drift")
    raw = M.flatten_two_current(tensor)
    correction = tuple(T.state_side(row) for row in raw)
    basis = M.canonical_row_basis(correction)
    coordinates = tuple(T.coordinates_in_rref(basis, row) for row in correction)
    labels = tuple((r0, r1) for r0 in range(P) for r1 in range(P))
    by_label = dict(zip(labels, coordinates))
    table_labels = tuple((t0, t1) for t0 in range(P - 1) for t1 in range(P))

    def next_vertical(label):
        return label[0], (label[1] + 1) % P

    def subtract(left, right):
        return tuple((a - b) % MOD for a, b in zip(left, right))

    response = {
        label: subtract(by_label[next_vertical(label)], by_label[label])
        for label in table_labels
    }
    zero = tuple(label for label in table_labels if response[label] == (0, 0))
    require(zero == EXPECTED_ZERO, ("response zero set", zero))

    amplitude = L.EXPECTED_AMPLITUDE
    inverse_amplitude = pow(amplitude, -1, MOD)

    def leak(vector):
        return (vector[1] - T.B * vector[0]) % MOD

    detector = {
        label: leak(by_label[label]) * inverse_amplitude % MOD
        for label in labels
    }
    require(frozenset(label for label in labels if detector[label])
            == L.X, "detector support")
    derivative = {
        label: (detector[next_vertical(label)] - detector[label]) % MOD
        for label in table_labels
    }
    transverse = tuple(label for label in table_labels if leak(response[label]))
    require(transverse == EXPECTED_TRANSVERSE, ("transverse cells", transverse))
    require(all(leak(response[label]) * inverse_amplitude % MOD
                == derivative[label] for label in table_labels),
            "projected response formula")

    generic = tuple(label for label in table_labels
                    if label not in transverse and label not in zero)
    require(len(generic) == 128, "generic response count")
    require(all(T.slope(response[label]) == T.B for label in generic),
            "generic response line")
    require((len(set(response[label] for label in generic)),
             len(set(response[label] for label in transverse)))
            == (44, 10), "response value split")
    require(len(set(response.values())) == 55, "response value count")
    transverse_classes = Counter(T.slope(response[label]) for label in transverse)
    require(tuple(sorted(transverse_classes.values(), reverse=True))
            == (2, 1, 1, 1, 1, 1, 1, 1, 1, 1),
            "transverse projective profile")
    repeated = tuple(label for label in transverse
                     if response[label] == response[(0, 10)])
    require(repeated == ((0, 10), (6, 7)), ("repeated transverse row", repeated))

    coordinate_matrices = tuple(
        tuple(tuple(response[(t0, t1)][coordinate] for t1 in range(P))
              for t0 in range(P - 1))
        for coordinate in range(2)
    )
    coordinate_ranks = tuple(M.rank_mod(matrix) for matrix in coordinate_matrices)
    derivative_matrix = tuple(
        tuple(derivative[(t0, t1)] for t1 in range(P))
        for t0 in range(P - 1)
    )
    derivative_rank = M.rank_mod(derivative_matrix)
    require((coordinate_ranks, derivative_rank) == ((11, 11), 4),
            ("response matrix ranks", coordinate_ranks, derivative_rank))

    # Exact reversal-odd projection.
    inverse_two = pow(2, -1, MOD)
    identity = ((1, 0), (0, 1))
    projector_minus = tuple(tuple(
        (identity[row][column] - T.REVERSAL_ACTION[row][column])
        * inverse_two % MOD
        for column in range(2)
    ) for row in range(2))

    def act(vector, matrix):
        return tuple(sum(vector[row] * matrix[row][column]
                         for row in range(2)) % MOD
                     for column in range(2))

    c = amplitude * pow(2 * T.B % MOD, -1, MOD) % MOD
    require(c == EXPECTED_C, ("odd projector scalar", c))
    minus_vector = (1, (-T.B) % MOD)
    require(all(
        act(by_label[label], projector_minus)
        == tuple((-c * detector[label] * value) % MOD for value in minus_vector)
        for label in labels
    ), "odd detector projector")
    require(all(
        act(response[label], projector_minus)
        == tuple((-c * derivative[label] * value) % MOD for value in minus_vector)
        for label in table_labels
    ), "odd response projector")

    signed_derivative = {
        label: ({0: 0, 1: 1, MOD - 1: -1, MOD - 2: -2}[derivative[label]])
        for label in table_labels
    }
    flux_support = tuple(
        (label, signed_derivative[label])
        for label in table_labels if signed_derivative[label]
    )
    require(tuple(label for label, _value in flux_support) == EXPECTED_TRANSVERSE,
            "flux support")

    # Constant-in-high-digit carry multiplicities annihilate the full vector
    # response by vertical telescoping, not only its quotient.
    natural_weight_sum = tuple(
        sum(P * (P - 1 - label[0]) * response[label][coordinate]
            for label in table_labels) % MOD
        for coordinate in range(2)
    )
    require(natural_weight_sum == (0, 0), "natural carry-weight cancellation")

    def edge_reversal(label):
        return (12 - label[0]) % P, (11 - label[1]) % P

    hostile_left = (1, 0)
    hostile_right = edge_reversal(hostile_left)
    require(hostile_right == (11, 11), "hostile reversal pair")
    require(response[hostile_left] == EXPECTED_HOSTILE_VECTOR,
            "hostile vector drift")
    require(response[hostile_right]
            == tuple((-value) % MOD for value in EXPECTED_HOSTILE_VECTOR),
            "hostile cancellation")
    require(hostile_left not in transverse and hostile_right not in transverse,
            "hostile pair touches exceptional flux")

    # Interior reversal-symmetric flux has only four independent weights.
    interior_representatives = ((3, 8), (3, 9), (6, 4), (6, 5))
    require(all(edge_reversal(edge_reversal(label)) == label
                for label in table_labels if label[0] != 0),
            "interior edge reversal")
    interior_flux_coefficients = tuple(
        signed_derivative[label]
        + signed_derivative[edge_reversal(label)]
        for label in interior_representatives
    )
    require(interior_flux_coefficients == (-2, 2, 2, -2),
            ("interior symmetric flux", interior_flux_coefficients))

    # The scalar projected carry response is not a split 2-cocycle.
    def scalar_gamma(left, right):
        if not carry(left, right):
            return 0
        return derivative[split_add(left, right)]

    witness = ((0, 1), (1, 0), (12, 9))
    x, y, z = witness
    associator = (
        scalar_gamma(y, z)
        - scalar_gamma(split_add(x, y), z)
        + scalar_gamma(x, split_add(y, z))
        - scalar_gamma(x, y)
    ) % MOD
    require(associator == MOD - 1, ("noncocycle witness", associator))

    semantic = digest_json((
        MOD, P, parent_hashes, grandparent_hashes, reconstruction,
        amplitude, T.B, c,
        tuple((label, response[label]) for label in table_labels),
        zero, generic, transverse, flux_support,
        tuple(coordinate_ranks), derivative_rank,
        projector_minus, natural_weight_sum,
        hostile_left, hostile_right, EXPECTED_HOSTILE_VECTOR,
        interior_flux_coefficients, witness, associator,
    ))
    require(semantic == EXPECTED_SEMANTIC_SHA256,
            ("semantic digest", semantic, EXPECTED_SEMANTIC_SHA256))

    source = Path(__file__).resolve()
    require(not any(isinstance(node, ast.Assert)
                    for node in ast.walk(ast.parse(source.read_text(encoding="utf-8")))),
            "Python assert node present")
    print("== THM-3662 LRC eleven-cell exceptional flux gate ==")
    print(f"field=p:{MOD};carry_table:12x13;parent_sha256_lf={parent_hashes}")
    print(f"grandparent_sha256_lf={grandparent_hashes};reconstruction_sha256={reconstruction}")
    print(f"response_partition=zero:17,generic_line:128,transverse:11;distinct_values=(1,44,10)")
    print(f"transverse_cells={transverse};flux_coefficients={flux_support}")
    print(f"response_matrix_ranks={coordinate_ranks};odd_projected_rank={derivative_rank}")
    print(f"odd_projector=Pi_-;scalar:c={c};Pi_-e=-c*g*(1,-b);Pi_-R=-c*Dg*(1,-b)")
    print(f"natural_carry_weights=13(12-t0);full_response_sum={natural_weight_sum}")
    print(f"positive_reversal_hostile={hostile_left, hostile_right};vector={EXPECTED_HOSTILE_VECTOR};sum:zero")
    print(f"interior_symmetric_flux_coefficients={interior_flux_coefficients};representatives={interior_representatives}")
    print(f"projected_response_2cocycle=REFUTED;witness={witness};associator=-1")
    print(f"response_sha256={digest_json(tuple((label, response[label]) for label in table_labels))}")
    print(f"semantic_sha256={semantic}")
    print(f"script_sha256_lf={lf_sha256(source)}")
    print("scope=static finite-field response/flux;not physical weights/chronology/positivity/LRC14")
    print("PASS")


if __name__ == "__main__":
    main()
