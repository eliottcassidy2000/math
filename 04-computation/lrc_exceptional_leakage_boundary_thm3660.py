#!/usr/bin/env python3
"""Exact exceptional detector and vertical boundary for THM-3657."""

from __future__ import annotations

import ast
from collections import Counter
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
PARENT_SCRIPT = ROOT / "04-computation/lrc_two_current_quotient_address_reversal_gate_thm3657.py"
PARENT_OUTPUT = ROOT / "05-knowledge/results/lrc_two_current_quotient_address_reversal_gate_thm3657.out"
EXPECTED_PARENT_HASHES = (
    "f0323550a039bd3c59bc3367a9a48503ff10db2b642e99e21d9499e984492ccd",
    "fdebbbba161149edec8c86b9f382ede76153e85aa01da635f0900075fe751628",
)
EXPECTED_SEMANTIC_SHA256 = "7811b56d46c9b6e0abdfcc86d6539b5b0169f5aaebc45b6e4431eed0b8c0b65e"

P = 13
MOD = 755373809845391722745761
EXPECTED_AMPLITUDE = 290125781062402394993132
X_PLUS = frozenset(((12, 0), (0, 11), (6, 5), (9, 3)))
X_MINUS = frozenset(((0, 12), (12, 1), (6, 7), (3, 9)))
X = X_PLUS | X_MINUS
EXPECTED_FULL_SUPPORT = (
    (0, 10), (0, 11), (0, 12),
    (3, 8), (3, 9),
    (6, 4), (6, 5), (6, 6), (6, 7),
    (9, 2), (9, 3),
    (12, 0), (12, 1), (12, 12),
)
EXPECTED_CARRY_SUPPORT = EXPECTED_FULL_SUPPORT[:-3]


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
    require(hashes == EXPECTED_PARENT_HASHES, ("THM-3657 parent hashes", hashes))
    spec = importlib.util.spec_from_file_location("thm3657_parent", PARENT_SCRIPT)
    require(spec is not None and spec.loader is not None, "THM-3657 loader")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module, hashes


def next_vertical(label):
    return label[0], (label[1] + 1) % P


def edge_reversal(label):
    # Reversal sends the oriented edge t -> t+(0,1) to
    # j(t+(0,1)) -> j(t), again in the positive vertical orientation.
    return (12 - label[0]) % P, (11 - label[1]) % P


def signed_lift(value):
    table = {0: 0, 1: 1, MOD - 1: -1, MOD - 2: -2}
    require(value in table, ("unexpected normalized boundary value", value))
    return table[value]


def main() -> None:
    T, parent_hashes = load_parent()
    M = T.M
    tensor, reconstruction = M.reconstruct_two_current()
    require(reconstruction == M.EXPECTED_PARENT_RECONSTRUCTION_SHA256[:2],
            "two-current reconstruction drift")
    raw = M.flatten_two_current(tensor)
    correction = tuple(T.state_side(row) for row in raw)
    basis = M.canonical_row_basis(correction)
    require(len(basis) == 2, "correction rank")
    require(M.rowspace_digest(basis) == T.G.EXPECTED_DIGESTS["state_side"],
            "correction plane digest")
    coordinates = tuple(T.coordinates_in_rref(basis, row) for row in correction)
    labels = tuple((r0, r1) for r0 in range(P) for r1 in range(P))
    by_label = dict(zip(labels, coordinates))

    def raw_leak(vector):
        return (vector[1] - T.B * vector[0]) % MOD

    amplitude = raw_leak(by_label[(12, 0)])
    require(amplitude == EXPECTED_AMPLITUDE, ("leak amplitude", amplitude))
    inverse_amplitude = pow(amplitude, -1, MOD)
    detector = {
        label: raw_leak(by_label[label]) * inverse_amplitude % MOD
        for label in labels
    }
    require(frozenset(label for label in labels if detector[label]) == X,
            "exceptional detector support")
    require(frozenset(label for label in labels if detector[label] == 1) == X_PLUS,
            "positive exceptional detector")
    require(frozenset(label for label in labels if detector[label] == MOD - 1)
            == X_MINUS, "negative exceptional detector")
    require(all(detector[T.reversal(label)] == (-detector[label]) % MOD
                for label in labels), "detector reversal oddness")
    require(all(raw_leak(T.act(vector)) == (-raw_leak(vector)) % MOD
                for vector in coordinates), "quotient action")

    derivative = {
        label: (detector[next_vertical(label)] - detector[label]) % MOD
        for label in labels
    }
    full_support = tuple(label for label in labels if derivative[label])
    incidence = tuple(
        label for label in labels
        if label in X or next_vertical(label) in X
    )
    require(full_support == EXPECTED_FULL_SUPPORT == incidence,
            ("full exceptional boundary", full_support, incidence))
    carry_labels = tuple(label for label in labels if label[0] != 12)
    carry_support = tuple(label for label in carry_labels if derivative[label])
    carry_incidence = tuple(
        label for label in carry_labels
        if label in X or next_vertical(label) in X
    )
    require(carry_support == EXPECTED_CARRY_SUPPORT == carry_incidence,
            ("carry-visible boundary", carry_support, carry_incidence))
    require(all(derivative[edge_reversal(label)] == derivative[label]
                for label in labels), "edge-reversal evenness")
    require(tuple(sum(derivative[(r0, r1)] for r1 in range(P)) % MOD
                  for r0 in range(P)) == (0,) * P,
            "vertical telescoping")

    signed = {label: signed_lift(value) for label, value in derivative.items()}
    full_counts = tuple(sorted(Counter(signed.values()).items()))
    carry_counts = tuple(sorted(Counter(signed[label]
                                        for label in carry_labels).items()))
    require(full_counts == ((-2, 2), (-1, 4), (0, 155), (1, 8)),
            ("full signed counts", full_counts))
    require(carry_counts == ((-2, 1), (-1, 4), (0, 145), (1, 6)),
            ("carry signed counts", carry_counts))
    full_l1 = sum(abs(value) for value in signed.values())
    full_energy = sum(value * value for value in signed.values())
    carry_l1 = sum(abs(signed[label]) for label in carry_labels)
    carry_energy = sum(signed[label] ** 2 for label in carry_labels)
    require((full_l1, full_energy, carry_l1, carry_energy) == (16, 20, 12, 14),
            "boundary norm ledger")

    weights = {label: 13 * (12 - label[0]) for label in carry_labels}
    weighted_nonzero = sum(weights[label] for label in carry_support)
    weighted_signed = sum(weights[label] * signed[label] for label in carry_labels)
    weighted_l1 = sum(weights[label] * abs(signed[label])
                      for label in carry_labels)
    weighted_energy = sum(weights[label] * signed[label] ** 2
                          for label in carry_labels)
    require((weighted_nonzero, weighted_signed, weighted_l1, weighted_energy)
            == (1092, 0, 1248, 1560), "carry multiplicity/energy ledger")

    semantic = digest_json((
        MOD, P, T.B, amplitude, parent_hashes, reconstruction,
        tuple((label, detector[label]) for label in labels),
        tuple((label, derivative[label]) for label in labels),
        tuple(sorted(X_PLUS)), tuple(sorted(X_MINUS)), full_support, carry_support,
        full_counts, carry_counts,
        full_l1, full_energy, carry_l1, carry_energy,
        weighted_nonzero, weighted_signed, weighted_l1, weighted_energy,
    ))
    require(semantic == EXPECTED_SEMANTIC_SHA256,
            ("semantic digest", semantic, EXPECTED_SEMANTIC_SHA256))

    source = Path(__file__).resolve()
    require(not any(isinstance(node, ast.Assert)
                    for node in ast.walk(ast.parse(source.read_text(encoding="utf-8")))),
            "Python assert node present")
    print("== THM-3660 LRC exceptional leakage functional and boundary ==")
    print(f"field=p:{MOD};labels:169;generic_slope:{T.B}")
    print(f"parent_sha256_lf={parent_hashes};reconstruction_sha256={reconstruction}")
    print(f"leak=lambda(x,y)=A^-1(y-bx);A={amplitude}")
    print(f"detector_support=X:8;X_plus={tuple(sorted(X_PLUS))};X_minus={tuple(sorted(X_MINUS))}")
    print("detector_values=+1 on X_plus;-1 on X_minus;0 elsewhere;reversal=odd")
    print(f"vertical_boundary=support:{len(full_support)};labels:{full_support};counts:{full_counts}")
    print(f"carry_visible_boundary=support:{len(carry_support)};labels:{carry_support};counts:{carry_counts}")
    print("boundary_reversal=even;column_telescoping=13/13")
    print(f"boundary_norms=(L1,E2):full:{(full_l1, full_energy)},carry:{(carry_l1, carry_energy)}")
    print(f"carry_pair_ledger=(nonzero,signed,L1,E2):{(weighted_nonzero, weighted_signed, weighted_l1, weighted_energy)}")
    print(f"detector_sha256={digest_json(tuple((label, detector[label]) for label in labels))}")
    print(f"boundary_sha256={digest_json(tuple((label, derivative[label]) for label in labels))}")
    print(f"semantic_sha256={semantic}")
    print(f"script_sha256_lf={lf_sha256(source)}")
    print("scope=static finite-field quotient detector/boundary;not current/chronology/positivity/LRC14")
    print("PASS")


if __name__ == "__main__":
    main()
