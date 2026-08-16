#!/usr/bin/env python3
"""Independent finite referee for the minimum U_full bridge coordinate.

This script does not import the primary probe.  It independently computes the
rank and complete kernel of the 2x2 marginal map over F_13, verifies that one
joint parity functional raises the rank from three to four, exhausts all
left-only and right-only parity functionals, and recomputes the two 13^3
character banks in the frozen split-prime chart.  It also audits the primary
record schema and confirms that no K4 computation occurs.
"""

from __future__ import annotations

import ast
from hashlib import sha256
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "04-computation/lrc_endpoint_ufull_minimal_joint_address_referee_20260816.py"
OUTPUT = "05-knowledge/results/lrc_endpoint_ufull_minimal_joint_address_referee_20260816.out"
PRIMARY = ROOT / "04-computation/lrc_endpoint_ufull_minimal_joint_address_gate_20260816.py"
PRIMARY_SHA256 = "fd5fcfc5f92385806f5ea5e77e854c823d2b3d6e34752c7aba777744201c8055"
PRIMARY_OUTPUT = ROOT / "05-knowledge/results/lrc_endpoint_ufull_minimal_joint_address_gate_20260816.out"
PRIMARY_OUTPUT_SHA256 = "3c68f4b3618abbcc2653eb1b2aee4729e59a57820749396ca757341ee7f91838"

P = 13
PRIME = 572252886246508880869
ROOT_OF_ORDER_N = 279936
N = 81750412320929840124
Q_H = (1, 0, 1)
Q_Q5 = (1, 0, 0)
BRIDGE = 389266878372286537904
EXPECTED_BANK_HASHES = (
    "b5ec7b38ed4c44abffa1e760b053f4275281017d2f7667fb44532f02c4969f7f",
    "d9c2b1e9d5cba6a635c4c75dccd0d0684b17ec3ee05718fabc808a902a95a86c",
)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_json(value: object) -> str:
    return sha256(
        json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    ).hexdigest()


def rank_mod(matrix: tuple[tuple[int, ...], ...], prime: int) -> int:
    rows = [list(value % prime for value in row) for row in matrix]
    if not rows:
        return 0
    column_count = len(rows[0])
    rank = 0
    for column in range(column_count):
        pivot = next(
            (index for index in range(rank, len(rows)) if rows[index][column]),
            None,
        )
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = pow(rows[rank][column], -1, prime)
        rows[rank] = [value * inverse % prime for value in rows[rank]]
        for index in range(len(rows)):
            if index == rank:
                continue
            multiplier = rows[index][column]
            if multiplier:
                rows[index] = [
                    (value - multiplier * pivot_value) % prime
                    for value, pivot_value in zip(rows[index], rows[rank])
                ]
        rank += 1
        if rank == len(rows):
            break
    return rank


def mat_vec(
    matrix: tuple[tuple[int, ...], ...], vector: tuple[int, ...], prime: int
) -> tuple[int, ...]:
    return tuple(
        sum(left * right for left, right in zip(row, vector)) % prime
        for row in matrix
    )


def dot(left: tuple[int, ...], right: tuple[int, ...], prime: int) -> int:
    return sum(a * b for a, b in zip(left, right)) % prime


def character_bank(
    matrix: tuple[tuple[int, int], tuple[int, int]], zeta: int
) -> tuple[int, ...]:
    rows = []
    for alpha in range(P):
        for beta in range(P):
            for tau in range(P):
                total = 0
                for i in range(2):
                    for j in range(2):
                        address = Q_H if (i ^ j) == 0 else Q_Q5
                        exponent = (
                            alpha * address[0]
                            + beta * address[1]
                            + tau * address[2]
                        ) % P
                        total = (
                            total + matrix[i][j] * pow(zeta, exponent, PRIME)
                        ) % PRIME
                rows.append(total)
    return tuple(rows)


def inverse_value(bank: tuple[int, ...], address: tuple[int, int, int], zeta: int) -> int:
    total = 0
    index = 0
    for alpha in range(P):
        for beta in range(P):
            for tau in range(P):
                exponent = -(
                    alpha * address[0]
                    + beta * address[1]
                    + tau * address[2]
                ) % P
                total = (total + bank[index] * pow(zeta, exponent, PRIME)) % PRIME
                index += 1
    return total * pow(P**3, -1, PRIME) % PRIME


def orthogonality(
    source: tuple[int, int, int], target: tuple[int, int, int], zeta: int
) -> int:
    total = 0
    for alpha in range(P):
        for beta in range(P):
            for tau in range(P):
                exponent = (
                    alpha * (source[0] - target[0])
                    + beta * (source[1] - target[1])
                    + tau * (source[2] - target[2])
                ) % P
                total = (total + pow(zeta, exponent, PRIME)) % PRIME
    return total


def primary_ast_audit() -> tuple[object, ...]:
    tree = ast.parse(PRIMARY.read_text(encoding="utf-8"))
    classes = {
        node.name: node for node in tree.body if isinstance(node, ast.ClassDef)
    }
    require("CompleteAtomRecord" in classes, "complete schema absent")
    schema = tuple(
        node.target.id
        for node in classes["CompleteAtomRecord"].body
        if isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name)
    )
    required = (
        "base",
        "root",
        "owner_sheet",
        "word_sheet",
        "source_sheet",
        "left_horizon",
        "right_horizon",
        "address",
    )
    require(all(name in schema for name in required), ("schema", schema))
    functions = {
        node.name: node for node in tree.body if isinstance(node, ast.FunctionDef)
    }
    address = functions["address_of"]
    args = tuple(argument.arg for argument in address.args.args)
    names = tuple(
        sorted(node.id for node in ast.walk(address) if isinstance(node, ast.Name))
    )
    require(args == ("record",) and "ell" not in names, "ell-dependent address")
    k4_calls = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        if isinstance(node.func, ast.Name):
            name = node.func.id
        elif isinstance(node.func, ast.Attribute):
            name = node.func.attr
        else:
            name = ""
        if "k4" in name.lower():
            k4_calls.append(name)
    require(not k4_calls, ("K4 call", k4_calls))
    return schema, args, names, tuple(k4_calls)


def security_certificate(path: Path) -> tuple[int, tuple[str, ...]]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    bad = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Assert):
            bad.append("Assert")
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
            if node.func.id in {"eval", "exec", "compile", "__import__"}:
                bad.append(node.func.id)
    require(not bad, ("security", bad))
    return len(tuple(ast.walk(tree))), tuple(bad)


def main() -> None:
    require(lf_sha256(PRIMARY) == PRIMARY_SHA256, "primary source hash drift")
    require(
        lf_sha256(PRIMARY_OUTPUT) == PRIMARY_OUTPUT_SHA256,
        "primary output hash drift",
    )
    output_text = PRIMARY_OUTPUT.read_text(encoding="utf-8")
    require("bridge_only_test=(positive=389266878372286537904,hostile=0" in output_text,
            "primary bridge line absent")
    require("no_K4_evaluated=True" in output_text, "primary K4 boundary absent")

    # Coordinates are x00,x01,x10,x11.  The four displayed rows are the two
    # row margins and two column margins.
    margin_map = (
        (1, 1, 0, 0),
        (0, 0, 1, 1),
        (1, 0, 1, 0),
        (0, 1, 0, 1),
    )
    checker = (1, -1, -1, 1)
    rank_margin = rank_mod(margin_map, P)
    kernel = tuple(
        vector
        for a in range(P)
        for b in range(P)
        for c in range(P)
        for d in range(P)
        if mat_vec(margin_map, (a, b, c, d), P) == (0, 0, 0, 0)
        for vector in ((a, b, c, d),)
    )
    expected_kernel = tuple(
        tuple(scale * value % P for value in checker) for scale in range(P)
    )
    require(rank_margin == 3, ("margin rank", rank_margin))
    require(set(kernel) == set(expected_kernel), ("kernel", kernel))
    pairing = dot(checker, checker, P)
    require(pairing == 4, ("pairing", pairing))
    rank_augmented = rank_mod(margin_map + (checker,), P)
    require(rank_augmented == 4, ("augmented rank", rank_augmented))

    left_only_values = {
        dot((a, a, b, b), checker, P)
        for a in range(P) for b in range(P)
    }
    right_only_values = {
        dot((a, b, a, b), checker, P)
        for a in range(P) for b in range(P)
    }
    require(left_only_values == right_only_values == {0}, "one-factor leak")

    require(BRIDGE % 4 == 0, "bridge quarter")
    lam = BRIDGE // 4
    positive = ((2 * lam, 0), (0, 2 * lam))
    hostile = ((lam, lam), (lam, lam))
    margins = lambda matrix: (
        (sum(matrix[0]), sum(matrix[1])),
        (matrix[0][0] + matrix[1][0], matrix[0][1] + matrix[1][1]),
    )
    require(margins(positive) == margins(hostile), "marginal hostile drift")

    zeta = pow(ROOT_OF_ORDER_N, N // P, PRIME)
    require(pow(zeta, P, PRIME) == 1 and zeta != 1, "order-13 root")
    orthogonality_table = tuple(
        (source, target, orthogonality(source, target, zeta))
        for source in (Q_H, Q_Q5)
        for target in (Q_H, Q_Q5)
    )
    require(
        all(value == (P**3 if source == target else 0)
            for source, target, value in orthogonality_table),
        ("orthogonality", orthogonality_table),
    )
    positive_bank = character_bank(positive, zeta)
    hostile_bank = character_bank(hostile, zeta)
    bank_hashes = (digest_json(positive_bank), digest_json(hostile_bank))
    require(bank_hashes == EXPECTED_BANK_HASHES, ("bank hashes", bank_hashes))
    values = tuple(
        tuple(inverse_value(bank, address, zeta) for address in (Q_H, Q_Q5))
        for bank in (positive_bank, hostile_bank)
    )
    require(values == ((BRIDGE, 0), (BRIDGE // 2, BRIDGE // 2)), values)
    bridges = tuple((row[0] - row[1]) % PRIME for row in values)
    require(bridges == (BRIDGE, 0), bridges)

    ast_audit = primary_ast_audit()
    security = security_certificate(ROOT / SCRIPT)
    semantic = (
        rank_margin,
        len(kernel),
        pairing,
        rank_augmented,
        tuple(sorted(left_only_values)),
        tuple(sorted(right_only_values)),
        orthogonality_table,
        bank_hashes,
        values,
        bridges,
        ast_audit,
        "minimum mixed-Haar coordinate only; no U_full ancestry realization",
    )

    print("LRC U_FULL MINIMUM JOINT-ADDRESS INDEPENDENT REFEREE")
    print("status=FINITE-EXACT independent API/kernel audit; LRC(14) OPEN")
    print(f"primary_hashes={(PRIMARY_SHA256, PRIMARY_OUTPUT_SHA256)}")
    print(f"margin_map=(rank={rank_margin},kernel_size={len(kernel)},kernel={kernel})")
    print(f"joint_coordinate=(checker={checker},pairing={pairing},augmented_rank={rank_augmented})")
    print(f"one_factor_exhaustion=(left={left_only_values},right={right_only_values})")
    print(f"orthogonality_table={orthogonality_table}")
    print(f"bank_hashes={bank_hashes}")
    print(f"inverse_values=(positive_H_q5={values[0]},hostile_H_q5={values[1]}); bridges={bridges}")
    print(f"primary_ast_audit={ast_audit}")
    print("verdict=one joint scalar is necessary and sufficient for a 2x2 lift; the calibrated positive recovers only the frozen bridge, while identical marginals admit the zero-bridge hostile")
    print("scope=no semantic map from actual endpoint cells to the joint coordinate is supplied; no K4/current/row/LRC consequence")
    print(f"security_ast={security}")
    print(f"semantic_sha256={digest_json(semantic)}")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
