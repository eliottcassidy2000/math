#!/usr/bin/env python3
"""Exact split/cyclic character bridge for the two mod-13 LRC digits."""

from __future__ import annotations

import ast
from hashlib import sha256
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PARENT_SCRIPT = ROOT / "04-computation/lrc_two_current_quotient_address_reversal_gate_thm3657.py"
PARENT_OUTPUT = ROOT / "05-knowledge/results/lrc_two_current_quotient_address_reversal_gate_thm3657.out"
EXPECTED_PARENT_HASHES = (
    "f0323550a039bd3c59bc3367a9a48503ff10db2b642e99e21d9499e984492ccd",
    "fdebbbba161149edec8c86b9f382ede76153e85aa01da635f0900075fe751628",
)
EXPECTED_SEMANTIC_SHA256 = "557033432ed61af275e58c73d3015346594f5d0dca67c1df3f2643e963e7ab3e"

P = 13
N = P * P
MOD = 755373809845391722745761
OMEGA = 298763986285447441216949
ZETA = 123453826432109539554819


def require(condition: bool, payload: object) -> None:
    if condition is not True:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_json(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()


def matrix_inverse(matrix):
    size = len(matrix)
    rows = [list(row) + [int(i == j) for j in range(size)]
            for i, row in enumerate(matrix)]
    determinant = 1
    for column in range(size):
        pivot = next((row for row in range(column, size)
                      if rows[row][column] % MOD), None)
        require(pivot is not None, ("singular block", column))
        if pivot != column:
            rows[column], rows[pivot] = rows[pivot], rows[column]
            determinant = -determinant % MOD
        pivot_value = rows[column][column] % MOD
        determinant = determinant * pivot_value % MOD
        inverse = pow(pivot_value, -1, MOD)
        rows[column] = [value * inverse % MOD for value in rows[column]]
        for row in range(size):
            if row == column:
                continue
            factor = rows[row][column] % MOD
            rows[row] = [(left - factor * right) % MOD
                         for left, right in zip(rows[row], rows[column])]
    require(all(rows[i][j] == int(i == j)
                for i in range(size) for j in range(size)),
            "left inverse reduction")
    return tuple(tuple(row[size:]) for row in rows), determinant


def multiply(left, right):
    return tuple(tuple(
        sum(left[row][middle] * right[middle][column]
            for middle in range(len(right))) % MOD
        for column in range(len(right[0])))
        for row in range(len(left)))


def coefficient(u: int, k: int) -> int:
    inverse_p = pow(P, -1, MOD)
    return inverse_p * sum(
        pow(OMEGA, k * r0, MOD) * pow(ZETA, -u * r0, MOD)
        for r0 in range(P)
    ) % MOD


def split_character(u: int, v: int, r0: int, r1: int) -> int:
    return pow(ZETA, u * r0 + v * r1, MOD)


def cyclic_character(k: int, r0: int, r1: int) -> int:
    return pow(OMEGA, k * (r0 + P * r1), MOD)


def main() -> None:
    parent_hashes = (lf_sha256(PARENT_SCRIPT), lf_sha256(PARENT_OUTPUT))
    require(parent_hashes == EXPECTED_PARENT_HASHES,
            ("THM-3657 parent hashes", parent_hashes))
    require(pow(OMEGA, N, MOD) == 1 and pow(OMEGA, P, MOD) != 1,
            "omega order is not 169")
    require(ZETA == pow(OMEGA, P, MOD), "zeta definition")
    require(pow(ZETA, P, MOD) == 1 and ZETA != 1, "zeta order is not 13")

    blocks = []
    inverses = []
    determinants = []
    nonzero_counts = []
    expansion_checks = 0
    inverse_expansion_checks = 0
    reversal_checks = 0
    identity = tuple(tuple(int(row == column) for column in range(P))
                     for row in range(P))
    for v in range(P):
        # Column h is psi_(v+13h), expanded in rows phi_(u,v).
        block = tuple(tuple(coefficient(u, v + P * h) for h in range(P))
                      for u in range(P))
        inverse, determinant = matrix_inverse(block)
        require(multiply(block, inverse) == identity
                and multiply(inverse, block) == identity,
                ("block inverse", v))
        blocks.append(block)
        inverses.append(inverse)
        determinants.append(determinant)
        nonzero_counts.append(sum(value != 0 for row in block for value in row))
        require(nonzero_counts[-1] == (P if v == 0 else P * P),
                ("block density", v, nonzero_counts[-1]))
        require(sum(value != 0 for row in inverse for value in row)
                == (P if v == 0 else P * P),
                ("inverse block density", v))

        for h in range(P):
            k = v + P * h
            for r0 in range(P):
                for r1 in range(P):
                    rebuilt = sum(
                        block[u][h] * split_character(u, v, r0, r1)
                        for u in range(P)
                    ) % MOD
                    require(rebuilt == cyclic_character(k, r0, r1),
                            ("cyclic expansion", v, h, r0, r1))
                    expansion_checks += 1
        for u in range(P):
            for r0 in range(P):
                for r1 in range(P):
                    rebuilt = sum(
                        inverse[h][u] * cyclic_character(v + P * h, r0, r1)
                        for h in range(P)
                    ) % MOD
                    require(rebuilt == split_character(u, v, r0, r1),
                            ("split expansion", v, u, r0, r1))
                    inverse_expansion_checks += 1

    blocks = tuple(blocks)
    inverses = tuple(inverses)
    require(nonzero_counts == [P] + [P * P] * (P - 1),
            "global block-density profile")

    # Reversal Rf(a)=f(-a-1) has the stated action in both bases.  Compare
    # coefficient vectors, not merely evaluations at sample addresses.
    for k in range(N):
        v = k % P
        target_k = (-k) % N
        target_v = target_k % P
        h_target = (target_k - target_v) // P
        transformed = [0] * P
        for u in range(P):
            target_u = (-u) % P
            transformed[target_u] = (
                blocks[v][u][(k - v) // P]
                * pow(ZETA, -(u + v), MOD)
            ) % MOD
        expected = [
            pow(OMEGA, -k, MOD) * blocks[target_v][u][h_target] % MOD
            for u in range(P)
        ]
        require(transformed == expected, ("reversal intertwining", k))
        reversal_checks += 1

    # Minimal carry hostile: identity on the 169 labels is not additive.
    split_left = ((12 + 1) % P, (0 + 0) % P)
    cyclic_sum = (12 + 1) % N
    cyclic_right = (cyclic_sum % P, cyclic_sum // P)
    require(split_left == (0, 0) and cyclic_right == (0, 1),
            "carry hostile drift")
    require(split_left != cyclic_right, "split/cyclic addition unexpectedly equal")

    block_digest = digest_json(blocks)
    inverse_digest = digest_json(inverses)
    semantic = digest_json((
        MOD, P, N, OMEGA, ZETA, parent_hashes,
        tuple(nonzero_counts), tuple(determinants),
        block_digest, inverse_digest,
        expansion_checks, inverse_expansion_checks, reversal_checks,
        split_left, cyclic_right,
    ))
    require(semantic == EXPECTED_SEMANTIC_SHA256,
            ("semantic digest", semantic, EXPECTED_SEMANTIC_SHA256))

    source = Path(__file__).resolve()
    require(not any(isinstance(node, ast.Assert)
                    for node in ast.walk(ast.parse(source.read_text(encoding="utf-8")))),
            "Python assert node present")
    print("== THM-3658 LRC mod-169 carry-Fourier block intertwiner ==")
    print(f"field=p:{MOD};digits={P}x{P};omega:{OMEGA};zeta:{ZETA}")
    print(f"parent_sha256_lf={parent_hashes}")
    print("bases=split:phi_(u,v)=zeta^(u*r0+v*r1);cyclic:psi_k=omega^(k*(r0+13*r1))")
    print(f"block_nonzero_counts={tuple(nonzero_counts)};block_ranks={(P,) * P}")
    print(f"block_determinants={tuple(determinants)}")
    print(f"block_sha256={block_digest};inverse_sha256={inverse_digest}")
    print(f"expansion_checks=(cyclic:{expansion_checks},split:{inverse_expansion_checks})")
    print(f"reversal=split:(u,v)->(-u,-v),phase:zeta^(-u-v);cyclic:k->-k,phase:omega^(-k);checks:{reversal_checks}")
    print("carry_hostile=(12,0)+(1,0):split=(0,0),cyclic_digits=(0,1)")
    print("scope=exact linear character-basis bridge;not label-group or convolution intertwiner")
    print(f"semantic_sha256={semantic}")
    print(f"script_sha256_lf={lf_sha256(source)}")
    print("PASS")


if __name__ == "__main__":
    main()
