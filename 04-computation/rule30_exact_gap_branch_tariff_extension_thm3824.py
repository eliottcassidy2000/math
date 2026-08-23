#!/usr/bin/env python3
"""Hostile audit of the exact-gap refinement of THM-3824.

For L=T_n^r Z^n and a fixed exponent e, THM-3824 defines

    K_e = {delta in L : 2^(-e) T_n delta lies in L}.

K_e is the largest translation lattice for the fixed linear branch, but it
need not preserve the ambient exact-gap shell

    X_e = {x in L : T_n x is in 2^e Z^n but not 2^(e+1) Z^n}.

This scratch companion constructs the actual image lattice L by column HNF
and checks that the maximal exact-gap-preserving subgroup of K_e is

    B_e = K_e intersect {delta : T_n delta is in 2^(e+1) Z^n}.

It also checks the closed formula for the elementary two-group K_e/B_e.
No Rule 30 prize consequence is claimed.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
import sys

from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


def require(condition: bool, label: object) -> None:
    if not condition:
        raise RuntimeError(f"requirement failed: {label}")


def matrix_t(n: int) -> Matrix:
    require(n >= 1, ("positive period", n))
    rows = [[0 for _ in range(n)] for _ in range(n)]
    for j in range(n):
        rows[j][(2 * j) % n] += 1
        rows[j][(2 * j + 1) % n] += 1
    return Matrix(rows)


def valuation_two(n: int) -> int:
    require(n >= 1, ("positive valuation input", n))
    answer = 0
    while n % 2 == 0:
        answer += 1
        n //= 2
    return answer


def image_basis(matrix: Matrix) -> Matrix:
    basis = hermite_normal_form(matrix)
    require(basis.rows == matrix.rows, "ambient HNF dimension")
    require(basis.cols == matrix.rank(), "HNF rank")
    return basis


def integral_matrix(matrix: Matrix, label: object) -> Matrix:
    require(all(value.q == 1 for value in matrix), (label, matrix))
    return matrix.applyfunc(int)


def coordinate_action(operator: Matrix, basis: Matrix) -> Matrix:
    rank = basis.cols
    pivot_rows = basis.T.rref()[1]
    require(len(pivot_rows) == rank, "basis pivot rows")
    square = basis.extract(pivot_rows, range(rank))
    target = (operator * basis).extract(pivot_rows, range(rank))
    action = integral_matrix(square.inv() * target, "integral action")
    require(basis * action == operator * basis, "action reconstruction")
    return action


def gf2_rank(matrix: Matrix) -> int:
    rows = [[int(value) & 1 for value in matrix.row(i)] for i in range(matrix.rows)]
    rank = 0
    column = 0
    while rank < len(rows) and column < matrix.cols:
        pivot = next((i for i in range(rank, len(rows)) if rows[i][column]), None)
        if pivot is None:
            column += 1
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for i in range(len(rows)):
            if i != rank and rows[i][column]:
                rows[i] = [a ^ b for a, b in zip(rows[i], rows[rank])]
        rank += 1
        column += 1
    return rank


def predicted_branch_dimension(d: int, s: int, exponent: int) -> int:
    """Dimension of K_e/B_e over F_2."""
    if d % 2 == 0:
        require(s == 0, ("even core is pre-core", d, s))
        return d // 2
    if s == 0:
        return d - 1 if exponent == 0 else d
    return d - 1


def dot_rows(rows: tuple[tuple[int, ...], ...], vector: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sum(a * b for a, b in zip(row, vector)) for row in rows)


def divisible(vector: tuple[int, ...], modulus: int) -> bool:
    return all(value % modulus == 0 for value in vector)


def main() -> None:
    records: list[dict[str, object]] = []
    lattice_cases = 0
    exponent_cases = 0
    residue_gates = 0
    preservation_gates = 0
    self_hostile_gates = 0
    adaptive_kernel_rank_gates = 0
    boundary_counts = {
        "pre_core_even": 0,
        "odd_core_boundary_e0": 0,
        "odd_core_boundary_positive_e": 0,
        "odd_post_core": 0,
        "zero_tariff": 0,
    }

    for n in range(1, 13):
        operator = matrix_t(n)
        a = valuation_two(n)
        for r in range(0, 8):
            basis = image_basis(operator**r)
            action = coordinate_action(operator, basis)
            raw_map = operator * basis
            d = n // (1 << min(a, r))
            s = max(0, r - a)
            require(basis.cols == d, ("effective rank", n, r, basis.cols, d))

            # Across all exact-gap branches at once, a translation preserves
            # the adaptive valuation function iff its response is literally
            # zero.  The coordinate stabilizer is therefore ker(action).
            adaptive_kernel_rank = d - action.rank()
            wanted_adaptive_rank = d // 2 if d % 2 == 0 else 0
            require(
                adaptive_kernel_rank == wanted_adaptive_rank,
                ("adaptive stabilizer rank", n, r, d, adaptive_kernel_rank),
            )
            adaptive_kernel_rank_gates += 1

            # The conceptual formula is the rank of T_d^s applied to the
            # normalized quotient image.  These rank controls independently
            # recover every piece of the closed form.
            effective = matrix_t(d)
            if d % 2 == 0:
                require(gf2_rank(effective) == d // 2, ("even rank", d))
            else:
                require(gf2_rank(effective) == d - 1, ("odd rank", d))
                if s >= 1:
                    require(gf2_rank(effective**s) == d - 1, ("stable odd rank", d, s))
            lattice_cases += 1

            action_rows = tuple(tuple(int(value) for value in action.row(i)) for i in range(d))
            raw_rows = tuple(tuple(int(value) for value in raw_map.row(i)) for i in range(n))

            for exponent in range(0, 5):
                modulus_e = 1 << exponent
                modulus_next = 1 << (exponent + 1)
                residue_universe = modulus_next**d
                if residue_universe > 131072:
                    continue

                k_count = 0
                branch_count = 0
                for vector in itertools.product(range(modulus_next), repeat=d):
                    action_value = dot_rows(action_rows, vector)
                    in_k = divisible(action_value, modulus_e)
                    if not in_k:
                        continue
                    k_count += 1
                    raw_value = dot_rows(raw_rows, vector)
                    in_b = divisible(raw_value, modulus_next)
                    if in_b:
                        # Adding this translation changes normalized raw
                        # output by an even vector, so it preserves the exact
                        # parity shell for every base point.
                        preservation_gates += 1
                    else:
                        # This translation is itself in X_e, while twice it
                        # leaves X_e.  Hence every missing K/B class is a
                        # direct maximality hostile.
                        require(divisible(raw_value, modulus_e), ("K implies D_e", n, r, exponent, vector))
                        doubled_raw = tuple(2 * value for value in raw_value)
                        require(divisible(doubled_raw, modulus_next), ("self-double hostile", n, r, exponent, vector))
                        branch_count += 1
                        self_hostile_gates += 1
                    residue_gates += 1

                b_count = k_count - branch_count
                require(b_count > 0 and k_count % b_count == 0, ("finite quotient count", n, r, exponent))
                index = k_count // b_count
                beta = predicted_branch_dimension(d, s, exponent)
                require(index == 1 << beta, ("branch tariff", n, r, exponent, d, s, index, beta))

                # 2K_e is contained in B_e, so the quotient of the checked
                # order is the elementary group (Z/2)^beta, not merely a
                # two-group of the same cardinality.
                require(index & (index - 1) == 0, ("two-power index", index))

                if d % 2 == 0:
                    boundary_counts["pre_core_even"] += 1
                elif s == 0 and exponent == 0:
                    boundary_counts["odd_core_boundary_e0"] += 1
                elif s == 0:
                    boundary_counts["odd_core_boundary_positive_e"] += 1
                else:
                    boundary_counts["odd_post_core"] += 1
                if beta == 0:
                    boundary_counts["zero_tariff"] += 1

                records.append(
                    {
                        "n": n,
                        "r": r,
                        "e": exponent,
                        "d": d,
                        "s": s,
                        "beta": beta,
                        "index": index,
                        "residues": residue_universe,
                        "adaptive_kernel_rank": adaptive_kernel_rank,
                    }
                )
                exponent_cases += 1

    # Canonical THM-3824 hostile: the unique nonzero K/B bit changes the
    # ambient exact common valuation from one to two.
    hostile_operator = matrix_t(2)
    hostile_x = Matrix([1, 1])
    hostile_delta = Matrix([2, 0])
    hostile_y = hostile_x + hostile_delta
    require(tuple(hostile_operator * hostile_x) == (2, 2), "canonical hostile left")
    require(tuple(hostile_operator * hostile_delta) == (2, 2), "canonical hostile tariff bit")
    require(tuple(hostile_operator * hostile_y) == (4, 4), "canonical hostile right")

    # The adaptive all-shell stabilizer is exactly the response kernel.  At
    # n=2,r=0 it is generated by (1,-1); at the one-point odd core it is zero,
    # and x=1 is already a self-doubling hostile.
    require(tuple(hostile_operator * Matrix([1, -1])) == (0, 0), "adaptive even kernel")
    one_operator = matrix_t(1)
    require(one_operator.rank() == 1, "one-core adaptive stabilizer zero")
    require(tuple(one_operator * Matrix([1])) == (2,), "one-core adaptive hostile left")
    require(tuple(one_operator * Matrix([2])) == (4,), "one-core adaptive hostile right")

    # Sharp dead boundary: after a power-of-two period has collapsed to the
    # one-point odd core and one further old-image step is retained, K_e/B_e
    # is trivial for every e.  Fixed-branch congruence already preserves the
    # ambient exact gap there.
    for exponent in range(0, 9):
        require(predicted_branch_dimension(1, 1, exponent) == 0, ("one-core post-core", exponent))

    semantic = hashlib.sha256(
        "\n".join(json.dumps(record, sort_keys=True) for record in records).encode("ascii")
    ).hexdigest()

    print("RULE30_EXACT_GAP_BRANCH_TARIFF_PROBE")
    print("status=PROVED_FORMULA_DERIVATION+FINITE_EXACT_HOSTILE_AUDIT;all_rule30_prizes_open")
    print("definitions=L=T_n^r*Z^n;K_e={delta_in_L:2^-e*T_n*delta_in_L}")
    print("exact_shell=X_e={x_in_L:T_n*x_in_2^e*Z^n_not_2^(e+1)*Z^n}")
    print("maximal_branch_subgroup=B_e=K_e_intersect_{T_n*delta_in_2^(e+1)*Z^n}")
    print("mechanism=normalized_response_parity_is_the_complete_exact-gap_sidecar")
    print("quotient=K_e/B_e=(Z/2)^beta")
    print("beta_pre_core_even=d/2")
    print("beta_odd_core_boundary=e==0?d-1:d")
    print("beta_odd_post_core=d-1")
    print("parameters=d=n/2^min(v2(n),r);s=max(0,r-v2(n))")
    print("maximality=every_delta_in_K_e_minus_B_e_is_a_self_hostile:x=delta_exact_then_x+delta=2delta_not_exact")
    print("canonical_hostile=n=2,r=0,e=1,x=(1,1),delta=(2,0),raw_outputs=(2,2)_then_(4,4)")
    print("sharp_zero_boundary=d=1,s>=1:beta=0_for_every_e")
    print("adaptive_all_shells=Stab_L(x_mapsto_v2(T_n*x))=ker(T_n|L)")
    print("adaptive_stabilizer_rank=d/2_pre_core_even;0_at_and_after_odd_core")
    print("adaptive_maximality=if_T_n*delta_nonzero_then_x=delta_has_gap_e_but_x+delta_has_gap_e+1")
    print("typed_analogy=THM3825_uses_valuation_quotient_remainder_to_decode_tags;here_one_normalized_parity_vector_preserves_a_branch_label_only")
    print(f"actual_image_lattice_cases={lattice_cases}")
    print(f"exhaustive_exponent_cases={exponent_cases}")
    print(f"K_residue_gates={residue_gates}")
    print(f"preserving_translation_gates={preservation_gates}")
    print(f"self_hostile_gates={self_hostile_gates}")
    print(f"adaptive_kernel_rank_gates={adaptive_kernel_rank_gates}")
    print(f"boundary_counts={json.dumps(boundary_counts, sort_keys=True)}")
    print("finite_universe=1<=n<=12;0<=r<=7;0<=e<=4;2^((e+1)*rank)<=131072")
    print(f"semantic_sha256={semantic}")
    print("PASS")


if __name__ == "__main__":
    main()
