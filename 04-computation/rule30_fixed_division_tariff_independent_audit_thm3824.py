#!/usr/bin/env python3
"""Independent full-lattice audit of the normalized Rule-30 division tariff.

This audit deliberately does not pass through a predicted effective period
matrix.  For each (n,r), it constructs the column-Hermite basis B of

    L = image(T_n**r) <= Z**n,

and solves T_n B = B M for the integral action M on that actual basis.  Smith
decompositions with unimodular transforms then construct, rather than merely
count, the two coordinate sublattices

    Kc = {z : M z       is 0 modulo 2**e},
    Jc = {z : T_n B z   is 0 modulo 2**e}.

Thus L/K, L/J, and J/K are obtained from direct lattice inclusions.  The last
quotient independently checks the optional post-integrality carry filtration.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
import sys

from sympy import Matrix, ZZ
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


def require(condition: bool, label: object) -> None:
    if not condition:
        raise RuntimeError(f"requirement failed: {label}")


def matrix_t(n: int) -> Matrix:
    require(n >= 1, ("positive period", n))
    entries = [[0 for _ in range(n)] for _ in range(n)]
    for row in range(n):
        entries[row][(2 * row) % n] += 1
        entries[row][(2 * row + 1) % n] += 1
    return Matrix(entries)


def valuation_two(n: int) -> int:
    require(n >= 1, ("valuation input", n))
    value = 0
    while n % 2 == 0:
        value += 1
        n //= 2
    return value


def integral_matrix(matrix: Matrix, label: object) -> Matrix:
    require(all(value.q == 1 for value in matrix), (label, matrix))
    return matrix.applyfunc(int)


def column_lattice_basis(generators: Matrix) -> Matrix:
    """Return the column-Hermite basis of the generated integer lattice."""
    basis = hermite_normal_form(generators)
    require(basis.rows == generators.rows, "HNF ambient rank")
    require(basis.cols == generators.rank(), "HNF lattice rank")
    return basis


def coordinate_action(operator: Matrix, basis: Matrix) -> Matrix:
    """Solve operator*basis = basis*M using independent pivot rows."""
    dimension = basis.cols
    pivot_rows = basis.T.rref()[1]
    require(len(pivot_rows) == dimension, "basis row pivots")
    square = basis.extract(pivot_rows, range(dimension))
    target = (operator * basis).extract(pivot_rows, range(dimension))
    action = integral_matrix(square.inv() * target, "integral coordinate action")
    require(basis * action == operator * basis, "coordinate action reconstruction")
    return action


def kernel_lattice_modulus(map_matrix: Matrix, modulus: int) -> tuple[Matrix, Matrix]:
    """Construct {z in Z^d : map_matrix*z = 0 mod modulus}.

    If D=U*A*V is Smith form, writing z=V*w makes the congruences diagonal.
    The returned columns are a genuine basis of the full-rank kernel lattice.
    """
    require(modulus >= 1, ("positive modulus", modulus))
    domain_matrix = DomainMatrix.from_Matrix(map_matrix).convert_to(ZZ)
    diagonal_dm, left_dm, right_dm = smith_normal_decomp(domain_matrix)
    diagonal = diagonal_dm.to_Matrix()
    left = left_dm.to_Matrix()
    right = right_dm.to_Matrix()
    require(left * map_matrix * right == diagonal, "Smith decomposition")
    require(abs(int(left.det())) == 1, "left Smith transform unimodular")
    require(abs(int(right.det())) == 1, "right Smith transform unimodular")

    dimension = map_matrix.cols
    multipliers = []
    for index in range(dimension):
        smith_value = (
            abs(int(diagonal[index, index]))
            if index < min(map_matrix.rows, dimension)
            else 0
        )
        multipliers.append(modulus // math.gcd(modulus, smith_value))
    kernel_basis = right * Matrix.diag(*multipliers)
    require(abs(int(kernel_basis.det())) >= 1, "full-rank congruence kernel")
    require(
        all(int(value) % modulus == 0 for value in map_matrix * kernel_basis),
        "kernel generators satisfy congruence",
    )
    return kernel_basis, diagonal


def quotient_orders(subgroup_basis: Matrix) -> tuple[int, ...]:
    """Nontrivial cyclic orders of Z^d/subgroup, largest first."""
    diagonal = smith_normal_form(subgroup_basis, domain=ZZ)
    orders = [
        abs(int(diagonal[index, index]))
        for index in range(subgroup_basis.rows)
        if abs(int(diagonal[index, index])) > 1
    ]
    return tuple(sorted(orders, reverse=True))


def predicted_total(d: int, exponent: int) -> tuple[int, ...]:
    if exponent == 0:
        return ()
    modulus = 1 << exponent
    if d % 2 == 0:
        return (modulus,) * (d // 2)
    return (modulus,) * (d - 1) + (() if exponent == 1 else (modulus // 2,))


def predicted_integrality(d: int, s: int, exponent: int) -> tuple[int, ...]:
    if exponent == 0:
        return ()
    modulus = 1 << exponent
    if d % 2 == 0:
        require(s == 0, ("even effective period", d, s))
        return (modulus,) * (d // 2)
    tail_exponent = max(exponent - s - 1, 0)
    return (modulus,) * (d - 1) + (
        () if tail_exponent == 0 else (1 << tail_exponent,)
    )


def predicted_carry(d: int, s: int, exponent: int) -> tuple[int, ...]:
    if exponent == 0 or d % 2 == 0:
        return ()
    carry_exponent = min(s, exponent - 1)
    return () if carry_exponent == 0 else (1 << carry_exponent,)


def membership_in_lattice(vector: Matrix, basis: Matrix) -> bool:
    """Test membership in a possibly non-full-rank column lattice."""
    pivot_rows = basis.T.rref()[1]
    square = basis.extract(pivot_rows, range(basis.cols))
    coordinates = square.inv() * vector.extract(pivot_rows, [0])
    if not all(value.q == 1 for value in coordinates):
        return False
    return basis * coordinates == vector


def exhaustive_definition_gate(
    operator: Matrix,
    basis: Matrix,
    action: Matrix,
    k_basis: Matrix,
    j_basis: Matrix,
    modulus: int,
) -> tuple[int, int, int]:
    """Compare constructed lattices with their ambient definitions mod q."""
    dimension = basis.cols
    residue_count = modulus**dimension
    require(residue_count <= 4096, ("bounded exhaustive gate", residue_count))
    k_inverse = k_basis.inv()
    j_inverse = j_basis.inv()
    k_count = 0
    j_count = 0
    carry_witnesses = 0
    for values in itertools.product(range(modulus), repeat=dimension):
        z = Matrix(values)
        delta = basis * z
        raw = operator * delta
        actual_j = all(int(value) % modulus == 0 for value in raw)
        actual_k = False
        if actual_j:
            divided = raw.applyfunc(lambda value: int(value) // modulus)
            actual_k = membership_in_lattice(divided, basis)
        coordinate_k = all(int(value) % modulus == 0 for value in action * z)
        constructed_k = all(value.q == 1 for value in k_inverse * z)
        constructed_j = all(value.q == 1 for value in j_inverse * z)
        require(
            actual_k == coordinate_k == constructed_k,
            ("K definition", values, actual_k, coordinate_k, constructed_k),
        )
        require(actual_j == constructed_j, ("J definition", values, actual_j, constructed_j))
        k_count += int(actual_k)
        j_count += int(actual_j)
        carry_witnesses += int(actual_j and not actual_k)
    return k_count, j_count, carry_witnesses


def main() -> None:
    records = []
    lattice_cases = 0
    tariff_cases = 0
    exhaustive_cases = 0
    exhaustive_residues = 0
    carry_witness_cases = 0
    boundary = {
        "e0": 0,
        "n1": 0,
        "r0": 0,
        "pre_core": 0,
        "core_boundary": 0,
        "post_core": 0,
    }

    for n in range(1, 19):
        operator = matrix_t(n)
        a = valuation_two(n)
        for r in range(0, 9):
            image_generators = operator**r
            basis = column_lattice_basis(image_generators)
            action = coordinate_action(operator, basis)
            dimension = basis.cols
            d = n // (1 << min(a, r))
            s = max(0, r - a)
            require(dimension == d, ("direct image rank", n, r, dimension, d))

            # This containment is the only fact needed to define the raw
            # quotient action; it is checked on the actual HNF basis.
            require(basis * action == operator * basis, ("T(L) subset L", n, r))
            lattice_cases += 1

            for exponent in range(0, 7):
                modulus = 1 << exponent
                k_basis, _ = kernel_lattice_modulus(action, modulus)
                j_basis, _ = kernel_lattice_modulus(operator * basis, modulus)

                # K <= J.  The coordinate matrix below presents J/K directly.
                inclusion = integral_matrix(j_basis.inv() * k_basis, "K inside J")
                require(j_basis * inclusion == k_basis, "carry inclusion reconstruction")

                total = quotient_orders(k_basis)
                integrality = quotient_orders(j_basis)
                carry = quotient_orders(inclusion)
                want_total = tuple(sorted(predicted_total(d, exponent), reverse=True))
                want_integrality = tuple(
                    sorted(predicted_integrality(d, s, exponent), reverse=True)
                )
                want_carry = tuple(sorted(predicted_carry(d, s, exponent), reverse=True))
                require(total == want_total, ("total tariff", n, r, exponent, total, want_total))
                require(
                    integrality == want_integrality,
                    ("integrality tariff", n, r, exponent, integrality, want_integrality),
                )
                require(carry == want_carry, ("carry tariff", n, r, exponent, carry, want_carry))
                require(
                    abs(int(k_basis.det()))
                    == abs(int(j_basis.det())) * abs(int(inclusion.det())),
                    "filtration index multiplication",
                )

                if exponent == 0:
                    boundary["e0"] += 1
                    require(total == integrality == carry == (), ("e=0", n, r))
                if n == 1:
                    boundary["n1"] += 1
                if r == 0:
                    boundary["r0"] += 1
                if r < a:
                    boundary["pre_core"] += 1
                    require(carry == (), ("pre-core has no carry", n, r, exponent))
                elif r == a:
                    boundary["core_boundary"] += 1
                    require(carry == (), ("core boundary has no carry", n, r, exponent))
                else:
                    boundary["post_core"] += 1

                if modulus**dimension <= 4096:
                    k_count, j_count, witnesses = exhaustive_definition_gate(
                        operator, basis, action, k_basis, j_basis, modulus
                    )
                    require(
                        modulus**dimension // k_count == math.prod(total),
                        ("exhaustive K index", n, r, exponent),
                    )
                    require(
                        modulus**dimension // j_count == math.prod(integrality),
                        ("exhaustive J index", n, r, exponent),
                    )
                    require(
                        j_count // k_count == math.prod(carry),
                        ("exhaustive carry index", n, r, exponent),
                    )
                    require(
                        (witnesses > 0) == (math.prod(carry) > 1),
                        ("carry witness equivalence", n, r, exponent, witnesses, carry),
                    )
                    exhaustive_cases += 1
                    exhaustive_residues += modulus**dimension
                    carry_witness_cases += int(witnesses > 0)

                records.append(
                    {
                        "n": n,
                        "r": r,
                        "e": exponent,
                        "rank": dimension,
                        "d": d,
                        "s": s,
                        "total": total,
                        "integrality": integrality,
                        "carry": carry,
                    }
                )
                tariff_cases += 1

    # Explicit exact-gap hostile: K is maximal for the fixed division N_e,
    # but need not preserve the stratum on which e is the exact valuation.
    operator = matrix_t(2)
    old_lattice_basis = Matrix.eye(2)  # r=0
    delta = Matrix([2, 0])
    x = Matrix([1, 1])
    x_prime = x + delta
    require(all(int(value) % 2 == 0 for value in operator * delta), "delta in K_2,0,1")
    require(
        min(valuation_two(abs(int(value))) for value in operator * x) == 1,
        "source exact exponent one",
    )
    require(
        min(valuation_two(abs(int(value))) for value in operator * x_prime) == 2,
        "translated source exact exponent two",
    )
    require(membership_in_lattice(delta, old_lattice_basis), "hostile delta in L")

    semantic = hashlib.sha256(
        "\n".join(json.dumps(record, sort_keys=True) for record in records).encode("ascii")
    ).hexdigest()
    print("RULE30_FIXED_DIVISION_TARIFF_INDEPENDENT_AUDIT_THM3824")
    print("status=PASS;formula_proved_subject_to_fixed-exponent_scope")
    print("method=full_Tn_power_column_HNF+coordinate_action+Smith_unimodular_kernel_lattices")
    print("definitions=K=L_intersect_Tn^(-1)(2^eL);J=L_intersect_Tn^(-1)(2^eZ^n)")
    print("e=0=all_three_quotients_trivial")
    print("e>=1,r<v2(n):L/K=(Z/2^e)^(d/2);J=K")
    print("e>=1,r>=v2(n):L/K=(Z/2^e)^(d-1)+Z/2^(e-1)")
    print("odd_integrality=L/J=(Z/2^e)^(d-1)+Z/2^max(e-s-1,0)")
    print("odd_carry=J/K=Z/2^min(s,e-1);Z/1_omitted")
    print(f"direct_image_lattice_cases={lattice_cases}")
    print(f"direct_tariff_filtration_cases={tariff_cases}")
    print(f"exhaustive_definition_cases={exhaustive_cases}")
    print(f"exhaustive_definition_residues={exhaustive_residues}")
    print(f"exhaustive_nontrivial_carry_cases={carry_witness_cases}")
    print(f"boundary_counts={json.dumps(boundary, sort_keys=True)}")
    print("universe=1<=n<=18;0<=r<=8;0<=e<=6;exhaustive_when_2^(e*rank)<=4096")
    print("scope_hostile=n=2,r=0,e=1;x=(1,1);delta=(2,0)_in_K;x+delta=(3,1)")
    print("scope_hostile_exact_gaps=1_then_2;K_does_not_preserve_exact-e_strata")
    print("scope_repair=K_is_maximal_for_fixed_N_e_on_the_2^e-divisibility_domain_only")
    print(f"semantic_sha256={semantic}")
    print("PASS")


if __name__ == "__main__":
    main()
