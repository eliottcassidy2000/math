#!/usr/bin/env python3
"""Finite-exact cubic target-projection no-go for the THM-3556 packet.

For the explicit packet Z=(L,T,U_*,S), take all 34 nonconstant target
monomials of total target degree at most three and all 561 pairwise
source Jacobians.  Their coefficient matrix represents the linear map from
an arbitrary bivector on that 34-dimensional nominal target space to a
source polynomial in Q[v,y].  Deleting the constant source-monomial row does
not lower its exact rank.  Consequently, any bivector that kills all
nonconstant source coefficients kills the constant coefficient as well.
This excludes a nonzero constant Jacobian even before imposing Pluecker
decomposability or polynomial-potential realizability.

We calculate with the integral diagonal rescaling

    (L,T,W,R)=(L,T,2U_*,2S).

Every target monomial, and hence every bracket column, is merely multiplied
by a nonzero power of two.  Thus this preserves the degree filtration,
ranks over Q and odd prime fields, and constant-polynomial reachability.

The script uses sparse integer polynomial arithmetic and python-flint's exact
integer-matrix rank.  All truth gates remain active under ``python -O``.
"""

from __future__ import annotations

from itertools import combinations, combinations_with_replacement

from flint import fmpz_mat, nmod_mat


Poly = dict[tuple[int, int], int]
PRIMES = (1000003, 1000033, 1000037, 2147483647)


def require(condition: bool, label: str) -> None:
    """Keep truth-bearing checks active under optimized execution."""
    if not bool(condition):
        raise ArithmeticError(f"FAILED: {label}")


def is_prime(integer: int) -> bool:
    """Deterministic trial-division check for the four named controls."""
    if integer < 2:
        return False
    if integer % 2 == 0:
        return integer == 2
    divisor = 3
    while divisor * divisor <= integer:
        if integer % divisor == 0:
            return False
        divisor += 2
    return True


def clean(poly: Poly) -> Poly:
    """Delete zero coefficients from a sparse bivariate polynomial."""
    return {monomial: coefficient
            for monomial, coefficient in poly.items() if coefficient}


def add(first: Poly, second: Poly, scale: int = 1) -> Poly:
    """Return first + scale*second."""
    result = dict(first)
    for monomial, coefficient in second.items():
        result[monomial] = result.get(monomial, 0) + scale * coefficient
    return clean(result)


def scalar_mul(scalar: int, poly: Poly) -> Poly:
    """Multiply a polynomial by an integer scalar."""
    return clean({monomial: scalar * coefficient
                  for monomial, coefficient in poly.items()})


def mul(first: Poly, second: Poly) -> Poly:
    """Multiply two sparse bivariate integer polynomials."""
    result: Poly = {}
    for (first_v, first_y), first_coefficient in first.items():
        for (second_v, second_y), second_coefficient in second.items():
            monomial = (first_v + second_v, first_y + second_y)
            result[monomial] = (
                result.get(monomial, 0)
                + first_coefficient * second_coefficient
            )
    return clean(result)


def derivative(poly: Poly, variable: int) -> Poly:
    """Differentiate in variable 0=v or 1=y."""
    result: Poly = {}
    for exponent, coefficient in poly.items():
        power = exponent[variable]
        if power:
            reduced = list(exponent)
            reduced[variable] -= 1
            result[tuple(reduced)] = coefficient * power
    return result


def product(polys: list[Poly]) -> Poly:
    """Multiply a finite list of sparse polynomials."""
    result: Poly = {(0, 0): 1}
    for poly in polys:
        result = mul(result, poly)
    return result


def jacobian(
    first_v: Poly,
    first_y: Poly,
    second_v: Poly,
    second_y: Poly,
) -> Poly:
    """Return first_v*second_y-first_y*second_v."""
    return add(
        mul(first_v, second_y),
        mul(first_y, second_v),
        scale=-1,
    )


def coefficient_rows(
    polys: list[Poly],
    force_constant: bool = False,
) -> tuple[list[tuple[int, int]], list[list[int]]]:
    """Build coefficient rows in deterministic graded-lex order."""
    support: set[tuple[int, int]] = set()
    for poly in polys:
        support.update(poly)
    if force_constant:
        support.add((0, 0))
    monomials = sorted(support, key=lambda monomial: (sum(monomial), monomial))
    rows = [
        [poly.get(monomial, 0) for poly in polys]
        for monomial in monomials
    ]
    return monomials, rows


def bracket_block(
    derivatives_v: list[Poly],
    derivatives_y: list[Poly],
    size: int,
) -> tuple[list[Poly], list[tuple[int, int]]]:
    """Return all pairwise brackets among the first ``size`` functions."""
    pairs = list(combinations(range(size), 2))
    brackets = [
        jacobian(
            derivatives_v[first],
            derivatives_y[first],
            derivatives_v[second],
            derivatives_y[second],
        )
        for first, second in pairs
    ]
    return brackets, pairs


def rank_pair(
    rows: list[list[int]],
    constant_index: int,
) -> tuple[int, int, list[list[int]]]:
    """Return exact ranks before/after deleting the constant row."""
    nonconstant_rows = [
        row for index, row in enumerate(rows) if index != constant_index
    ]
    return (
        fmpz_mat(rows).rank(),
        fmpz_mat(nonconstant_rows).rank(),
        nonconstant_rows,
    )


def modular_rank_pairs(
    rows: list[list[int]],
    nonconstant_rows: list[list[int]],
) -> list[tuple[int, int, int]]:
    """Return full/deleted ranks over each named prime field."""
    return [
        (
            prime,
            nmod_mat(rows, prime).rank(),
            nmod_mat(nonconstant_rows, prime).rank(),
        )
        for prime in PRIMES
    ]


def main() -> None:
    require(all(is_prime(prime) for prime in PRIMES), "named moduli are prime")

    v: Poly = {(1, 0): 1}
    y: Poly = {(0, 1): 1}
    y_squared = mul(y, y)
    y_cubed = mul(y_squared, y)
    v_squared = mul(v, v)

    # W=2U_*=2+2y-y^2-3vy^2+9vy.
    W: Poly = {
        (0, 0): 2,
        (0, 1): 2,
        (0, 2): -1,
        (1, 1): 9,
        (1, 2): -3,
    }
    # T=y^2-3vW, R=2S=2y^3-9vWy, L=v^2(4vW-y^2).
    T = add(y_squared, mul(v, W), scale=-3)
    R = add(scalar_mul(2, y_cubed), mul(mul(v, W), y), scale=-9)
    L = mul(v_squared, add(scalar_mul(4, mul(v, W)), y_squared, scale=-1))
    packet = [L, T, W, R]
    names = ["L", "T", "W=2U", "R=2S"]

    # Rescaled cusp-square identity.  This supplies one explicit relation.
    relation = add(mul(R, R), mul(mul(T, T), T), scale=-4)
    relation = add(relation, mul(L, mul(W, W)), scale=-27)
    require(relation == {}, "R^2-4T^3-27LW^2 identity")

    # Exactly the 4+10+20=34 nonconstant monomials of total target degree <=3.
    basis: list[Poly] = []
    labels: list[str] = []
    target_multiindices: list[tuple[int, int, int, int]] = []
    for degree in range(1, 4):
        for factors in combinations_with_replacement(range(4), degree):
            multiindex = tuple(factors.count(index) for index in range(4))
            target_multiindices.append(multiindex)
            basis.append(product([packet[index] for index in factors]))
            labels.append("*".join(names[index] for index in factors))
    require(len(basis) == 34, "34 nominal cubic target monomials")
    require(len(set(target_multiindices)) == 34, "target monomials are distinct")

    composition_monomials, composition_rows = coefficient_rows(basis)
    composition_matrix = fmpz_mat(composition_rows)
    composition_rank = composition_matrix.rank()
    require(composition_matrix.nrows() == 111, "pullback coefficient row count")
    require(composition_matrix.ncols() == 34, "pullback coefficient column count")
    require(composition_rank == 33, "unique target pullback relation")

    derivatives_v = [derivative(poly, 0) for poly in basis]
    derivatives_y = [derivative(poly, 1) for poly in basis]

    # Independent inherited control: the degree-<=2 prefix must reproduce the
    # quadratic agent's exact 139x91 matrix and ranks 67/67.
    quadratic_brackets, quadratic_pairs = bracket_block(
        derivatives_v, derivatives_y, 14
    )
    require(len(quadratic_pairs) == 91, "quadratic wedge column count")
    quadratic_monomials, quadratic_rows = coefficient_rows(
        quadratic_brackets, force_constant=True
    )
    quadratic_constant_index = quadratic_monomials.index((0, 0))
    (
        quadratic_full_rank,
        quadratic_deleted_rank,
        quadratic_nonconstant_rows,
    ) = rank_pair(quadratic_rows, quadratic_constant_index)
    require(len(quadratic_rows) == 139, "quadratic coefficient row count")
    require(len(quadratic_rows[0]) == 91, "quadratic coefficient column count")
    require(
        (quadratic_full_rank, quadratic_deleted_rank) == (67, 67),
        "quadratic exact ranks",
    )
    quadratic_modular = modular_rank_pairs(
        quadratic_rows, quadratic_nonconstant_rows
    )
    require(
        all((full, deleted) == (67, 67)
            for _, full, deleted in quadratic_modular),
        "quadratic modular rank controls",
    )

    # Cubic universe requested in the theorem scout: retain all 34 nominal
    # functions and all C(34,2)=561 columns despite the unique relation.
    brackets, pairs = bracket_block(derivatives_v, derivatives_y, 34)
    require(len(pairs) == 561, "cubic wedge column count")
    bracket_monomials, full_rows = coefficient_rows(
        brackets, force_constant=True
    )
    constant_index = bracket_monomials.index((0, 0))
    exact_full_rank, exact_deleted_rank, nonconstant_rows = rank_pair(
        full_rows, constant_index
    )
    full_matrix = fmpz_mat(full_rows)
    require(full_matrix.nrows() == 336, "cubic coefficient row count")
    require(full_matrix.ncols() == 561, "cubic coefficient column count")
    require(
        (exact_full_rank, exact_deleted_rank) == (187, 187),
        "cubic exact full/deleted ranks",
    )

    modular_results = modular_rank_pairs(full_rows, nonconstant_rows)
    require(
        all((full, deleted) == (187, 187)
            for _, full, deleted in modular_results),
        "cubic modular rank controls",
    )

    # Linear-algebra consequence: equality says the constant coefficient row
    # lies in the span of the nonconstant rows.  Hence a column coefficient
    # vector annihilated by every nonconstant row is annihilated by the
    # constant row too.  This applies to arbitrary bivectors, a relaxation of
    # the decomposable coefficient vectors produced by dA wedge dB.
    constant_excluded = exact_full_rank == exact_deleted_rank
    require(constant_excluded, "constant row is redundant")

    print("THM-3556 CUBIC TARGET-PROJECTION BIVECTOR AUDIT")
    print("base_field=Q; extension_scope=every characteristic-zero field")
    print("packet_rescaling=(L,T,W,R)=(L,T,2U_*,2S)")
    print("rescaling_effect=nonzero diagonal scaling of target functions and bracket columns")
    print("target_universe=all nonconstant monomials of total target degree<=3 in four packet coordinates")
    print(f"target_function_count={len(basis)}")
    print(f"target_pullback_coefficient_shape={composition_matrix.nrows()}x{composition_matrix.ncols()}")
    print(f"target_pullback_exact_rank={composition_rank}")
    print("unique_target_relation=R^2-4*T^3-27*L*W^2=0")
    print(f"bivector_columns=binomial(34,2)={len(brackets)}")
    print(
        "source_coefficient_rows="
        f"{len(bracket_monomials)};total_degree_range="
        f"{min(map(sum, bracket_monomials))}..{max(map(sum, bracket_monomials))}"
    )
    print(f"coefficient_matrix_shape={full_matrix.nrows()}x{full_matrix.ncols()}")
    print(
        f"constant_row_index={constant_index};"
        f"deleted_shape={len(nonconstant_rows)}x{len(nonconstant_rows[0])}"
    )
    print(f"exact_ranks_full_deleted={exact_full_rank},{exact_deleted_rank}")
    print(f"modular_ranks_prime_full_deleted={modular_results}")
    print(
        "quadratic_positive_control="
        f"shape=139x91;exact={quadratic_full_rank},{quadratic_deleted_rank};"
        f"modular={quadratic_modular}"
    )
    print(f"constant_row_in_nonconstant_row_span={constant_excluded}")
    print("VERDICT: no arbitrary bivector on the degree-<=3 target space yields a nonzero constant.")
    print("COROLLARY: no pair of target polynomials of total target degree<=3 has nonzero constant source Jacobian.")
    print("SCOPE: this is the fixed explicit packet; target degree>=4 remains open.")


if __name__ == "__main__":
    main()
