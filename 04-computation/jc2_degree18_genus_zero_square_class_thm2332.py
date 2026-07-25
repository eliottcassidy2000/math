#!/usr/bin/env python3
"""Exact referee for THM-2332's genus-zero square-class reduction.

The script verifies the depressed-cubic Mordell form, the uniform
irreducibility coefficient mechanism on every support of size at least
three, and the complete connected degree-three ramification-signature
ledger.  The order/maximal-order discriminant relation and
Riemann--Hurwitz implication are mathematical steps in the theorem.
Every executable check remains active under optimized Python.
"""

from __future__ import annotations

from itertools import combinations, product

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def compose(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    """Composition left after right."""
    return tuple(left[right[index]] for index in range(3))


def inverse(permutation: tuple[int, ...]) -> tuple[int, ...]:
    result = [0, 0, 0]
    for source, target in enumerate(permutation):
        result[target] = source
    return tuple(result)


def ramification_contribution(permutation: tuple[int, ...]) -> int:
    """Three minus the number of cycles."""
    seen: set[int] = set()
    cycles = 0
    for start in range(3):
        if start in seen:
            continue
        cycles += 1
        current = start
        while current not in seen:
            seen.add(current)
            current = permutation[current]
    return 3 - cycles


def transitive(generators: tuple[tuple[int, ...], ...]) -> bool:
    """Check transitivity of the subgroup generated on three sheets."""
    identity = (0, 1, 2)
    subgroup = {identity}
    frontier = [identity]
    steps = tuple(generators) + tuple(inverse(item) for item in generators)
    while frontier:
        current = frontier.pop()
        for generator in steps:
            candidate = compose(generator, current)
            if candidate not in subgroup:
                subgroup.add(candidate)
                frontier.append(candidate)
    return {permutation[0] for permutation in subgroup} == {0, 1, 2}


def main() -> None:
    u, y, v = sp.symbols("u y v")
    bvar, cvar, dvar, wvar = sp.symbols("B C D W")

    spectral = (
        -5878656 * wvar * y
        - 26040609 * u**3
        + 49601160 * bvar * u**2
        + 1607445 * u**2 * y**2
        - 20995200 * bvar**2 * u
        - 2857680 * bvar * u * y**2
        - 52907904 * dvar * u
        - 138915 * u * y**4
        + 777600 * bvar**2 * y**2
        + 33592320 * bvar * dvar
        - 5598720 * bvar * cvar * y
        + 78120 * bvar * y**4
        + 1959552 * dvar * y**2
        - 435456 * cvar * y**3
        + 1127 * y**6
    )

    polynomial_u = sp.Poly(spectral, u)
    a3 = polynomial_u.coeff_monomial(u**3)
    a2 = polynomial_u.coeff_monomial(u**2)
    a1 = polynomial_u.coeff_monomial(u)
    a0 = polynomial_u.coeff_monomial(1)
    depressed_p = sp.factor((3 * a3 * a1 - a2**2) / (3 * a3**2))
    depressed_q = sp.factor(
        (
            27 * a3**2 * a0
            - 9 * a3 * a2 * a1
            + 2 * a2**3
        )
        / (27 * a3**3)
    )

    quartic_p = (
        245 * y**4
        + 1890 * bvar * y**2
        - 24300 * bvar**2
        + 122472 * dvar
    )
    sextic_q = (
        539 * y**6
        + 11340 * bvar * y**4
        + 183708 * cvar * y**3
        + (72900 * bvar**2 - 367416 * dvar) * y**2
        + (2361960 * bvar * cvar + 2480058 * wvar) * y
    )
    require(
        sp.cancel(depressed_p - sp.Rational(16, 964467) * quartic_p)
        == 0,
        "depressed quartic coefficient changed",
    )
    require(
        sp.cancel(depressed_q - sp.Rational(64, 703096443) * sextic_q)
        == 0,
        "depressed sextic coefficient changed",
    )

    branch = sp.discriminant(spectral, u)
    mordell = 4 * quartic_p**3 + 49 * sextic_q**2
    branch_scalar = -2099434729254912
    require(
        sp.expand(branch - branch_scalar * mordell) == 0,
        "spectral discriminant is not the stated Mordell polynomial",
    )
    require(
        sp.factorint(abs(branch_scalar)) == {2: 12, 3: 21, 7: 2},
        "discriminant scalar factorization changed",
    )
    require(
        sp.Poly(mordell, y).degree() == 12
        and sp.Poly(mordell, y).LC() == 73060029,
        "Mordell polynomial lost degree twelve",
    )

    infinity_polynomial = (
        1127 - 138915 * v + 1607445 * v**2 - 26040609 * v**3
    )
    infinity_discriminant = sp.discriminant(infinity_polynomial, v)
    require(
        infinity_discriminant == -153384762202971019112448,
        "weighted infinity cubic is not separable",
    )
    require(
        sp.Poly(branch, y).LC() == infinity_discriminant,
        "finite branch polynomial leading coefficient changed",
    )

    # Uniform irreducibility on every support with C or W nonzero.
    aa, bb, cc = sp.symbols("aa bb cc")
    root_substitution = sp.Poly(
        sp.expand(spectral.subs(u, aa * y**2 + bb * y + cc)),
        y,
    )
    require(
        sp.expand(
            root_substitution.coeff_monomial(y**6)
            - infinity_polynomial.subs(v, aa)
        )
        == 0,
        "putative-root y^6 coefficient changed",
    )
    require(
        sp.expand(
            root_substitution.coeff_monomial(y**5)
            - bb * sp.diff(infinity_polynomial, v).subs(v, aa)
        )
        == 0,
        "putative-root y^5 coefficient changed",
    )
    require(
        sp.expand(
            root_substitution.coeff_monomial(y**3).subs(bb, 0)
            + 435456 * cvar
        )
        == 0,
        "putative-root y^3 coefficient changed",
    )
    require(
        sp.expand(
            root_substitution.coeff_monomial(y).subs({bb: 0, cvar: 0})
            + 5878656 * wvar
        )
        == 0,
        "putative-root y coefficient changed",
    )
    labels = ("B", "C", "D", "W")
    high_supports = [
        support
        for size in (3, 4)
        for support in combinations(labels, size)
    ]
    require(
        len(high_supports) == 5
        and all("C" in support or "W" in support for support in high_supports),
        "a high-support stratum escaped the irreducibility coordinates",
    )

    # Connected tame degree-three covers of P1 with genus zero have total
    # ramification four. Enumerate all branch cycles in S3 to retain the
    # transitivity/product-one sidecars, not only the integer partition.
    permutations = tuple(product(range(3), repeat=3))
    permutations = tuple(
        item for item in permutations if sorted(item) == [0, 1, 2]
    )
    identity = (0, 1, 2)
    nonidentity = tuple(item for item in permutations if item != identity)
    signatures: set[tuple[int, int]] = set()
    signature_witnesses: dict[tuple[int, int], tuple[tuple[int, ...], ...]] = {}
    for length in range(2, 5):
        for branch_cycles in product(nonidentity, repeat=length):
            total = identity
            for cycle in branch_cycles:
                total = compose(cycle, total)
            if total != identity or not transitive(branch_cycles):
                continue
            contributions = tuple(
                ramification_contribution(cycle) for cycle in branch_cycles
            )
            if sum(contributions) != 4:
                continue
            transpositions = contributions.count(1)
            three_cycles = contributions.count(2)
            signature = (transpositions, three_cycles)
            signatures.add(signature)
            signature_witnesses.setdefault(signature, branch_cycles)
    require(
        signatures == {(0, 2), (2, 1), (4, 0)},
        "connected genus-zero degree-three branch signatures changed",
    )
    require(
        all(
            transpositions + 2 * three_cycles == 4
            for transpositions, three_cycles in signatures
        ),
        "Riemann--Hurwitz signature sum changed",
    )

    odd_square_class_degrees = tuple(
        sorted(transpositions for transpositions, _ in signatures)
    )
    require(
        odd_square_class_degrees == (0, 2, 4),
        "odd discriminant square-class degrees changed",
    )
    square_factor_degrees = tuple(
        (12 - degree) // 2 for degree in odd_square_class_degrees
    )
    require(
        square_factor_degrees == (6, 5, 4),
        "square-factor degrees changed",
    )

    print(f"depressed_p=16/964467*({quartic_p})")
    print(f"depressed_q=64/703096443*({sextic_q})")
    print(f"branch_identity=Delta={branch_scalar}*(4*P^3+49*Q^2)")
    print("branch_scalar_factorization=-2^12*3^21*7^2")
    print("mordell_degree=12")
    print("mordell_leading_coefficient=73060029")
    print(f"infinity_cubic={infinity_polynomial}")
    print(f"infinity_discriminant={infinity_discriminant}")
    print("high_support_irreducibility=PASS_C_or_W_coefficient_obstruction")
    print("connected_degree3_genus0_signatures=(3,3),(3,2,2),(2,2,2,2)")
    print("transposition_counts=0,2,4")
    print("odd_square_class_degrees=0,2,4")
    print("square_factor_degrees=6,5,4")
    print("index_discriminant_and_Riemann_Hurwitz=MATHEMATICAL_PROOF_REQUIRED")
    print("status=THM2332_GENUS_ZERO_SQUARE_CLASS_TRAP_EXACT_REFEREE")


if __name__ == "__main__":
    main()
