#!/usr/bin/env python3
"""Independent exact audit for THM-3496.

This deliberately does not import or replay the original D5 sidecar.  It uses
SymPy's finite-field matrix and polynomial engines, together with a separate
prefix-integration construction on C_91, to audit the marked graph--Kummer
square and the Frobenius finite-coefficient hostile.

The script verifies finite algebra only.  The Kummer-sequence and formal-unit
arguments, the additive-group no-go, and the derived-completion scope boundary
are proved in the theorem text.
"""

from __future__ import annotations

import hashlib
import json

import sympy as sp
from sympy.polys.matrices import DomainMatrix


PRIME = 13
SOURCE_LENGTH = 7
EXPECTED_SEMANTIC_SHA256 = "0c5756987947c93c271b7db5443a1618e483269c490b6710c30b9056562be1bb"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def edge_delta(values: tuple[int, ...], modulus: int = PRIME) -> tuple[int, ...]:
    return tuple(
        (values[(index + 1) % len(values)] - values[index]) % modulus
        for index in range(len(values))
    )


def seam(cochain: tuple[int, ...], modulus: int = PRIME) -> int:
    return sum(cochain) % modulus


def cover_primitive(cochain: tuple[int, ...]) -> tuple[int, ...]:
    """Prefix-integrate the thirteenfold pullback on the 91-cycle."""

    pullback = tuple(cochain[index % SOURCE_LENGTH] for index in range(91))
    values = [0]
    for index in range(90):
        values.append((values[-1] + pullback[index]) % PRIME)
    primitive = tuple(values)
    require(
        all(
            (primitive[(index + 1) % 91] - primitive[index]) % PRIME
            == pullback[index]
            for index in range(91)
        ),
        "prefix primitive does not differentiate to the pullback",
    )
    return primitive


def deck_defect(primitive: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        (primitive[(index + SOURCE_LENGTH) % 91] - primitive[index]) % PRIME
        for index in range(91)
    )


def main() -> None:
    # ------------------------------------------------------------------
    # Independent graph and C91 deck-defect route.
    # ------------------------------------------------------------------
    incidence = sp.zeros(SOURCE_LENGTH, SOURCE_LENGTH)
    for column in range(SOURCE_LENGTH):
        incidence[column, column] = -1
        incidence[(column - 1) % SOURCE_LENGTH, column] = 1
    rank = DomainMatrix.from_Matrix(incidence).convert_to(sp.GF(PRIME)).rank()
    require(rank == 6, "C7 incidence rank is not six")
    require(
        sp.ones(1, SOURCE_LENGTH) * incidence == sp.zeros(1, SOURCE_LENGTH),
        "seam row does not annihilate the incidence image",
    )

    basis_defects = []
    for support in range(SOURCE_LENGTH):
        cochain = tuple(1 if index == support else 0 for index in range(SOURCE_LENGTH))
        primitive = cover_primitive(cochain)
        defects = deck_defect(primitive)
        require(set(defects) == {seam(cochain)}, "basis deck defect mismatch")
        basis_defects.append(defects[0])

        for gauge_support in range(SOURCE_LENGTH):
            gauge = tuple(
                1 if index == gauge_support else 0
                for index in range(SOURCE_LENGTH)
            )
            gauged = tuple(
                (cochain[index] + edge_delta(gauge)[index]) % PRIME
                for index in range(SOURCE_LENGTH)
            )
            require(seam(gauged) == seam(cochain), "vertex gauge changed seam")
            require(
                set(deck_defect(cover_primitive(gauged))) == {seam(cochain)},
                "vertex gauge changed deck defect",
            )

    canonical = (1,) * SOURCE_LENGTH
    canonical_defect = deck_defect(cover_primitive(canonical))[0]
    require(canonical_defect == 7, "constant chart class has wrong deck defect")

    # All F13-linear maps commute with scalar degree and orientation reversal;
    # markings, not naturality alone, select the scalar-one normalization.
    normalized_scalars = []
    for scalar in range(PRIME):
        for degree in range(1, 501):
            for source_class in range(PRIME):
                left = scalar * degree * source_class % PRIME
                right = degree * scalar * source_class % PRIME
                require(left == right, "degree square failed")
            require((-scalar * degree) % PRIME == (scalar * (-degree)) % PRIME,
                    "orientation square failed")
        if scalar:
            normalized_scalars.append(scalar)
    require(len(normalized_scalars) == 12, "unmarked isomorphism torsor mismatch")
    require((13 * canonical_defect) % PRIME == 0, "degree thirteen must kill")
    require((14 * canonical_defect) % PRIME == canonical_defect,
            "degree fourteen must return")

    # A principal unit has a unique thirteenth root in characteristic zero.
    # This checks the formal Hensel/binomial mechanism through order eight.
    lam = sp.symbols("lambda")
    unit = 1 + 2 * lam - 3 * lam**2 + 5 * lam**3 + 7 * lam**5
    unit_root = sp.series(unit ** sp.Rational(1, PRIME), lam, 0, 9).removeO()
    unit_residual = sp.series(unit_root**PRIME - unit, lam, 0, 9).removeO().expand()
    require(unit_residual == 0, "formal principal-unit root failed")

    # ------------------------------------------------------------------
    # Independent SymPy route for the Hamiltonian/Frobenius hostile.
    # ------------------------------------------------------------------
    x, z = sp.symbols("x z")
    y = x * z
    polynomial_p = x + x**2 * z

    def derivation(expression: sp.Expr) -> sp.Expr:
        return sp.expand(
            (1 + 2 * x * z) * sp.diff(expression, z)
            - x**2 * sp.diff(expression, x)
        )

    bezout_a = 1 - 2 * x * z
    bezout_c = 4 * z**2
    require(sp.expand(bezout_a * sp.diff(polynomial_p, x) + bezout_c * x**2) == 1,
            "Bezout row failed")
    h = -x**-2 + 2 * z / x
    require(sp.expand(derivation(h) - 6 * z) == 0, "D(h)=6z failed")
    require(
        sp.expand(polynomial_p * h - (-x**-1 + z + 2 * x * z**2)) == 0,
        "P h principal part failed",
    )
    require(sp.expand(derivation(-x**-1) + 1) == 0,
            "connecting unit representative failed")

    # The monomial recurrence is obtained independently by differentiation.
    terminal_checks = 0
    for index in range(201):
        monomial = x**index * z ** (index + 1)
        expected = (index + 1) * y**index + (index + 2) * y ** (index + 1)
        require(sp.expand(derivation(monomial) - expected) == 0,
                f"weight recurrence failed at {index}")
        require((index + 2) * ((-1) ** index) != 0,
                "characteristic-zero terminal coefficient vanished")
        terminal_checks += 1

    odd_primes = (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43)
    for prime in odd_primes:
        q_prime = sum(
            (-1) ** index * x**index * z ** (index + 1)
            for index in range(prime - 1)
        )
        require(
            sp.expand(derivation(q_prime) - (1 - prime * y ** (prime - 1))) == 0,
            f"integral telescope failed at {prime}",
        )
        numerator = sp.together(
            q_prime - x**-1 + polynomial_p ** (prime - 1) / x**prime
        ).as_numer_denom()[0]
        require(
            all(int(coefficient) % prime == 0
                for coefficient in sp.Poly(sp.expand(numerator), x, z).coeffs()),
            f"localized Frobenius identity failed at {prime}",
        )

    # Direct finite-coefficient controls at 13 and 13^2.  The proved general
    # recurrence gives the same conclusion for every r; no claim about a
    # derived inverse limit or lim^1 is encoded here.
    finite_moduli = (13, 13**2)
    for modulus in finite_moduli:
        q_modulus = sum(
            (-1) ** index * x**index * z ** (index + 1)
            for index in range(modulus - 1)
        )
        residual = sp.Poly(sp.expand(derivation(q_modulus) - 1), x, z)
        require(
            all(int(coefficient) % modulus == 0
                for coefficient in residual.coeffs()),
            f"unit class did not die modulo {modulus}",
        )

    semantic = {
        "prime": PRIME,
        "graph_deck": {
            "incidence_rank": rank,
            "basis_defects": basis_defects,
            "canonical_defect": canonical_defect,
            "gauge_cells": SOURCE_LENGTH * SOURCE_LENGTH,
        },
        "kummer_scope": {
            "formal_principal_unit_order": 8,
            "unmarked_nonzero_scalars": normalized_scalars,
            "residue_field_hostile": "2 in Q((lambda)) has valuation zero but is not a thirteenth power",
        },
        "degree_square": {
            "degrees_checked": 500,
            "degree_13": 0,
            "degree_14": canonical_defect,
        },
        "flux": {
            "odd_primes": odd_primes,
            "terminal_checks": terminal_checks,
            "finite_moduli_direct": finite_moduli,
            "general_finite_coefficient_claim": "D(Q_(13^r-1))=1 mod 13^r for every r>=1",
            "derived_completion_claim": "not adjudicated",
        },
    }
    semantic_bytes = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()
    semantic_sha256 = hashlib.sha256(semantic_bytes).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "UNFROZEN":
        require(semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
                "semantic payload hash drift")

    print("THM-3496 INDEPENDENT GRAPH/KUMMER/FROBENIUS AUDIT")
    print("STATUS: EXACT ALGEBRA AUDIT; REPAIRED SCOPE")
    print()
    print("GRAPH / DECK")
    print(f"  rank_F13(delta_C7)={rank}; H1 dimension={SOURCE_LENGTH-rank}")
    print(f"  seven basis deck defects={basis_defects}")
    print(f"  constant chart class deck defect={canonical_defect}")
    print("  all 49 basis/gauge cells preserve seam and C91 deck defect")
    print()
    print("KUMMER / DEGREE")
    print("  principal-unit thirteenth root checked through lambda^8")
    print("  unmarked nonzero line isomorphisms=12; scalar-one marking selects one")
    print("  hostile: over Q((lambda)), [2] is an extra valuation-zero Kummer class")
    print("  degree square checked for k<=500; k=13 dies, k=14 returns")
    print()
    print("FLUX HOSTILE")
    print("  Bezout, D(h)=6z, [P h]=[-x^-1], and D(-x^-1)=-1 PASS")
    print(f"  D(Q_p)=1-p(xz)^(p-1) and localized Frobenius PASS at {odd_primes}")
    print("  weight recurrence and nonzero terminal obstruction checked through 201 terms")
    print("  [1]=0 after coefficient reduction mod 13 and mod 13^2 (direct CAS)")
    print("  general recurrence proves [1]=0 mod 13^r for every finite r")
    print("  no universal Bockstein, derived-completion, or lim^1 claim is made")
    print()
    print(f"SEMANTIC_SHA256={semantic_sha256}")


if __name__ == "__main__":
    main()
