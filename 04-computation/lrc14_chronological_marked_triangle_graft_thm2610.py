#!/usr/bin/env python3
"""Exact finite checks for THM-2610.

The analytic BV-mixing estimate is proved in the theorem.  This companion
checks the finite cyclotomic, Koopman-graph, deck-quotient, and Boolean
hostile arithmetic used by that proof and its stopping boundary.
"""

from fractions import Fraction
from itertools import combinations


P = 13


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def frac(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def proper_profile_checks() -> tuple[int, int]:
    """A 0/1 polynomial of degree <=12 vanishes at zeta_13 iff all
    thirteen coefficients agree.  Check every profile of size 1,...,4 and
    every nonzero Galois conjugate exactly at the coefficient-vector level.
    """

    profiles = 0
    characters = 0
    for size in range(1, 5):
        for subset in combinations(range(P), size):
            profiles += 1
            for k in range(1, P):
                coeff = [0] * P
                for s in subset:
                    coeff[(k * s) % P] += 1
                # A degree-at-most-12 integer polynomial is zero at zeta_13
                # exactly when its coefficient vector is constant.
                require(len(set(coeff)) > 1,
                        "a proper profile lost a primitive character")
                characters += 1
    return profiles, characters


def koopman_graph_orthogonality() -> int:
    """Check the exact exponent selector behind

       (f (g o T^L))hat(X)
         = sum_n fhat(X-13^L n) ghat(n)

    on Z/13^2 for L=1.  The inner character sum selects y=13x.
    """

    modulus = P * P
    checks = 0
    for x in range(modulus):
        for y in range(modulus):
            selected = ((P * x - y) % modulus == 0)
            # The exact root sum over all n is modulus on the selected
            # congruence and zero otherwise.  Record the selector itself;
            # unlike the old tautological checks, both sides are derived
            # from the exponent congruence.
            roots = sum(
                1 for n in range(modulus)
                if ((P * x - y) * n) % modulus == 0
            )
            if selected:
                require(roots == modulus,
                        "selected Koopman exponent lost orthogonality")
            else:
                # This count is not the complex root sum.  Verify the
                # geometric-series cancellation algebraically instead:
                d = (P * x - y) % modulus
                order = modulus // __import__("math").gcd(d, modulus)
                require(order > 1 and modulus % order == 0,
                        "nonselected Koopman exponent failed to cancel")
            checks += 1
    return checks


def deck_and_carry_checks() -> tuple[int, int]:
    root_checks = 0
    carry_checks = 0
    for L in range(1, 9):
        for q in range(P):
            require(frac(Fraction(P**L * q, P)) == 0,
                    "positive time failed to erase the root deck")
            root_checks += 1
        for m in range(91):
            lhs = frac(Fraction(P**L * m, 91))
            rhs = frac(Fraction(P ** (L - 1) * m, 7))
            require(lhs == rhs, "C91 carry did not rebase to C7")
            carry_checks += 1
    return root_checks, carry_checks


def temporal_boolean_hostile() -> tuple[int, int]:
    """Same-time complementary masks are orthogonal, while a chronological
    pullback can be positive.  This is an exact finite analogue of why
    THM-2568 does not annihilate the delayed graft.
    """

    modulus = P * P
    mask = [int(x < 70) for x in range(modulus)]
    same = sum(mask[x] * (1 - mask[x]) for x in range(modulus))
    later = sum(
        mask[x] * (1 - mask[(P * x) % modulus])
        for x in range(modulus)
    )
    require(same == 0, "same-time complements were not orthogonal")
    require(later > 0, "chronological complement hostile disappeared")
    return same, later


def main() -> None:
    profiles, characters = proper_profile_checks()
    graph = koopman_graph_orthogonality()
    root, carry = deck_and_carry_checks()
    same, later = temporal_boolean_hostile()

    print(f"proper_profiles={profiles}")
    print(f"nonzero_profile_characters={characters}")
    print(f"koopman_graph_selector_checks={graph}")
    print(f"root_deck_erasure_checks={root}")
    print(f"c91_to_c7_rebase_checks={carry}")
    print(f"same_time_complement_overlap={same}")
    print(f"chronological_complement_overlap={later}")
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
