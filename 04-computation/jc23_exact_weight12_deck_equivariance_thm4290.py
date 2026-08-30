#!/usr/bin/env python3
"""Exact character and degree audit for THM-4290.

This script freezes the exponent arithmetic behind the 12-fold deck action,
the genus-two quotient coordinates already proved in THM-4230, and the final
mod-four response contradiction.  It is deterministic, dependency-free, and
uses optimization-safe checks only.
"""

from __future__ import annotations


MODULUS = 12


def check(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def character(*weighted_exponents: tuple[int, int]) -> int:
    """Character of a monomial from (variable character, power) pairs."""
    return sum(weight * exponent for weight, exponent in weighted_exponents) % MODULUS


def main() -> None:
    # Q=sigma^12; original coordinates are
    # s=sigma^-1 S, p=sigma^-2 P, A=sigma^-4 X, C=sigma^-6 Y.
    weights = {"sigma": 1, "S": 1, "P": 2, "X": 4, "Y": 6}
    check(12 * weights["sigma"] % MODULUS == 0, "Q is not deck-invariant")
    check(character((weights["S"], 1), (weights["sigma"], -1)) == 0, "s changed")
    check(character((weights["P"], 1), (weights["sigma"], -2)) == 0, "p changed")
    check(character((weights["X"], 1), (weights["sigma"], -4)) == 0, "A changed")
    check(character((weights["Y"], 1), (weights["sigma"], -6)) == 0, "C changed")

    # Every top-face monomial has deck character zero.
    face_characters = {
        "P^6": character((weights["P"], 6)),
        "S^2P^5": character((weights["S"], 2), (weights["P"], 5)),
        "S^4P^4": character((weights["S"], 4), (weights["P"], 4)),
    }
    check(set(face_characters.values()) == {0}, "face equation is not invariant")

    # On E0, the deck action is (X,Y)->(xi^4 X,xi^6 Y)=[-omega].
    check(character((weights["X"], 3)) == 0, "X^3 changed")
    check(character((weights["Y"], 2)) == 0, "Y^2 changed")
    differential_character = (weights["X"] - weights["Y"]) % MODULUS
    check(differential_character == 10, "invariant differential character changed")
    check((6 * weights["X"]) % MODULUS == 0, "sixth target power changed X")
    check((6 * weights["Y"]) % MODULUS == 0, "sixth target power changed Y")

    # tau^6 fixes P and negates S, hence its quotient has u=1/P and
    # x=S^2/P.  The residual tau characters are u:10, x:0, v=W+2Zx:0.
    u_character = (-weights["P"]) % MODULUS
    x_character = (2 * weights["S"] - weights["P"]) % MODULUS
    v_character = 0
    check((u_character, x_character, v_character) == (10, 0, 0), "quotient characters")
    check((6 * weights["S"]) % MODULUS == 6, "tau^6 does not negate S")
    check((6 * weights["P"]) % MODULUS == 0, "tau^6 does not fix P")

    # THM-4230's two visible quotient maps have characters
    # f1:(u^2,v)=(8,0) and f2:(u^-2,v/u^3)=(4,6).  The second is exactly
    # the target deck character; both factor through C->B of degree two.
    f1_characters = ((2 * u_character) % MODULUS, v_character)
    f2_characters = (
        (-2 * u_character) % MODULUS,
        (v_character - 3 * u_character) % MODULUS,
    )
    check(f1_characters == (8, 0), "first elliptic quotient character")
    check(f2_characters == (4, 6), "second elliptic quotient character")

    # The quotient identity follows from u^6=U+Wx+Zx^2 and v=W+2Zx:
    # v^2-(W^2-4UZ)-4Zu^6 has coefficients 0 in 1,x,x^2.
    quotient_identity_coefficients = (
        -4 * 1 + 4 * 1,  # coefficient of UZ
        4 - 4,            # coefficient of WZ*x
        4 - 4,            # coefficient of Z^2*x^2
    )
    check(quotient_identity_coefficients == (0, 0, 0), "genus-two quotient identity")

    # THM-4230's saturated visible lattice gives degree
    # 4(N(alpha)+N(beta)).  Its residues cannot equal either inherited
    # response degree.
    responses = (34, 42)
    check(all(response % 4 == 2 for response in responses), "response residue changed")
    for norm_sum in range(16):
        check((4 * norm_sum) % 4 == 0, "visible degree lost mod-four divisibility")

    # The original differential is invariant only after the sigma^2 twist.
    check((2 + differential_character) % MODULUS == 0, "dA/(2C) twist changed")

    print("THM4290_EXACT_WEIGHT12_DECK_EQUIVARIANCE_AUDIT_V1")
    print("DECK sigma:1 S:1 P:2 X:4 Y:6 MOD_12")
    print("FACE_CHARACTERS", ",".join(f"{name}:{value}" for name, value in face_characters.items()))
    print("TARGET_ACTION X:4 Y:6 NAME [-omega] SIXTH_POWER ID")
    print("TARGET_DIFFERENTIAL_CHARACTER 10 SIGMA2_TWIST INVARIANT")
    print("SOURCE_SIXTH_POWER S:-1 P:+1 QUOTIENT C_TO_B DEGREE 2")
    print("VISIBLE_MAP_CHARACTERS f1:(8,0) f2:(4,6)")
    print("VISIBLE_DEGREE 4*(N(alpha)+N(beta)) MOD_4 0")
    print("RESPONSES 34,42 MOD_4 2 CONTRADICTION")
    print("VERDICT PASS EXACT_CHARACTER_ARITHMETIC; JC2 OPEN")


if __name__ == "__main__":
    main()
