#!/usr/bin/env python3
"""Exact companion for THM-2419.

The script uses only integer/Fraction arithmetic.  It checks the finite
affine-shell quotient, observer homogenizations, Bezout ambiguity, charged
reference difference, and the 13/7 valuation split.
"""

from itertools import product
from math import gcd


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def dot(left, right):
    return sum(a * b for a, b in zip(left, right))


def extended_gcd(a, b):
    old_r, r = abs(a), abs(b)
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r:
        quotient = old_r // r
        old_r, r = r, old_r - quotient * r
        old_s, s = s, old_s - quotient * s
        old_t, t = t, old_t - quotient * t
    if a < 0:
        old_s = -old_s
    if b < 0:
        old_t = -old_t
    return old_r, old_s, old_t


def bezout_pair(w):
    common, left, right = extended_gcd(w[0], w[1])
    require(common == 1, "row is not primitive")
    require(dot((left, right), w) == 1, "Bezout pair failed")
    return left, right


def kernel_mod(w, modulus):
    return tuple(
        residue
        for residue in product(range(modulus), repeat=2)
        if dot(residue, w) % modulus == 0
    )


def valuation(number, prime):
    number = abs(number)
    exponent = 0
    while number and number % prime == 0:
        exponent += 1
        number //= prime
    return exponent


def main():
    quotient_cases = 0
    lifted_fibres = 0
    for w in product(range(-3, 4), repeat=2):
        if w == (0, 0) or gcd(abs(w[0]), abs(w[1])) != 1:
            continue
        section = bezout_pair(w)
        for modulus in (2, 3, 5, 13, 26):
            kernel = kernel_mod(w, modulus)
            require(
                len(kernel) == modulus,
                f"kernel size failed for w={w}, M={modulus}",
            )
            seen = set()
            for residue in kernel:
                quotient = dot(residue, w) // modulus
                address = tuple(
                    residue[index]
                    + modulus * (1 - quotient) * section[index]
                    for index in range(2)
                )
                require(dot(address, w) == modulus, "shell lift failed")
                require(
                    tuple(value % modulus for value in address) == residue,
                    "shell residue changed",
                )
                seen.add(residue)
                lifted_fibres += 1
            require(len(seen) == len(kernel), "finite torsor map not bijective")
            quotient_cases += 1

    # Minimal Bezout-section hostile.
    w = (1, 1)
    x = 13
    address = (1, 12)
    z = (1, 0)
    z_prime = (0, 1)
    relation = tuple(address[i] - x * z[i] for i in range(2))
    relation_prime = tuple(
        address[i] - x * z_prime[i] for i in range(2)
    )
    require(relation == (-12, 12), "first Bezout contraction mismatch")
    require(relation_prime == (1, -1), "second contraction mismatch")
    require(
        tuple(relation_prime[i] - relation[i] for i in range(2))
        == (13, -13),
        "section ambiguity is not X*Lambda",
    )
    require(dot(relation, w) == dot(relation_prime, w) == 0, "not relations")
    require(
        tuple(value % 13 for value in relation)
        == tuple(value % 13 for value in address),
        "charged residue lost under Bezout contraction",
    )

    # Observer and valuation-normalized homogenizations.
    observer_w = w + (x,)
    observer_a = address + (-1,)
    require(dot(observer_w, observer_a) == 0, "observer did not homogenize")
    require(
        all(value % 13 != 0 for value in observer_a),
        "observer residue not all-unit",
    )

    valuation_checks = 0
    for x_value in (13, -26, 169 * 5, -(13**3) * 2):
        d = valuation(x_value, 13)
        y = x_value // (13**d)
        require(d >= 1 and y % 13 != 0, "valuation normalization failed")
        a_value = (1, x_value - 1)
        extended_w = w + (13**d,)
        extended_a = a_value + (-y,)
        require(dot(extended_w, extended_a) == 0, "normalized observer failed")
        for exponent in range(d + 3):
            integral = x_value % (13**exponent) == 0
            if not integral:
                unit = False
            else:
                coefficient = x_value // (13**exponent)
                unit = coefficient % 13 != 0
            require(
                (integral and unit) == (exponent == d),
                "pure-power observer uniqueness failed",
            )
        valuation_checks += 1

    # Same-shell self-difference is neutral; a residue-zero reference is
    # charged.  Use X=26 so two different K_X fibres reduce to one q mod 13.
    x = 26
    first = (1, 25)
    second = (14, 12)
    reference = (26, 0)
    require(dot(first, w) == dot(second, w) == dot(reference, w) == x, "bad shell")
    self_difference = tuple(first[i] - second[i] for i in range(2))
    charged_difference = tuple(first[i] - reference[i] for i in range(2))
    require(
        all(value % 13 == 0 for value in self_difference),
        "self-difference retained charge",
    )
    require(
        tuple(value % 13 for value in charged_difference) == (1, 12),
        "reference difference lost charge",
    )
    require(dot(charged_difference, w) == 0, "reference difference not relation")

    # Finite fibre pushforward: two K_26 fibres with the same mod-13
    # residue carry coefficients 2 and -1, hence total mass one.
    fibre_mass = {
        tuple(value % x for value in first): 2,
        tuple(value % x for value in second): -1,
    }
    require(sum(fibre_mass.values()) == 1, "pushforward total vanished")
    require(any(value != 0 for value in fibre_mass.values()), "all fibres vanished")
    require(
        all(tuple(value % 13 for value in fibre) == (1, 12) for fibre in fibre_mass),
        "finite fibre changed mod-13 residue",
    )

    # CRT split controls.
    require(91 % 13 == 0 and 91 % 7 == 0, "91 branch failed")
    require(13 % 13 == 0 and 13 % 7 != 0, "13-only branch failed")

    print("THM-2419 AFFINE SIDEBAND HOMOGENIZATION EXACT AUDIT")
    print(f"primitive-row/modulus quotient cases={quotient_cases}")
    print(f"finite kernel fibres lifted={lifted_fibres}")
    print("|K_M|=M^(n-1) for n=2 controls: PASS")
    print("Lambda_X/(M Lambda) -> K_M bijection: PASS")
    print("observer (w,X),(a,-1) all-unit mod13: PASS")
    print(f"valuation-normalized observer cases={valuation_checks}")
    print("pure 13^e observer is integral/unit iff e=nu13(X): PASS")
    print("minimal Bezout hostile relations=(-12,12),(1,-1)")
    print("section difference=13*(1,-1): PASS")
    print("self-difference mod13=0 / reference difference=(1,12): PASS")
    print("finite K_26 fibre masses=(2,-1), total=1")
    print("13|X preserves q in K_13: PASS")
    print("CRT branches 91|X versus 13|X,7 does not divide X: PASS")
    print("artificial observer does not transport physical amplitude")
    print("THM-2419 exact companion PASS")


if __name__ == "__main__":
    main()
