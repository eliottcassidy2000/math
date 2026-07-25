#!/usr/bin/env python3
"""Exact checks for THM-2334.

The proof is symbolic.  This companion checks:

* primitive relation-lattice reduction modulo N;
* a finite-support instance of the exact-address pushforward;
* the character-twist formula, representative invariance, and inversion;
* the 13^2 target-character orthogonality kernel;
* target neutrality modulo 13 and activity modulo 7 of the word dilation;
* an exact positive-factor hostile with a zero exact-address orbit sum.
"""

from fractions import Fraction
from itertools import product
from math import gcd


def dot(x, y):
    return sum(a * b for a, b in zip(x, y))


def add_poly(a, b):
    return tuple(x + y for x, y in zip(a, b))


def mul_poly(a, b):
    n = len(a)
    out = [Fraction(0) for _ in range(n)]
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[(i + j) % n] += x * y
    return tuple(out)


def shift_poly(a, shift):
    n = len(a)
    out = [Fraction(0) for _ in range(n)]
    for i, x in enumerate(a):
        out[(i + shift) % n] += x
    return tuple(out)


def conj_poly(a):
    n = len(a)
    out = [Fraction(0) for _ in range(n)]
    for i, x in enumerate(a):
        out[(-i) % n] += x
    return tuple(out)


def monomial_poly(n, exponent, coefficient):
    out = [Fraction(0) for _ in range(n)]
    out[exponent % n] = coefficient
    return tuple(out)


def reduce_prime_cyclotomic(a):
    """Coordinates in Q[zeta_p] on 1,zeta,...,zeta^(p-2)."""
    last = a[-1]
    return tuple(x - last for x in a[:-1])


def relation_quotient_checks():
    w = (2, 3, 5)
    z = (-1, 1, 0)
    assert dot(z, w) == 1
    total = 0
    for n in range(2, 14):
        kernel = [
            y
            for y in product(range(n), repeat=3)
            if dot(y, w) % n == 0
        ]
        assert len(kernel) == n**2
        for y in kernel:
            h = dot(y, w) // n
            r = tuple(y[i] - n * h * z[i] for i in range(3))
            assert dot(r, w) == 0
            assert all((r[i] - y[i]) % n == 0 for i in range(3))
        total += len(kernel)
    return total


def finite_pushforward_check():
    w = (1, 2, 3)
    dilation = 4
    x_frequency = 0
    deepest_coordinate = 2
    multiplier = 1
    y_frequency = x_frequency + multiplier * w[deepest_coordinate]
    deep_coefficient = Fraction(2, 5)

    endpoint = (
        {-1: Fraction(1, 3), 0: Fraction(1), 1: Fraction(-1, 4)},
        {-1: Fraction(-2, 5), 0: Fraction(3, 2), 1: Fraction(1, 5)},
        {-1: Fraction(1, 6), 0: Fraction(4, 3), 1: Fraction(2, 7)},
    )
    word = (
        {-1: Fraction(1, 2), 0: Fraction(2), 1: Fraction(-1, 3)},
        {-1: Fraction(1, 7), 0: Fraction(5, 4), 1: Fraction(3, 8)},
        {-1: Fraction(-1, 5), 0: Fraction(7, 6), 1: Fraction(1, 9)},
    )

    support = (-1, 0, 1)
    left_terms = []
    for u in product(support, repeat=3):
        au = endpoint[0][u[0]] * endpoint[1][u[1]] * endpoint[2][u[2]]
        for beta in product(support, repeat=3):
            combined = tuple(u[i] + dilation * beta[i] for i in range(3))
            if dot(combined, w) != x_frequency:
                continue
            bb = word[0][beta[0]] * word[1][beta[1]] * word[2][beta[2]]
            left_terms.append((u, beta, combined, au * bb))

    right_terms = []
    for v in product(support, repeat=3):
        if dot(v, w) != y_frequency:
            continue
        av = endpoint[0][v[0]] * endpoint[1][v[1]] * endpoint[2][v[2]]
        right_terms.append((v, av))

    left_total = sum((term[3] for term in left_terms), Fraction(0))
    right_total = sum((term[1] for term in right_terms), Fraction(0))
    current = deep_coefficient * left_total * right_total
    assert current != 0

    address_weights = {}
    term_count = 0
    for u, beta, combined, left_weight in left_terms:
        for v, right_weight in right_terms:
            r = list(combined)
            r[deepest_coordinate] += multiplier
            r = tuple(r[i] - v[i] for i in range(3))
            assert dot(r, w) == 0
            address_weights[r] = (
                address_weights.get(r, Fraction(0))
                + deep_coefficient * left_weight * right_weight
            )
            term_count += 1
    assert sum(address_weights.values(), Fraction(0)) == current

    modulus = 5
    zero_poly = tuple(Fraction(0) for _ in range(modulus))
    direct_transform = {}
    factor_transform = {}
    for ell in product(range(modulus), repeat=3):
        direct = zero_poly
        for r, coefficient in address_weights.items():
            direct = add_poly(
                direct,
                monomial_poly(modulus, dot(ell, r), coefficient),
            )
        direct_transform[ell] = direct

        left_twist = zero_poly
        for _, _, combined, coefficient in left_terms:
            left_twist = add_poly(
                left_twist,
                monomial_poly(modulus, dot(ell, combined), coefficient),
            )
        right_twist = zero_poly
        for v, coefficient in right_terms:
            right_twist = add_poly(
                right_twist,
                monomial_poly(modulus, dot(ell, v), coefficient),
            )
        deep_phase = monomial_poly(
            modulus,
            multiplier * ell[deepest_coordinate],
            deep_coefficient,
        )
        factorized = mul_poly(
            mul_poly(deep_phase, left_twist),
            conj_poly(right_twist),
        )
        factor_transform[ell] = factorized
        assert direct == factorized

        for scalar in range(modulus):
            shifted = tuple(
                (ell[i] + scalar * w[i]) % modulus for i in range(3)
            )
            if shifted in direct_transform:
                assert direct_transform[shifted] == direct

    residue_weights = {}
    for r, coefficient in address_weights.items():
        residue = tuple(x % modulus for x in r)
        residue_weights[residue] = residue_weights.get(residue, Fraction(0)) + coefficient

    for y in product(range(modulus), repeat=3):
        inverse_poly = zero_poly
        for ell, transformed in direct_transform.items():
            inverse_poly = add_poly(
                inverse_poly,
                shift_poly(transformed, -dot(ell, y)),
            )
        reduced = reduce_prime_cyclotomic(inverse_poly)
        expected = modulus**3 * residue_weights.get(y, Fraction(0))
        assert reduced[0] == expected
        assert all(value == 0 for value in reduced[1:])

    return term_count, len(address_weights), current, len(direct_transform)


def target_character_kernel_check():
    prime = 13
    group = list(product(range(prime), repeat=2))
    checked = 0
    for delta in group:
        histogram = [0 for _ in range(prime)]
        for ell in group:
            histogram[dot(ell, delta) % prime] += 1
        if delta == (0, 0):
            assert histogram[0] == prime**2
            assert sum(histogram[1:]) == 0
        else:
            assert histogram == [prime for _ in range(prime)]
        checked += 1
    return len(group), checked


def word_charge_check():
    residues_mod_7 = set()
    for shallow_depth in range(1, 13):
        dilation = 13 ** (shallow_depth + 1)
        assert dilation % 13 == 0
        assert gcd(dilation, 7) == 1
        residues_mod_7.add(dilation % 7)
    assert residues_mod_7 == {1, 6}
    return sorted(residues_mod_7)


def positive_factor_hostile():
    # f_1 has Fourier coefficients below and is strictly positive because
    # its constant coefficient dominates twice the sum of the off-zero
    # absolute coefficients.  f_2 is the Poisson kernel with q=9/10.
    a1 = {
        -3: Fraction(10, 81),
        -2: Fraction(-10, 81),
        -1: Fraction(1, 90),
        0: Fraction(1),
        1: Fraction(1, 90),
        2: Fraction(-10, 81),
        3: Fraction(10, 81),
    }
    q = Fraction(9, 10)
    lower_bound = Fraction(1) - 2 * (
        abs(a1[1]) + abs(a1[2]) + abs(a1[3])
    )
    assert lower_bound == Fraction(196, 405)
    assert lower_bound > 0

    a = {k: value * q ** abs(k) for k, value in a1.items()}
    expected = {
        -3: Fraction(9, 100),
        -2: Fraction(-1, 10),
        -1: Fraction(1, 100),
        0: Fraction(1),
        1: Fraction(1, 100),
        2: Fraction(-1, 10),
        3: Fraction(9, 100),
    }
    assert a == expected
    endpoint_coefficient = sum(a.values(), Fraction(0))
    assert endpoint_coefficient == 1

    terms = [
        a[k] * a[k - 1]
        for k in sorted(a)
        if k - 1 in a
    ]
    assert len(terms) == 6
    assert all(term != 0 for term in terms)
    address_coefficient = sum(terms, Fraction(0))
    assert address_coefficient == 0
    full_current = endpoint_coefficient**2
    assert full_current == 1
    return lower_bound, len(terms), address_coefficient, full_current


def main():
    quotient_points = relation_quotient_checks()
    term_count, address_count, current, twist_count = finite_pushforward_check()
    target_characters, target_kernels = target_character_kernel_check()
    mod_7_residues = word_charge_check()
    lower, hostile_terms, hostile_coefficient, hostile_current = positive_factor_hostile()

    print("THM-2334 exact companion")
    print(f"relation_quotient_points_checked={quotient_points}")
    print(f"finite_pushforward_terms={term_count}")
    print(f"finite_exact_addresses={address_count}")
    print(f"finite_current={current}")
    print(f"mod5_character_twists={twist_count}")
    print("character_factorization=PASS")
    print("representative_invariance=PASS")
    print("residue_inversion=PASS")
    print(f"target_character_count={target_characters}")
    print(f"target_orthogonality_kernels={target_kernels}")
    print("word_charge_mod13=0")
    print(f"word_charge_mod7_units={mod_7_residues}")
    print(f"hostile_f1_lower_bound={lower}")
    print(f"hostile_full_current={hostile_current}")
    print(f"hostile_address_terms={hostile_terms}")
    print(f"hostile_address_coefficient={hostile_coefficient}")
    print("all checks passed")


if __name__ == "__main__":
    main()
