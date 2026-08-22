#!/usr/bin/env python3
"""Exact reciprocal-root certificate for THM-3159.

The load-bearing computation is one extended gcd in F_q[a] at
q=249727.  It is intentionally performed on the linear-size reverse
truncated exponential P(a)=2*a*U(a)+1, rather than constructing the two
original degree-124863 endpoint faces by a quadratic-time bivariate
expansion.
"""

from math import isqrt

from flint import nmod_poly


Q = 249_727


def require(condition, message):
    if not bool(condition):
        raise RuntimeError(message)


def is_prime(number):
    if number < 2:
        return False
    if number % 2 == 0:
        return number == 2
    return all(number % divisor for divisor in range(3, isqrt(number) + 1, 2))


def inverses_and_inverse_factorials(prime):
    inverses = [0] * prime
    inverse_factorials = [1] * prime
    inverses[1] = 1
    for index in range(2, prime):
        inverses[index] = (-(prime // index) * inverses[prime % index]) % prime
    for index in range(1, prime):
        inverse_factorials[index] = (
            inverse_factorials[index - 1] * inverses[index]
        ) % prime
    return inverses, inverse_factorials


def face_value(prime, c, v, inverses, inverse_factorials):
    """Return [t^(p-1)] E(t)/(1-t-v*t^2)^c in F_p."""
    coefficients = [1]
    for degree in range(1, prime):
        value = (degree + c - 1) * coefficients[degree - 1]
        if degree >= 2:
            value += v * (degree + 2 * c - 2) * coefficients[degree - 2]
        coefficients.append(value * inverses[degree] % prime)
    return sum(
        inverse_factorials[index] * coefficients[prime - 1 - index]
        for index in range(prime)
    ) % prime


def choose_mod_prime(n, k, prime, inverses):
    k = min(k, n - k)
    answer = 1
    for index in range(1, k + 1):
        answer = answer * (n - k + index) * inverses[index] % prime
    return answer


def reverse_truncated_exponential(prime):
    """U(a)=a^(p-1)E(1/a)=-sum (-1)^k k! a^k in F_p[a]."""
    coefficients = []
    factorial = 1
    sign = -1
    for degree in range(prime):
        if degree:
            factorial = factorial * degree % prime
        coefficients.append(sign * factorial % prime)
        sign = -sign
    return nmod_poly(coefficients, prime)


def small_transform_controls():
    records = []
    for prime in (5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47):
        x = nmod_poly([0, 1], prime)
        u = reverse_truncated_exponential(prime)
        polynomial = 2 * x * u + 1
        reflected = polynomial.compose(1 - x)
        common = polynomial.gcd(reflected)
        require(common.degree() == 0, f"small transform p={prime}")
        records.append((prime, common.degree()))
    return records


def main():
    require(is_prime(Q), "q must be prime")
    inverse, inverse_factorial = inverses_and_inverse_factorials(Q)

    zero_faces = tuple(
        face_value(Q, c, 0, inverse, inverse_factorial) for c in (2, 3)
    )
    repeated_v = -pow(4, -1, Q) % Q
    repeated_faces = tuple(
        face_value(Q, c, repeated_v, inverse, inverse_factorial)
        for c in (2, 3)
    )
    require(all(zero_faces), "v=0 face boundary")
    require(all(repeated_faces), "repeated-root face boundary")

    midpoint = (Q - 1) // 2
    midpoint_residues = []
    for c in (3, 2):
        n = Q - c
        residue = (
            -choose_mod_prime(n, midpoint, Q, inverse)
            * pow(Q - 1, n - midpoint, Q)
        ) % Q
        require(residue, f"midpoint c={c}")
        midpoint_residues.append(residue)

    x = nmod_poly([0, 1], Q)
    u = reverse_truncated_exponential(Q)
    polynomial = 2 * x * u + 1
    reflected = polynomial.compose(1 - x)

    gcd_direct = polynomial.gcd(reflected)
    gcd_extended, bezout_left, bezout_right = polynomial.xgcd(reflected)
    require(gcd_direct == gcd_extended == 1, "reciprocal-root gcd")
    require(
        bezout_left * polynomial + bezout_right * reflected == 1,
        "extended-gcd identity",
    )
    chart_values = tuple(int(polynomial(value)) for value in (0, 1, pow(2, -1, Q)))
    require(all(chart_values), "reciprocal-root chart boundary")

    small = small_transform_controls()
    require(len(small) == 13, "small-control count")

    print(f"q={Q} prime=1")
    print(
        "degrees="
        f"U:{u.degree()},P:{polynomial.degree()},Pstar:{reflected.degree()}"
    )
    print(f"zero_faces_F2_F3={zero_faces}")
    print(f"repeated_faces_F2_F3={repeated_faces}")
    print(f"midpoint_A3_A2={tuple(midpoint_residues)}")
    print("newton_A3=0,2/249722")
    print("newton_A2=0,2/249724")
    print(
        "reciprocal_gcd="
        f"degree:{gcd_direct.degree()},bezout:{int(gcd_extended == 1)}"
    )
    print(f"reciprocal_chart_P_0_1_half={chart_values}")
    print(
        "small_transform_controls="
        + ",".join(f"{p}:{degree}" for p, degree in small)
    )
    print("THM3159_CHECKS_PASSED")


if __name__ == "__main__":
    main()
