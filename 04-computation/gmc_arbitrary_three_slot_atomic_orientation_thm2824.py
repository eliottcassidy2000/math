#!/usr/bin/env python3
"""Exact companion for THM-2824.

The script verifies the binary quadratic/cubic reduction for normalized
factorial monomials f_n=s^n/n!, the strict cubic orientation t_111>0, the
Wronskian formula, and the universal atomic proof certificate:

* the hockey-stick K/H expansion;
* positivity of the three-step q weight;
* strict separation of the tilted q and interval means;
* the exact likelihood-ratio recursion for successive rows;
* strict monotonicity of R_n; and
* the final weighted-secant factorization of every atom.

The proof is algebraic in the theorem.  This finite exact program audits every
formula and equality boundary through c<=50.  It uses only integers and
fractions, has explicit exception gates, and has no scratch dependency.
"""

from fractions import Fraction
from functools import lru_cache
from math import comb, factorial


MAX_C = 50


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


@lru_cache(maxsize=None)
def moment(indices):
    total = sum(indices)
    denominator = 1
    for index in indices:
        denominator *= factorial(index)
    return factorial(total) // denominator


def bilinear(left, right):
    return sum(
        a_value * b_value * moment(tuple(sorted((a_index, b_index))))
        for a_index, a_value in left.items()
        for b_index, b_value in right.items()
    )


def trilinear(left, middle, right):
    return sum(
        a_value
        * b_value
        * c_value
        * moment(tuple(sorted((a_index, b_index, c_index))))
        for a_index, a_value in left.items()
        for b_index, b_value in middle.items()
        for c_index, c_value in right.items()
    )


def add_sparse(*polynomials):
    out = {}
    for polynomial in polynomials:
        for index, value in polynomial.items():
            out[index] = out.get(index, 0) + value
            if out[index] == 0:
                del out[index]
    return out


def scale_sparse(polynomial, scalar):
    return {
        index: scalar * value
        for index, value in polynomial.items()
        if scalar * value
    }


def multiply_sparse(left, right):
    out = {}
    for i, a in left.items():
        for j, b in right.items():
            out[i + j] = out.get(i + j, Fraction(0)) + a * b
    return {index: value for index, value in out.items() if value}


def derivative_sparse(polynomial):
    return {
        index - 1: index * value
        for index, value in polynomial.items()
        if index
    }


def normalized_monomial(index):
    return {index: Fraction(1, factorial(index))}


def factorial_ratio(high, low):
    require(high >= low, "factorial ratio order drift")
    return factorial(high) // factorial(low)


def invariants(a, b, c):
    f_a = {a: 1}
    f_b = {b: 1}
    f_c = {c: 1}
    u = add_sparse(f_b, scale_sparse(f_a, -1))
    v = add_sparse(f_c, scale_sparse(f_b, -1))

    g11 = bilinear(u, u)
    g12 = bilinear(u, v)
    g22 = bilinear(v, v)
    t111 = trilinear(u, u, u)
    t112 = trilinear(u, u, v)
    t122 = trilinear(u, v, v)
    t222 = trilinear(v, v, v)

    i1 = (
        3 * t112 * g11 * g22
        - t222 * g11 * g11
        - 2 * t111 * g12 * g22
    )
    i2 = (
        3 * t122 * g11 * g22
        - 2 * t222 * g12 * g11
        - t111 * g22 * g22
    )
    d_total = 2 * t222 * g12 - 3 * t122 * g22
    return (
        u,
        v,
        g11,
        g12,
        g22,
        t111,
        t112,
        t122,
        t222,
        i1,
        i2,
        d_total,
    )


def atomic_d(i, b, c):
    d_i = {i + 1: 1, i: -1}
    v = {c: 1, b: -1}
    g22 = bilinear(v, v)
    t222 = trilinear(v, v, v)
    return (
        2 * t222 * bilinear(d_i, v)
        - 3 * trilinear(d_i, v, v) * g22
    )


def hockey_h(n, j):
    require(n >= 0 and j >= 0, "hockey H index drift")
    return comb(n + j, j)


def hockey_k(n, m):
    require(n >= 0 and m >= 1, "hockey K index drift")
    return comb(n + m, m - 1)


def verify_tilted_ratio_certificate(b, c):
    """Audit the universal R_n monotonicity certificate for one (b,c)."""

    require(1 <= b < c, "tilted-ratio support order drift")
    gap = c - b
    alpha = comb(2 * b, b)
    beta = comb(b + c, b)
    gamma = comb(2 * c, c)

    gamma_over_beta_product = Fraction(1)
    for k in range(1, gap + 1):
        gamma_over_beta_product *= Fraction(2 * b + gap + k, b + k)
    require(
        Fraction(gamma, beta) == gamma_over_beta_product,
        f"gamma/beta product drift at {(b, c)}",
    )
    require(
        gamma_over_beta_product >= 2**gap,
        f"positive middle q-block floor drift at {(b, c)}",
    )

    gamma_over_alpha_product = Fraction(1)
    for t in range(1, gap + 1):
        gamma_over_alpha_product *= Fraction(4 * (b + t) - 2, b + t)
    require(
        Fraction(gamma, alpha) == gamma_over_alpha_product,
        f"gamma/alpha product drift at {(b, c)}",
    )
    require(
        gamma_over_alpha_product >= 3**gap > 2 * gap - 1,
        f"centered q-mean floor drift at {(b, c)}",
    )

    first_weight = alpha - 2 * beta + gamma
    middle_weight = gamma - 2 * beta
    require(first_weight > 0, f"first q block lost positivity at {(b, c)}")
    require(middle_weight >= 0, f"middle q block lost positivity at {(b, c)}")
    require(
        (middle_weight == 0) == (gap == 1),
        f"middle q equality boundary drift at {(b, c)}",
    )

    q = tuple(
        first_weight
        if j < 2 * b
        else middle_weight
        if j < b + c
        else gamma
        for j in range(2 * c)
    )
    require(q[0] > 0 and q[1] > 0, f"q variance support drift at {(b, c)}")

    centered_direct = sum((j - (c - 1)) * q[j] for j in range(2 * c))
    centered_closed = (
        c * gamma
        + (gap - 1) * (b + c) * beta
        - (2 * gap - 1) * b * alpha
    )
    require(
        centered_direct == centered_closed,
        f"centered q-mean identity drift at {(b, c)}",
    )
    require(centered_closed > 0, f"centered q-mean sign drift at {(b, c)}")

    a_values = []
    b_values = []
    r_values = []
    previous_mu_q = None
    previous_mu_a = None
    previous_variance_q = None
    previous_r = None
    tilt_checks = 0

    for n in range(MAX_C + 1):
        for m in (b, c, 2 * b, b + c, 2 * c):
            require(
                hockey_k(n, m) == sum(hockey_h(n, j) for j in range(m)),
                f"hockey-stick identity drift at {(n, m)}",
            )

        a_n = hockey_k(n, c) - hockey_k(n, b)
        b_n = (
            alpha * hockey_k(n, 2 * b)
            - 2 * beta * hockey_k(n, b + c)
            + gamma * hockey_k(n, 2 * c)
        )
        require(
            a_n == sum(hockey_h(n, j) for j in range(b, c)),
            f"A hockey expansion drift at {(n, b, c)}",
        )
        require(
            b_n == sum(q[j] * hockey_h(n, j) for j in range(2 * c)),
            f"B hockey expansion drift at {(n, b, c)}",
        )
        require(a_n > 0 and b_n > 0, f"positive transformed mass drift at {(n, b, c)}")

        q_first = sum(j * q[j] * hockey_h(n, j) for j in range(2 * c))
        q_second = sum(j * j * q[j] * hockey_h(n, j) for j in range(2 * c))
        a_first = sum(j * hockey_h(n, j) for j in range(b, c))
        mu_q = Fraction(q_first, b_n)
        mu_a = Fraction(a_first, a_n)
        variance_q = Fraction(q_second, b_n) - mu_q * mu_q
        require(mu_q > c - 1, f"q tilted-mean separation drift at {(n, b, c)}")
        require(mu_a <= c - 1, f"interval mean support drift at {(n, b, c)}")
        require(variance_q > 0, f"q tilted variance drift at {(n, b, c)}")

        r_n = Fraction(b_n, a_n)
        if n:
            require(
                mu_q - previous_mu_q
                == previous_variance_q / (n + previous_mu_q),
                f"q mean-tilt recursion drift at {(n, b, c)}",
            )
            require(
                r_n / previous_r
                == (n + previous_mu_q) / (n + previous_mu_a),
                f"R likelihood-ratio recursion drift at {(n, b, c)}",
            )
            require(r_n > previous_r, f"R monotonicity drift at {(n, b, c)}")
            tilt_checks += 1

        a_values.append(a_n)
        b_values.append(b_n)
        r_values.append(r_n)
        previous_mu_q = mu_q
        previous_mu_a = mu_a
        previous_variance_q = variance_q
        previous_r = r_n

    return tuple(a_values), tuple(b_values), tuple(r_values), tilt_checks


def verify_wronskian(a, b, c):
    p = b - a
    r = c - b
    a_scale = factorial_ratio(b, a)
    b_scale = factorial_ratio(c, b)

    f_a = normalized_monomial(a)
    x = {p: Fraction(1, a_scale)}
    y = {r: Fraction(1, b_scale)}
    one = {0: Fraction(1)}
    u = multiply_sparse(f_a, add_sparse(x, scale_sparse(one, -1)))
    v = multiply_sparse(
        multiply_sparse(f_a, x),
        add_sparse(y, scale_sparse(one, -1)),
    )
    actual = add_sparse(
        multiply_sparse(u, derivative_sparse(v)),
        scale_sparse(multiply_sparse(derivative_sparse(u), v), -1),
    )

    prefactor = {
        2 * a + p - 1: Fraction(1, factorial(a) ** 2 * a_scale)
    }
    xy = multiply_sparse(x, y)
    bracket = add_sparse(
        {0: Fraction(p)},
        scale_sparse(xy, r),
        scale_sparse(y, -(p + r)),
    )
    expected = multiply_sparse(prefactor, bracket)
    require(actual == expected, f"Wronskian identity drift at {(a, b, c)}")

    # B^(1/r)>A^(1/p), written without radicals.
    require(
        b_scale**p > a_scale**r,
        f"root-mean orientation drift at {(a, b, c)}",
    )


def verify_t111_ratio(a, b, t111):
    r = b - a
    f_values = []
    for k in range(4):
        numerator = factorial(3 * a + k * r)
        denominator = factorial(a) ** 3 * factorial_ratio(b, a) ** k
        f_values.append(Fraction(numerator, denominator))
    reconstructed = f_values[3] - 3 * f_values[2] + 3 * f_values[1] - f_values[0]
    require(reconstructed == t111, f"t111 reconstruction drift at {(a, b)}")
    ratios = tuple(f_values[k + 1] / f_values[k] for k in range(3))
    products = tuple(
        Fraction(
            factorial(3 * a + (k + 1) * r),
            factorial(3 * a + k * r) * factorial_ratio(b, a),
        )
        for k in range(3)
    )
    require(ratios == products, f"ratio product drift at {(a, b)}")
    if a == 0:
        require(ratios[0] == 1, f"a=0 rho0 drift at {(a, b)}")
        require(ratios[1] >= 2 and ratios[2] >= 3, "a=0 ratio floor drift")
    else:
        require(
            ratios[0] > 1 and ratios[1] > 2**r and ratios[2] >= 3**r,
            f"a>0 ratio floor drift at {(a, b)}",
        )
    require(t111 > 0, f"strict cubic orientation failed at {(a, b)}")


def main():
    support_count = 0
    atomic_count = 0
    atomic_positive = 0
    atomic_zero = []
    d_total_zero = []
    i1_negative = 0
    i2_negative = 0
    divisibility_identity_checks = 0
    telescoping_checks = 0
    wronskian_checks = 0
    ratio_checks = 0
    tilted_pair_checks = 0
    tilted_row_checks = 0
    atomic_factorization_checks = 0

    smallest_positive_atomic = None
    smallest_positive_atomic_index = None

    for c in range(2, MAX_C + 1):
        for b in range(1, c):
            a_values, b_values, r_values, tilt_checks = verify_tilted_ratio_certificate(
                b, c
            )
            tilted_pair_checks += 1
            tilted_row_checks += tilt_checks
            secant_a = sum(a_values[n] for n in range(b - 1, c - 1))
            secant_b = sum(b_values[n] for n in range(b - 1, c - 1))
            require(secant_a > 0 and secant_b > 0, "long secant mass drift")

            for i in range(b):
                value = atomic_d(i, b, c)
                factorized = 6 * (
                    a_values[i] * secant_b - b_values[i] * secant_a
                )
                require(
                    value == factorized,
                    f"atomic weighted-secant factorization drift at {(i, b, c)}",
                )
                require(
                    (value == 0) == (i == b - 1 and c == b + 1),
                    f"atomic universal equality boundary drift at {(i, b, c)}",
                )
                require(
                    all(r_values[n + 1] > r_values[n] for n in range(c - 2)),
                    f"needed R prefix monotonicity drift at {(i, b, c)}",
                )
                atomic_factorization_checks += 1
                require(value >= 0, f"atomic orientation failed at {(i, b, c)}")
                atomic_count += 1
                if value == 0:
                    atomic_zero.append((i, b, c))
                else:
                    atomic_positive += 1
                    if smallest_positive_atomic is None or value < smallest_positive_atomic:
                        smallest_positive_atomic = value
                        smallest_positive_atomic_index = (i, b, c)

            for a in range(b):
                (
                    _u,
                    _v,
                    g11,
                    g12,
                    g22,
                    t111,
                    t112,
                    t122,
                    t222,
                    i1,
                    i2,
                    d_total,
                ) = invariants(a, b, c)
                require(
                    g11 > 0 and g22 > 0 and g11 * g22 - g12 * g12 > 0,
                    f"Gram positivity drift at {(a, b, c)}",
                )
                verify_t111_ratio(a, b, t111)
                ratio_checks += 1
                verify_wronskian(a, b, c)
                wronskian_checks += 1

                require(
                    i2 == -g11 * d_total - t111 * g22 * g22,
                    f"I2/D identity drift at {(a, b, c)}",
                )
                divisibility_identity_checks += 1

                atomic_sum = sum(atomic_d(i, b, c) for i in range(a, b))
                require(
                    d_total == atomic_sum,
                    f"atomic telescoping drift at {(a, b, c)}",
                )
                telescoping_checks += 1
                require(d_total >= 0, f"finite D total failed at {(a, b, c)}")
                require(i2 < 0, f"finite common-null obstruction failed at {(a, b, c)}")
                if i1 < 0:
                    i1_negative += 1
                if i2 < 0:
                    i2_negative += 1
                if d_total == 0:
                    d_total_zero.append((a, b, c))
                support_count += 1

    expected_count = (MAX_C + 1) * MAX_C * (MAX_C - 1) // 6
    expected_zero = tuple((i, i + 1, i + 2) for i in range(MAX_C - 1))
    require(support_count == expected_count, "support triple count drift")
    require(atomic_count == expected_count, "atomic triple count drift")
    require(tuple(atomic_zero) == expected_zero, "atomic equality boundary drift")
    require(tuple(d_total_zero) == expected_zero, "total equality boundary drift")
    require(i1_negative == support_count, "finite I1 sign drift")
    require(i2_negative == support_count, "finite I2 sign drift")

    print("THM-2824 ARBITRARY THREE-SLOT FACTORIAL ATOMIC ORIENTATION")
    print(f"scan_max_c={MAX_C}; support_triples={support_count}; atomic_triples={atomic_count}")
    print(
        "gram_positive=all; t111_positive=all; "
        f"ratio_proof_checks={ratio_checks}; wronskian_checks={wronskian_checks}"
    )
    print(
        f"tilted_pair_certificates={tilted_pair_checks}; "
        f"tilted_row_recursions={tilted_row_checks}; "
        f"atomic_secant_factorizations={atomic_factorization_checks}"
    )
    print(
        f"divisibility_identity_checks={divisibility_identity_checks}; "
        f"atomic_telescoping_checks={telescoping_checks}"
    )
    print(
        f"atomic_nonnegative={atomic_count}/{atomic_count}; "
        f"atomic_positive={atomic_positive}; atomic_zero={len(atomic_zero)}"
    )
    print(
        "atomic_zero_boundary=(i,b,c)=(j,j+1,j+2); "
        f"j=0..{MAX_C - 2}"
    )
    print(
        f"smallest_positive_atomic={smallest_positive_atomic}; "
        f"index={smallest_positive_atomic_index}"
    )
    print(
        f"I1_negative={i1_negative}/{support_count}; "
        f"I2_negative={i2_negative}/{support_count}; "
        f"D_total_zero={len(d_total_zero)}"
    )
    print(
        "scope=finite exact audit of the universal hockey-stick/tilted-mean "
        "proof; arbitrary three-slot factorial detection and the lowest-slot-"
        "positive two-charge Gaussian moment-six consequence are proof-complete "
        "candidates"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
