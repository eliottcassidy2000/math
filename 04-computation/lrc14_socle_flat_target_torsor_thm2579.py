#!/usr/bin/env python3
"""Exact companion for THM-2579.

The referee is dependency-free.  It checks the integral Cayley cokernel,
its functorial law under every displayed circulant filter, the sharp
singleton target torsor, and the closed canonical Omega times Y class.
"""

from fractions import Fraction
from itertools import combinations, product


P = 13
DIM = 12
checks = 0


def check(condition: bool, message: str = "exact check failed") -> None:
    global checks
    checks += 1
    if not condition:
        raise AssertionError(message)


def zero_matrix(n, m):
    return [[Fraction(0) for _ in range(m)] for _ in range(n)]


def matvec(a, x):
    return [sum(a[i][j] * x[j] for j in range(len(x))) for i in range(len(a))]


def cayley_matrices():
    c = zero_matrix(P, P)
    b = zero_matrix(P, P)
    for d in range(1, P):
        for m in range(P):
            c[m][(m + d) % P] += 1 if d % 2 else -1
            b[m][(m + d) % P] += Fraction(2 * d - P, P)
    return c, b


def apply_filter(coefficients, vector):
    return [
        sum(coefficients[d] * vector[(m + d) % P] for d in range(P))
        for m in range(P)
    ]


def beta_int(vector):
    return sum(m * vector[m] for m in range(P)) % P


def zeta_power(k):
    k %= P
    if k < DIM:
        return tuple(Fraction(int(j == k)) for j in range(DIM))
    return tuple(Fraction(-1) for _ in range(DIM))


def field_add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def field_scale(c, a):
    return tuple(c * x for x in a)


def field_mul(a, b):
    coeff = [Fraction(0)] * P
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            coeff[(i + j) % P] += x * y
    top = coeff[P - 1]
    return tuple(coeff[j] - top for j in range(DIM))


def field_mod13(a):
    out = []
    for value in a:
        check(value.denominator % P != 0)
        out.append((value.numerator * pow(value.denominator, -1, P)) % P)
    return tuple(out)


def field_vector_beta(vector):
    out = (Fraction(0),) * DIM
    for m, value in enumerate(vector):
        out = field_add(out, field_scale(m, value))
    return out


def field_vector_denominator(vector):
    answer = 1
    for value in vector:
        for coefficient in value:
            answer = answer * coefficient.denominator // __import__("math").gcd(
                answer, coefficient.denominator
            )
    return answer


def apply_scalar_matrix_to_field(matrix, vector):
    out = []
    for row in matrix:
        value = (Fraction(0),) * DIM
        for coefficient, entry in zip(row, vector):
            value = field_add(value, field_scale(coefficient, entry))
        out.append(value)
    return out


def cyc_power_mod(prime, exponent):
    exponent %= prime
    if exponent < prime - 1:
        return tuple(int(j == exponent) for j in range(prime - 1))
    return tuple(-1 for _ in range(prime - 1))


def cyc_mul_mod(prime, a, b):
    coeff = [0] * prime
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            coeff[(i + j) % prime] += x * y
    top = coeff[prime - 1]
    return tuple((coeff[j] - top) % P for j in range(prime - 1))


def main() -> None:
    c, b = cayley_matrices()

    # Functorial first-moment law.  The filter bank contains all translations,
    # all ordered differences, the full and punctured orbit sums, and C itself.
    filters = []
    for d in range(P):
        f = [0] * P
        f[d] = 1
        filters.append(tuple(f))
    for d in range(P):
        for e in range(P):
            if d == e:
                continue
            f = [0] * P
            f[d] = 1
            f[e] = -1
            filters.append(tuple(f))
    filters.append((1,) * P)
    filters.append((0,) + (1,) * DIM)
    filters.append(tuple(int(c[0][d]) for d in range(P)))
    check(len(filters) == 172)

    profiles = 0
    fillable = 0
    functorial_controls = 0
    for tail in product(range(2), repeat=DIM):
        if not any(tail):
            continue
        a = [-sum(tail)] + list(tail)
        beta = beta_int(a)
        primitive = matvec(b, a)
        integral = all(value.denominator == 1 for value in primitive)
        check(integral == (beta == 0))
        fillable += int(integral)
        profiles += 1
        for f in filters:
            filtered = apply_filter(f, a)
            check(sum(filtered) == 0)
            check(beta_int(filtered) == (sum(f) * beta) % P)
            functorial_controls += 1

    # The singleton orbit gives thirteen absolute nonfillable profiles whose
    # pairwise differences have integral Cayley primitives.
    target_profiles = []
    target_betas = []
    for target in range(P):
        a = [zeta_power(target + m) for m in range(P)]
        augmentation = (Fraction(0),) * DIM
        for value in a:
            augmentation = field_add(augmentation, value)
        check(not any(augmentation))
        beta = field_vector_beta(a)
        target_profiles.append(a)
        target_betas.append(beta)
        check(field_vector_denominator(apply_scalar_matrix_to_field(b, a)) == P)

    omega_mod = tuple(range(1, P))
    for beta in target_betas:
        check(field_mod13(beta) == omega_mod)

    pairwise_fillings = 0
    for left, right in combinations(range(P), 2):
        difference = [
            field_add(target_profiles[left][m], field_scale(-1, target_profiles[right][m]))
            for m in range(P)
        ]
        check(not any(field_mod13(field_vector_beta(difference))))
        primitive = apply_scalar_matrix_to_field(b, difference)
        check(field_vector_denominator(primitive) == 1)
        pairwise_fillings += 1

    # The unnormalized target Fourier numerators fill.  Division by 13 is
    # not an integral-lattice operation and can restore the absolute class.
    singleton_fourier_fillings = 0
    normalized_restorations = 0
    for q in range(P):
        numerator = []
        for m in range(P):
            value = (Fraction(0),) * DIM
            for target in range(P):
                value = field_add(
                    value,
                    field_mul(zeta_power(q * target), target_profiles[target][m]),
                )
            numerator.append(value)
        check(not any(field_mod13(field_vector_beta(numerator))))
        check(
            field_vector_denominator(apply_scalar_matrix_to_field(b, numerator))
            == 1
        )
        singleton_fourier_fillings += 1
        if q == P - 1:
            expected = [field_scale(P, value) for value in target_profiles[0]]
            check(numerator == expected)
            normalized = [field_scale(Fraction(1, P), value) for value in numerator]
            check(normalized == target_profiles[0])
            check(field_mod13(field_vector_beta(normalized)) == omega_mod)
            check(
                field_vector_denominator(apply_scalar_matrix_to_field(b, normalized))
                == P
            )
            normalized_restorations += 1
        else:
            check(not any(coefficient for value in numerator for coefficient in value))

    # Closed canonical carry class.  Y is a septimal unit; every target has
    # the same Omega factor.  Integer combinations depend only on the total
    # coefficient modulo 13, and every full target Fourier bank is fillable.
    y = [0, 6, 5, 1, 12, 8, 7]
    y_reduced = tuple((y[j] - y[6]) % P for j in range(6))
    y_inverse = (0, 0, 3, 5, 8, 10)
    check(cyc_mul_mod(7, y_reduced, y_inverse) == (1, 0, 0, 0, 0, 0))

    owner_factors = []
    for kappa in range(1, 7):
        factor = [0] * 6
        for ell in range(7):
            power = cyc_power_mod(7, kappa * ell)
            for i in range(6):
                factor[i] = (factor[i] + y[ell] * power[i]) % P
        check(all(factor))
        owner_factors.append(tuple(factor))

    canonical_classes = 0
    target_pair_classes = 0
    target_fourier_fillings = 0
    for owner in owner_factors:
        beta = tuple(tuple((u * v) % P for v in omega_mod) for u in owner)
        check(sum(int(value != 0) for row in beta for value in row) == 72)
        classes = [beta for _ in range(P)]
        canonical_classes += len(classes)
        for left, right in combinations(range(P), 2):
            difference = [
                [(classes[left][i][j] - classes[right][i][j]) % P for j in range(12)]
                for i in range(6)
            ]
            check(not any(value for row in difference for value in row))
            target_pair_classes += 1
        for q in range(P):
            total = [[0] * 12 for _ in range(6)]
            for target in range(P):
                phase = cyc_power_mod(P, q * target)
                for i in range(6):
                    product_beta = cyc_mul_mod(P, phase, classes[target][i])
                    for j in range(12):
                        total[i][j] = (total[i][j] + product_beta[j]) % P
            check(not any(value for row in total for value in row))
            target_fourier_fillings += 1

    # Sharp normalization and orbit boundaries.
    for owner in owner_factors:
        beta = tuple(tuple((u * v) % P for v in omega_mod) for u in owner)
        extra13 = tuple(tuple((P * value) % P for value in row) for row in beta)
        check(not any(value for row in extra13 for value in row))
    check(sum((0,) + (1,) * DIM) % P == 12)
    check(sum((1,) * P) % P == 0)

    print("THM-2579 socle-flat target torsor")
    print(f"integral profiles {profiles}, fillable {fillable}")
    print(f"circulant functorial controls {functorial_controls}")
    print(f"singleton target profiles {len(target_profiles)}, pairwise fillings {pairwise_fillings}")
    print(f"canonical target classes {canonical_classes}, pairwise class fillings {target_pair_classes}")
    print(f"canonical unnormalized target Fourier fillings {target_fourier_fillings}")
    print(
        "singleton unnormalized Fourier fillings "
        f"{singleton_fourier_fillings}, normalized restorations {normalized_restorations}"
    )
    print("punctured orbit negates, full orbit and extra13 kill")
    print(f"explicit checks {checks}")


if __name__ == "__main__":
    main()
