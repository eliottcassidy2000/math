#!/usr/bin/env python3
"""Exact companion for THM-2585.

The canonical carrier is reconstructed from the fixed THM-2550(B) interval
module.  No THM-2571 or THM-2579 output is imported.  All arithmetic is over
integers, Fraction, and exact cyclotomic power bases.
"""

from collections import Counter
from fractions import Fraction
import hashlib
import importlib.util
from math import gcd
from pathlib import Path


P = 13
OWNER_P = 7
checks = 0


def check(condition: bool, message: str = "exact check failed") -> None:
    global checks
    checks += 1
    if not condition:
        raise AssertionError(message)


def intersect_interval(intervals, left, right):
    out = []
    for a, b in intervals:
        lo = max(a, left)
        hi = min(b, right)
        if lo < hi:
            out.append((lo, hi))
    return out


def interval_length(intervals):
    return sum(b - a for a, b in intervals)


def cyc_power_int(prime, exponent):
    """Power of zeta_prime in the basis 1,zeta,...,zeta^(prime-2)."""
    exponent %= prime
    if exponent < prime - 1:
        return tuple(int(j == exponent) for j in range(prime - 1))
    return tuple(-1 for _ in range(prime - 1))


def cyc_power_mod(prime, exponent, modulus=P):
    return tuple(x % modulus for x in cyc_power_int(prime, exponent))


def cyc_mul_mod(prime, a, b, modulus=P):
    coeff = [0] * prime
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            coeff[(i + j) % prime] += x * y
    top = coeff[prime - 1]
    return tuple((coeff[j] - top) % modulus for j in range(prime - 1))


def vec_add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def vec_scale(c, a):
    return tuple(c * x for x in a)


def determinant_mod(matrix, modulus=P):
    a = [[x % modulus for x in row] for row in matrix]
    det = 1
    for col in range(len(a)):
        pivot = next((r for r in range(col, len(a)) if a[r][col]), None)
        if pivot is None:
            return 0
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
            det = -det
        value = a[col][col] % modulus
        det = det * value % modulus
        inverse = pow(value, -1, modulus)
        for j in range(col, len(a)):
            a[col][j] = a[col][j] * inverse % modulus
        for r in range(col + 1, len(a)):
            factor = a[r][col]
            if factor:
                for j in range(col, len(a)):
                    a[r][j] = (a[r][j] - factor * a[col][j]) % modulus
    return det % modulus


def multiplication_determinant_7(poly):
    columns = []
    for j in range(OWNER_P - 1):
        basis = tuple(int(i == j) for i in range(OWNER_P - 1))
        columns.append(cyc_mul_mod(OWNER_P, poly, basis))
    matrix = [[columns[j][i] for j in range(OWNER_P - 1)] for i in range(OWNER_P - 1)]
    return determinant_mod(matrix)


def reconstruct_canonical_carrier():
    """Reconstruct the fixed THM-2550(B)/THM-2309(25) primitive tensor."""
    module_path = Path(__file__).with_name(
        "lrc14_replica_dichotomy_typed_row_opus_20260727.py"
    )
    spec = importlib.util.spec_from_file_location("thm2550b_exact_2585", module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    t_den = module.T
    clock = P**6
    check(t_den == 297836897838480)
    check(clock == 4826809)
    check(clock * t_den == 1437601819018855810320)

    word = module.build_word_Ta()
    word_digit = [[None] * P for _ in range(OWNER_P)]
    for ell in range(OWNER_P):
        q_ell = module.subtract_comb(
            word, module.C1, 182, 26 * ell - 13, 26 * ell + 13
        )
        for h in range(P):
            word_digit[ell][h] = intersect_interval(
                q_ell, h * t_den // P, (h + 1) * t_den // P
            )

    x_fraction = [
        [[Fraction(0) for _ in range(P)] for _ in range(P)]
        for _ in range(OWNER_P)
    ]
    for ell in range(OWNER_P):
        q_prefix = [module.make_prefix(word_digit[ell][h]) for h in range(P)]
        q_length = [interval_length(word_digit[ell][h]) for h in range(P)]
        for s in range(P):
            f_ls = module.build_F(ell, s)
            for r in range(P):
                e_lsr = module.intersect_comb(
                    f_ls, module.C3, 182, 14 * r - 13, 14 * r + 13
                )
                value = Fraction(0)
                for h in range(P):
                    e_h = intersect_interval(
                        e_lsr, h * t_den // P, (h + 1) * t_den // P
                    )
                    if not e_h or not word_digit[ell][h]:
                        continue
                    starts, lens, pref = q_prefix[h]
                    acc_r, acc_p = module.sweep_acc(
                        e_h, clock % t_den, starts, lens, pref
                    )
                    value += module.IR_from_acc(
                        clock,
                        interval_length(e_h),
                        q_length[h],
                        acc_r,
                        acc_p,
                    )
                x_fraction[ell][s][r] = value

    raw = [[[0 for _ in range(P)] for _ in range(P)] for _ in range(OWNER_P)]
    nonzero = []
    positive = 0
    for ell in range(OWNER_P):
        for s in range(P):
            for r in range(P):
                scaled = x_fraction[ell][s][r] * clock * t_den
                check(scaled.denominator == 1)
                value = scaled.numerator
                raw[ell][s][r] = value
                check((value > 0) == (r != 0))
                if value:
                    positive += 1
                    nonzero.append(abs(value))

    raw_gcd = 0
    for value in nonzero:
        raw_gcd = gcd(raw_gcd, value)
    check(positive == OWNER_P * P * (P - 1))
    check(raw_gcd == P)

    primitive = [
        [[raw[ell][s][r] // P for r in range(P)] for s in range(P)]
        for ell in range(OWNER_P)
    ]
    denominator = clock * t_den // P
    check(denominator == 110584755309142754640)
    serialization = ";".join(
        ",".join(str(primitive[ell][s][r]) for r in range(P))
        for ell in range(OWNER_P)
        for s in range(P)
    )
    digest = hashlib.sha256(serialization.encode("ascii")).hexdigest()
    check(digest == "a66ba96d31a33354468392b1dabc19865e6e925158efdab059fad9a98d4390f4")
    return primitive, positive, raw_gcd, denominator, digest


def first_moment_vectors():
    beta = []
    for r in range(P):
        value = [0] * (P - 1)
        for m in range(P):
            power = cyc_power_mod(P, m * r)
            for j in range(P - 1):
                value[j] = (value[j] + m * power[j]) % P
        beta.append(tuple(value))
    omega = tuple(range(1, P))
    check(beta[0] == (0,) * (P - 1))
    for r in range(1, P):
        inv_r = pow(r, -1, P)
        check(beta[r] == tuple(inv_r * x % P for x in omega))
    for phase in range(P):
        check(cyc_mul_mod(P, cyc_power_mod(P, phase), omega) == omega)
    return beta, omega


def owner_factor(y, kappa):
    factor = [0] * (OWNER_P - 1)
    for ell, coefficient in enumerate(y):
        power = cyc_power_mod(OWNER_P, kappa * ell)
        for i in range(OWNER_P - 1):
            factor[i] = (factor[i] + coefficient * power[i]) % P
    return tuple(factor)


def cayley_apply(vector):
    out = [0] * P
    for d in range(1, P):
        sign = 1 if d % 2 else -1
        for m in range(P):
            out[m] += sign * vector[(m + d) % P]
    return out


def exact_hostiles(omega):
    # Positive signal: one common Fourier carrier at target shift s=1 and
    # deep root r=1.  The normalized q=-1 section returns its nonzero class.
    singleton_nonzero = 0
    for q in range(P):
        divided = []
        for m in range(P):
            numerator = (0,) * (P - 1)
            for b in range(P):
                numerator = vec_add(
                    numerator, cyc_power_int(P, (q + 1) * b + m)
                )
            expected = (
                vec_scale(P, cyc_power_int(P, m))
                if q == P - 1
                else (0,) * (P - 1)
            )
            check(numerator == expected)
            check(all(value % P == 0 for value in numerator))
            divided.append(tuple(value // P for value in numerator))
        beta = [0] * (P - 1)
        for m, value in enumerate(divided):
            for j in range(P - 1):
                beta[j] = (beta[j] + m * value[j]) % P
        check((tuple(beta) == omega) == (q == P - 1))
        singleton_nonzero += int(any(beta))
    check(singleton_nonzero == 1)

    # Same-coset hostile: every profile represents the same nonzero Cayley
    # class, and every unnormalized DFT numerator has zero Bockstein, yet no
    # numerator has coefficientwise content 13.  Thus normalization is not an
    # integral operation on arbitrary independent representatives.
    a = [-1, 1] + [0] * (P - 2)
    boundary = cayley_apply(a)
    check(sum(a) == 0)
    check(sum(boundary) == 0)
    check(sum(m * a[m] for m in range(P)) % P == 1)
    check(sum(m * boundary[m] for m in range(P)) % P == 0)
    content = 0
    for value in boundary:
        content = gcd(content, abs(value))
    check(content == 1)

    profiles = []
    for b in range(P):
        profiles.append([a[m] + (boundary[m] if b == 0 else 0) for m in range(P)])
        check(sum(profiles[-1]) == 0)
        check(sum(m * profiles[-1][m] for m in range(P)) % P == 1)

    fillable_numerators = 0
    integral_normalizations = 0
    for q in range(P):
        numerator = []
        for m in range(P):
            value = (0,) * (P - 1)
            for b in range(P):
                value = vec_add(
                    value,
                    vec_scale(profiles[b][m], cyc_power_int(P, q * b)),
                )
            numerator.append(value)
        augmentation = [sum(numerator[m][j] for m in range(P)) for j in range(P - 1)]
        check(not any(augmentation))
        beta = [0] * (P - 1)
        for m, value in enumerate(numerator):
            for j in range(P - 1):
                beta[j] = (beta[j] + m * value[j]) % P
        check(not any(beta))
        fillable_numerators += 1
        integral = all(value % P == 0 for row in numerator for value in row)
        integral_normalizations += int(integral)
        check(not integral)
    check(fillable_numerators == P)
    check(integral_normalizations == 0)
    return singleton_nonzero, fillable_numerators, integral_normalizations, content


def main() -> None:
    primitive, positive, raw_gcd, denominator, digest = reconstruct_canonical_carrier()
    beta13, omega = first_moment_vectors()

    # Exact inverse-target orthogonality.  On a common Fourier carrier the
    # normalized numerator is therefore a literal target-shift slice.
    orthogonality = 0
    for q in range(P):
        for s in range(P):
            total = (0,) * (P - 1)
            for b in range(P):
                total = vec_add(total, cyc_power_int(P, b * (q + s)))
            expected = (
                vec_scale(P, cyc_power_int(P, 0))
                if (q + s) % P == 0
                else (0,) * (P - 1)
            )
            check(total == expected)
            orthogonality += 1

    yq = []
    for q in range(P):
        s = (-q) % P
        y = []
        for ell in range(OWNER_P):
            value = 0
            for r in range(1, P):
                value += primitive[ell][s][r] * pow(r, -1, P)
            y.append(value % P)
        yq.append(tuple(y))

    expected_yq = [
        (0, 9, 9, 0, 0, 4, 4),
        (5, 9, 7, 11, 11, 4, 9),
        (3, 11, 7, 10, 10, 11, 1),
        (9, 11, 4, 2, 10, 5, 10),
        (11, 5, 4, 10, 12, 11, 8),
        (11, 1, 11, 2, 8, 1, 1),
        (6, 9, 12, 8, 4, 4, 7),
        (7, 6, 9, 9, 5, 1, 4),
        (2, 12, 12, 5, 11, 2, 12),
        (2, 5, 2, 1, 3, 9, 8),
        (4, 3, 8, 3, 11, 9, 2),
        (10, 12, 2, 3, 3, 6, 2),
        (8, 4, 9, 2, 2, 6, 4),
    ]
    check(yq == expected_yq)

    determinants = []
    support_histogram = Counter()
    full_support = []
    beta_profiles = 0
    for q, y in enumerate(yq):
        reduced = tuple((y[j] - y[OWNER_P - 1]) % P for j in range(OWNER_P - 1))
        determinant = multiplication_determinant_7(reduced)
        determinants.append(determinant)
        check(determinant != 0)
        full_q = []
        s = (-q) % P
        for kappa in range(1, OWNER_P):
            factor = owner_factor(y, kappa)
            check(any(factor))

            direct = [[0] * (P - 1) for _ in range(OWNER_P - 1)]
            for ell in range(OWNER_P):
                z7 = cyc_power_mod(OWNER_P, kappa * ell)
                for r in range(1, P):
                    coefficient = primitive[ell][s][r] % P
                    for i in range(OWNER_P - 1):
                        for j in range(P - 1):
                            direct[i][j] = (
                                direct[i][j]
                                + coefficient * z7[i] * beta13[r][j]
                            ) % P
            expected = [
                [factor[i] * omega[j] % P for j in range(P - 1)]
                for i in range(OWNER_P - 1)
            ]
            check(direct == expected)
            support = sum(int(value != 0) for row in direct for value in row)
            check(support in {48, 60, 72})
            support_histogram[support] += 1
            if support == 72:
                full_q.append(kappa)
            beta_profiles += 1
        full_support.append(tuple(full_q))

    check(determinants == [7, 6, 5, 2, 7, 10, 7, 7, 10, 7, 2, 5, 6])
    check(support_histogram == Counter({72: 38, 60: 32, 48: 8}))
    expected_full = [
        (),
        (3, 4),
        (1, 3),
        (2, 3, 4, 6),
        (1, 2, 3, 5, 6),
        (2, 5),
        (1, 2, 3, 6),
        (1, 4, 5, 6),
        (2, 5),
        (1, 2, 4, 5, 6),
        (1, 3, 4, 5),
        (4, 6),
        (3, 4),
    ]
    check(full_support == expected_full)
    check(beta_profiles == 78)

    global_y = tuple(sum(yq[q][ell] for q in range(P)) % P for ell in range(OWNER_P))
    check(global_y == (0, 6, 5, 1, 12, 8, 7))
    for kappa in range(1, OWNER_P):
        check(any(owner_factor(global_y, kappa)))

    # Translation T_a x(ell,s,r)=x(ell,s-a,r) permutes q to q+a.
    translation_checks = 0
    for a in range(P):
        check(len({(q + a) % P for q in range(P)}) == P)
        for q in range(P):
            for ell in range(OWNER_P):
                shifted_s = (-q - a) % P
                value = sum(
                    primitive[ell][shifted_s][r] * pow(r, -1, P)
                    for r in range(1, P)
                ) % P
                check(value == yq[(q + a) % P][ell])
            translation_checks += 1

    # Primitive reduction does not commute with the mod-13 Bockstein:
    # beta(13 D_q)=0, whereas beta(D_q) is nonzero in all 78 profiles.
    for q in range(P):
        for kappa in range(1, OWNER_P):
            factor = owner_factor(yq[q], kappa)
            check(not any((P * value) % P for value in factor))
            check(any(factor))

    singleton_nonzero, hostile_fillable, hostile_integral, hostile_content = exact_hostiles(omega)

    full_text = []
    for q, owners in enumerate(full_support):
        label = "-" if not owners else ",".join(str(x) for x in owners)
        full_text.append(f"q{q}:{label}")

    print("THM-2585 saturated normalized target projector")
    print(
        f"canonical carry cells {positive}/1183, raw gcd {raw_gcd}, "
        f"primitive denominator {denominator}"
    )
    print(f"canonical carry digest {digest}")
    print(f"target orthogonality identities {orthogonality}/169")
    print("slice Y determinants " + ",".join(str(x) for x in determinants))
    print("slice Y sum " + ",".join(str(x) for x in global_y))
    print(
        "slice Bocksteins 78/78, support histogram "
        f"48:{support_histogram[48]},60:{support_histogram[60]},72:{support_histogram[72]}"
    )
    print("full-support owner sets " + ";".join(full_text))
    print(f"translation-permuted section checks {translation_checks}/169")
    print(
        "singleton normalized sections "
        f"nonzero {singleton_nonzero}/13 at q=12"
    )
    print(
        "same-coset hostile unnormalized fillable "
        f"{hostile_fillable}/13, normalized integral {hostile_integral}/13, "
        f"boundary content {hostile_content}"
    )
    print(f"explicit checks {checks}")


if __name__ == "__main__":
    main()
