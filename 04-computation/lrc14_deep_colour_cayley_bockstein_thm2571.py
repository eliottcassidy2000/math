#!/usr/bin/env python3
"""Exact companion for THM-2571.

Only Python integers and Fraction arithmetic are used.  Cyclotomic elements
are represented in the power basis 1,zeta,...,zeta^11 of Q(zeta_13).
"""

from fractions import Fraction
import hashlib
import importlib.util
from itertools import combinations, product
from math import gcd
from pathlib import Path


P = 13
DIM = P - 1
checks = 0


def check(condition: bool, message: str = "exact check failed") -> None:
    global checks
    checks += 1
    if not condition:
        raise AssertionError(message)


def lcm(a: int, b: int) -> int:
    return abs(a * b) // gcd(a, b) if a and b else 0


def zero_matrix(n, m):
    return [[Fraction(0) for _ in range(m)] for _ in range(n)]


def matmul(a, b):
    out = zero_matrix(len(a), len(b[0]))
    for i in range(len(a)):
        for k in range(len(b)):
            if a[i][k] == 0:
                continue
            for j in range(len(b[0])):
                out[i][j] += a[i][k] * b[k][j]
    return out


def matvec(a, x):
    return [sum(a[i][j] * x[j] for j in range(len(x))) for i in range(len(a))]


def determinant(a):
    a = [row[:] for row in a]
    det = Fraction(1)
    for col in range(len(a)):
        pivot = next((r for r in range(col, len(a)) if a[r][col]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
            det = -det
        p = a[col][col]
        det *= p
        for j in range(col, len(a)):
            a[col][j] /= p
        for r in range(col + 1, len(a)):
            q = a[r][col]
            if q:
                for j in range(col, len(a)):
                    a[r][j] -= q * a[col][j]
    return det


def rank_mod(a, prime):
    a = [[int(x) % prime for x in row] for row in a]
    rank = 0
    for col in range(len(a[0])):
        pivot = next((r for r in range(rank, len(a)) if a[r][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        inv = pow(a[rank][col], -1, prime)
        a[rank] = [(inv * x) % prime for x in a[rank]]
        for r in range(len(a)):
            if r == rank:
                continue
            q = a[r][col]
            if q:
                a[r] = [(x - q * y) % prime for x, y in zip(a[r], a[rank])]
        rank += 1
    return rank


def shift_matrix(tau=1):
    out = zero_matrix(P, P)
    for v in range(P):
        out[v][(v + tau) % P] = 1
    return out


def cayley_matrices(tau=1):
    c = zero_matrix(P, P)
    b = zero_matrix(P, P)
    for d in range(1, P):
        sign = 1 if d % 2 else -1
        for v in range(P):
            c[v][(v + d * tau) % P] += sign
            b[v][(v + d * tau) % P] += Fraction(2 * d - P, P)
    return c, b


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
    coeff = [Fraction(0) for _ in range(P)]
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            coeff[(i + j) % P] += x * y
    top = coeff[P - 1]
    return tuple(coeff[j] - top for j in range(DIM))


def field_conj(a):
    out = (Fraction(0),) * DIM
    for j, value in enumerate(a):
        out = field_add(out, field_scale(value, zeta_power(-j)))
    return out


def field_is_zero(a):
    return not any(a)


def apply_scalar_matrix_to_field(matrix, vector):
    out = []
    for i in range(len(matrix)):
        value = (Fraction(0),) * DIM
        for j, scalar in enumerate(matrix[i]):
            value = field_add(value, field_scale(scalar, vector[j]))
        out.append(value)
    return out


def field_denominator(vector):
    ans = 1
    for value in vector:
        for coefficient in value:
            ans = lcm(ans, coefficient.denominator)
    return ans


def cyc_power_mod(prime, exponent):
    exponent %= prime
    if exponent < prime - 1:
        return tuple(int(j == exponent) for j in range(prime - 1))
    return tuple(-1 for _ in range(prime - 1))


def cyc_mul_mod(prime, a, b, modulus=13):
    coeff = [0] * prime
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            coeff[(i + j) % prime] += x * y
    top = coeff[prime - 1]
    return tuple((coeff[j] - top) % modulus for j in range(prime - 1))


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


def canonical_carry_bockstein():
    """Reconstruct the canonical THM-2550(B) old/future diagonal carrier."""
    module_path = Path(__file__).with_name(
        "lrc14_replica_dichotomy_typed_row_opus_20260727.py"
    )
    spec = importlib.util.spec_from_file_location("thm2550b_exact", module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    t_den = module.T
    clock = P**6
    check(t_den == 297836897838480)
    check(clock == 4826809)
    check(clock * t_den == 1437601819018855810320)

    word = module.build_word_Ta()
    word_phase = []
    word_digit = [[None] * P for _ in range(7)]
    for ell in range(7):
        q_ell = module.subtract_comb(
            word, module.C1, 182, 26 * ell - 13, 26 * ell + 13
        )
        word_phase.append(q_ell)
        for h in range(P):
            word_digit[ell][h] = intersect_interval(
                q_ell, h * t_den // P, (h + 1) * t_den // P
            )

    x_fraction = [[[Fraction(0) for _ in range(P)] for _ in range(P)] for _ in range(7)]
    for ell in range(7):
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

    raw = [[[0 for _ in range(P)] for _ in range(P)] for _ in range(7)]
    nonzero = []
    positive = 0
    for ell in range(7):
        for s in range(P):
            for r in range(P):
                scaled = x_fraction[ell][s][r] * clock * t_den
                check(scaled.denominator == 1)
                raw[ell][s][r] = scaled.numerator
                if scaled:
                    positive += 1
                    nonzero.append(abs(scaled.numerator))
                    check(r != 0)
                else:
                    check(r == 0)
    raw_gcd = 0
    for value in nonzero:
        raw_gcd = gcd(raw_gcd, value)
    check(positive == 1092)
    check(raw_gcd == P)

    def v13(value):
        valuation = 0
        while value and value % P == 0:
            value //= P
            valuation += 1
        return valuation

    valuations = [v13(value) for value in nonzero]
    check(min(valuations) == 1)
    check(max(valuations) == 3)
    primitive = [
        [[raw[ell][s][r] // P for r in range(P)] for s in range(P)]
        for ell in range(7)
    ]
    primitive_denominator = clock * t_den // P
    check(primitive_denominator == 110584755309142754640)
    check(sum(primitive[0][0]) == 36582596186548944)

    serialization = ";".join(
        ",".join(str(primitive[ell][s][r]) for r in range(P))
        for ell in range(7)
        for s in range(P)
    )
    digest = hashlib.sha256(serialization.encode("ascii")).hexdigest()
    check(digest == "a66ba96d31a33354468392b1dabc19865e6e925158efdab059fad9a98d4390f4")

    beta13 = []
    for r in range(P):
        value = [0] * 12
        for m in range(P):
            power = cyc_power_mod(P, m * r)
            for j in range(12):
                value[j] = (value[j] + m * power[j]) % P
        beta13.append(tuple(value))

    # Closed socle certificate.  In O_13/13 the first-moment vector is
    # r^(-1) Omega, while every target modulation fixes Omega.  Contracting
    # the primitive carrier over (s,r) therefore leaves one septimal
    # polynomial Y, independent of the target character.
    omega = tuple(range(1, P))
    for r in range(1, P):
        inverse_r = pow(r, -1, P)
        check(beta13[r] == tuple((inverse_r * value) % P for value in omega))
    for phase in range(P):
        check(cyc_mul_mod(P, cyc_power_mod(P, phase), omega) == omega)

    y = []
    for ell in range(7):
        value = 0
        for s in range(P):
            for r in range(1, P):
                value += primitive[ell][s][r] * pow(r, -1, P)
        y.append(value % P)
    check(y == [0, 6, 5, 1, 12, 8, 7])

    # Reduce Y modulo Phi_7 and exhibit its inverse.  This is the structural
    # reason all six nonzero owner characters survive, rather than a
    # seventy-eight-case nonvanishing census.
    y_reduced = tuple((y[j] - y[6]) % P for j in range(6))
    y_inverse = (0, 0, 3, 5, 8, 10)
    check(cyc_mul_mod(7, y_reduced, y_inverse) == (1, 0, 0, 0, 0, 0))

    pre_zero = []
    for ell in range(7):
        for s in range(P):
            value = [0] * 12
            for r in range(P):
                for j in range(12):
                    value[j] = (
                        value[j] + primitive[ell][s][r] * beta13[r][j]
                    ) % P
            if not any(value):
                pre_zero.append((ell, s))
    check(pre_zero == [(0, 0), (3, 0), (4, 0)])

    beta_nonzero = 0
    samples = {}
    for kappa in range(1, 7):
        owner_factor = [0] * 6
        for ell in range(7):
            z7 = cyc_power_mod(7, kappa * ell)
            for i in range(6):
                owner_factor[i] = (owner_factor[i] + y[ell] * z7[i]) % P
        check(all(owner_factor))
        for target in range(P):
            beta = [[0] * 12 for _ in range(6)]
            for ell in range(7):
                z7 = cyc_power_mod(7, kappa * ell)
                for s in range(P):
                    phase13 = cyc_power_mod(P, target * s)
                    for r in range(1, P):
                        phase_beta = cyc_mul_mod(P, phase13, beta13[r])
                        coefficient = primitive[ell][s][r] % P
                        if not coefficient:
                            continue
                        for i in range(6):
                            if not z7[i]:
                                continue
                            for j in range(12):
                                beta[i][j] = (
                                    beta[i][j]
                                    + coefficient * z7[i] * phase_beta[j]
                                ) % P
            support = sum(int(value != 0) for row in beta for value in row)
            expected = [
                [(owner_factor[i] * omega[j]) % P for j in range(12)]
                for i in range(6)
            ]
            check(beta == expected)
            check(support == 72)
            beta_nonzero += 1
            if (kappa, target) in {(1, 0), (1, 1), (6, 12)}:
                samples[(kappa, target)] = beta

    check(beta_nonzero == 78)
    for key in samples:
        check(sum(int(value != 0) for row in samples[key] for value in row) == 72)
    check(samples[(1, 0)][0] == [6, 12, 5, 11, 4, 10, 3, 9, 2, 8, 1, 7])
    check(samples[(1, 1)][0] == [6, 12, 5, 11, 4, 10, 3, 9, 2, 8, 1, 7])
    check(samples[(6, 12)][0] == [7, 1, 8, 2, 9, 3, 10, 4, 11, 5, 12, 6])

    return {
        "positive": positive,
        "raw_gcd": raw_gcd,
        "denominator": primitive_denominator,
        "digest": digest,
        "pre_nonzero": 91 - len(pre_zero),
        "beta_nonzero": beta_nonzero,
        "y": y,
    }


def main() -> None:
    c, b = cayley_matrices(1)
    cb = matmul(c, b)
    for i in range(P):
        for j in range(P):
            check(cb[i][j] == Fraction(int(i == j), 1) - Fraction(1, P))

    # Restrict C to the augmentation lattice in the basis e_j-e_12.
    restriction = zero_matrix(DIM, DIM)
    for j in range(DIM):
        x = [Fraction(0)] * P
        x[j] = 1
        x[P - 1] = -1
        y = matvec(c, x)
        check(sum(y) == 0)
        check(sum(r * y[r] for r in range(P)) % P == 0)
        for i in range(DIM):
            restriction[i][j] = y[i]
    det = determinant(restriction)
    check(abs(det) == P)
    check(rank_mod(restriction, P) == DIM - 1)

    # Exhaust a 4,095-profile integral slice.  B gives an integral primitive
    # exactly when the cyclic first moment vanishes mod 13.
    integral_profiles = 0
    obstructed_profiles = 0
    for tail in product(range(2), repeat=DIM):
        if not any(tail):
            continue
        a = [-sum(tail)] + list(tail)
        moment = sum(m * a[m] for m in range(P)) % P
        primitive = matvec(b, a)
        integral = all(x.denominator == 1 for x in primitive)
        check(integral == (moment == 0))
        check(matvec(c, primitive) == [Fraction(x) for x in a])
        integral_profiles += int(integral)
        obstructed_profiles += int(not integral)

    # Minimal lawful singleton H=delta_(1,1,0), at q1=q2=0.
    singleton = [zeta_power(m) for m in range(P)]
    singleton_sum = (Fraction(0),) * DIM
    for value in singleton:
        singleton_sum = field_add(singleton_sum, value)
    check(field_is_zero(singleton_sum))
    singleton_beta = (Fraction(0),) * DIM
    for m, value in enumerate(singleton):
        singleton_beta = field_add(singleton_beta, field_scale(m, value))
    check(not field_is_zero(singleton_beta))
    check(any(int(x) % P for x in singleton_beta))
    zeta_minus_one = field_add(zeta_power(1), field_scale(-1, zeta_power(0)))
    check(field_mul(zeta_minus_one, singleton_beta) == field_scale(P, zeta_power(0)))
    singleton_primitive = apply_scalar_matrix_to_field(b, singleton)
    check(field_denominator(singleton_primitive) == P)
    check(apply_scalar_matrix_to_field(c, singleton_primitive) == singleton)

    # Uniform primitive hostile: a_0=144, a_m=-12.  It has positive norm
    # but vanishing first Bockstein and an integral Cayley primitive.
    uniform = [144] + [-12] * DIM
    check(sum(uniform) == 0)
    uniform_beta = sum(m * uniform[m] for m in range(P))
    check(uniform_beta % P == 0)
    uniform_primitive = matvec(b, uniform)
    check(all(x.denominator == 1 for x in uniform_primitive))
    check(matvec(c, uniform_primitive) == [Fraction(x) for x in uniform])

    # Physical displacement formula on every diagonal-free singleton and a
    # deterministic target phase.  Multiplication by zeta^u-1 recovers 13.
    displacement_controls = 0
    for r in range(P):
        for s in range(P):
            for t in range(P):
                if r == t:
                    continue
                u = (r - t) % P
                q1 = (r + 2 * s + t) % P
                q2 = (2 * r + s + 3 * t) % P
                phase = zeta_power((q1 * s + q2 * t) % P)
                beta_u = (Fraction(0),) * DIM
                for m in range(P):
                    beta_u = field_add(beta_u, field_scale(m, zeta_power(m * u)))
                lhs = field_mul(
                    field_add(zeta_power(u), field_scale(-1, zeta_power(0))),
                    beta_u,
                )
                check(lhs == field_scale(P, zeta_power(0)))
                phased = field_mul(phase, beta_u)
                direct = (Fraction(0),) * DIM
                for m in range(P):
                    direct = field_add(
                        direct,
                        field_scale(m, field_mul(phase, zeta_power(m * u))),
                    )
                check(direct == phased)
                check(not field_is_zero(phased))
                displacement_controls += 1

    # The deep x gain/replica tensor is plaquette-flat.  Use the singleton
    # Galois orbit and all 18 rows from the three duty classes x six replicas.
    row_scales = [d for d in (229692, 440232, 440244) for _ in range(6)]
    plaquettes = 0
    for m, n in combinations(range(1, P), 2):
        for i, j in combinations(range(len(row_scales)), 2):
            left = field_mul(
                field_scale(-row_scales[i], singleton[m]),
                field_scale(-row_scales[j], singleton[n]),
            )
            right = field_mul(
                field_scale(-row_scales[j], singleton[m]),
                field_scale(-row_scales[i], singleton[n]),
            )
            check(left == right)
            plaquettes += 1

    # Galois conjugacy and positive norm pairing on every nonempty Boolean G.
    galois_profiles = 0
    for tail in product(range(2), repeat=DIM):
        if not any(tail):
            continue
        anchors = []
        for m in range(1, P):
            value = (Fraction(0),) * DIM
            for u, weight in enumerate((0,) + tail):
                value = field_add(value, field_scale(weight, zeta_power(m * u)))
            anchors.append(value)
            check(not field_is_zero(value))
        for m in range(1, 7):
            check(anchors[P - 1 - m] == field_conj(anchors[m - 1]))
        galois_profiles += 1

    carry = canonical_carry_bockstein()

    print("THM-2571 deep-colour Cayley filling and Bockstein")
    print(f"Cayley augmentation determinant {abs(det)}, mod13 rank {DIM - 1}")
    print(
        "integral slice profiles "
        f"{integral_profiles + obstructed_profiles}, fillable {integral_profiles}, "
        f"obstructed {obstructed_profiles}"
    )
    print(f"singleton primitive denominator {field_denominator(singleton_primitive)}")
    print(f"physical displacement controls {displacement_controls}")
    print(f"rank-one plaquette controls {plaquettes}")
    print(f"nonzero Galois-conjugacy profiles {galois_profiles}")
    print("uniform hostile beta 0, singleton hostile beta nonzero")
    print(
        "canonical carry cells "
        f"{carry['positive']}/1183, primitive denominator {carry['denominator']}"
    )
    print(
        "canonical carry Bockstein "
        f"{carry['beta_nonzero']}/78, precontraction {carry['pre_nonzero']}/91"
    )
    print("canonical carry Y " + ",".join(str(value) for value in carry["y"]))
    print(f"canonical carry digest {carry['digest']}")
    print(f"explicit checks {checks}")


if __name__ == "__main__":
    main()
