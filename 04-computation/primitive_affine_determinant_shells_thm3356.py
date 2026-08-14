#!/usr/bin/env python3
"""Exact companion for THM-3356 primitive affine determinant shells.

The exact universe has four parts.

1.  The 22,890 disconnected-low formal affine rows are divided by their
    common affine content and compared with the 14,168 primitive rows.
2.  Every primitive row is checked as a determinant shell, parabolic bouquet,
    Pythagorean Lorentz shell, and fixed U-spine resultant fingerprint.
3.  General unramified sum-of-two-squares root atlases are checked on selected
    primitive anchors, shell charges, and split-prime-power moduli.
4.  The coherent rational Gaussian-toggle integrality gate is checked for all
    split primes below 100 on a bounded family of kernel anchors and shells.

The finite sweeps support the displayed identities and sharp finite ledgers;
the universal theorem rests on the on-page algebra.
"""

from collections import Counter
from fractions import Fraction
from hashlib import sha256
from math import gcd, isqrt


DMAX = 8
SLOPE_EXCESS = 7
CMAX = 46
LRC_DMAX = Fraction(186_636_088_362, 11_773_143_757_375)
EXPECTED_PRIMITIVE = 14_168
EXPECTED_RANK_DISTRIBUTION = ((0, 25), (1, 229), (2, 4214),
                              (3, 7168), (4, 2371), (5, 161))
EXPECTED_FIBRE_DISTRIBUTION = ((1, 12_236), (2, 1512), (3, 112), (4, 182),
                               (5, 28), (6, 14), (7, 14), (8, 70))
EXPECTED_COHERENT_GATES = ((5, 909), (13, 135), (17, 68),
                           (29, 24), (37, 14), (41, 16))
EXPECTED_SEMANTIC = "d356a53ce4b237c7a809d207f2d6ed5cbadfd82228568620f062af60c4d2f6af"


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def egcd(a, b):
    old_r, r = a, b
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t
    return old_r, old_s, old_t


def det(x, y):
    return x[0] * y[1] - x[1] * y[0]


def dot(x, y):
    return x[0] * y[0] + x[1] * y[1]


def norm(x):
    return dot(x, x)


def content(x):
    return gcd(abs(x[0]), abs(x[1]))


def add(x, y):
    return x[0] + y[0], x[1] + y[1]


def scale(k, x):
    return k * x[0], k * x[1]


def gauss_mul(z, w):
    return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]


def gauss_conj(z):
    return z[0], -z[1]


def gauss_div_exact(z, w):
    den = norm(w)
    real_num = z[0] * w[0] + z[1] * w[1]
    imag_num = z[1] * w[0] - z[0] * w[1]
    if real_num % den or imag_num % den:
        return None
    return real_num // den, imag_num // den


def mat_vec(M, x):
    return M[0][0] * x[0] + M[0][1] * x[1], M[1][0] * x[0] + M[1][1] * x[1]


def mat_det(M):
    return M[0][0] * M[1][1] - M[0][1] * M[1][0]


def phi(x):
    m, n = x
    return m * m - n * n, 2 * m * n, m * m + n * n


def lorentz(X, Y):
    return X[2] * Y[2] - X[0] * Y[0] - X[1] * Y[1]


def factorint(n):
    n = abs(n)
    factors = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            factors[p] = factors.get(p, 0) + 1
            n //= p
        p = 3 if p == 2 else p + 2
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors


def valuation(n, p):
    value = 0
    while n % p == 0:
        value += 1
        n //= p
    return value


def is_prime(n):
    if n < 2:
        return False
    return factorint(n) == {n: 1}


def split_prime_rep(p):
    for b in range(1, isqrt(p) + 1):
        a2 = p - b * b
        a = isqrt(a2)
        if a > b and a * a == a2:
            return a, b
    raise RuntimeError(("no split representation", p))


def bezout_v(u):
    """Return v with det(v,u)=1 for primitive u."""
    D, d = u
    g, r, s = egcd(d, D)
    require(g == 1, (u, g))
    v = r, -s
    require(det(v, u) == 1, (u, v, det(v, u)))
    return v


def parabolic(u, x):
    """T_u(x)=x+det(x,u)u."""
    return add(x, scale(det(x, u), u))


def original_rows():
    rows = []
    for d in range(1, DMAX + 1):
        for a in range(SLOPE_EXCESS * d + 1):
            for c in tuple(range(-CMAX, 0)) + tuple(range(1, CMAX + 1)):
                if (a == 0 and c < 0) or (a == SLOPE_EXCESS * d and c > 0):
                    continue
                for p0 in range(1, d + 1):
                    if (a * p0 + c) % d == 0:
                        q0 = p0 + (a * p0 + c) // d
                        rows.append((d, a, c, p0, q0))
    require(len(rows) == 22_890, len(rows))
    return tuple(rows)


def primitive_row(row):
    d, a, c, p0, _ = row
    g = gcd(gcd(d, a), abs(c))
    dd, aa, cc = d // g, a // g, c // g
    pp = next(x for x in range(1, dd + 1) if (aa * x + cc) % dd == 0)
    qq = pp + (aa * pp + cc) // dd
    require((p0 - pp) % dd == 0, (row, (dd, aa, cc, pp, qq)))
    return dd, aa, cc, pp, qq


def shell_and_fingerprint_audit(row):
    d, a, c, p0, q0 = row
    D = d + a
    u = D, d
    x0 = q0, p0
    require(gcd(D, d) == 1, row)
    require(det(x0, u) == c, (row, det(x0, u)))

    v = bezout_v(u)
    k0 = det(v, x0)
    require(x0 == add(scale(c, v), scale(k0, u)), (row, v, k0))

    A, h, C0 = norm(u), dot(x0, u), norm(x0)
    require(A * C0 - h * h == c * c, (row, A, h, C0))
    require((2 * h) ** 2 - 4 * A * C0 == -4 * c * c, row)

    R_minus = c * c + (a - c) * (a - c)
    R_plus = (c + d) * (c + d) + (D - c) * (D - c)
    require(R_plus - R_minus == 2 * d * D, (row, R_minus, R_plus))
    minus_factors = factorint(R_minus)
    plus_factors = factorint(R_plus)
    split_support = {
        prime for prime in set(minus_factors) | set(plus_factors)
        if prime % 4 == 1
    }

    checks = 0
    for n in range(40):
        p = p0 + d * n
        q = q0 + D * n
        x = q, p
        require(x == add(x0, scale(n, u)), (row, n, x))
        require(det(x, u) == c, (row, n, det(x, u)))
        require(norm(x) == A * n * n + 2 * h * n + C0, (row, n))
        require(lorentz(phi(x), phi(u)) == 2 * c * c, (row, n))
        require(content(x) == gcd(abs(c), abs(k0 + n)), (row, n, x, k0))
        require(parabolic(u, x) == add(x, scale(c, u)), (row, n))

        g = content(x)
        y = x[0] // g, x[1] // g
        eps = 2 if y[0] % 2 and y[1] % 2 else 1
        require(content(phi(x)) == g * g * eps, (row, n, x, y))
        require(lorentz(phi(y), phi(u)) == 2 * (c // g) ** 2, (row, n, y))

        C_p = 2 * p * p + 2 * p + 1
        C_q = 2 * q * q + 2 * q + 1
        d_minus = gcd(C_p, q - p)
        d_plus = gcd(C_p, p + q + 1)
        require(gcd(C_p, C_q) == d_minus * d_plus, (row, n, C_p, C_q))
        require(gcd(d_minus, d_plus) == 1, (row, n, d_minus, d_plus))
        require(R_minus % d_minus == 0, (row, n, d_minus, R_minus))
        require(R_plus % d_plus == 0, (row, n, d_plus, R_plus))

        L_minus = q - p
        L_plus = p + q + 1
        require(
            a * a * C_p
            == R_minus + 2 * d * L_minus * (d * L_minus + a - 2 * c),
            ("minus resultant", row, n),
        )
        b = 2 * d + a
        require(
            b * b * C_p
            == R_plus + 2 * d * L_plus * (d * L_plus + a - 2 * c),
            ("plus resultant", row, n),
        )
        channels = (
            (a, c, L_minus, a * (p + 1) - c, R_minus, d_minus, minus_factors),
            (b, c + d, L_plus, b * (p + 1) - (c + d),
             R_plus, d_plus, plus_factors),
        )
        for coefficient, constant, H_value, J_value, resultant, channel, factors in channels:
            require(coefficient * p + constant == d * H_value,
                    ("linear branch", row, n, coefficient, constant, H_value))
            require(gcd(C_p, resultant) == gcd(C_p, d * H_value * J_value),
                    ("resultant gcd", row, n, resultant, H_value, J_value))
            for prime in factors:
                if C_p % prime or prime == 2:
                    continue
                regular = (d * gcd(coefficient, constant)) % prime != 0
                if not regular:
                    continue
                H_hit = (d * H_value) % prime == 0
                J_hit = J_value % prime == 0
                require(H_hit != J_hit,
                        ("branch selector", row, n, prime, H_value, J_value))
                if H_hit:
                    require(valuation(channel, prime)
                            == min(valuation(C_p, prime), valuation(resultant, prime)),
                            ("matching valuation", row, n, prime, channel, C_p, resultant))
                else:
                    require(channel % prime != 0,
                            ("conjugate overcount", row, n, prime, channel))
        checks += 1

    return len(split_support), checks, R_minus, R_plus


def roots_mod_quadratic(u, x0, c, modulus):
    A, h, C0 = norm(u), dot(x0, u), norm(x0)
    require(gcd(modulus, 2 * A * c) == 1, (u, x0, c, modulus))
    roots = tuple(n for n in range(modulus) if (A * n * n + 2 * h * n + C0) % modulus == 0)
    return roots


def root_atlas_audit():
    moduli = (5, 13, 17, 25, 29, 37, 41, 65, 85, 125, 145, 169, 185, 221, 325)
    atlas_checks = 0
    pair_checks = 0
    allocation_checks = 0
    for D in range(1, 9):
        for d in range(0, D + 1):
            if gcd(D, d) != 1:
                continue
            u = D, d
            v = bezout_v(u)
            for c in (-7, -3, -1, 1, 2, 5, 8):
                for k in (-3, 0, 4):
                    x0 = add(scale(c, v), scale(k, u))
                    A, h = norm(u), dot(x0, u)
                    require(A * norm(x0) - h * h == c * c, (u, x0, c))
                    for M in moduli:
                        if gcd(M, 2 * A * c) != 1:
                            continue
                        factors = factorint(M)
                        if any(p % 4 != 1 for p in factors):
                            continue
                        roots = roots_mod_quadratic(u, x0, c, M)
                        require(len(roots) == 2 ** len(factors), (u, x0, c, M, roots))
                        reflected = {(-2 * h * pow(A, -1, M) - n) % M for n in roots}
                        require(reflected == set(roots), (u, x0, c, M, roots, reflected))
                        for n in roots:
                            x = add(x0, scale(n, u))
                            require(norm(x) % M == 0, (u, x0, c, M, n))
                            require(gcd(gcd(abs(x[0]), abs(x[1])), M) == 1,
                                    ("nonzero isotropic", u, x0, c, M, n, x))
                            L = A * n + h, c
                            require(norm(L) == A * norm(x), (u, x0, c, M, n, L, x))
                            allocation = (1, 0)
                            for prime, exponent in factors.items():
                                aa, bb = split_prime_rep(prime)
                                pi = aa, bb
                                pi_bar = aa, -bb
                                current = L
                                plus_depth = 0
                                minus_depth = 0
                                while (quotient := gauss_div_exact(current, pi)) is not None:
                                    plus_depth += 1
                                    current = quotient
                                current = L
                                while (quotient := gauss_div_exact(current, pi_bar)) is not None:
                                    minus_depth += 1
                                    current = quotient
                                require((plus_depth >= exponent) != (minus_depth >= exponent),
                                        ("allocation branch", u, x0, c, M, n,
                                         prime, exponent, plus_depth, minus_depth, L))
                                chosen = pi if plus_depth >= exponent else pi_bar
                                for _ in range(exponent):
                                    allocation = gauss_mul(allocation, chosen)
                            require(norm(allocation) == M, (u, x0, c, M, n, allocation))
                            allocation_checks += 1
                        for r in roots:
                            for s in roots:
                                delta_minus = gcd(M, s - r)
                                delta_plus = gcd(M, A * (r + s) + 2 * h)
                                require(delta_minus * delta_plus == M,
                                        (u, x0, c, M, r, s, delta_minus, delta_plus))
                                require(gcd(delta_minus, delta_plus) == 1,
                                        (u, x0, c, M, r, s, delta_minus, delta_plus))
                                z_r = add(x0, scale(r, u))
                                z_s = add(x0, scale(s, u))
                                same_product = gauss_mul(z_s, gauss_conj(z_r))
                                opposite_product = gauss_mul(z_r, z_s)
                                require(gcd(M, content(same_product)) == delta_minus,
                                        ("same channel", u, x0, c, M, r, s,
                                         same_product, delta_minus))
                                require(gcd(M, content(opposite_product)) == delta_plus,
                                        ("opposite channel", u, x0, c, M, r, s,
                                         opposite_product, delta_plus))
                                pair_checks += 1
                        atlas_checks += 1
    return atlas_checks, allocation_checks, pair_checks


def gaussian_gate_audit():
    gate_checks = 0
    primes = tuple(p for p in range(5, 100) if p % 4 == 1 and is_prime(p))
    for p in primes:
        a, b = split_prime_rep(p)
        d0, e0 = a * a - b * b, 2 * a * b
        matrices = (
            ((d0, e0), (-e0, d0)),
            ((d0, -e0), (e0, d0)),
        )
        for H in matrices:
            require(mat_det(H) == p * p, (p, H, mat_det(H)))
            require(gcd(gcd(abs(H[0][0]), abs(H[0][1])),
                        gcd(abs(H[1][0]), abs(H[1][1]))) == 1, (p, H))
            anchors = []
            for D in range(1, 3 * p + 1):
                for d in range(0, D + 1):
                    u = D, d
                    if gcd(D, d) == 1 and all(z % p == 0 for z in mat_vec(H, u)):
                        anchors.append(u)
                        if len(anchors) == 4:
                            break
                if len(anchors) == 4:
                    break
            require(anchors, ("no kernel anchor", p, H))
            for u in anchors:
                U0 = mat_vec(H, u)
                U = U0[0] // p, U0[1] // p
                require(norm(U) == norm(u), (p, H, u, U))
                v = bezout_v(u)
                for c in range(-2 * p, 2 * p + 1):
                    if c == 0:
                        continue
                    for k in (-p - 1, -1, 0, 1, p + 2):
                        x = add(scale(c, v), scale(k, u))
                        Hx = mat_vec(H, x)
                        integral = Hx[0] % p == 0 and Hx[1] % p == 0
                        require(integral == (c % p == 0),
                                ("toggle gate", p, H, u, c, k, x, Hx))
                        if integral:
                            y = Hx[0] // p, Hx[1] // p
                            require(det(y, U) == c, (p, H, u, c, k, y, U))
                            require(norm(y) == norm(x), (p, H, u, c, k, y))
                        gate_checks += 1
    return len(primes), gate_checks


def main():
    original = original_rows()
    primitive = tuple(sorted({primitive_row(row) for row in original}))
    require(len(primitive) == EXPECTED_PRIMITIVE, len(primitive))
    by_d = tuple(sum(row[0] == d for row in primitive) for d in range(1, DMAX + 1))
    expected_by_d = tuple(2 * CMAX * SLOPE_EXCESS * sum(
        1 for a in range(1, d + 1) if gcd(a, d) == 1
    ) for d in range(1, DMAX + 1))
    require(by_d == expected_by_d, (by_d, expected_by_d))
    require(sum(by_d) == 2 * CMAX * SLOPE_EXCESS * 22, by_d)

    rank_distribution = Counter()
    fibre_distribution = Counter()
    coherent_gate_distribution = Counter()
    coherent_gate_simple = 0
    coherent_gate_deep = 0
    shell_checks = 0
    low_unit = 0
    a_one = 0
    digest = sha256()
    for row in primitive:
        rank, checks, R_minus, R_plus = shell_and_fingerprint_audit(row)
        rank_distribution[rank] += 1
        shell_checks += checks
        d, a, c, _, _ = row
        D = d + a
        M = min(DMAX // d, CMAX // abs(c))
        require(M >= 1, (row, M))
        fibre_distribution[M] += 1
        anchor_factors = factorint(D * D + d * d)
        row_gates = tuple(prime for prime in anchor_factors
                          if prime % 4 == 1 and c % prime == 0)
        require(len(row_gates) <= 1, ("two coherent gates", row, row_gates))
        for prime in row_gates:
            coherent_gate_distribution[prime] += 1
            if anchor_factors[prime] == 1:
                coherent_gate_simple += 1
            else:
                coherent_gate_deep += 1
        low_unit += int(d + D < 8 and abs(c) == 1)
        a_one += int(a == 1)
        digest.update((repr((row, R_minus, R_plus, rank)) + "\n").encode())
    require(tuple(sorted(rank_distribution.items())) == EXPECTED_RANK_DISTRIBUTION,
            rank_distribution)
    require(tuple(sorted(fibre_distribution.items())) == EXPECTED_FIBRE_DISTRIBUTION,
            fibre_distribution)
    require(sum(count * (M * (M + 1) // 2)
                for M, count in fibre_distribution.items()) == len(original),
            fibre_distribution)
    require(tuple(sorted(coherent_gate_distribution.items())) == EXPECTED_COHERENT_GATES,
            coherent_gate_distribution)
    require((coherent_gate_simple, coherent_gate_deep) == (978, 188),
            (coherent_gate_simple, coherent_gate_deep))
    require(low_unit == 17, low_unit)
    require(a_one == 736, a_one)

    p, q = 14_426_006, 28_851_968
    require(q - 2 * p == -44, (p, q))
    require(gcd(p, q) == 2 and p >= 264 and p < q < 8 * p, (p, q, gcd(p, q)))
    C_p, C_q = 2 * p * p + 2 * p + 1, 2 * q * q + 2 * q + 1
    d_minus = gcd(C_p, q - p)
    d_plus = gcd(C_p, p + q + 1)
    require((d_minus, d_plus) == (3961, 3965), (d_minus, d_plus))
    require(factorint(d_minus) == {17: 1, 233: 1}, factorint(d_minus))
    require(factorint(d_plus) == {5: 1, 13: 1, 61: 1}, factorint(d_plus))
    require(gcd(C_p, C_q) == 15_705_365, gcd(C_p, C_q))
    require((C_p, C_q) == (416_219_327_076_085, 1_664_872_172_649_985),
            (C_p, C_q))
    require((q - p, p + q + 1) == (14_425_962, 43_277_975),
            (q - p, p + q + 1))
    require(p % 3961 == 44 and (3 * p - 43) % 3965 == 0, p)

    overcount_hostiles = (
        (1, 2, -1, 1, 2, "minus"),
        (1, 6, -5, 1, 2, "plus"),
        (5, 1, -46, 51, 52, "d-exception"),
    )
    for d, a, c, pp, qq, kind in overcount_hostiles:
        Cp = 2 * pp * pp + 2 * pp + 1
        Rm = c * c + (a - c) ** 2
        Rp = (c + d) ** 2 + (d + a - c) ** 2
        gm = gcd(Cp, qq - pp)
        gp = gcd(Cp, pp + qq + 1)
        if kind == "minus":
            require((Cp, Rm, gm, a * (pp + 1) - c) == (5, 10, 1, 5),
                    (kind, Cp, Rm, gm))
        elif kind == "plus":
            require((Cp, Rp, gp, (2 * d + a) * (pp + 1) - (c + d))
                    == (5, 160, 1, 20), (kind, Cp, Rp, gp))
        else:
            require((gcd(Cp, Rm), gcd(Cp, Rp), gm, gp) == (5, 5, 1, 1),
                    (kind, Cp, Rm, Rp, gm, gp))

    carrier_cutoffs = []
    for d in (6, 2):
        floor = Fraction(1, 105)
        coefficient = Fraction(12 * d, 35) + 10
        ratio = coefficient / (floor - LRC_DMAX / 5)
        carrier_cutoffs.append(ratio.numerator // ratio.denominator + 1)
    require(tuple(carrier_cutoffs) == (1898, 1682), carrier_cutoffs)

    atlas_checks, allocation_checks, atlas_pairs = root_atlas_audit()
    split_primes, gate_checks = gaussian_gate_audit()
    semantic = digest.hexdigest()
    require(semantic == EXPECTED_SEMANTIC, semantic)

    print("THM-3356 PRIMITIVE AFFINE DETERMINANT-SHELL AUDIT")
    print("formal_rows", len(original), "primitive_rows", len(primitive), "by_d", by_d)
    print("shell_lorentz_gcd_resultant_checks", shell_checks)
    print("low_unit_farey_orbits", low_unit, "u_spine_parallel_rows", a_one)
    print("potential_split_rank_distribution", tuple(sorted(rank_distribution.items())))
    print("normalization_fibre_depth_distribution", tuple(sorted(fibre_distribution.items())))
    print("coherent_shell_gate_rows", sum(coherent_gate_distribution.values()),
          "by_prime", tuple(sorted(coherent_gate_distribution.items())),
          "simple_deep", (coherent_gate_simple, coherent_gate_deep))
    print("rank5_witness", (p, q), "channels", (d_minus, d_plus),
          "common", gcd(C_p, C_q))
    print("norm85_toggle_carrier_cutoffs", tuple(carrier_cutoffs))
    print("unramified_root_atlases", atlas_checks, "allocation_checks", allocation_checks,
          "ordered_root_pairs", atlas_pairs)
    print("coherent_gaussian_gate_primes", split_primes, "checks", gate_checks)
    print("semantic_sha256", semantic)
    print("status=ALL EXACT CHECKS PASSED; no LRC physical-tail closure claimed")


if __name__ == "__main__":
    main()
