#!/usr/bin/env python3
"""Exact controls for THM-3370's Berggren/two-cube norm collision."""

from hashlib import sha256
from math import gcd, isqrt


BOX = 5_000
EXPECTED_MOD7 = (2, 3, 4)
EXPECTED_MOD9 = (1, 2, 4, 6, 7)
EXPECTED_MOD63 = (2, 4, 10, 11, 16, 24, 25, 31, 37, 38, 46, 51, 52, 58, 60)
EXPECTED_HITS = (
    (13_712_211, 107, 232, 3_703, 1_851),
    (27_127_737_027, 360, 3_003, 164_705, 82_352),
    (48_789_299_691, 1_907, 3_472, 220_883, 110_441),
    (122_025_161_043, 107, 4_960, 349_321, 174_660),
    (181_178_773_803, 4_403, 4_576, 425_651, 212_825),
)
PELL_D = 33
PELL_U = 23
PELL_V = 4
EXPECTED_PELL_FIRST = (
    2_835,
    987,
    499,
    -488,
    1_417,
    8_037_227,
)
EXPECTED_NONPRIMITIVE_SHARP = (591, 1_183, 120, -69, 27_441, 3_049)


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def u_spine(r):
    return (2 * r + 1, 2 * r * (r + 1), 2 * r * (r + 1) + 1)


def scalar(r):
    return 4 * r * (r + 1) + 3


def matvec(matrix, vector):
    return tuple(sum(row[j] * vector[j] for j in range(3)) for row in matrix)


def cube_sum_residues(modulus):
    cubes = {pow(x, 3, modulus) for x in range(modulus)}
    return {(x + y) % modulus for x in cubes for y in cubes}


def allowed_depths(modulus):
    sums = cube_sum_residues(modulus)
    return tuple(
        r for r in range(modulus) if scalar(r) % modulus in sums
    )


def factor(n):
    result = []
    p = 2
    while p * p <= n:
        if n % p:
            p = 3 if p == 2 else p + 2
            continue
        exponent = 0
        while n % p == 0:
            n //= p
            exponent += 1
        result.append((p, exponent))
        p = 3 if p == 2 else p + 2
    if n > 1:
        result.append((n, 1))
    return tuple(result)


def cube_root_floor(n):
    lo, hi = 0, 1
    while hi**3 <= n:
        hi *= 2
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if mid**3 <= n:
            lo = mid
        else:
            hi = mid
    return lo


def scan_box(bound):
    hits = []
    for y in range(2, bound + 1):
        y3 = y**3
        for x in range(1, y):
            candidate_square = x**3 + y3 - 2
            a = isqrt(candidate_square)
            if a * a != candidate_square or a % 2 == 0:
                continue
            r = (a - 1) // 2
            hits.append((candidate_square + 2, x, y, a, r))
    return tuple(sorted(hits))


def main():
    U = ((1, -2, 2), (2, -1, 2), (2, -2, 3))
    for r in range(1, BOX + 1):
        a, b, c = u_spine(r)
        require(a * a + b * b == c * c, ("Pythagorean", r))
        require(c - b == 1, ("spine plane", r))
        require(matvec(U, (a, b, c)) == u_spine(r + 1), ("U action", r))
        require(scalar(r) == a * a + 2 == 2 * c + 1 == 2 * b + 3, ("Q", r))
        require(c == (r + 1) ** 2 + r**2, ("Gaussian norm", r))
        require(gcd(c, scalar(r)) == 1, ("coprime carriers", r))

    mod7 = allowed_depths(7)
    mod9 = allowed_depths(9)
    mod63 = allowed_depths(63)
    require(mod7 == EXPECTED_MOD7, ("mod7", mod7))
    require(mod9 == EXPECTED_MOD9, ("mod9", mod9))
    require(mod63 == EXPECTED_MOD63, ("mod63", mod63))

    hits = scan_box(BOX)
    require(hits == EXPECTED_HITS, ("box hits", hits))
    details = []
    for n, x, y, a, r in hits:
        require(n == x**3 + y**3 == a * a + 2 == scalar(r), ("identity", x, y))
        g = gcd(x, y)
        X, Y = x // g, y // g
        d = x + y
        q = x * x - x * y + y * y
        q0 = X * X - X * Y + Y * Y
        require(n == d * q and q == g * g * q0, ("cube factorization", x, y))

        numerator = d**3 - n
        require(numerator % (3 * d) == 0, ("good divisor integrality", x, y))
        s = numerator // (3 * d)
        delta_numerator = 4 * n - d**3
        require(delta_numerator % (3 * d) == 0, ("discriminant integrality", x, y))
        delta = delta_numerator // (3 * d)
        require(s == x * y and delta == (y - x) ** 2, ("good divisor", x, y))
        require(0 < y - x < d, ("positive distinct", x, y))

        q0_factors = factor(q0)
        v3 = dict(q0_factors).get(3, 0)
        require(v3 <= 1, ("v3 primitive cofactor", x, y, q0_factors))
        for prime, _exponent in q0_factors:
            if prime == 3:
                continue
            require(prime % 24 in (1, 19), ("mod24 cofactor", x, y, prime))

        details.append(
            (x, y, a, r, n, g, d, q0, factor(n), factor(d), q0_factors)
        )

    first_n, first_x, first_y, first_a, first_r = hits[0]
    max_y_below = cube_root_floor(first_n - 1)
    require(max_y_below == 239, max_y_below)
    for y in range(2, max_y_below + 1):
        for x in range(1, y):
            n = x**3 + y**3
            if n >= first_n:
                continue
            a = isqrt(n - 2)
            require(a * a != n - 2 or a % 2 == 0, ("smaller hit", x, y, n))
    require((first_x, first_y, first_a, first_r) == (107, 232, 3703, 1851), "first")

    require(PELL_U * PELL_U - PELL_D * PELL_V * PELL_V == 1, "Pell unit")
    pell_a, pell_e, pell_d = 63, 21, 11
    require((pell_d, PELL_D) == (11, 3 * pell_d), "Pell fibre")
    pell_rows = []
    for index in range(6):
        require(pell_a % 2 == pell_e % 2 == pell_d % 2 == 1, ("Pell parity", index))
        pell_x = (pell_d + pell_e) // 2
        pell_y = (pell_d - pell_e) // 2
        pell_n = pell_a * pell_a + 2
        pell_r = (pell_a - 1) // 2
        require(pell_x**3 + pell_y**3 == pell_n == scalar(pell_r), ("Pell identity", index))
        require((2 * pell_a) ** 2 - PELL_D * pell_e**2 == pell_d**3 - 8, ("Pell norm", index))
        pell_rows.append((pell_a, pell_e, pell_x, pell_y, pell_r, pell_n))
        next_a = PELL_U * pell_a + (PELL_D * PELL_V // 2) * pell_e
        next_e = 2 * PELL_V * pell_a + PELL_U * pell_e
        require(next_a > pell_a and next_e > pell_e, ("Pell growth", index))
        pell_a, pell_e = next_a, next_e
    require(pell_rows[1] == EXPECTED_PELL_FIRST, ("Pell first", pell_rows[1]))

    sharp_r, sharp_a, sharp_x, sharp_y, sharp_q, sharp_q0 = EXPECTED_NONPRIMITIVE_SHARP
    require(scalar(sharp_r) == sharp_a * sharp_a + 2 == sharp_x**3 + sharp_y**3, "sharp identity")
    sharp_g = gcd(sharp_x, sharp_y)
    require(sharp_g == 3, sharp_g)
    require(sharp_q == sharp_x * sharp_x - sharp_x * sharp_y + sharp_y * sharp_y, "sharp q")
    require(sharp_q == sharp_g * sharp_g * sharp_q0, "sharp primitive reduction")
    require(factor(sharp_q) == ((3, 2), (3_049, 1)), factor(sharp_q))
    require(factor(sharp_q0) == ((3_049, 1),) and 3_049 % 24 == 1, "sharp q0")

    semantic = sha256()
    semantic.update(
        repr(
            (
                mod7,
                mod9,
                mod63,
                tuple(details),
                tuple(pell_rows),
                EXPECTED_NONPRIMITIVE_SHARP,
            )
        ).encode()
    )

    print("BERGGREN U-SPINE / TWO-DISTINCT-CUBES NORM COLLISION")
    print("status=PROVED filters/Pell family + FINITE-EXACT positive box/minimum")
    print("Q_r=4r(r+1)+3=(2r+1)^2+2=2C_r+1;C_r=(r+1)^2+r^2")
    print("gcd(C_r,Q_r)=1;primitive_Eisenstein_cofactor_primes_mod24=(1,19)")
    print(f"allowed_r_mod7={mod7};allowed_r_mod9={mod9};allowed_r_mod63={mod63}")
    print(f"box=1<=x<y<={BOX};pairs={BOX * (BOX - 1) // 2};hits={len(hits)}")
    for detail in details:
        print(f"hit={detail}")
    print(
        f"global_minimum={first_n}={first_x}^3+{first_y}^3={first_a}^2+2="
        f"Q_{first_r};all_smaller_solutions_have_y<={max_y_below}"
    )
    print(f"pell_unit=({PELL_U},{PELL_V});pell_unit_norm=1;fixed_cube_sum={pell_d}")
    print(f"pell_signed_seed={pell_rows[0]}")
    print(f"pell_first_signed_iterate={pell_rows[1]}")
    print(f"nonprimitive_sharp={EXPECTED_NONPRIMITIVE_SHARP}")
    print(f"semantic_sha256={semantic.hexdigest()}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
