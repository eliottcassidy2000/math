#!/usr/bin/env python3
"""Finite hostile audit of the proved predecessor digital-skeleton lemma."""

from math import comb


LIMIT = 1000
DIRECT_UNIT_LIMIT = 120


def require(condition, data):
    if not condition:
        raise RuntimeError(data)


def prime_divisors(n):
    ans = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            ans.append(p)
            while n % p == 0:
                n //= p
        p += 1 if p == 2 else 2
    if n > 1:
        ans.append(n)
    return tuple(ans)


def vp_integer(n, p):
    require(n > 0, (n, p, "positive valuation input"))
    answer = 0
    while n % p == 0:
        n //= p
        answer += 1
    return answer


def vp_factorial(n, p):
    answer = 0
    while n:
        n //= p
        answer += n
    return answer


def weight(n, j, p):
    return (
        vp_factorial(n, p)
        - vp_factorial(j, p)
        - vp_factorial(n - j, p)
        + vp_factorial(2 * j, p)
    )


def hull(points):
    answer = []
    for point in points:
        x2, y2 = point
        while len(answer) >= 2:
            x0, y0 = answer[-2]
            x1, y1 = answer[-1]
            if (y1 - y0) * (x2 - x1) >= (y2 - y1) * (x1 - x0):
                answer.pop()
            else:
                break
        answer.append(point)
    return tuple(answer)


def slope_numerator_at_scale(p, power):
    return vp_factorial(2 * power, p)


def z_mod_p(n, j, p):
    m = n - j
    answer = 0
    rising = 1
    for ell in range(m + 1):
        if ell:
            rising = rising * (2 * j + ell) % p
        term = comb(m, ell) * pow(n + 1, m - ell, p) * rising
        answer = (answer - term if ell % 2 else answer + term) % p
    return answer


def main():
    pair_count = 0
    point_count = 0
    vertex_count = 0
    direct_unit_checks = 0
    binary_slope_controls = []
    odd_slope_controls = []
    for n in range(2, LIMIT + 1):
        for p in prime_divisors(n):
            pair_count += 1
            m = n // p
            half = (p - 1) // 2
            points = []
            for j in range(n + 1):
                q, a = divmod(j, p)
                actual = weight(n, j, p)
                h_q = 2 * q + weight(m, q, p)
                if a == 0:
                    predicted = h_q
                elif p == 2:
                    predicted = h_q + 2 + vp_integer(m - q, 2)
                elif a <= half:
                    predicted = h_q + 1 + vp_integer(m - q, p)
                else:
                    predicted = 2 * (q + 1) + weight(m, q + 1, p)
                require(actual == predicted, (n, p, j, actual, predicted))
                points.append((j, actual))
                point_count += 1

            vertices = hull(points)
            residues = {j % p for j, _ in vertices}
            allowed = {0} if p == 2 else {0, half}
            require(residues <= allowed, (n, p, residues, allowed))
            vertex_count += len(vertices)

            # The residue restriction is exactly the unit mechanism:
            # j=0 mod p uses the binomial/rising-factorial split, while
            # j=(p-1)/2 mod p makes the first rising factor vanish.
            require(all(j % p in allowed for j, _ in vertices),
                    (n, p, vertices, allowed))
            if n <= DIRECT_UNIT_LIMIT:
                for j, _ in vertices:
                    require(z_mod_p(n, j, p) == 1, (n, p, j))
                    direct_unit_checks += 1

    for exponent in range(1, 11):
        power = 2**exponent
        numerator = slope_numerator_at_scale(2, power)
        require(numerator == 2 * power - 1,
                (2, exponent, numerator, power))
        binary_slope_controls.append((exponent, numerator, power))
    for p in (3, 5, 7, 11, 13):
        for exponent in range(1, 6):
            power = p**exponent
            numerator = slope_numerator_at_scale(p, power)
            require(numerator == 2 * (power - 1) // (p - 1),
                    (p, exponent, numerator, power))
            odd_slope_controls.append((p, exponent, numerator, power))

    print("THM-3161 PREDECESSOR DIGITAL SKELETON HOSTILE AUDIT")
    print("universe=2<=N<=%d, every prime p|N" % LIMIT)
    print("pairs=%d points=%d lower_hull_vertices=%d" % (pair_count, point_count, vertex_count))
    print("one-digit recurrence failures=0")
    print("vertex residue failures=0")
    print("direct Z_j unit checks through N=%d: %d failures=0" % (DIRECT_UNIT_LIMIT, direct_unit_checks))
    print("binary slope controls=%s" % (tuple(binary_slope_controls),))
    print("odd slope control count=%d" % len(odd_slope_controls))
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
