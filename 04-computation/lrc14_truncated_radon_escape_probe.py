#!/usr/bin/env python3
"""Exact prototype for the punctured-stalk truncated Radon observable.

This is deliberately independent of the THM-2436/2506 atlas programs.  It
checks the linear-algebra kernel, constructs an integral sharpness witness for
the seven-of-twelve slope count, replays the two-row THM-2506 hostile, and
exhibits the ordered-section carry under a physical CRT translation.
"""

from __future__ import annotations

from itertools import combinations


P = 13
R = 7
ZERO = (0,) * 12
ONE = (1,) + (0,) * 11


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def add(a: tuple[int, ...], b: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(x + y for x, y in zip(a, b))


def neg(a: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(-x for x in a)


def mul(a: tuple[int, ...], b: tuple[int, ...]) -> tuple[int, ...]:
    """Multiply in Z[z]/(1+z+...+z^12)."""
    raw = [0] * 23
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            raw[i + j] += x * y
    for e in range(22, 11, -1):
        c = raw[e]
        if c:
            # z^e = -(z^(e-12)+...+z^(e-1)).
            for j in range(12):
                raw[e - 12 + j] -= c
            raw[e] = 0
    return tuple(raw[:12])


def zeta_pow(e: int) -> tuple[int, ...]:
    e %= P
    if e < 12:
        ans = [0] * 12
        ans[e] = 1
        return tuple(ans)
    return (-1,) * 12


def dft(vec: list[int], alpha: int) -> tuple[int, ...]:
    ans = ZERO
    for v, coefficient in enumerate(vec):
        if coefficient:
            term = tuple(coefficient * x for x in zeta_pow(-alpha * v))
            ans = add(ans, term)
    return ans


def poly_mul_linear(
    coefficients: list[tuple[int, ...]], root: tuple[int, ...]
) -> list[tuple[int, ...]]:
    out = [ZERO] * (len(coefficients) + 1)
    for i, coefficient in enumerate(coefficients):
        out[i] = add(out[i], neg(mul(coefficient, root)))
        out[i + 1] = add(out[i + 1], coefficient)
    return out


def coefficient_to_column(a: tuple[int, ...]) -> list[int]:
    """Return q with sum_h q[h] z^(-h) = a and q[1]=0."""
    q = [0] * P
    q[0] = a[0]
    for exponent in range(1, 12):
        q[P - exponent] = a[exponent]
    require(dft(q, 1) == a, "cyclotomic coefficient lift failed")
    return q


def radon(d: list[list[int]], tau: int) -> list[int]:
    return [sum(d[(v - tau * r) % P][r] for r in range(R)) for v in range(P)]


def l1_matrix(d: list[list[int]]) -> int:
    return sum(abs(x) for row in d for x in row)


def rank_mod(matrix: list[list[int]], prime: int) -> int:
    a = [[x % prime for x in row] for row in matrix]
    rows = len(a)
    cols = len(a[0])
    rank = 0
    for col in range(cols):
        pivot = next((i for i in range(rank, rows) if a[i][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        inv = pow(a[rank][col], -1, prime)
        a[rank] = [(x * inv) % prime for x in a[rank]]
        for i in range(rows):
            if i != rank and a[i][col]:
                c = a[i][col]
                a[i] = [(x - c * y) % prime for x, y in zip(a[i], a[rank])]
        rank += 1
        if rank == cols:
            break
    return rank


def row_zero_radon_matrix(slopes: range) -> list[list[int]]:
    """Use d(h,6)=-sum_{r<6}d(h,r), so the domain has dimension 78."""
    matrix: list[list[int]] = []
    for tau in slopes:
        for v in range(P):
            row = [0] * (P * (R - 1))
            for h in range(P):
                for r in range(R - 1):
                    index = h * (R - 1) + r
                    if h == (v - tau * r) % P:
                        row[index] += 1
                    if h == (v - tau * (R - 1)) % P:
                        row[index] -= 1
            matrix.append(row)
    return matrix


def det_bareiss(matrix: list[list[int]]) -> int:
    a = [row[:] for row in matrix]
    n = len(a)
    sign = 1
    previous = 1
    for k in range(n - 1):
        pivot = next((i for i in range(k, n) if a[i][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            sign *= -1
        p = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = a[i][j] * p - a[i][k] * a[k][j]
                require(numerator % previous == 0, "Bareiss division was not exact")
                a[i][j] = numerator // previous
        previous = p
        for i in range(k + 1, n):
            a[i][k] = 0
        for j in range(k + 1, n):
            a[k][j] = 0
    return sign * a[-1][-1]


def norm(a: tuple[int, ...]) -> int:
    columns = []
    for j in range(12):
        columns.append(mul(a, zeta_pow(j)))
    multiplication_matrix = [
        [columns[j][i] for j in range(12)] for i in range(12)
    ]
    return det_bareiss(multiplication_matrix)


def sharp_seven_of_twelve_witness() -> list[list[int]]:
    # P_1(X)=(X-1) product_{tau=1}^5 (X-zeta^(-tau)).
    coefficients = [ONE]
    for tau in range(6):
        coefficients = poly_mul_linear(coefficients, zeta_pow(-tau))
    require(len(coefficients) == R, "sharp polynomial has wrong degree")

    columns = [coefficient_to_column(a) for a in coefficients]
    d = [[columns[r][h] for r in range(R)] for h in range(P)]
    require(all(sum(row) == 0 for row in d), "sharp witness is not row-zero")
    require(
        any(d[h] != d[0] for h in range(1, P)),
        "sharp witness is h-independent",
    )
    for alpha in range(1, P):
        bad = [tau for tau in range(1, P) if dft(radon(d, tau), alpha) == ZERO]
        require(bad == [1, 2, 3, 4, 5], "sharp bad-slope set drifted")
    return d


def thm2506_two_row_hostile() -> list[list[int]]:
    d = [[0] * R for _ in range(P)]
    d[0][5] = 1
    d[0][3] = -1
    d[1][5] = 1
    d[1][4] = -1
    require(all(sum(row) == 0 for row in d), "two-row hostile is not row-zero")
    return d


def translate_crt(d: list[list[int]], kappa: int) -> list[list[int]]:
    a = kappa % P
    c = kappa % R
    return [
        [d[(h - a) % P][(r - c) % R] for r in range(R)] for h in range(P)
    ]


def main() -> None:
    full_matrix = row_zero_radon_matrix(range(1, P))
    ranks = {prime: rank_mod(full_matrix, prime) for prime in (101, 103)}
    require(ranks == {101: 72, 103: 72}, "Radon rank drifted")
    print("row-zero domain dimension: 78")
    print("all-nonzero-slope ranks mod 101,103:", ranks[101], ranks[103])
    print("kernel dimension: 6 (the h-independent zero-sum columns)")

    six_slope_banks = 0
    for slopes in combinations(range(1, P), 6):
        matrix = row_zero_radon_matrix(slopes)
        require(rank_mod(matrix, 101) == 72, f"six-slope rank failed at {slopes}")
        six_slope_banks += 1
    require(six_slope_banks == 924, "wrong number of six-slope banks")
    print("all 924 six-nonzero-slope banks have exact rational rank: 72")

    sharp = sharp_seven_of_twelve_witness()
    print("sharp integral witness L1:", l1_matrix(sharp))
    print("sharp bad nonzero slopes for every alpha: 1 2 3 4 5")
    print("sharp good nonzero slopes for every alpha: 7")

    hostile = thm2506_two_row_hostile()
    require(l1_matrix(hostile) == 4, "two-row hostile L1 drifted")
    good = [tau for tau in range(1, P) if radon(hostile, tau) != [0] * P]
    require(good == list(range(1, P)), "two-row hostile lost a slope")
    norms = []
    for tau in good:
        rv = radon(hostile, tau)
        require(sum(rv) == 0, "Radon output is not zero-sum")
        require(sum(abs(x) for x in rv) <= 4, "Radon L1 contraction failed")
        colours = [dft(rv, alpha) for alpha in range(1, P)]
        require(all(theta != ZERO for theta in colours), "a C13 colour vanished")
        n = norm(colours[0])
        require(n != 0, "a nonzero colour has zero norm")
        norms.append(abs(n))
    print("THM-2506 two-row hostile L1 / good slopes:", 4, len(good))
    print("THM-2506 hostile absolute C13 norm range:", min(norms), max(norms))

    moved = translate_crt(hostile, 2)
    before = radon(hostile, 1)
    after = radon(moved, 1)
    require(
        not any(after == before[s:] + before[:s] for s in range(P)),
        "CRT carry hostile became a pure translation",
    )
    print("CRT-translation carry witness kappa/tau: 2 1")
    print("before:", " ".join(map(str, before)))
    print("after: ", " ".join(map(str, after)))
    print("after is not a cyclic translate of before: PASS")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
