#!/usr/bin/env python3
"""Exact companion for THM-3335.

This script checks consequence objects, not a fitted recurrence.  It combines
an independent direct square-triangular scan, a Pell/Markov state compiler,
primitive Euclid-parameter classification, two exact LRC maximin engines, and
the tournament/FC/JC hostile boundaries used by the theorem.
"""

from __future__ import annotations

from fractions import Fraction
from functools import reduce
from hashlib import sha256
from itertools import combinations
from math import gcd, isqrt


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def triangular(n: int) -> int:
    return n * (n + 1) // 2


def square_root(n: int) -> int | None:
    if n < 0:
        return None
    r = isqrt(n)
    return r if r * r == n else None


def phi(m: int, n: int) -> tuple[int, int, int]:
    return m * m - n * n, 2 * m * n, m * m + n * n


def lorentz(x: tuple[int, int, int], y: tuple[int, int, int]) -> int:
    return x[2] * y[2] - x[0] * y[0] - x[1] * y[1]


def matrix_mul(
    a: tuple[tuple[int, int], tuple[int, int]],
    b: tuple[tuple[int, int], tuple[int, int]],
) -> tuple[tuple[int, int], tuple[int, int]]:
    return (
        (a[0][0] * b[0][0] + a[0][1] * b[1][0],
         a[0][0] * b[0][1] + a[0][1] * b[1][1]),
        (a[1][0] * b[0][0] + a[1][1] * b[1][0],
         a[1][0] * b[0][1] + a[1][1] * b[1][1]),
    )


def matrix_pow(
    a: tuple[tuple[int, int], tuple[int, int]], exponent: int
) -> tuple[tuple[int, int], tuple[int, int]]:
    result = ((1, 0), (0, 1))
    base = a
    while exponent:
        if exponent & 1:
            result = matrix_mul(result, base)
        base = matrix_mul(base, base)
        exponent //= 2
    return result


PELL_STEP = ((3, 8), (1, 3))


def selector_state(k: int) -> tuple[int, int]:
    """Return (x_k,q_k), with x^2-8q^2=1 and (x_0,q_0)=(1,0)."""
    p = matrix_pow(PELL_STEP, k)
    return p[0][0], p[1][0]


def rational_series(
    numerator: tuple[int, ...], denominator: tuple[int, ...], count: int
) -> tuple[Fraction, ...]:
    require(denominator[0] == 1, "series denominator must have unit constant")
    values: list[Fraction] = []
    for n in range(count):
        rhs = Fraction(numerator[n] if n < len(numerator) else 0)
        for j in range(1, min(n + 1, len(denominator))):
            rhs -= denominator[j] * values[n - j]
        values.append(rhs)
    return tuple(values)


def circle_distance(speed: int, t: Fraction) -> Fraction:
    residue = (speed * t.numerator) % t.denominator
    residue = min(residue, t.denominator - residue)
    return Fraction(residue, t.denominator)


def margin(speeds: tuple[int, ...], t: Fraction) -> Fraction:
    return min(circle_distance(speed, t) for speed in speeds)


def pair_candidates(speeds: tuple[int, ...]) -> set[Fraction]:
    """THM-757/THM-1002 peak candidates on [0,1/2]."""
    speeds = tuple(sorted(set(speeds)))
    candidates: set[Fraction] = {Fraction(0), Fraction(1, 2)}
    denominators = {2 * speed for speed in speeds}
    for a, b in combinations(speeds, 2):
        denominators.add(a + b)
        denominators.add(abs(a - b))
    denominators.discard(0)
    for denominator in denominators:
        for numerator in range(1, denominator // 2 + 1):
            candidates.add(Fraction(numerator, denominator))
    return candidates


def exact_max_pair(speeds: tuple[int, ...]) -> tuple[Fraction, Fraction]:
    best = Fraction(-1)
    best_t = Fraction(0)
    for t in pair_candidates(speeds):
        value = margin(speeds, t)
        if value > best or (value == best and t < best_t):
            best, best_t = value, t
    return best, best_t


def affine_piece(speed: int, midpoint: Fraction) -> tuple[int, int]:
    x = speed * midpoint
    floor_x = x.numerator // x.denominator
    local = x - floor_x
    return (speed, -floor_x) if local <= Fraction(1, 2) else (-speed, floor_x + 1)


def exact_max_cells(speeds: tuple[int, ...]) -> tuple[Fraction, Fraction]:
    """Independent lower-envelope enumeration on [0,1]."""
    breaks = {Fraction(0), Fraction(1)}
    for speed in speeds:
        for numerator in range(2 * speed + 1):
            breaks.add(Fraction(numerator, 2 * speed))
    ordered = sorted(breaks)
    candidates = set(ordered)
    for left, right in zip(ordered, ordered[1:]):
        midpoint = (left + right) / 2
        lines = [affine_piece(speed, midpoint) for speed in speeds]
        for (slope_a, shift_a), (slope_b, shift_b) in combinations(lines, 2):
            if slope_a == slope_b:
                continue
            t = Fraction(shift_b - shift_a, slope_a - slope_b)
            if left <= t <= right:
                candidates.add(t)
    best = max((margin(speeds, t), -t, t) for t in candidates)
    return best[0], best[2]


def clock_row(j: int, n: int) -> tuple[int, ...]:
    a = n + 1
    return tuple(sorted([a * r for r in range(1, 14) if r != j]
                        + [(j + 1) * a - 1]))


def clock_numerator(j: int) -> int:
    if j % 2 == 1:
        return 2
    if j in (8, 10):
        return 3
    if j == 12:
        return 5
    raise ValueError(j)


def qsqrt_add(a: tuple[Fraction, Fraction], b: tuple[Fraction, Fraction]) -> tuple[Fraction, Fraction]:
    return a[0] + b[0], a[1] + b[1]


def qsqrt_mul(a: tuple[Fraction, Fraction], b: tuple[Fraction, Fraction]) -> tuple[Fraction, Fraction]:
    return a[0] * b[0] + 2 * a[1] * b[1], a[0] * b[1] + a[1] * b[0]


def qsqrt_pow(a: tuple[Fraction, Fraction], exponent: int) -> tuple[Fraction, Fraction]:
    result = (Fraction(1), Fraction(0))
    base = a
    while exponent:
        if exponent & 1:
            result = qsqrt_mul(result, base)
        base = qsqrt_mul(base, base)
        exponent //= 2
    return result


def qsqrt_half(a: tuple[Fraction, Fraction]) -> tuple[Fraction, Fraction]:
    return a[0] / 2, a[1] / 2


def pell_numbers(count: int) -> list[int]:
    values = [0, 1]
    while len(values) < count:
        values.append(2 * values[-1] + values[-2])
    return values


def main() -> None:
    digest_rows: list[str] = []
    print("THM-3335 SQUARE-TRIANGULAR / PELL / MARKOV / PYTHAGOREAN EXACT AUDIT")
    print("universe: direct triangular scan n<=100000; primitive Euclid m<=500;")
    print("          selector depths k<=30; LRC ladder j=7..13,n<=200;")
    print("          independent exact max replay j=2..13 at n=1")

    # 1. Square-triangular selector and sequence compiler.
    states = [selector_state(k) for k in range(31)]
    rows = []
    for k, (x, q) in enumerate(states):
        n = (x - 1) // 2
        y = q * q
        b_leg = 4 * y
        c = b_leg + 1
        ew_order = b_leg + 2
        require(x * x - 8 * q * q == 1, f"Pell equation k={k}")
        require(triangular(n) == y, f"square-triangular k={k}")
        require(phi(n + 1, n) == (x, b_leg, c), f"consecutive phi k={k}")
        require(x * x + b_leg * b_leg == c * c, f"Pythagoras k={k}")
        require(2 * ew_order - 3 == x * x, f"EW gate k={k}")
        rows.append((n, q, x, 2 * q, b_leg, c, ew_order))
    for k in range(29):
        x0, q0 = states[k]
        x1, q1 = states[k + 1]
        require((x1, q1) == (3 * x0 + 8 * q0, x0 + 3 * q0), "matrix recurrence")
    for sequence, shift in (
        ([row[0] for row in rows], 2),
        ([row[4] for row in rows], 8),
        ([row[5] for row in rows], -24),
        ([row[6] for row in rows], -56),
    ):
        coefficient = 6 if shift == 2 else 34
        for k in range(29):
            require(sequence[k + 2] == coefficient * sequence[k + 1] - sequence[k] + shift,
                    f"affine recurrence {coefficient},{shift},k={k}")

    direct_hits = []
    for n in range(100001):
        q = square_root(triangular(n))
        if q is not None:
            direct_hits.append((n, q))
    recurrence_hits = [(row[0], row[1]) for row in rows if row[0] <= 100000]
    require(direct_hits == recurrence_hits, "direct scan versus Pell recurrence")

    d6 = (1, -6, 1)
    d6_affine = (1, -7, 7, -1)
    d34_affine = (1, -35, 35, -1)
    require(rational_series((1, -3), d6, 25) == tuple(Fraction(row[2]) for row in rows[:25]), "X(z)")
    require(rational_series((0, 1), d6, 25) == tuple(Fraction(row[1]) for row in rows[:25]), "Q(z)")
    require(rational_series((0, 1, 1), d6_affine, 25) == tuple(Fraction(row[0]) for row in rows[:25]), "N(z)")
    require(rational_series((0, 1, 1), d34_affine, 25) == tuple(Fraction(row[1] ** 2) for row in rows[:25]), "Y(z)")
    require(rational_series((1, -30, 5), d34_affine, 25) == tuple(Fraction(row[5]) for row in rows[:25]), "C(z)")
    require(rational_series((2, -64, 6), d34_affine, 25) == tuple(Fraction(row[6]) for row in rows[:25]), "W(z)")

    print("\n[1] selector/compiler")
    print("  first positive rows (n,q,x,even_root,even_leg,hypotenuse,EW_candidate):")
    for k, row in enumerate(rows[1:6], start=1):
        print(f"    k={k}: {row}")
    print(f"  direct square-triangular hits through 100000: {direct_hits}")
    print("  2x2 matrix, affine recurrences, and six rational generating functions: PASS")
    digest_rows.append(f"hits={direct_hits}")

    # 2. Ambient primitive square-even-leg classification.
    square_even_count = 0
    consecutive_square_even = []
    for m in range(2, 501):
        for n in range(1, m):
            if gcd(m, n) != 1 or (m - n) % 2 == 0:
                continue
            even_leg = 2 * m * n
            if square_root(even_leg) is None:
                continue
            square_even_count += 1
            sheet_a = square_root(m) is not None and n % 2 == 0 and square_root(n // 2) is not None
            sheet_b = m % 2 == 0 and square_root(m // 2) is not None and square_root(n) is not None
            require(sheet_a ^ sheet_b, f"valuation sheet m={m},n={n}")
            if m - n == 1:
                consecutive_square_even.append((m, n))
    require(square_even_count == 134, "primitive square-even count")
    require(consecutive_square_even == [(2, 1), (9, 8), (50, 49), (289, 288)],
            "consecutive square-even intersection")
    hostile_gap = phi(3, 2)       # c-b=1, but b=12 is not square.
    hostile_square = phi(8, 1)    # b=16 square, but c-b=49.
    require(hostile_gap == (5, 12, 13), "gap-one hostile")
    require(hostile_square == (63, 16, 65), "square-leg hostile")
    require(15 * 15 + 36 * 36 == 39 * 39 and gcd(gcd(15, 36), 39) == 3,
            "nonprimitive scale hostile")
    print("\n[2] primitive square-even-leg boundary")
    print(f"  primitive Euclid hits m<=500: {square_even_count}; valuation sheets: PASS")
    print(f"  consecutive intersection: {consecutive_square_even}")
    print("  hostiles: (5,12,13) gap-one/non-square; (63,16,65) square/non-gap-one;")
    print("            (15,36,39) shows scaling can create a square leg but loses unit gap")
    digest_rows.append(f"square_even={square_even_count}:{consecutive_square_even}")

    # 3. Fixed-two Markov branch and corrected Chebyshev--Pell evaluator.
    pell = pell_numbers(2 * 31 + 2)
    markov_rows = []
    for k in range(1, 15):
        x, q = states[k]
        s = 2 * q
        a, b = x - s, x + s
        require((a, b) == (pell[2 * k - 1], pell[2 * k + 1]), f"Pell odd coordinates k={k}")
        require(a * a + b * b + 4 == 6 * a * b, f"fixed-two Markov k={k}")
        require(a * b == s * s + 1, f"Cassini/product k={k}")
        require((x, s * s, a * b) == phi((x + 1) // 2, (x - 1) // 2), f"Markov-to-Phi k={k}")
        if k < 14:
            next_x, next_q = states[k + 1]
            require((b, 6 * b - a) == (next_x - 2 * next_q, next_x + 2 * next_q),
                    f"Markov mutation k={k}")
            require((3 * x + 4 * s, 2 * x + 3 * s) == (next_x, 2 * next_q),
                    f"half-Hadamard conjugacy k={k}")
        markov_rows.append((a, b, x, s, a * b))

    for k in range(11):
        plus = qsqrt_pow((Fraction(1), Fraction(1)), 2 * k)
        minus = qsqrt_pow((Fraction(-1), Fraction(1)), 2 * k)
        e_even = qsqrt_half(qsqrt_add(plus, minus))
        o_even = qsqrt_half(qsqrt_add(plus, (-minus[0], -minus[1])))
        require(e_even == (Fraction(states[k][0]), Fraction(0)), f"E_2k evaluator k={k}")
        require(o_even == (Fraction(0), Fraction(2 * states[k][1])), f"O_2k evaluator k={k}")
        plus_odd = qsqrt_pow((Fraction(1), Fraction(1)), 2 * k + 1)
        minus_odd = qsqrt_pow((Fraction(-1), Fraction(1)), 2 * k + 1)
        e_odd = qsqrt_half(qsqrt_add(plus_odd, minus_odd))
        require(e_odd == (Fraction(0), Fraction(pell[2 * k + 1])), f"E_2k+1 evaluator k={k}")

    # Corrected THM-1880 recursion checked in Q(sqrt(2)).
    e_prev, o_prev = (Fraction(1), Fraction(0)), (Fraction(0), Fraction(0))
    root2 = (Fraction(0), Fraction(1))
    for degree in range(1, 16):
        e_now = qsqrt_add(qsqrt_mul(root2, e_prev), o_prev)
        o_now = qsqrt_add(e_prev, qsqrt_mul(root2, o_prev))
        direct_plus = qsqrt_pow((Fraction(1), Fraction(1)), degree)
        direct_minus = qsqrt_pow((Fraction(-1), Fraction(1)), degree)
        direct_e = qsqrt_half(qsqrt_add(direct_plus, direct_minus))
        direct_o = qsqrt_half(qsqrt_add(direct_plus, (-direct_minus[0], -direct_minus[1])))
        require((e_now, o_now) == (direct_e, direct_o), f"corrected E/O recursion n={degree}")
        e_prev, o_prev = e_now, o_now

    # Full silver sequence alternates Cassini; even decimation has constant sign.
    for n in range(1, 30):
        require(pell[n - 1] * pell[n + 1] - pell[n] * pell[n] == (-1) ** n,
                f"Pell Cassini n={n}")
    for k in range(2, 29):
        q0, q1, q2 = states[k - 1][1], states[k][1], states[k + 1][1]
        require(q0 * q2 - q1 * q1 == -1, f"decimated Cassini k={k}")

    print("\n[3] fixed-two Markov and parity-typed evaluator")
    print(f"  first branches (a,b,x,s,ab): {markov_rows[:4]}")
    print("  Markov mutation <-> Pell-unit matrix; corrected E/O evaluator: PASS")
    print("  full silver Cassini alternates; even decimation has constant -1: PASS")
    digest_rows.append(f"markov={markov_rows[:5]}")

    # 4. Cannonball/even-odd square split and bounded rare intersections.
    cannonball_hits = []
    for cap in range(1, 100001):
        total = cap * (cap + 1) * (2 * cap + 1) // 6
        root = square_root(total)
        if root is not None:
            cannonball_hits.append((cap, root))
    require(cannonball_hits == [(1, 1), (24, 70)], "cannonball scan")
    even_sum = sum((2 * r) ** 2 for r in range(1, 13))
    odd_sum = sum((2 * r - 1) ** 2 for r in range(1, 13))
    require((even_sum, odd_sum, even_sum + odd_sum) == (2600, 2300, 4900), "even/odd split")
    require(markov_rows[2] == (29, 169, 99, 70, 4901), "cannonball Markov address")
    guy_hits = []
    for k in range(1, 31):
        q = states[k][1]
        j = (isqrt(8 * q + 1) - 1) // 2
        if triangular(j) == q:
            guy_hits.append((k, q, q * q, j))
    require(guy_hits == [(1, 1, 1, 1), (2, 6, 36, 3)], "bounded Guy intersection")
    print("\n[4] exceptional sums-of-squares intersections")
    print(f"  cannonball hits N<=100000: {cannonball_hits}")
    print("  70^2 = sum_{r<=12}(2r)^2 + sum_{r<=12}(2r-1)^2 = 2600+2300")
    print("  (29,169) -> (99,4900,4901); bounded Guy-root hits k<=30: " + str(guy_hits))
    digest_rows.append(f"cannon={cannonball_hits}:{even_sum}:{odd_sum}:{guy_hits}")

    # 5. Tournament order maps and their sharp separation.
    ew_candidates = [row[6] for row in rows]
    require(ew_candidates[:5] == [2, 6, 146, 4902, 166466], "EW candidate prefix")
    square_arc_ew_hits = []
    for n in range(1, 100001):
        q = square_root(triangular(n))
        r = square_root(2 * n - 1)
        if q is not None and r is not None:
            square_arc_ew_hits.append((n + 1, n, q, r))
    require(square_arc_ew_hits == [(2, 1, 1, 1)], "square-arc/EW intersection scan")
    for odd_x in range(1, 100, 2):
        order = (odd_x * odd_x + 3) // 2
        triple = (odd_x, order - 2, order - 1)
        require(triple == phi((odd_x + 1) // 2, (odd_x - 1) // 2), "general EW spine")
        require(triple[0] ** 2 + triple[1] ** 2 == triple[2] ** 2, "general EW Pythagoras")
    require((13, 84, 85) == phi(7, 6), "N=86 arithmetic hostile")
    print("\n[5] tournament gates")
    print(f"  square-leg skew-EW arithmetic candidates: {ew_candidates[:5]} ...")
    print(f"  square universal-arc-count order AND skew-EW gate through 100000: {square_arc_ew_hits}")
    print("  proof factorization leaves only order 2 globally; N=86 -> (13,84,85) is the")
    print("  necessity-only control; its restricted symmetric-ansatz failure does not decide sufficiency")
    digest_rows.append(f"ew={ew_candidates[:8]}:arc={square_arc_ew_hits}")

    # 6. Seven-clock LRC ladder and exact failure boundary.
    exact_ladder = []
    for j in range(7, 14):
        kappa = clock_numerator(j)
        require(gcd(kappa, j) == 1 and 2 <= kappa <= j // 2, f"clock numerator j={j}")
        for n in range(1, 201):
            a = n + 1
            row = clock_row(j, n)
            phase = Fraction(kappa, j * a)
            require(len(row) == 13 and reduce(gcd, row) == 1, f"primitive clock row j={j},n={n}")
            require(margin(row, phase) == Fraction(1, j), f"clock margin j={j},n={n}")
        pair_value, pair_time = exact_max_pair(clock_row(j, 1))
        cell_value, cell_time = exact_max_cells(clock_row(j, 1))
        require(pair_value == cell_value == Fraction(1, j), f"independent exact M j={j}")
        exact_ladder.append((j, pair_value, pair_time, cell_time))

    boundary_values = {}
    expected_boundary = {
        2: Fraction(2, 17), 3: Fraction(2, 17),
        4: Fraction(2, 19), 5: Fraction(2, 19), 6: Fraction(2, 23),
    }
    for j in range(2, 7):
        pair_value, pair_time = exact_max_pair(clock_row(j, 1))
        cell_value, _ = exact_max_cells(clock_row(j, 1))
        require(pair_value == cell_value == expected_boundary[j], f"lower-rung hostile j={j}")
        boundary_values[j] = (pair_value, pair_time)

    def owner_data(j: int, n: int) -> tuple[int, int, str, int]:
        c_norm = 2 * n * n + 2 * n + 1
        determinants = [((n + 1) * (1 if i == j else 0) - n * i, i)
                        for i in range(1, 14)]
        d_max = max(abs(value) for value, _ in determinants)
        active = [(value, i) for value, i in determinants if abs(value) == d_max]
        require(len(active) == 1, f"unique determinant owner j={j},n={n}")
        value, index = active[0]
        owner = ("+" if value > 0 else "-") + f"c_{index}"
        return c_norm, d_max, owner, c_norm - 91 * d_max

    same_address = {owner_data(j, 8) for j in range(7, 13)}
    require(len(same_address) == 1, "same compressed Gaussian/Kelvin scalar tuple j=7..12")
    require([Fraction(1, j) for j in range(7, 13)]
            == [exact_max_pair(clock_row(j, 8))[0] for j in range(7, 13)],
            "different exact phase maxima at same address")
    cutoff_13n = next(n for n in range(1, 2000) if owner_data(12, n)[3] >= 0)
    cutoff_12n = next(n for n in range(1, 2000) if owner_data(13, n)[3] >= 0)
    require((cutoff_13n, cutoff_12n) == (591, 545), "Kelvin cutoffs")
    residual_pell_13n = [row[0] for row in rows[1:] if row[0] < cutoff_13n]
    residual_pell_12n = [row[0] for row in rows[1:] if row[0] < cutoff_12n]
    require(residual_pell_13n == residual_pell_12n == [1, 8, 49, 288], "four Pell false positives")

    print("\n[6] seven-clock LRC benchmark")
    print("  exact ladder (j,M,pair_arg,cell_arg):")
    for item in exact_ladder:
        print(f"    {item}")
    print(f"  lower-rung n=1 hostiles: {boundary_values}")
    print(f"  j=7..12 same (C,D,owner,F) at n=8: {next(iter(same_address))}")
    print("  but exact M values are 1/7,...,1/12; labelled clock placement is essential")
    print("  Kelvin safe cutoffs j<=12/j=13: 591/545; residual Pell depths: 1,8,49,288")
    digest_rows.append(f"ladder={exact_ladder}:boundary={boundary_values}:address={same_address}")

    # 7. Quotient and lacunarity stopping controls.
    u, up, v = (8, 1), (7, 4), (1, 0)
    require(u[0] * u[0] + u[1] * u[1] == up[0] * up[0] + up[1] * up[1] == 65,
            "norm-65 collision")
    require((abs(u[0] * v[1] - u[1] * v[0]),
             abs(up[0] * v[1] - up[1] * v[0])) == (1, 4), "represented norm determinants")
    require((lorentz(phi(*u), phi(*v)), lorentz(phi(*up), phi(*v))) == (2, 32),
            "represented norm Lorentz shells")
    for coordinate in range(7):
        sequence = [row[coordinate] for row in rows[1:14]]
        if coordinate == 5 or coordinate == 6:
            # Hypotenuse/order also grow fast after the first positive row.
            pass
        require(all(3 * sequence[k + 1] >= 7 * sequence[k] for k in range(len(sequence) - 1)),
                f"7/3 lacunarity coordinate={coordinate}")
    for coordinate in (0, 1):
        sequence = [row[coordinate] for row in markov_rows]
        require(all(3 * sequence[k + 1] >= 7 * sequence[k] for k in range(len(sequence) - 1)),
                f"7/3 lacunarity Markov coordinate={coordinate}")
    print("\n[7] stopping boundaries")
    print("  norm 65 collision: represented determinants 1/4, Lorentz shells 2/32;")
    print("  a torus-invariant norm quotient cannot retain the Pythagorean/Farey selector")
    print("  every positive selector coordinate is >=7/3-lacunary: 13-term packets")
    print("  from one selector stream are LRC-safe by THM-928(A'); hard use must be mixed")
    digest_rows.append("norm65=1,4:2,32;lacunary=PASS")

    semantic_digest = sha256("\n".join(digest_rows).encode()).hexdigest()
    print(f"\nsemantic_digest_sha256={semantic_digest}")
    print("ALL EXACT CHECKS PASS")


if __name__ == "__main__":
    main()
