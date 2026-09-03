#!/usr/bin/env python3
"""Dependency-free exact checks for the binomial-diagonal hierarchy.

This file intentionally rebuilds only from the source-normal monomial law

  x^a u^b p^c y^e
    = x^(a+2b+e) t^(b+c+2e) (1+x^2 t)^(c+e).

It does not import any repository computation.
"""

from math import comb, factorial
import sys


sys.stdout.reconfigure(newline="\n")

CHECKS = 0


def check(condition, label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError("check failed: " + label)


def ceil_div(a, b):
    return -((-a) // b)


def polynomial_binom(z, q):
    """Polynomial continuation z(z-1)...(z-q+1)/q! for integer z."""
    if q < 0:
        return 0
    value = 1
    for j in range(q):
        value *= z - j
    return value // factorial(q)


def weights(ell, q, m):
    s = ceil_div(ell, 2)
    return [(-1) ** (n - s) * comb(m + q - n, q)
            for n in range(s, m + 1)]


def generator_value(ell, q, m, a, b, c, e):
    """Evaluate L_(m,ell,q) on one spanning-law generator."""
    s = ceil_div(ell, 2)
    n0 = b + c + 2 * e
    r0 = a + 2 * b + e
    degree = c + e
    value = 0
    for k in range(degree + 1):
        n = n0 + k
        r = r0 + 2 * k
        if s <= n <= m and r == 2 * n - ell:
            value += ((-1) ** (n - s)
                      * comb(m + q - n, q)
                      * comb(degree, k))
    return value


def finite_difference_checks():
    count = 0
    for degree in range(0, 17):
        for q in range(0, 10):
            for x in range(-10, 31):
                lhs = sum(((-1) ** k) * comb(degree, k)
                          * polynomial_binom(x - k, q)
                          for k in range(degree + 1))
                rhs = (polynomial_binom(x - degree, q - degree)
                       if degree <= q else 0)
                check(lhs == rhs,
                      f"finite difference N={degree},q={q},X={x}")
                count += 1
    return count


def partial_boundary_checks():
    count = 0
    for degree in range(1, 24):
        for q in range(degree):
            for missing in range(q + 1, degree + 1):
                upto = degree - missing
                lhs = sum(((-1) ** k) * comb(degree, k)
                          * comb(degree + q - missing - k, q)
                          for k in range(upto + 1))
                rhs = ((-1) ** upto
                       * comb(degree - q - 1, missing - q - 1))
                check(lhs == rhs,
                      f"partial identity N={degree},q={q},R={missing}")
                count += 1
    return count


def check_all_visible_generators_zero(ell, q, m, depth):
    count = 0
    for a in range(depth + 1):
        for b in range(depth - a + 1):
            # Only n0=b+c+2e <= m can meet the retained rows.
            for e in range((m - b) // 2 + 1):
                for c in range(m - b - 2 * e + 1):
                    value = generator_value(ell, q, m, a, b, c, e)
                    check(value == 0,
                          (f"annihilation ell={ell},q={q},m={m},d={depth},"
                           f"g=({a},{b},{c},{e}),v={value}"))
                    count += 1
    return count


def unit_hostile(ell, q, m):
    threshold = m + q - ell
    depth = threshold + 1
    terminal = m + q + 1
    e = terminal % 2
    c = (terminal - 3 * e) // 2
    a = depth
    b = 0
    return (a, b, c, e)


def main_hierarchy_checks():
    generator_checks = 0
    hostile_checks = 0
    cases = 0
    for ell in range(1, 25):
        qmax = ceil_div(ell, 3) - 1
        for q in range(qmax + 1):
            for extra_rows in range(6):
                m = ell - q + extra_rows
                threshold = m + q - ell
                generator_checks += check_all_visible_generators_zero(
                    ell, q, m, threshold)
                hostile = unit_hostile(ell, q, m)
                value = generator_value(ell, q, m, *hostile)
                expected = (-1) ** (m - ceil_div(ell, 2))
                check(sum(hostile[:2]) == threshold + 1,
                      f"hostile depth ell={ell},q={q},m={m}")
                check(value == expected and abs(value) == 1,
                      f"unit hostile ell={ell},q={q},m={m},v={value}")
                hostile_checks += 2
                cases += 1
    return cases, generator_checks, hostile_checks


def p0_annihilates(ell, q, m):
    for e in range(m // 2 + 1):
        for c in range(m - 2 * e + 1):
            if generator_value(ell, q, m, 0, 0, c, e) != 0:
                return False
    return True


def p0_phase_checks():
    count = 0
    cases = 0
    # For ell >= 2 and m >= ceil(ell/2), the two inequalities are not
    # merely sufficient: together they exactly characterize P_0 vanishing.
    for ell in range(2, 25):
        s = ceil_div(ell, 2)
        for q in range(0, ell + 5):
            for m in range(s, ell + 6):
                actual = p0_annihilates(ell, q, m)
                predicted = q < ceil_div(ell, 3) and m >= ell - q
                check(actual == predicted,
                      f"P0 phase ell={ell},q={q},m={m}")
                count += 1
                cases += 1
    return cases, count


def sharp_boundary_checks():
    count = 0
    for ell in range(2, 41):
        s = ceil_div(ell, 2)

        # First forbidden finite-difference order.  A minimal-degree exact
        # diagonal generator starts at m=ell-q and gives a unit value.
        q = ceil_div(ell, 3)
        m = ell - q
        degree = q
        e = ell - 2 * degree
        c = 3 * degree - ell
        value = generator_value(ell, q, m, 0, 0, c, e)
        check(c >= 0 and e >= 0 and c + e == degree,
              f"order-bound representation ell={ell}")
        check(value == (-1) ** (m - s) and abs(value) == 1,
              f"order-bound unit ell={ell},v={value}")
        count += 2

        # One row below the allowed range.  The maximal-degree exact
        # diagonal generator starts at s and terminates at m+q+1.
        for q in range(ceil_div(ell, 3)):
            m = ell - q - 1
            e = ell % 2
            c = (ell - 3 * e) // 2
            value = generator_value(ell, q, m, 0, 0, c, e)
            check(m >= s, f"row-bound retained start ell={ell},q={q}")
            check(value == (-1) ** (m - s) and abs(value) == 1,
                  f"row-bound unit ell={ell},q={q},v={value}")
            count += 2
    return count


def p0_parameter_bound_examples():
    """Concrete unit hostiles for the two sharp P_0 parameter boundaries."""
    ell = 8
    s = ceil_div(ell, 2)

    # First forbidden order q=ceil(ell/3), at the last-row normalization
    # m=ell-q.  Here p*y^2 has (a,b,c,e)=(0,0,1,2).
    order_q = ceil_div(ell, 3)
    order_m = ell - order_q
    order_hostile = (0, 0, 1, 2)
    order_value = generator_value(
        ell, order_q, order_m, *order_hostile)
    check(order_q == ceil_div(ell, 3) and order_m == ell - order_q,
          "P0 order-bound example parameters")
    check(sum(order_hostile[:2]) == 0,
          "P0 order-bound example depth")
    check(order_value == (-1) ** (order_m - s) and abs(order_value) == 1,
          "P0 order-bound example unit")

    # First forbidden row m=ell-q-1 at the largest allowed q.  Here p^4
    # has (a,b,c,e)=(0,0,4,0).
    row_q = ceil_div(ell, 3) - 1
    row_m = ell - row_q - 1
    row_hostile = (0, 0, 4, 0)
    row_value = generator_value(ell, row_q, row_m, *row_hostile)
    check(row_q < ceil_div(ell, 3) and row_m == ell - row_q - 1,
          "P0 row-bound example parameters")
    check(sum(row_hostile[:2]) == 0,
          "P0 row-bound example depth")
    check(row_value == (-1) ** (row_m - s) and abs(row_value) == 1,
          "P0 row-bound example unit")

    return ((ell, order_q, order_m, order_hostile, order_value),
            (ell, row_q, row_m, row_hostile, row_value))


def specialization_checks():
    expected = {
        (8, 2, 8): [15, -10, 6, -3, 1],
        (10, 3, 9): [35, -20, 10, -4, 1],
        (10, 3, 10): [56, -35, 20, -10, 4, -1],
    }
    labels = {
        (8, 2, 8): "row8_P2_triangular",
        (10, 3, 9): "row9_P2_tetrahedral",
        (10, 3, 10): "row10_P3_tetrahedral",
    }
    lines = []
    for key, want in expected.items():
        got = weights(*key)
        check(got == want, f"specialization {key}: {got}")
        hostile = unit_hostile(*key)
        value = generator_value(*key, *hostile)
        lines.append((labels[key], got, hostile, value))
    return lines


def main():
    fd = finite_difference_checks()
    partial = partial_boundary_checks()
    cases, generators, hostiles = main_hierarchy_checks()
    p0_cases, p0 = p0_phase_checks()
    sharp = sharp_boundary_checks()
    order_example, row_example = p0_parameter_bound_examples()
    specs = specialization_checks()

    print("GENERAL BINOMIAL-DIAGONAL HIERARCHY EXACT VERIFIER")
    print(f"finite_difference_identity_checks={fd}")
    print(f"partial_boundary_identity_checks={partial}")
    print(f"main_parameter_cases={cases}")
    print(f"visible_generator_annihilation_checks={generators}")
    print(f"unit_hostile_checks={hostiles}")
    print(f"P0_phase_cases={p0_cases}")
    print(f"P0_phase_checks={p0}")
    print(f"sharp_boundary_checks={sharp}")
    ell, q, m, hostile, value = order_example
    print("P0_order_bound_unit_hostile="
          f"(ell,q,m;a,b,c,e)=({ell},{q},{m};"
          f"{hostile[0]},{hostile[1]},{hostile[2]},{hostile[3]}); "
          f"value={value}")
    ell, q, m, hostile, value = row_example
    print("P0_row_bound_unit_hostile="
          f"(ell,q,m;a,b,c,e)=({ell},{q},{m};"
          f"{hostile[0]},{hostile[1]},{hostile[2]},{hostile[3]}); "
          f"value={value}")
    for label, stencil, hostile, value in specs:
        print(f"{label}_stencil={stencil}")
        print(f"{label}_unit_hostile={hostile}; value={value}")
    print(f"checks={CHECKS}")
    print("PASS")


if __name__ == "__main__":
    main()

