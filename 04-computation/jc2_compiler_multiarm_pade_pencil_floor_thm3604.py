#!/usr/bin/env python3
"""Exact controls for provisional THM-3604.

The universal pencil floors are proof-driven consequences of THM-3589,
THM-3600, and THM-3550.  This companion checks the integer invoices, sharp
global-function controls, nodal target changes, and the A13 specialization
with exact SymPy arithmetic.  Finite rows are controls, not extrapolations.
"""

from itertools import product

import sympy as sp


CHECKS = 0


def require(label, condition):
    """Record one explicit gate and raise on failure."""
    global CHECKS
    if condition != True:
        raise RuntimeError(f"FAIL {label}: {condition}")
    CHECKS += 1


def degree_in(expr, variable):
    """Exact polynomial degree in one selected variable."""
    return sp.Poly(sp.expand(expr), variable).degree()


def c_order(expr, c):
    """Lowest c exponent of a nonzero polynomial."""
    polynomial = sp.Poly(sp.expand(expr), c)
    return min(monomial[0] for monomial, coefficient in polynomial.terms() if coefficient)


def is_composite(integer):
    """Elementary exact composite predicate for the small control range."""
    if integer < 4:
        return False
    return any(integer % divisor == 0 for divisor in range(2, int(integer**0.5) + 1))


print("THM-3604 exact companion -- provisional compiler multiarm pencil floor")
print("scope=finite exact controls; every-line, Pade, and planar composite inputs are proof-driven")


print("SECTION integer multiarm and leading-c invoices")
invoice_rows = 0
for n in range(2, 7):
    for h in range(2, 9):
        member_a_floor = 2 * h
        basis_a_floor = 3 * h
        member_order_floor = 2 * n * (h - 1)
        basis_order_floor = 3 * n * (h - 1)
        total_unit = h + n * (h - 1)

        require(
            f"member exact tariff n={n} h={h}",
            n * (member_a_floor - member_a_floor // h) == member_order_floor,
        )
        require(
            f"basis exact tariff n={n} h={h}",
            n * (basis_a_floor - basis_a_floor // h) == basis_order_floor,
        )
        require(
            f"member total connection n={n} h={h}",
            member_a_floor + member_order_floor == 2 * total_unit,
        )
        require(
            f"basis total connection n={n} h={h}",
            basis_a_floor + basis_order_floor == 3 * total_unit,
        )
        require(
            f"member requested shadow n={n} h={h}",
            2 * total_unit >= max(2 * h, 6),
        )
        require(
            f"basis requested shadow n={n} h={h}",
            3 * total_unit >= max(3 * h, 8),
        )

        for d in range(member_a_floor, 4 * h + 1):
            tariff = n * (d - d // h)
            require(
                f"member monotone tariff n={n} h={h} d={d}",
                tariff >= member_order_floor,
            )

        for d in range(basis_a_floor, 5 * h + 1):
            tariff = n * (d - d // h)
            require(
                f"basis monotone tariff n={n} h={h} d={d}",
                tariff >= basis_order_floor,
            )
        invoice_rows += 1

print(f"PASS invoice_rows={invoice_rows} n=2..6 h=2..8")


print("SECTION sharp e^r Hermite-Pade controls")
b, c, a = sp.symbols("b c a")
power_rows = 0
for n in range(2, 5):
    for h in range(2, 7):
        E = sp.expand(sp.prod(b - root for root in range(1, h)))
        E0 = E.subs(b, 0)
        e_chart = sp.expand(a * E.subs(b, c**n * a))

        require(f"nonzero central slope n={n} h={h}", E0 != 0)
        require(f"e chart degree n={n} h={h}", degree_in(e_chart, a) == h)
        require(
            f"e central arm degree n={n} h={h}",
            degree_in(e_chart.subs(c, 0), a) == 1,
        )

        for r in range(1, 5):
            expression = sp.expand(e_chart**r)
            d_a = degree_in(expression, a)
            arm_degree = degree_in(expression.subs(c, 0), a)
            leading_a = sp.Poly(expression, a).LC()
            order = c_order(leading_a, c)
            tariff = n * (d_a - d_a // h)

            require(f"power a-degree n={n} h={h} r={r}", d_a == r * h)
            require(f"power arm degree n={n} h={h} r={r}", arm_degree == r)
            require(
                f"power leading order n={n} h={h} r={r}",
                order == n * r * (h - 1),
            )
            require(f"power exact tariff n={n} h={h} r={r}", order == tariff)
            require(
                f"power exact total n={n} h={h} r={r}",
                sp.Poly(expression, a, c).total_degree()
                == r * (h + n * (h - 1)),
            )
            power_rows += 1

print(f"PASS sharp_power_rows={power_rows} equality=d_a/(arm degree)=h")


print("SECTION nodal (2,3) arm boundary and target-GL2 hostiles")
s = sp.symbols("s")
nodal_two = s**2 - 1
nodal_three = s * (s**2 - 1)
require("nodal first collision", nodal_two.subs(s, 1) == nodal_two.subs(s, -1))
require("nodal second collision", nodal_three.subs(s, 1) == nodal_three.subs(s, -1))
require("nodal collision parameters distinct", sp.Integer(1) != sp.Integer(-1))
require("nodal degree two", degree_in(nodal_two, s) == 2)
require("nodal degree three", degree_in(nodal_three, s) == 3)
require(
    "nodal immersion",
    sp.degree(sp.gcd(sp.diff(nodal_two, s), sp.diff(nodal_three, s)), s) == 0,
)

gl2_rows = 0
for alpha, beta, gamma, delta in product((-1, 0, 1), repeat=4):
    determinant = alpha * delta - beta * gamma
    if determinant == 0:
        continue
    first = sp.expand(alpha * nodal_two + beta * nodal_three)
    second = sp.expand(gamma * nodal_two + delta * nodal_three)
    first_degree = degree_in(first, s)
    second_degree = degree_in(second, s)

    require(
        f"GL2 every-line floor matrix={alpha,beta,gamma,delta}",
        min(first_degree, second_degree) >= 2,
    )
    require(
        f"GL2 basis height matrix={alpha,beta,gamma,delta}",
        max(first_degree, second_degree) >= 3,
    )
    require(
        f"GL2 immersion matrix={alpha,beta,gamma,delta}",
        sp.degree(sp.gcd(sp.diff(first, s), sp.diff(second, s)), s) == 0,
    )
    gl2_rows += 1

print(f"PASS nodal_GL2_rows={gl2_rows} arm_degrees=2,3 immersed_noninjective")


print("SECTION global nodal lifts and zero-Jacobian warning")
nodal_chart_rows = 0
for n in range(2, 5):
    for h in range(2, 6):
        E = sp.expand(sp.prod(b - root for root in range(1, h)))
        E0 = E.subs(b, 0)
        e_chart = sp.expand(a * E.subs(b, c**n * a))
        first = sp.expand(e_chart**2 - 1)
        second = sp.expand(e_chart**3 - e_chart)
        first_arm = sp.expand(first.subs(c, 0))
        second_arm = sp.expand(second.subs(c, 0))

        require(f"lift first degree n={n} h={h}", degree_in(first, a) == 2 * h)
        require(f"lift second degree n={n} h={h}", degree_in(second, a) == 3 * h)
        require(
            f"lift first order n={n} h={h}",
            c_order(sp.Poly(first, a).LC(), c) == 2 * n * (h - 1),
        )
        require(
            f"lift second order n={n} h={h}",
            c_order(sp.Poly(second, a).LC(), c) == 3 * n * (h - 1),
        )
        require(
            f"lift first total n={n} h={h}",
            sp.Poly(first, a, c).total_degree() == 2 * (h + n * (h - 1)),
        )
        require(
            f"lift second total n={n} h={h}",
            sp.Poly(second, a, c).total_degree() == 3 * (h + n * (h - 1)),
        )
        require(f"lift first arm degree n={n} h={h}", degree_in(first_arm, a) == 2)
        require(f"lift second arm degree n={n} h={h}", degree_in(second_arm, a) == 3)
        require(
            f"lift arm immersion n={n} h={h}",
            sp.degree(sp.gcd(sp.diff(first_arm, a), sp.diff(second_arm, a)), a) == 0,
        )
        require(
            f"lift collision plus n={n} h={h}",
            first_arm.subs(a, 1 / E0) == 0 and second_arm.subs(a, 1 / E0) == 0,
        )
        require(
            f"lift collision minus n={n} h={h}",
            first_arm.subs(a, -1 / E0) == 0 and second_arm.subs(a, -1 / E0) == 0,
        )
        jacobian = sp.expand(
            sp.diff(first, a) * sp.diff(second, c)
            - sp.diff(first, c) * sp.diff(second, a)
        )
        require(f"lift zero Jacobian n={n} h={h}", jacobian == 0)
        nodal_chart_rows += 1

print(
    f"PASS global_nodal_rows={nodal_chart_rows} exact_floors=2h,3h "
    "warning=Jacobian_zero"
)


print("SECTION total-degree/composite arithmetic separation")
composite_rows = 0
for n in range(2, 7):
    for h in range(2, 13):
        total_unit = h + n * (h - 1)
        member_floor = 2 * total_unit
        height_floor = 3 * total_unit
        require(
            f"member strong floor n={n} h={h}",
            member_floor >= max(2 * h, 6),
        )
        require(f"member floor composite n={n} h={h}", is_composite(member_floor))
        require(
            f"height strong floor n={n} h={h}",
            height_floor >= max(3 * h, 8),
        )
        require(f"height floor composite n={n} h={h}", is_composite(height_floor))
        composite_rows += 1

print(f"PASS composite_floor_rows={composite_rows} n=2..6 h=2..12")


print("SECTION A13 specialization")
kappa_13 = sp.Rational(72, 91) ** 3
e_13 = sp.expand(a * (c**6 * a**2 + kappa_13**2))
first_13 = sp.expand(e_13**2 - 1)
second_13 = sp.expand(e_13**3 - e_13)
require("A13 e degree", degree_in(e_13, a) == 3)
require("A13 member a floor", degree_in(first_13, a) == 6)
require("A13 basis a floor", degree_in(second_13, a) == 9)
require("A13 member leading order", c_order(sp.Poly(first_13, a).LC(), c) == 12)
require("A13 basis leading order", c_order(sp.Poly(second_13, a).LC(), c) == 18)
require("A13 member total degree", sp.Poly(first_13, a, c).total_degree() == 18)
require("A13 basis total degree", sp.Poly(second_13, a, c).total_degree() == 27)
require("A13 member requested shadow", 18 >= max(2 * 3, 6))
require("A13 height requested shadow", 27 >= max(3 * 3, 8))
require(
    "A13 nodal warning",
    sp.expand(
        sp.diff(first_13, a) * sp.diff(second_13, c)
        - sp.diff(first_13, c) * sp.diff(second_13, a)
    )
    == 0,
)
print("PASS A13 h=3 a_floors=6,9 leading_c_orders=12,18 total_floors=18,27")


print("SECTION consequence ledger")
print("every-line=each completed target direction has central-arm degree>=2; each basis max>=3")
print("multiarm=Pade descent multiplies arm degree by h; the leading-c tariff adds to total degree")
print("planar=central chart is a nonautomorphic Keller map, so total degrees are composite>=6")
print("sharp=e^r and nodal (2,3) attain invoices but the displayed nodal lift has Jacobian zero")
print("nonconsequence=no global Darboux pair and no JC(2) counterexample is constructed")
print(f"CHECKS={CHECKS}")
print("FINAL PASS -- THM-3604 finite exact controls")
