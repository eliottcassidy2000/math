#!/usr/bin/env python3
"""Exact combinatorial companion for THM-3941's all-degree carrier theorem.

Reproduction:
  python3 04-computation/jc2_all_degree_centered_pole_carrier_routing_thm3941.py
  python3 -O 04-computation/jc2_all_degree_centered_pole_carrier_routing_thm3941.py
"""

from __future__ import annotations

import hashlib
import json

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def retained(exponent: int, e: int) -> bool:
    return exponent % e == 0


def infinity_allowed(m: int) -> bool:
    return m == 0 or m % 3 != 0


def c3_order_allowed(r: int) -> bool:
    return r >= 1 and r % 3 != 0


def c2_order_allowed(r: int) -> bool:
    return r >= 1 and r % 2 == 1


def carrier_rows(n: int):
    rows = []
    for m in range(n + 1):
        if not infinity_allowed(m):
            continue
        finite_degree = n - m
        if finite_degree == 0:
            rows.append(("root_regular", m, ()))
            continue
        if c3_order_allowed(finite_degree):
            rows.append(("C3", m, (finite_degree,)))
        if c2_order_allowed(finite_degree):
            rows.append(("C2", m, (finite_degree,)))
        for r in range(1, finite_degree):
            s = finite_degree - r
            if r <= s and c2_order_allowed(r) and c2_order_allowed(s):
                rows.append(("C2xC2", m, (r, s)))
    return rows


# Local character gates: only nontrivial inertia characters survive trace.
for exponent in range(-18, 0):
    gate(retained(exponent, 2) == (exponent % 2 == 0),
         f"C2 retained-exponent rule at {exponent}")
    gate(retained(exponent, 3) == (exponent % 3 == 0),
         f"C3 retained-exponent rule at {exponent}")

gate(all(not c2_order_allowed(2*j) for j in range(1, 8)),
     "C2 finite pole orders must be odd")
gate(all(not c3_order_allowed(3*j) for j in range(1, 8)),
     "C3 finite pole orders cannot be divisible by three")
gate(all(not infinity_allowed(3*j) for j in range(1, 8)),
     "positive infinity orders cannot be divisible by three")


# The three canonical cubic polynomial projections.
u = sp.symbols("u")
A_C3 = u**3
A_C2 = u**3 + u**2
A_C2xC2 = u**3 - 3*u
gate(sp.factor(sp.diff(A_C3, u)) == 3*u**2,
     "C3 carrier has one finite e=3 point")
gate(sp.factor(sp.diff(A_C2, u)) == u*(3*u+2),
     "one-C2 carrier has the selected e=2 point at zero")
gate(A_C2.subs(u, 0) == 0 and A_C2.subs(u, -sp.Rational(2, 3)) == sp.Rational(4, 27),
     "one-C2 carrier has distinct finite critical values")
gate(sp.expand(sp.diff(A_C2xC2, u) - 3*(u-1)*(u+1)) == 0,
     "two-C2 carrier has critical points plus and minus one")
gate(A_C2xC2.subs(u, 1) == -2 and A_C2xC2.subs(u, -1) == 2,
     "two-C2 carrier has distinct critical values")


# Shared homogeneous infinity-root address: multiplicity forces a=C=0.
w, a, color, c, d = sp.symbols("w a C c d")
phi_infinity_chart = a + color*w + c*w**2 + d*w**3
gate(phi_infinity_chart.subs(w, 0) == a,
     "homogeneous infinity-root value is a")
gate(sp.diff(phi_infinity_chart, w).subs(w, 0) == color,
     "homogeneous infinity-root derivative is C")


rows3 = carrier_rows(3)
rows4 = carrier_rows(4)
rows5 = carrier_rows(5)
expected3 = [
    ("C2", 0, (3,)),
    ("C3", 1, (2,)),
    ("C2xC2", 1, (1, 1)),
    ("C3", 2, (1,)),
    ("C2", 2, (1,)),
]
expected4 = [
    ("C3", 0, (4,)),
    ("C2xC2", 0, (1, 3)),
    ("C2", 1, (3,)),
    ("C3", 2, (2,)),
    ("C2xC2", 2, (1, 1)),
    ("root_regular", 4, ()),
]
expected5 = [
    ("C3", 0, (5,)),
    ("C2", 0, (5,)),
    ("C3", 1, (4,)),
    ("C2xC2", 1, (1, 3)),
    ("C2", 2, (3,)),
    ("C3", 4, (1,)),
    ("C2", 4, (1,)),
    ("root_regular", 5, ()),
]
gate(rows3 == expected3, "degree three recovers THM-3933/3936's five pole rows")
gate(rows4 == expected4, "degree four recovers THM-3938 plus its root-regular exit")
gate(rows5 == expected5, "degree five has exactly seven nonregular rows plus root-regular")
gate(sum(row[0] != "root_regular" for row in rows3) == 5,
     "degree three nonregular row count")
gate(sum(row[0] != "root_regular" for row in rows4) == 5,
     "degree four nonregular row count")
gate(sum(row[0] != "root_regular" for row in rows5) == 7,
     "degree five nonregular row count")


# A second, generating-function enumeration of the all-degree carrier count.
x = sp.symbols("x")
M = 1 + (x + x**2)/(1 - x**3)       # allowed infinity orders
R3 = (x + x**2)/(1 - x**3)          # one pole at the C3 point
R2 = x/(1 - x**2)                    # one pole at a C2 point
P2 = x**2/((1 - x**2)*(1 - x**4))   # unordered pair of odd C2 orders
nonregular_rational = M * (R3 + R2 + P2)
nonregular_series = sp.series(nonregular_rational, x, 0, 31).removeO().expand()
for n in range(1, 31):
    enumerated = sum(row[0] != "root_regular" for row in carrier_rows(n))
    gate(nonregular_series.coeff(x, n) == enumerated,
         f"degree {n} generating-function row count")


# The rational generating function has an exact period-twelve quadratic
# coefficient law.  This is an ordinal invoice, not an asymptotic fit.
# For N=12*j+r, the triple below is (coefficient of j^2,j,constant).
quasipolynomial_coefficients = [
    (6, 16, 0),
    (6, 10, 2),
    (6, 14, 4),
    (6, 16, 5),
    (6, 16, 5),
    (6, 14, 7),
    (6, 22, 10),
    (6, 16, 9),
    (6, 20, 12),
    (6, 22, 14),
    (6, 22, 15),
    (6, 20, 16),
]
y = x**12
quasipolynomial_numerator = sp.expand(sum(
    x**residue * (
        aa * y * (1 + y)
        + bb * y * (1 - y)
        + cc * (1 - y)**2
    )
    for residue, (aa, bb, cc) in enumerate(quasipolynomial_coefficients)
))
nonregular_numerator, nonregular_denominator = sp.together(
    nonregular_rational
).as_numer_denom()
gate(
    sp.expand(
        nonregular_numerator * (1 - y)**3
        - nonregular_denominator * quasipolynomial_numerator
    ) == 0,
    "period-twelve quadratic carrier-count identity",
)

# Within a fixed degree, the existing deterministic row order gives a
# zero-based natural-number address for every generated color-division task.
for n in range(1, 31):
    nonregular_rows = [row for row in carrier_rows(n) if row[0] != "root_regular"]
    ordinal = {row: index for index, row in enumerate(nonregular_rows)}
    gate(
        len(ordinal) == len(nonregular_rows)
        and all(nonregular_rows[ordinal[row]] == row for row in nonregular_rows),
        f"degree {n} carrier-to-natural ordinal bijection",
    )


# Hostile controls: forbidden residues never leak into generated rows.
for n in range(1, 31):
    for carrier, m, orders in carrier_rows(n):
        gate(infinity_allowed(m), f"degree {n} infinity residue control")
        if carrier == "C3":
            gate(len(orders) == 1 and c3_order_allowed(orders[0]),
                 f"degree {n} C3 residue control")
        elif carrier == "C2":
            gate(len(orders) == 1 and c2_order_allowed(orders[0]),
                 f"degree {n} C2 residue control")
        elif carrier == "C2xC2":
            gate(len(orders) == 2 and orders[0] <= orders[1]
                 and all(c2_order_allowed(order) for order in orders),
                 f"degree {n} two-C2 parity/orientation control")
        else:
            gate(carrier == "root_regular" and not orders and m == n,
                 f"degree {n} root-regular control")
        gate(m + sum(orders) == n, f"degree {n} pole-divisor degree control")


summary = {
    "checks": CHECKS,
    "scope": "all N;centered trace-zero;degA=3;carrier routing only",
    "carriers": "C3 one pole;C2 one pole;C2xC2 two unordered poles",
    "characters": "C2 odd orders;C3 nonzero mod3 orders;infinity nonzero mod3",
    "exits": "shared A-address non-unibranch;no finite pole root-regular",
    "degree3": "5 nonregular rows",
    "degree4": "5 nonregular rows",
    "degree5": "7 nonregular rows",
    "generating_function": "M*(R3+R2+P2);M=1+R3;P2=x2/((1-x2)(1-x4))",
    "quasipolynomial": "period12;(a,b,c)=6:16:0,6:10:2,6:14:4,6:16:5,6:16:5,6:14:7,6:22:10,6:16:9,6:20:12,6:22:14,6:22:15,6:20:16",
    "ordinal": "zero-based deterministic carrier-row rank at fixed N",
    "boundary": "no all-degree color-divisibility claim",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3941 all-degree centered pole-carrier exact companion")
print(f"CHECKS={CHECKS}")
print("SCOPE=all N;centered trace-zero;degA=3;carrier routing only")
print("CARRIERS=C3 one pole;C2 one pole;C2xC2 two unordered poles")
print("CHARACTERS=C2 odd orders;C3 nonzero-mod3;infinity nonzero-mod3")
print("EXITS=shared A-address non-unibranch;no finite pole root-regular")
print("DEGREE3=5 nonregular rows")
print("DEGREE4=5 nonregular rows")
print("DEGREE5=7 nonregular rows")
print("GF=M*(R3+R2+P2), M=1+R3, P2=x^2/((1-x^2)(1-x^4))")
print("QUASIPOLY_PERIOD12=(6,16,0);(6,10,2);(6,14,4);(6,16,5);(6,16,5);(6,14,7);(6,22,10);(6,16,9);(6,20,12);(6,22,14);(6,22,15);(6,20,16)")
print("ORDINAL=zero-based deterministic carrier-row rank at fixed N")
print("BOUNDARY=no all-degree color-divisibility claim")
print(f"SEMANTIC_SHA256={semantic}")
