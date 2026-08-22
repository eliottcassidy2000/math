"""Exact symbolic companion for THM-3684."""

import ast
import hashlib
from pathlib import Path

import sympy as sp


GATES = 0


def require(condition, label):
    global GATES
    if not condition:
        raise RuntimeError(label)
    GATES += 1


x, q, d = sp.symbols("x q d")
D = 1 + x**2 * q

B = sp.expand((d - 1) * (d + 2) ** 2)
C = sp.expand(d * (d + 2))
E = sp.expand((d - 1) * (d + 3))

b = sp.expand(B.subs(d, D))
c = sp.expand(x * C.subs(d, D))
e = sp.expand(x ** -2 * E.subs(d, D))


def jacobian_xq(left, right):
    return sp.expand(
        sp.diff(left, x) * sp.diff(right, q)
        - sp.diff(left, q) * sp.diff(right, x)
    )


def jacobian_xd(left, right):
    return sp.expand(
        sp.diff(left, x) * sp.diff(right, d)
        - sp.diff(left, d) * sp.diff(right, x)
    )


print("THM-3684 exact companion -- compiler-generator Laurent-weight no-mate")
print("status=PROVED VERIFIED-EXACT")

require(sp.expand(c**2 * e - b * (b + 4)) == 0, "compiler surface relation")
require(sp.expand(e - q * (D + 3)) == 0, "e is polynomial in source coordinates")
require(sp.expand(B.subs(d, D) - (D - 1) * (D + 2) ** 2) == 0, "b Laurent form")
require(sp.expand(c - x * D * (D + 2)) == 0, "c Laurent form")
require(sp.expand(e - x ** -2 * (D - 1) * (D + 3)) == 0, "e Laurent form")

# The change (x,q) -> (x,d=1+x^2 q) has determinant x^2.
for i in range(7):
    for j in range(7):
        left = x**i * q**j
        left_laurent = sp.expand(x ** (i - 2 * j) * (d - 1) ** j)
        require(
            sp.simplify(left.subs(q, (d - 1) / x**2) - left_laurent) == 0,
            f"source monomial Laurent conversion {i,j}",
        )

# Master pure-weight identity.  If F=x^a S(d) and G=x^m K(d), then
# Jac_(x,q)(F,G)=x^(a+m+1)(a S K'-m S'K).
families = (
    ("b", 0, B, 3),
    ("c", 1, C, 2),
    ("e", -2, E, 2),
)
for name, a, S, degree_S in families:
    require(sp.Poly(S, d).degree() == degree_S, f"{name} D-degree")
    for m in range(-7, 8):
        for n in range(8):
            K = (d + 5) ** n
            actual = sp.expand(x**2 * jacobian_xd(x**a * S, x**m * K))
            expected = sp.expand(
                x ** (a + m + 1)
                * (a * S * sp.diff(K, d) - m * sp.diff(S, d) * K)
            )
            require(sp.expand(actual - expected) == 0, f"master identity {name,m,n}")

# A constant can occur only in Laurent weight zero, hence only for
# m=-a-1.  The resulting one-variable operator has the displayed leading
# coefficient on every nonzero degree-n polynomial K.
for name, a, S, degree_S in families:
    m = -a - 1
    for n in range(13):
        K = d**n + sum(sp.Rational(j + 1, n + 2) * d**j for j in range(n))
        operator = sp.Poly(
            sp.expand(a * S * sp.diff(K, d) - m * sp.diff(S, d) * K), d
        )
        expected_degree = degree_S + n - 1
        expected_lead = a * n + (a + 1) * degree_S
        require(operator.degree() == expected_degree, f"{name} ODE degree n={n}")
        require(operator.LC() == expected_lead, f"{name} ODE leading coefficient n={n}")
        require(expected_degree > 0, f"{name} ODE positive degree n={n}")
        require(expected_lead != 0, f"{name} ODE nonzero lead n={n}")

# Direct source-coordinate controls for a rectangular monomial bank.
for name, source_generator in (("b", b), ("c", c), ("e", e)):
    for i in range(6):
        for j in range(6):
            source_monomial = x**i * q**j
            laurent_generator = {
                "b": B,
                "c": x * C,
                "e": x ** -2 * E,
            }[name]
            laurent_monomial = x ** (i - 2 * j) * (d - 1) ** j
            transformed = sp.expand(
                jacobian_xq(source_generator, source_monomial).subs(q, (d - 1) / x**2)
            )
            chart = sp.expand(x**2 * jacobian_xd(laurent_generator, laurent_monomial))
            require(sp.simplify(transformed - chart) == 0, f"direct chart control {name,i,j}")

# The companion itself contains no Python assert statements, so -O cannot
# remove a mathematical gate.
tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "assert-free AST")

semantic = "|".join(
    f"{name}:a={a}:S={sp.sstr(S)}:lead={a}*n+{a + 1}*{degree_S}"
    for name, a, S, degree_S in families
)
print(f"semantic_sha256={hashlib.sha256(semantic.encode()).hexdigest()}")
print(f"CHECKS={GATES}")
print("RESULT=PASS")
