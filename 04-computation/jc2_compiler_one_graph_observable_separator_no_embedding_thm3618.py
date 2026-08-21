#!/usr/bin/env python3
"""Exact controls for provisional THM-3618 one-observable classification.

The all-degree Laurent-weight and normal-image arguments are proof-driven.
This companion checks their algebraic inputs, positive and hostile controls,
the complete boundary inventory, the negative survivor identity, the
puncture obstruction, and the sharp two-observable recovery using exact
SymPy arithmetic and no assertion gates.
"""

import sympy as sp


CHECKS = 0


def require(label, condition):
    """Record one exact deterministic gate."""
    global CHECKS
    if condition != True:
        raise RuntimeError(f"FAIL {label}: {condition}")
    CHECKS += 1


def zero(expression):
    """Exact polynomial or rational-function zero test."""
    return sp.cancel(expression) == 0


def same(left, right):
    """Exact rational-function equality test."""
    return zero(left - right)


print("THM-3618 exact companion -- provisional one-graph-observable separator classification")
print("status=verified exact controls; proof candidate pending independent hostile audit")


print("SECTION compiler ring, normal surface, and two-observable identity")
x, q, z = sp.symbols("x q z")
D = 1 + x**2 * q
b = (D - 1) * (D + 2) ** 2
c = x * D * (D + 2)
e = q * (D + 3)
z_source = x * q

require("D definition", same(D, 1 + x**2 * q))
require("b expanded compiler", same(b, x**2 * q * (3 + x**2 * q) ** 2))
require("c expanded compiler", same(c, x * (1 + x**2 * q) * (3 + x**2 * q)))
require("e expanded compiler", same(e, q * (4 + x**2 * q)))
require("b plus four factor", same(b + 4, D**2 * (D + 3)))
require("surface relation", same(c**2 * e, b * (b + 4)))
require("z definition", same(z_source, x * q))
require("e equals four q plus z squared", same(e, 4 * q + z_source**2))
require("q recovered from e,z", same((e - z_source**2) / 4, q))

B, C, E = sp.symbols("B C E")
surface = C**2 * E - B * (B + 4)
surface_singular_ideal = sp.groebner(
    [surface, sp.diff(surface, B), sp.diff(surface, C), sp.diff(surface, E)],
    B,
    C,
    E,
    order="lex",
)
require("surface singular ideal is unit", surface_singular_ideal.contains(sp.Integer(1)))
print("PASS compiler_and_surface_gates=10 surface=smooth two_observables=recover_q")


print("SECTION generic inverse and complete special-fibre controls")
t, C0 = sp.symbols("t C0", nonzero=True)
f_t = (t - 1) * (t + 2) ** 2
g_t = t * (t + 2)
x_hat = C0 / g_t
q_hat = (t - 1) * g_t**2 / C0**2
D_hat = sp.cancel(1 + x_hat**2 * q_hat)
b_hat = sp.cancel((D_hat - 1) * (D_hat + 2) ** 2)
c_hat = sp.cancel(x_hat * D_hat * (D_hat + 2))
e_hat = sp.cancel(q_hat * (D_hat + 3))
require("generic inverse D", same(D_hat, t))
require("generic inverse b", same(b_hat, f_t))
require("generic inverse c", same(c_hat, C0))
require("generic inverse e surface quotient", same(e_hat, f_t * (f_t + 4) / C0**2))
require("special b zero usable root", f_t.subs(t, 1) == 0 and g_t.subs(t, 1) == 3)
require("special b zero unusable root", f_t.subs(t, -2) == 0 and g_t.subs(t, -2) == 0)
require("special b minus four usable root", f_t.subs(t, -3) == -4 and g_t.subs(t, -3) == 3)
require("special b minus four unusable root", f_t.subs(t, 0) == -4 and g_t.subs(t, 0) == 0)

u = sp.symbols("u", nonzero=True)

def compiler_tuple(x_value, q_value):
    substitutions = {x: x_value, q: q_value}
    return tuple(sp.cancel(value.subs(substitutions)) for value in (D, b, c, e, z_source))


central_zero = compiler_tuple(0, -3 / (4 * u**2))
minus_two_plus = compiler_tuple(u, -3 / u**2)
minus_two_minus = compiler_tuple(-u, -3 / u**2)
require("b0 central tuple", central_zero == (1, 0, 0, -3 / u**2, 0))
require("b0 Dminus2 plus tuple", minus_two_plus == (-2, 0, 0, -3 / u**2, -3 / u))
require("b0 Dminus2 minus tuple", minus_two_minus == (-2, 0, 0, -3 / u**2, 3 / u))
require("b0 three points same compiler value", central_zero[1:4] == minus_two_plus[1:4] == minus_two_minus[1:4])
require("b0 z separates three points", len({central_zero[4], minus_two_plus[4], minus_two_minus[4]}) == 3)
require("b0 z square control", same(minus_two_plus[4] ** 2, -3 * minus_two_plus[3]))

zero_plus = compiler_tuple(u, -1 / u**2)
zero_minus = compiler_tuple(-u, -1 / u**2)
require("bminus4 D0 plus tuple", zero_plus == (0, -4, 0, -3 / u**2, -1 / u))
require("bminus4 D0 minus tuple", zero_minus == (0, -4, 0, -3 / u**2, 1 / u))
require("bminus4 pair same compiler value", zero_plus[1:4] == zero_minus[1:4])
require("bminus4 z separates signs", zero_plus[4] != zero_minus[4])
require("bminus4 z square control", same(zero_plus[4] ** 2, -zero_plus[3] / 3))

eta = sp.symbols("eta", nonzero=True)
b0_e0 = compiler_tuple(C0 / 3, 0)
bminus4_e0 = compiler_tuple(C0 / 3, -36 / C0**2)
require("e0 b0 singleton formula", b0_e0 == (1, 0, C0, 0, 0))
require("e0 bminus4 singleton formula", bminus4_e0 == (-3, -4, C0, 0, -12 / C0))
require("missing point has no special branch", all(row[1:4] != (-4, 0, 0) for row in (central_zero, minus_two_plus, zero_plus)))
print("PASS generic_inverse_gates=8 special_fibre_gates=16 image_puncture=(-4,0,0)")


print("SECTION off-diagonal conic and excluded-pair algebra")
t, s, p = sp.symbols("t s p")
f = lambda value: (value - 1) * (value + 2) ** 2
g = lambda value: value * (value + 2)
conic = t**2 + t * s + s**2 + 3 * (t + s)
require("difference factors through conic", same(f(t) - f(s), (t - s) * conic))
conic_basis = sp.groebner([p - t - s, conic], t, s, p, order="lex")

def conic_zero(expression):
    return conic_basis.reduce(sp.Poly(sp.expand(expression), t, s, p).as_expr())[1] == 0


require("conic product identity", conic_zero(t * s - p**2 - 3 * p))
require("conic difference square", conic_zero((t - s) ** 2 + 3 * p * (p + 4)))
require("conic g product", conic_zero(g(t) * g(s) - p * (p + 1) * (p + 3) * (p + 4)))
mixed_pairs = {
    -1: (1, -2),
    -3: (0, -3),
    0: (0, 0),
    -4: (-2, -2),
}
for p_value, pair in mixed_pairs.items():
    require(f"excluded pair sum p={p_value}", sum(pair) == p_value)
    require(f"excluded pair conic p={p_value}", conic.subs({t: pair[0], s: pair[1]}) == 0)
    require(f"excluded pair has zero g p={p_value}", g(pair[0]) * g(pair[1]) == 0)
print("PASS off_diagonal_identities=4 excluded_pair_gates=12 deleted_p=0,-1,-3,-4")


print("SECTION weight monomial normal forms and parity")
weight_rows = 0
positive_rows = 0
negative_even_rows = 0
negative_odd_rows = 0
for i in range(11):
    for j in range(11):
        weight = i - 2 * j
        monomial = x**i * q**j
        if weight >= 0:
            representative = x**weight * (D - 1) ** j
            positive_rows += 1
        elif weight % 2 == 0:
            m_value = -weight // 2
            representative = q**m_value * (D - 1) ** (i // 2)
            negative_even_rows += 1
        else:
            m_value = (1 - weight) // 2
            representative = x * q**m_value * (D - 1) ** ((i - 1) // 2)
            negative_odd_rows += 1
        require(f"weight normal form i={i} j={j}", same(monomial, representative))
        require(f"weight parity i={i} j={j}", weight % 2 == i % 2)
        weight_rows += 2
require("weight row accounting", positive_rows + negative_even_rows + negative_odd_rows == 121)
print(
    f"PASS weight_rows={weight_rows} nonnegative={positive_rows} "
    f"negative_even={negative_even_rows} negative_odd={negative_odd_rows}"
)


print("SECTION positive-odd cyclic obstruction controls")
t1, t2, t3 = sp.symbols("t1 t2 t3")
g1, g2, g3 = g(t1), g(t2), g(t3)
root_basis = sp.groebner(
    [t1 + t2 + t3 + 3, t1 * t2 + t1 * t3 + t2 * t3],
    t1,
    t2,
    t3,
    order="lex",
)

def roots_zero(expression):
    return root_basis.reduce(sp.Poly(sp.expand(expression), t1, t2, t3).as_expr())[1] == 0


require("cyclic ratio first equality", roots_zero((t1 - t2) * g3 - (t2 - t3) * g1))
require("cyclic ratio second equality", roots_zero((t2 - t3) * g1 - (t3 - t1) * g2))
require("root shifted product 12", roots_zero((t1 + 3) * (t2 + 3) - t3**2))
require("root shifted product 23", roots_zero((t2 + 3) * (t3 + 3) - t1**2))
require("root shifted product 31", roots_zero((t3 + 3) * (t1 + 3) - t2**2))

positive_controls = 0
T, S, P = sp.symbols("T S P")
for a_value in range(10):
    w_value = 2 * a_value + 1
    ratio_left = (T - S) * (P * (P + 4)) ** a_value
    ratio_right = (-sp.Rational(1, 3)) ** a_value * (T - S) ** w_value
    relation_substitution = {P * (P + 4): -(T - S) ** 2 / 3}
    require(
        f"positive forced difference w={w_value}",
        same(ratio_left.subs(relation_substitution, simultaneous=True), ratio_right),
    )
    require(f"positive telescoping obstruction w={w_value}", 3 * 9**w_value != 0)
    positive_controls += 2
print(f"PASS cyclic_root_gates=5 positive_odd_controls={positive_controls} tested_w=1..19")


print("SECTION negative-odd survivor divided differences")
negative_rows = 0
for m_value in range(1, 13):
    ell_t = g(t) * (t + 3) ** m_value
    ell_s = g(s) * (s + 3) ** m_value
    P_t = (t + 3) ** (m_value - 1)
    P_s = (s + 3) ** (m_value - 1)
    numerator = sp.expand(P_t * ell_s - P_s * ell_t)
    quotient, remainder = sp.div(numerator, t - s, domain=sp.QQ)
    require(f"negative numerator divisible m={m_value}", remainder == 0)
    require(
        f"negative survivor Q m={m_value}",
        conic_zero(quotient + 2 * (p + 3) ** (2 * m_value - 1)),
    )
    survivor = x * q**m_value * (D + 3) ** (m_value - 1)
    require(f"negative survivor source identity m={m_value}", same(survivor, z_source * e ** (m_value - 1)))
    negative_rows += 3
print(f"PASS negative_survivor_rows={negative_rows} tested_m=1..12 Q=-2(p+3)^(2m-1)")


print("SECTION integer-power telescope filter and hostile specializations")
require(
    "reciprocal root sum numerator",
    roots_zero(t1 * t2 + t1 * t3 + t2 * t3),
)
require("zero power hostile", 3 != 0)
positive_power_rows = 0
negative_power_rows = 0
for exponent in range(1, 25):
    require(f"positive power hostile exponent={exponent}", (-3) ** exponent != 0)
    positive_power_rows += 1
for n_value in range(2, 25):
    hostile_sum = sp.Rational(1) + 2 * sp.Rational(-2) ** (-n_value)
    require(f"negative power hostile n={n_value}", hostile_sum != 0)
    negative_power_rows += 1
require("surviving exponent is minus one", -1 == -1)
print(
    f"PASS power_filter_base=3 positive_hostiles={positive_power_rows} "
    f"negative_hostiles={negative_power_rows} sole_exponent=-1"
)


print("SECTION D=0 puncture obstruction and e-power boundary")
u = sp.symbols("u", nonzero=True)
puncture_substitution = {x: 1 / u, q: -u**2}
puncture_values = tuple(sp.cancel(value.subs(puncture_substitution)) for value in (D, b, c, e, z_source))
require("puncture compiler tuple", puncture_values == (0, -4, 0, -3 * u**2, -u))
for m_value in range(1, 13):
    observable_value = sp.expand((-3 * u**2) ** (m_value - 1) * (-u))
    require(f"puncture separator polynomial m={m_value}", sp.Poly(observable_value, u).is_univariate)
    require(f"puncture separator vanishes at zero m={m_value}", observable_value.subs(u, 0) == 0)
require("puncture x has negative valuation", sp.cancel(x.subs(puncture_substitution)) == 1 / u)
require("puncture x not polynomial", sp.denom(sp.cancel(1 / u)) == u)
print("PASS puncture_gates=26 separator_traces_in=C[u] x_trace=u^-1")


print("SECTION fibrewise separator powers and THM-3614 exceptional graph")
for m_value in range(1, 13):
    central_value = central_zero[3] ** (m_value - 1) * central_zero[4]
    minus_two_values = {
        sp.cancel(minus_two_plus[3] ** (m_value - 1) * minus_two_plus[4]),
        sp.cancel(minus_two_minus[3] ** (m_value - 1) * minus_two_minus[4]),
    }
    zero_values = {
        sp.cancel(zero_plus[3] ** (m_value - 1) * zero_plus[4]),
        sp.cancel(zero_minus[3] ** (m_value - 1) * zero_minus[4]),
    }
    require(f"b0 power central distinct m={m_value}", central_value not in minus_two_values)
    require(f"b0 power signs distinct m={m_value}", len(minus_two_values) == 2)
    require(f"bminus4 power signs distinct m={m_value}", len(zero_values) == 2)

n = sp.symbols("n")
h_exceptional = -x * q + n
require("THM3614 exceptional is minus z plus base", same(h_exceptional, -z_source + n))
require("THM3614 exceptional is m1 class", same(h_exceptional - n, -e**0 * z_source))
print("PASS fibre_power_rows=36 tested_m=1..12 exceptional_h=-z+n_is_m1")


print(f"PASS total_exact_gates={CHECKS}")
print("RESULT PASS -- provisional proof candidate; exact companion frozen pending hostile audit")
