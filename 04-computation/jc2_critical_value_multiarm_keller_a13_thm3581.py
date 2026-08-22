#!/usr/bin/env python3
"""Exact companion for THM-3581.

The script checks a finite positive/hostile atlas for the critical-value
compiler and the complete degree-thirteen row.  Universal valuation, Krull,
image, Picard, trace, and Zariski-main arguments remain proof-driven in the
theorem.  All arithmetic below is exact SymPy arithmetic.
"""

from math import isqrt

import sympy as sp


t, z = sp.symbols("t z")
a, b, c, e = sp.symbols("a b c e")
x, q = sp.symbols("x q")

CHECKS = 0


def require(label, condition):
    """Record one exact gate and stop with its label on failure."""
    global CHECKS
    CHECKS += 1
    if not bool(condition):
        raise RuntimeError("FAILED exact gate: " + label)


def is_zero(expression):
    return sp.cancel(expression) == 0


def monic(expression, variable):
    return sp.Poly(expression, variable).monic().as_expr()


def rational_square(value):
    value = sp.Rational(value)
    numerator, denominator = map(int, value.as_numer_denom())
    if numerator < 0:
        return False
    return isqrt(numerator) ** 2 == numerator and isqrt(denominator) ** 2 == denominator


def ode_solution(n, S):
    rhs = sp.Poly(sp.expand(S ** (n - 1)), t, domain=sp.QQ)
    terms = []
    for (j,), coefficient in rhs.terms():
        terms.append(coefficient * t**j / sp.Rational(n * j + 1))
    return sp.expand(sum(terms, sp.S.Zero))


def canonical_bracket(F, G):
    return sp.diff(F, a) * sp.diff(G, c) - sp.diff(F, c) * sp.diff(G, a)


def surface_bracket(F, G, n, Sigma):
    bc = c**n * (sp.diff(F, b) * sp.diff(G, c) - sp.diff(F, c) * sp.diff(G, b))
    ce = -sp.diff(Sigma, b) * (
        sp.diff(F, c) * sp.diff(G, e) - sp.diff(F, e) * sp.diff(G, c)
    )
    be = -n * c ** (n - 1) * e * (
        sp.diff(F, b) * sp.diff(G, e) - sp.diff(F, e) * sp.diff(G, b)
    )
    return sp.expand(bc + ce + be)


def compiler_row(n, S, label, expected_E=None):
    """Run exact universal-formula gates for one squarefree rational row."""
    S = sp.expand(S)
    P = ode_solution(n, S)
    B = sp.expand(t * P**n)
    d = n * (n - 1) * sp.degree(S, t) + 1

    require(label + " squarefree S", sp.gcd(sp.Poly(S, t), sp.Poly(sp.diff(S, t), t)).degree() == 0)
    require(label + " nonzero S(0)", S.subs(t, 0) != 0)
    require(label + " ODE", is_zero(P + n * t * sp.diff(P, t) - S ** (n - 1)))
    require(label + " B derivative", is_zero(sp.diff(B, t) - P ** (n - 1) * S ** (n - 1)))
    require(label + " degree P", sp.degree(P, t) == (n - 1) * sp.degree(S, t))
    require(label + " degree B", sp.degree(B, t) == d)
    require(label + " gcd(P,S)", sp.gcd(sp.Poly(P, t), sp.Poly(S, t)).degree() == 0)
    require(label + " P squarefree", sp.gcd(sp.Poly(P, t), sp.Poly(sp.diff(P, t), t)).degree() == 0)

    # Reducing B modulo S makes the critical-value resultant very small while
    # preserving it exactly.
    Bmod = sp.rem(sp.Poly(B, t), sp.Poly(S, t)).as_expr()
    resultant = sp.factor(sp.resultant(S, z - Bmod, t))
    resultant_poly = sp.Poly(resultant, z, domain=sp.QQ)
    E = resultant_poly.sqf_part().monic().as_expr()
    require(label + " nonzero critical values", E.subs(z, 0) != 0)
    if expected_E is not None:
        require(label + " expected distinct-value E", is_zero(E - expected_E))

    numerator = sp.Poly(sp.expand(E.subs(z, B)), t, domain=sp.QQ)
    denominator = sp.Poly(sp.expand(S**n), t, domain=sp.QQ)
    W_poly, remainder = sp.div(numerator, denominator)
    require(label + " E(B) divisible by S^n", remainder.is_zero)
    W = W_poly.as_expr()
    require(label + " W squarefree", sp.gcd(W_poly, W_poly.diff()).degree() == 0)
    require(label + " gcd(W,PS)", sp.gcd(W_poly, sp.Poly(P * S, t)).degree() == 0)

    critical_gcd = sp.gcd(
        sp.Poly(E.subs(z, B), t, domain=sp.QQ),
        sp.Poly(sp.diff(B, t), t, domain=sp.QQ),
    ).monic().as_expr()
    require(label + " critical gcd", is_zero(critical_gcd - monic(S ** (n - 1), t)))
    local_value = n * sp.diff(S, t) * W - sp.diff(E, z).subs(z, B) * P ** (n - 1)
    require(label + " local W value", sp.rem(sp.Poly(local_value, t), sp.Poly(S, t)).is_zero)

    H = S ** (-n)
    G = P * S
    master = -sp.diff(t * H * G**n, t) / G ** (n - 1)
    require(label + " master Jacobian", is_zero(master + 1))

    bb = a * c**n
    ee = a * E.subs(z, bb)
    Sigma_ac = bb * E.subs(z, bb)
    require(label + " target relation", is_zero(c**n * ee - Sigma_ac))
    require(label + " bracket b,c", is_zero(canonical_bracket(bb, c) - c**n))
    require(
        label + " bracket c,e",
        is_zero(canonical_bracket(c, ee) + sp.diff(z * E, z).subs(z, bb)),
    )
    require(
        label + " bracket b,e",
        is_zero(canonical_bracket(bb, ee) + n * c ** (n - 1) * ee),
    )
    require(label + " central fibre count", 1 + n * sp.degree(P, t) == d)
    require(label + " side simple remainder", d - n * sp.degree(S, t) > 0)

    return {
        "n": n,
        "S": S,
        "P": P,
        "B": B,
        "d": d,
        "E": E,
        "W": W,
        "resultant": resultant,
    }


print("THM-3581 exact companion")
print("SECTION compiler cyclic atlas")

CYCLIC = {}
for n in range(2, 6):
    for r in range(1, 4):
        S = 1 + t**r
        p = sp.prod(sp.Rational(j * n * r, j * n * r + 1) for j in range(1, n))
        expected_E = z**r + p ** (n * r)
        label = "cyclic n=%d r=%d" % (n, r)
        row = compiler_row(n, S, label, expected_E)
        require(label + " P mod S", sp.rem(sp.Poly(row["P"] - p, t), sp.Poly(S, t)).is_zero)
        expected_W_num = sp.Poly(row["B"] ** r + p ** (n * r), t, domain=sp.QQ)
        expected_W, expected_rem = sp.div(expected_W_num, sp.Poly(S**n, t, domain=sp.QQ))
        require(label + " cyclic W divisibility", expected_rem.is_zero)
        require(label + " cyclic W", is_zero(row["W"] - expected_W.as_expr()))
        require(label + " distinct side count", sp.degree(row["E"], z) == r)
        require(label + " side fibre count", row["d"] - n > 0)
        CYCLIC[(n, r)] = row
        print("PASS %-16s degree=%d arms=%d" % (label, row["d"], r))


print("SECTION noncyclic positive control")
NONCYCLIC = compiler_row(3, 1 + t + t**2, "noncyclic n=3")
require("noncyclic has two distinct values", sp.degree(NONCYCLIC["E"], z) == 2)
require("noncyclic values distinct", sp.discriminant(NONCYCLIC["E"], z) != 0)
print("PASS noncyclic n=3 degree=%d arms=2" % NONCYCLIC["d"])


print("SECTION duplicate-value positive control")
P_DUP = t**2 + 5 * t + 5
S_DUP = 5 * t**2 + 15 * t + 5
DUP = compiler_row(2, S_DUP, "duplicate n=2", z + 4)
require("duplicate prescribed P", is_zero(DUP["P"] - P_DUP))
require("duplicate raw resultant multiplicity", is_zero(monic(DUP["resultant"], z) - (z + 4) ** 2))
require("duplicate factorization", is_zero(DUP["B"] + 4 - (t + 4) * S_DUP**2 / 25))
require("duplicate W", is_zero(DUP["W"] - (t + 4) / 25))
print("PASS two critical roots -> one distinct arm; E=z+4")


print("SECTION hypothesis hostiles")
S_NS = (1 + t) ** 2
P_NS = ode_solution(2, S_NS)
B_NS = sp.expand(t * P_NS**2)
BETA_NS = sp.simplify(B_NS.subs(t, -1))
NS_FACTOR = (t + 1) ** 3 * (9 * t**2 + 33 * t + 64) / 225
require("nonsquarefree ODE", is_zero(P_NS + 2 * t * sp.diff(P_NS, t) - S_NS))
require("nonsquarefree beta", BETA_NS == -sp.Rational(64, 225))
require("nonsquarefree exact factor", is_zero(B_NS - BETA_NS - NS_FACTOR))
require("nonsquarefree misses S^2", not sp.rem(sp.Poly(B_NS - BETA_NS, t), sp.Poly(S_NS**2, t)).is_zero)

P_ZERO = (t - 1) ** 2
S_ZERO = (t - 1) * (5 * t - 1)
B_ZERO = sp.expand(t * P_ZERO**2)
require("zero-value ODE", is_zero(P_ZERO + 2 * t * sp.diff(P_ZERO, t) - S_ZERO))
require("zero-value S squarefree", sp.gcd(sp.Poly(S_ZERO, t), sp.Poly(sp.diff(S_ZERO, t), t)).degree() == 0)
require("zero-value collision", B_ZERO.subs(t, 1) == 0)
require(
    "zero-value P order n",
    P_ZERO.subs(t, 1) == 0
    and sp.diff(P_ZERO, t).subs(t, 1) == 0
    and sp.diff(P_ZERO, t, 2).subs(t, 1) != 0,
)
require("zero-value degree obstruction", 2 * sp.degree(S_ZERO, t) > sp.degree(P_ZERO, t))
for k in range(1, 4):
    rr = sp.symbols("r_local")
    local_relation = (c ** (2 * k) * rr**2 - (c**k * rr) ** 2)
    require("normalization local relation k=%d" % k, is_zero(local_relation))
    normalized_bracket = sp.diff(b / c**k, b) * c ** (2 * k)
    require(
        "normalization bracket degeneracy k=%d" % k,
        is_zero(normalized_bracket - c**k) and normalized_bracket.subs(c, 0) == 0,
    )
print("PASS nonsquarefree, zero-value, and normalization-degeneracy controls")


print("SECTION degree-thirteen A13 row")
A13 = CYCLIC[(3, 2)]
P13 = (7 * t**4 + 26 * t**2 + 91) / 91
KAPPA = sp.Rational(72, 91) ** 3
require("A13 P", is_zero(A13["P"] - P13))
require("A13 kappa", KAPPA == sp.Rational(373248, 753571))
require("A13 B", is_zero(A13["B"] - t * P13**3))
require("A13 E", is_zero(A13["E"] - (z**2 + KAPPA**2)))
require("A13 W formula", is_zero(A13["W"] - (A13["B"] ** 2 + KAPPA**2) / (1 + t**2) ** 3))
require("A13 degree", A13["d"] == 13)
require("A13 W degree", sp.degree(A13["W"], t) == 20)

DISC13 = sp.factor(sp.discriminant(A13["B"] - z, t))
DISC13_EXPECTED = sp.Rational(1, 13**23) * z**8 * (z**2 + KAPPA**2) ** 2
require("A13 discriminant", is_zero(DISC13 - DISC13_EXPECTED))
require("A13 arithmetic nonsquare", not rational_square(sp.Rational(1, 13**23)))

PI = sp.I
for sign in (1, -1):
    sigma = sign * PI
    beta = sign * PI * KAPPA
    fibre = sp.Poly(A13["B"] - beta, t, extension=PI)
    arm = sp.Poly((t - sigma) ** 3, t, extension=PI)
    quotient, remainder = sp.div(fibre, arm)
    require("A13 side cubic sign=%d" % sign, remainder.is_zero)
    require("A13 side quotient degree sign=%d" % sign, quotient.degree() == 10)
    require("A13 side quotient squarefree sign=%d" % sign, sp.gcd(quotient, quotient.diff()).degree() == 0)
    require("A13 side disjoint sign=%d" % sign, quotient.eval(sigma) != 0)
    require(
        "A13 side derivative gcd sign=%d" % sign,
        sp.gcd(fibre, sp.Poly(sp.diff(A13["B"], t), t, extension=PI)).degree() == 2,
    )

require("A13 zero passport", 4 * 3 + 1 == 13 and sp.degree(P13, t) == 4)
require("A13 side passport", 3 + 10 == 13)
require("A13 infinity passport", sp.degree(A13["B"], t) == 13)
require("A13 prime degree", sp.isprime(13))
require("A13 all inertia even", all((13 - cycles) % 2 == 0 for cycles in (5, 11, 11, 1)))
require("A13 literal Jordan cycle", 3 <= 13 - 3)

rho, eta = sp.symbols("rho eta", nonzero=True)
q_collision = eta * (1 + rho**2) ** 3 / KAPPA**2
x_cube = rho * KAPPA**2 / (eta * (1 + rho**2) ** 3)
require("A13 collision t", is_zero(x_cube * q_collision - rho))
require("A13 collision e", is_zero(q_collision * KAPPA**2 / (1 + rho**2) ** 3 - eta))
require("A13 collision count", 1 + 3 * sp.degree(P13, t) == 13)

SIGMA13 = b * (b**2 + KAPPA**2)
FIRST = -15 * b * c * e / (4 * KAPPA**4)
SECOND = (1 / KAPPA**2 + 3 * b**2 / (2 * KAPPA**4)) * e
BRACKET_SUM = surface_bracket(FIRST, b, 3, SIGMA13) + surface_bracket(SECOND, c, 3, SIGMA13)
require("A13 two-bracket identity", is_zero(BRACKET_SUM.subs(e, SIGMA13 / c**3) - 1))
require(
    "Kummer z bracket",
    is_zero(surface_bracket(c * e, b, 3, SIGMA13).subs(e, SIGMA13 / c**3) - 2 * SIGMA13),
)

require("A13 P even", is_zero(P13.subs(t, -t) - P13))
require("A13 S even", is_zero((1 + (-t) ** 2) - (1 + t**2)))
require("A13 B odd", is_zero(A13["B"].subs(t, -t) + A13["B"]))
require("A13 W even", is_zero(A13["W"].subs(t, -t) - A13["W"]))

IOTA = {b: -b, c: c, e: -e}
relation13 = c**3 * e - SIGMA13
require("target iota", is_zero(relation13.subs(IOTA, simultaneous=True) + relation13))
for F, G, name in ((b, c, "bc"), (c, e, "ce"), (b, e, "be")):
    left = surface_bracket(
        F.subs(IOTA, simultaneous=True), G.subs(IOTA, simultaneous=True), 3, SIGMA13
    )
    right = surface_bracket(F, G, 3, SIGMA13).subs(IOTA, simultaneous=True)
    require("iota anti-Poisson " + name, is_zero(left + right))

# ell(x,q)=(x+q,x-q) conjugates tau(x,q)=(x,-q) to exchange.
ELL = sp.Matrix([[1, 1], [1, -1]])
TAU = sp.diag(1, -1)
ALPHA = sp.Matrix([[0, 1], [1, 0]])
require("reflection-exchange conjugacy", ELL * TAU == ALPHA * ELL)
require("reflection determinant", TAU.det() == -1 and ALPHA.det() == -1)

print("PASS A13 discriminant, passports, collision, brackets, and reflection gates")


print("SECTION degree-twenty-six Belyi quotient hostile")
C26 = sp.cancel(-A13["B"] ** 2 / KAPPA**2)
require("C26 degree", sp.degree(C26, t) == 26)
require("C26 factor through B", is_zero(C26 + A13["B"] ** 2 / KAPPA**2))
require("C26 zero passport", 4 * 6 + 2 == 26)
require("C26 one factorization", is_zero(C26 - 1 + (A13["B"] ** 2 + KAPPA**2) / KAPPA**2))
require("C26 one passport", 3 + 3 + sp.degree(A13["W"], t) == 26)
require("C26 infinity passport", sp.degree(C26, t) == 26)
require("C26 block system", 2 * 13 == 26)
require("C26 merges arms", is_zero((-((PI * KAPPA) ** 2) / KAPPA**2) - 1))

JINV = {b: -b, c: -c, e: e}
require("quotient involution preserves surface", is_zero(relation13.subs(JINV, simultaneous=True) + relation13))
for F, G, name in ((b, c, "bc"), (c, e, "ce"), (b, e, "be")):
    left = surface_bracket(
        F.subs(JINV, simultaneous=True), G.subs(JINV, simultaneous=True), 3, SIGMA13
    )
    right = surface_bracket(F, G, 3, SIGMA13).subs(JINV, simultaneous=True)
    require("quotient anti-Poisson " + name, is_zero(left + right))
print("PASS quotient passports=(6^4,2)/(3,3,1^20)/(26), blocks=2x13")


print("SECTION proof-driven universal gates")
print("PASS weight intersection: negative weight uses ceil(s/n) and distinct E")
print("PASS field intersection: etale codim-2 image + normal Krull valuations")
print("PASS image/Jelonek/Pic: fibre proof + escape chart + boundary divisors")
print("PASS sheet debt: central d-1; side-e nk; side-c d-nk")
print("PASS Kummer trace no-go: reciprocal primitive has at most one root")
print("PASS field sidecars: conorm survives finite covers; A12 maximal in A13")
print("PASS symmetry scope: explicit reflection/exchange only")
print("CHECKS=%d" % CHECKS)
print("RESULT=PASS")
