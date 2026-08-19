"""Exact probes for the JC(2) / dephasing / square-row synthesis.

This companion checks finite algebraic identities used in the accompanying
research reflection.  It uses explicit require/raise gates so optimized
Python executes every truth-bearing check.
"""

from math import isqrt

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def jac_xy(f, g, x, y):
    return sp.expand(sp.diff(f, x) * sp.diff(g, y) - sp.diff(f, y) * sp.diff(g, x))


def det2(a, b):
    return a[0] * b[1] - a[1] * b[0]


def norm2(a):
    return a[0] ** 2 + a[1] ** 2


def three_square_forbidden(number):
    while number % 4 == 0:
        number //= 4
    return number % 8 == 7


print("SECTION A. Equal-sum rows, square norms, and the Legendre twin.")
for n in range(1, 13):
    pairs = [(x, n + 1 - x) for x in range(1, n + 1) if x < n + 1 - x]
    require(len(pairs) == n // 2, f"strict-pair count failed at row {n}")
    completed_length = 2 * len(pairs) + (n % 2)
    require(completed_length == n, f"oriented completion failed at row {n}")
print("   strict count=floor(n/2); two orientations plus the odd diagonal give n: PASS")

fixed_y_rows = [
    [x * x + (n + 1) ** 2 for x in range(1, n + 1)] for n in range(1, 6)
]
expected_rows = [[5], [10, 13], [17, 20, 25], [26, 29, 34, 41], [37, 40, 45, 52, 61]]
require(fixed_y_rows == expected_rows, "fixed-y triangular array mismatch")
print(f"   fixed-y rows 1..5: {fixed_y_rows}: PASS")

fillers = [(2 * t + 1) ** 2 + 2 for t in range(6)]
for t, value in enumerate(fillers):
    n = 2 * t + 1
    require(value == n * n + 1 + 1, "three-square filler identity failed")
    require(value % 8 == 3 and value % 4 == 3, "filler residue failed")
    require(not three_square_forbidden(value), "filler incorrectly forbidden")
    twin = value + 4
    require(twin % 8 == 7 and three_square_forbidden(twin), "Legendre twin failed")
    for exponent in range(5):
        require(three_square_forbidden((4**exponent) * twin), "4-adic obstruction tower failed")
require(fillers == [3, 11, 27, 51, 83, 123], "user filler sequence mismatch")
print(f"   user fillers Q_t=(2t+1)^2+2: {fillers}: PASS")
print("   Q_t is 3 mod 8 and is three squares; Q_t+4 is 7 mod 8 and starts a 4^a(8b+7) tower: PASS")

z = sp.symbols("z")
pair_count_gf = z**2 / ((1 - z) * (1 - z**2))
completed_gf = sp.simplify(2 * pair_count_gf + z / (1 - z**2))
require(sp.simplify(completed_gf - z / (1 - z) ** 2) == 0, "row-count GF failed")
filler_gf = z * (3 + 2 * z**2 + 3 * z**4) / (1 - z**2) ** 3
series_fillers = sum(((2 * t + 1) ** 2 + 2) * z ** (2 * t + 1) for t in range(8))
require(
    sp.series(filler_gf, z, 0, 16).removeO().expand() == series_fillers.expand(),
    "filler GF failed",
)
print("   2R(z)+z/(1-z^2)=z/(1-z)^2 and the rational filler GF: PASS")

for t in range(1, 30):
    n = 2 * t + 1
    row_sum = sum(x * x + (n + 1 - x) ** 2 for x in range(1, t + 1))
    formula = t * (t + 1) * (8 * t + 7) // 3
    require(row_sum == formula, "odd-row norm sum failed")
require(three_square_forbidden(124), "S_3 control should be forbidden")
require(not three_square_forbidden(260), "S_4 control should not be forbidden")
print("   odd-row norm sum S_t=t(t+1)(8t+7)/3; S_3=124 forbidden, S_4=260 not: PASS")

require(4 * 4 + 7 * 7 == 65 and 1 + 8 * 8 == 65, "first duplicate failed")
print("   fixed-y rows enumerate representations, not distinct values: 65=(4,7)=(1,8): PASS")


print("\nSECTION B. Strong dephasing: conductances lead, triangle flux returns next.")
a, b, d, p1, p2, p3 = sp.symbols("a b d p1 p2 p3", real=True)
I = sp.I
H = sp.Matrix([[0, a, -I * d], [a, 0, b], [I * d, b, 0]])
Pdiag = sp.diag(p1, p2, p3)


def commutator_generator(matrix):
    return -I * (H * matrix - matrix * H)


K2 = commutator_generator(commutator_generator(Pdiag))
K3 = commutator_generator(K2)
diag2 = [sp.factor(K2[index, index]) for index in range(3)]
diag3 = [sp.factor(K3[index, index]) for index in range(3)]
expected2 = [
    2 * (a**2 * (p2 - p1) + d**2 * (p3 - p1)),
    2 * (a**2 * (p1 - p2) + b**2 * (p3 - p2)),
    2 * (b**2 * (p2 - p3) + d**2 * (p1 - p3)),
]
expected3 = [6 * a * b * d * (p3 - p2), 6 * a * b * d * (p1 - p3), 6 * a * b * d * (p2 - p1)]
require(all(sp.expand(left - right) == 0 for left, right in zip(diag2, expected2)), "K^2 diagonal failed")
require(all(sp.expand(left - right) == 0 for left, right in zip(diag3, expected3)), "K^3 triangle failed")
print("   Pi K^2 Pi=2L_c with c=|H_xy|^2: PASS")
print("   Pi K^3 Pi is the 6 Im(H12 H23 H31) circulation: PASS")

x, y, qx, qy = sp.symbols("x y qx qy")
Pplus = x + sp.Rational(1, 2) * x**2 * y**2 + sp.Rational(1, 3) * x**3 * y + sp.Rational(1, 3) * x * y**3
Pminus = x + sp.Rational(1, 2) * x**2 * y**2 + sp.Rational(1, 3) * x**3 * y - sp.Rational(1, 3) * x * y**3
Qlinear = qx * x + qy * y


def response_block(Psource):
    response = sp.Poly(jac_xy(Psource, Qlinear, x, y), x, y)
    row1 = response.coeff_monomial(x**2 * y)
    row2 = response.coeff_monomial(x * y**2)
    return sp.Matrix(
        [
            [sp.expand(row1).coeff(qx), sp.expand(row1).coeff(qy)],
            [sp.expand(row2).coeff(qx), sp.expand(row2).coeff(qy)],
        ]
    )


Aplus = response_block(Pplus)
Aminus = response_block(Pminus)
require(Aplus == sp.Matrix([[-1, 1], [-1, 1]]), "plus response block failed")
require(Aminus == sp.Matrix([[-1, 1], [1, 1]]), "minus response block failed")
require(Aplus.applyfunc(lambda entry: entry**2) == Aminus.applyfunc(lambda entry: entry**2), "conductance shadows differ")
require(Aplus.rank() == 1 and Aminus.rank() == 2, "rank hostile failed")
print(f"   same unit conductances, different response ranks: {Aplus.tolist()} vs {Aminus.tolist()}: PASS")


print("\nSECTION C. Reflection-equivariant Keller PDE: three low transverse boxes close.")
u, w, lead = sp.symbols("u w lead", nonzero=True)


def reflection_jacobian(Amap, Bmap):
    return sp.expand(
        sp.diff(Amap, u) * (Bmap + 2 * w * sp.diff(Bmap, w))
        - 2 * w * sp.diff(Amap, w) * sp.diff(Bmap, u)
    )


for degree in range(1, 9):
    aa = lead * u**degree

    bb11 = -sp.diff(aa, u) / 3
    residual11 = sp.Poly(reflection_jacobian(u + aa * w, 1 + bb11 * w) - 1, w).coeff_monomial(w**2)
    predicted11 = -lead**2 * degree * (degree + 2) * u ** (2 * degree - 2) / 3
    require(sp.simplify(residual11 - predicted11) == 0 and predicted11 != 0, "(1,1) box failed")

    bb12 = -sp.diff(aa, u) / 3
    dd12 = (3 * sp.diff(aa, u) ** 2 - 2 * aa * sp.diff(aa, u, 2)) / 15
    residual12 = sp.Poly(reflection_jacobian(u + aa * w, 1 + bb12 * w + dd12 * w**2) - 1, w).coeff_monomial(w**3)
    predicted12 = lead**3 * degree * (degree + 2) * (degree + 4) * u ** (3 * degree - 3) / 15
    require(sp.simplify(residual12 - predicted12) == 0 and predicted12 != 0, "(1,2) box failed")

    bb21 = -sp.diff(aa, u) / 3
    dd21_prime = sp.diff(aa, u) ** 2 - sp.Rational(2, 3) * aa * sp.diff(aa, u, 2)
    dd21 = sp.integrate(dd21_prime, u)
    residual21 = sp.Poly(reflection_jacobian(u + aa * w + dd21 * w**2, 1 + bb21 * w) - 1, w).coeff_monomial(w**3)
    predicted21 = -lead**3 * degree**2 * (degree + 2) * (2 * degree + 1) * u ** (3 * degree - 3) / (9 * (2 * degree - 1))
    require(sp.simplify(residual21 - predicted21) == 0 and predicted21 != 0, "(2,1) box failed")
print("   leading obstructions for boxes (1,1), (1,2), (2,1) are nonzero for degrees 1..8: PASS")
print("   formulas contain N(N+2), N(N+2)(N+4), and N^2(N+2)(2N+1): PASS")


print("\nSECTION D. Exact homogeneous-jet closure of the (deg P,deg Q)=(6,4), K=2 stratum.")
r0, r1, r2, r3, Acoef, Bcoef, Ccoef, shear, alpha, beta = sp.symbols(
    "r0 r1 r2 r3 Acoef Bcoef Ccoef shear alpha beta"
)
Hxy = x * y
R3 = r0 * y**3 + r1 * x * y**2 + r2 * x**2 * y + r3 * x**3
P5 = sp.Rational(3, 2) * Hxy * R3
require(jac_xy(Hxy**3, R3, x, y) + jac_xy(P5, Hxy**2, x, y) == 0, "E7 identity failed")
cross53 = sp.factor(jac_xy(P5, R3, x, y))
require(sp.expand(cross53.subs(x, 0) - sp.Rational(9, 2) * r0**2 * y**6) == 0, "x-axis E6 failed")
require(sp.expand(cross53.subs(y, 0) + sp.Rational(9, 2) * r3**2 * x**6) == 0, "y-axis E6 failed")

Q2 = Acoef * x**2 + Bcoef * Hxy + Ccoef * y**2
P4 = sp.Rational(3, 2) * Hxy * (Acoef * x**2 + Ccoef * y**2) + shear * Hxy**2
require(jac_xy(Hxy**3, Q2, x, y) + jac_xy(P4, Hxy**2, x, y) == 0, "E6 layer failed")
cross42 = sp.factor(jac_xy(P4, Q2, x, y))
require(sp.expand(cross42.subs(x, 0) - 3 * Ccoef**2 * y**4) == 0, "x-axis E4 failed")
require(sp.expand(cross42.subs(y, 0) + 3 * Acoef**2 * x**4) == 0, "y-axis E4 failed")

Q1 = alpha * x + beta * y
P3 = sp.Rational(3, 2) * Hxy * Q1
require(jac_xy(Hxy**3, Q1, x, y) + jac_xy(P3, Hxy**2, x, y) == 0, "E5 layer failed")
last_product = sp.factor(jac_xy(P3, Q1, x, y))
expected_last = -sp.Rational(3, 2) * (alpha * x - beta * y) * (alpha * x + beta * y)
require(sp.expand(last_product - expected_last) == 0, "E2 contradiction failed")
print("   E7 proportionality, E6/E4 axis squares, and E2 final product: PASS")
print("   the K=2 (6,4) stratum forces Q_1=0, contradicting the constant Jacobian row: PASS")


print("\nSECTION E. Infinite norm-balanced cancellation gadgets and their cross leaks.")
base_u = (1, 8)
base_up = (4, 7)
for parameter in range(3, 25):
    v = (parameter + 2, 3 * parameter + 1)
    vp = (parameter - 1, 3 * parameter + 2)
    require(norm2(base_u) == norm2(base_up) == 65, "base norm balance failed")
    require(norm2(v) == norm2(vp) == 10 * parameter**2 + 10 * parameter + 5, "moving norm balance failed")
    require((base_u[0] + v[0], base_u[1] + v[1]) == (base_up[0] + vp[0], base_up[1] + vp[1]), "equal-sum fibre failed")
    require(det2(base_u, v) == -5 * (parameter + 3), "negative determinant arm failed")
    require(det2(base_up, vp) == 5 * (parameter + 3), "positive determinant arm failed")
    require(det2(base_u, vp) == -5 * (parameter - 2), "first cross leak failed")
    require(det2(base_up, v) == 5 * (parameter - 2), "second cross leak failed")
print("   norm-balanced parallelogram family and both singleton cross leaks for 3<=r<25: PASS")
print("   r=4 has norms 65/205, cancelling determinants +/-35, cross leaks +/-10: PASS")


print("\nSECTION F. Cusp-square parametrization: all six natural output pairs are obstructed.")
v, yy = sp.symbols("v yy")
U = sp.Function("U")(v, yy)
T = yy**2 - 6 * v * U
S = yy**3 - 9 * v * U * yy
L = v**2 * (8 * v * U - yy**2)


def jac_vy(first, second):
    return sp.factor(
        sp.diff(first, v) * sp.diff(second, yy)
        - sp.diff(first, yy) * sp.diff(second, v)
    )


require(sp.factor(S**2 - T**3 - 27 * L * U**2) == 0, "cusp-square identity failed")
require(
    sp.factor(jac_vy(T, S) - 54 * v * U * (U + v * sp.diff(U, v))) == 0,
    "Jac(T,S) factor failed",
)
for first, second, label in [(L, T, "L,T"), (L, U, "L,U"), (L, S, "L,S")]:
    require(sp.factor(jac_vy(first, second)).subs(v, 0) == 0, f"Jac({label}) lacks v factor")
require(
    sp.factor(jac_vy(T, U) + 2 * (yy * sp.diff(U, v) + 3 * U * sp.diff(U, yy))) == 0,
    "Jac(T,U) transport equation failed",
)
expected_us_axis = 3 * yy * (
    yy * sp.Subs(sp.diff(U, v), v, 0)
    + 3 * U.subs(v, 0) * sp.diff(U.subs(v, 0), yy)
)
require(
    sp.simplify(jac_vy(U, S).subs(v, 0) - expected_us_axis) == 0,
    "Jac(U,S) axis factor failed",
)
print("   S^2=T^3+27LU^2 and Jac(T,S)=54vU(U+vU_v): PASS")
print("   Jac(L,T), Jac(L,U), Jac(L,S) have factor v; Jac(U,S)|v=0 has factor yy: PASS")
print("   Jac(T,U)=-2(yy U_v+3U U_y), the remaining polynomial transport obstruction: PASS")


print("\nSECTION G. An everywhere-immersive four-coordinate cusp packet.")
Upacket = 1 + yy - sp.Rational(1, 2) * yy**2 - sp.Rational(3, 2) * v * yy * (yy - 3)
Tpacket = yy**2 - 6 * v * Upacket
Spacket = yy**3 - 9 * v * Upacket * yy
Lpacket = v**2 * (8 * v * Upacket - yy**2)
packet = [Lpacket, Tpacket, Upacket, Spacket]
minors = [
    sp.expand(jac_vy(packet[first], packet[second]))
    for first in range(4)
    for second in range(first + 1, 4)
]
require(sp.factor(Spacket**2 - Tpacket**3 - 27 * Lpacket * Upacket**2) == 0, "immersive packet identity failed")
minor_ideal = sp.groebner(minors, v, yy, order="lex")
require(minor_ideal.reduce(sp.Integer(1))[1] == 0, "packet minors have a common affine zero")
weights = sp.symbols("weight0:6")
constant_linear = sp.Poly(sum(weight * minor for weight, minor in zip(weights, minors)) - 1, v, yy)
linear_solution = sp.linsolve([coefficient for _, coefficient in constant_linear.terms()], weights)
require(linear_solution is sp.EmptySet or linear_solution == sp.EmptySet, "constant-linear Keller projection unexpectedly exists")
print("   the six 2x2 minors generate the unit ideal: the A^2->A^4 cusp packet is immersive: PASS")
print("   no constant linear combination of its six minors equals 1: PASS")


print("\nSECTION H. Root-multiplicity divisor squeeze in the first Euclidean chamber.")
for g_degree in range(2, 8):
    for multiplicity in range(1, g_degree):
        Hmonomial = x**multiplicity * y ** (g_degree - multiplicity)
        for b_exponent in range(2, 6):
            r_degree = g_degree * b_exponent - 1
            for u_valuation in range(r_degree + 1):
                Rmonomial = x**u_valuation * y ** (r_degree - u_valuation)
                coefficient = g_degree * (multiplicity * b_exponent - u_valuation) - multiplicity
                require(coefficient != 0, "valuation resonance unexpectedly occurred")
                expected_bracket = (
                    coefficient
                    * x ** (multiplicity + u_valuation - 1)
                    * y ** (g_degree * (b_exponent + 1) - multiplicity - u_valuation - 2)
                )
                require(
                    sp.expand(jac_xy(Hmonomial, Rmonomial, x, y) - expected_bracket) == 0,
                    "root-multiplicity valuation formula failed",
                )
print("   v_L(Jac(H,R))=e+v_L(R)-1; its coefficient is g(eb-u)-e !=0: PASS")


def forced_exponent(root_multiplicity, chamber_k):
    numerator = root_multiplicity * (chamber_k - 1) + 1
    return (numerator + 1) // 2


require(forced_exponent(1, 1) == 1, "(3,2) squarefree exponent failed")
require(forced_exponent(1, 3) == 2, "(5,4) squarefree exponent failed")
require(forced_exponent(2, 3) == 3 and forced_exponent(3, 3) == 4, "repeated-root exponent failed")
print("   D_k(H)=product L^ceil((e(k-1)+1)/2), k=2b-a: exact sample exponents PASS")
print("   squarefree (a,b)=(5,4) forces H^2|Q_(m-1); multiplicities (2,3) force (3,4): PASS")

print("\nALL CHECKS PASSED")
