#!/usr/bin/env python3
"""Exact controls for THM-3290's Archimedes-flatness GMC(3)/GVC(3) family.

Two independent routes are compared throughout:

  (R1) direct three-variable polynomial algebra -- build the polynomial, apply
       the operator Delta = 4 d_x d_y + d_t^2 the required number of times, and
       read the constant;
  (R2) the sphere reduction proved in the theorem -- a single generalized
       binomial coefficient times a Beta integral.

Exact integer and rational arithmetic only; no floating point, no randomness,
no imported executable.  Every gate is an explicit ``require`` so that ordinary
and ``-O`` replay are byte-identical.
"""

from fractions import Fraction
from math import comb, factorial


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


# ------------------------------------------------ three-variable polynomials
# A monomial x^i y^j t^k is the key (i, j, k); values are exact integers.


def pmul(a, b):
    out = {}
    for (ia, ja, ka), va in a.items():
        for (ib, jb, kb), vb in b.items():
            key = (ia + ib, ja + jb, ka + kb)
            out[key] = out.get(key, 0) + va * vb
    return {k: v for k, v in out.items() if v}


def padd(*polys):
    out = {}
    for p in polys:
        for k, v in p.items():
            out[k] = out.get(k, 0) + v
    return {k: v for k, v in out.items() if v}


def pscale(p, c):
    return {k: v * c for k, v in p.items() if v * c}


def ppow(p, n):
    result = {(0, 0, 0): 1}
    for _ in range(n):
        result = pmul(result, p)
    return result


def laplace(p):
    """Delta = 4 d_x d_y + d_t^2."""
    out = {}
    for (i, j, k), v in p.items():
        if i >= 1 and j >= 1:
            key = (i - 1, j - 1, k)
            out[key] = out.get(key, 0) + 4 * i * j * v
        if k >= 2:
            key = (i, j, k - 2)
            out[key] = out.get(key, 0) + k * (k - 1) * v
    return {k: v for k, v in out.items() if v}


def divide_by_x(p):
    """Exact division by x; fails loudly if the quotient is not polynomial."""
    out = {}
    for (i, j, k), v in p.items():
        require(i >= 1, "divisibility by x")
        out[(i - 1, j, k)] = v
    return out


def degree(p):
    return max(sum(k) for k in p)


def homogeneous(p):
    return len({sum(k) for k in p}) == 1


def double_factorial(n):
    result = 1
    while n > 1:
        result *= n
        n -= 2
    return result


X = {(1, 0, 0): 1}
Y = {(0, 1, 0): 1}
T = {(0, 0, 1): 1}
RHO = padd(ppow(T, 2), pmul(X, Y))                     # rho = t^2 + xy
A_FORM = padd(RHO, ppow(X, 2))                         # A = rho + x^2


def c_poly(nu):
    """C_nu = (rho^(2nu+1) - t^2 A^(2nu)) / x."""
    numerator = padd(ppow(RHO, 2 * nu + 1),
                     pscale(pmul(ppow(T, 2), ppow(A_FORM, 2 * nu)), -1))
    return divide_by_x(numerator)


def r_poly(nu):
    """R_nu = A * C_nu^2."""
    return pmul(A_FORM, ppow(c_poly(nu), 2))


# ------------------------------------- 1.  the supplied object, verbatim

C_SUPPLIED = padd(pmul(Y, ppow(RHO, 2)),
                  pscale(pmul(pmul(X, ppow(T, 2)), RHO), -2),
                  pscale(pmul(ppow(X, 3), ppow(T, 2)), -1))
require(C_SUPPLIED == c_poly(1),
        "supplied C equals (rho^3 - t^2 A^2)/x")
require(pmul(X, C_SUPPLIED)
        == padd(ppow(RHO, 3),
                pscale(pmul(ppow(T, 2), ppow(A_FORM, 2)), -1)),
        "key identity x C = rho^3 - t^2 A^2")

P_SUPPLIED = pmul(A_FORM, ppow(C_SUPPLIED, 2))
require(degree(P_SUPPLIED) == 12 and len(P_SUPPLIED) == 23,
        "supplied P has degree 12 and 23 terms")
require(homogeneous(A_FORM) and homogeneous(C_SUPPLIED)
        and homogeneous(P_SUPPLIED),
        "A, C, P are homogeneous")
require(laplace(ppow(X, 2)) == {},
        "x^2 is Delta-harmonic: x is an isotropic linear form")
require(laplace(RHO) == {(0, 0, 0): 6},
        "Delta rho = 2n = 6, so rho is the quadratic form of Delta")

Q_SUPPLIED = ppow(X, 2)
SUPPLIED_ROWS = []
SUPPLIED_MAX = 5
for m in range(1, SUPPLIED_MAX + 1):
    power = ppow(P_SUPPLIED, m)
    vanishing = power
    for _ in range(6 * m):
        vanishing = laplace(vanishing)
    require(vanishing == {}, "Delta^(6m) (P^m) = 0")
    witness = pmul(Q_SUPPLIED, power)
    for _ in range(6 * m + 1):
        witness = laplace(witness)
    value = witness.get((0, 0, 0), 0)
    predicted = (2 ** (8 * m + 1) * factorial(6 * m + 1) * factorial(2 * m)
                 * double_factorial(12 * m + 3) // double_factorial(4 * m + 1))
    require(value == predicted,
            "Delta^(6m+1)(Q P^m) = 2^(8m+1)(6m+1)!(2m)!(12m+3)!!/(4m+1)!!")
    SUPPLIED_ROWS.append((m, value))
require(len(SUPPLIED_ROWS) == SUPPLIED_MAX, "supplied bank is non-vacuous")

# The stated closed form is an m>=1 statement: at m=0 it disagrees.
require(laplace(Q_SUPPLIED) == {}, "Delta(Q) = 0")
require(2 ** 1 * factorial(1) * factorial(0)
        * double_factorial(3) // double_factorial(1) == 6,
        "the m=0 instance of the closed form would read 6, not 0")


# --------------------------------- 2.  Gaussian-moment dictionary and route R2


def moment(p):
    """L(f) = (exp(Delta/2) f)(0); for f homogeneous of degree 2k this is
    Delta^k f / (2^k k!), and equals E[f(G)] for a standard Gaussian G in the
    coordinates that diagonalize rho."""
    total = Fraction(0)
    current = dict(p)
    order = 0
    while current:
        constant = current.get((0, 0, 0), 0)
        if constant:
            total += Fraction(constant, 2 ** order * factorial(order))
        current = laplace(current)
        order += 1
    return total


for k in range(6):
    require(moment(ppow(RHO, k)) == double_factorial(2 * k + 1),
            "L(rho^k) = (2k+1)!! -- the n=3 radial moments")

# Route R2: the closed form proved in the theorem.
BETA = {}


def beta_value(k):
    """F_(2k)(1) = int_0^1 (1-u^2)^(2k) du = 2^(2k)(2k)!/(4k+1)!!."""
    if k not in BETA:
        BETA[k] = Fraction(2 ** (2 * k) * factorial(2 * k),
                           double_factorial(4 * k + 1))
    return BETA[k]


def generalized_binomial(top, low):
    """C(top, low) for any integer top and low >= 0."""
    if low < 0:
        return 0
    numerator = 1
    for step in range(low):
        numerator *= top - step
    return Fraction(numerator, factorial(low))


def sphere_average_R2(nu, k, delta):
    """<x^(2delta) R_nu^k> over the unit sphere, by the theorem."""
    return generalized_binomial(k - nu, k - delta) * beta_value(k)


def moment_R2(nu, k, delta):
    """L(x^(2delta) R_nu^k) = (2K+1)!! <...>, 2K the total degree."""
    total_degree = 2 * delta + (8 * nu + 4) * k
    require(total_degree % 2 == 0, "even total degree")
    return (double_factorial(total_degree + 1)
            * sphere_average_R2(nu, k, delta))


# Route R1 versus route R2 on a deterministic grid.
CROSS_ROWS = []
for nu in (1, 2, 3):
    r_nu = r_poly(nu)
    require(degree(r_nu) == 8 * nu + 4, "deg R_nu = 8nu+4")
    require(homogeneous(r_nu), "R_nu homogeneous")
    for k in (1, 2):
        if degree(r_nu) * k > 60:
            continue
        power = ppow(r_nu, k)
        for delta in (0, nu):
            direct = moment(pmul(ppow(X, 2 * delta), power))
            closed = moment_R2(nu, k, delta)
            require(direct == closed,
                    "route R1 (polynomial) equals route R2 (closed form)")
            CROSS_ROWS.append((nu, k, delta, direct))
require(len(CROSS_ROWS) >= 10, "cross-route bank is non-vacuous")


# ------------------------------------- 3.  the sharp threshold and the family

THRESHOLD_ROWS = []
for nu in (1, 2, 3, 4, 5):
    zeros = [k for k in range(1, 13) if sphere_average_R2(nu, k, 0) == 0]
    nonzeros = [k for k in range(1, 13) if sphere_average_R2(nu, k, 0) != 0]
    require(zeros == list(range(nu, 13)) and nonzeros == list(range(1, nu)),
            "<R_nu^k> vanishes exactly for k >= nu")
    for k in range(1, nu):
        require(sphere_average_R2(nu, k, 0)
                == (-1) ** k * comb(nu - 1, k) * beta_value(k),
                "sub-threshold closed form (-1)^k C(nu-1,k) F_(2k)(1)")
    THRESHOLD_ROWS.append((nu, nu))

# P_nu = R_nu^nu is a genuine counterexample: every exponent is >= nu.
FAMILY_ROWS = []
for nu in (1, 2, 3, 4, 5):
    for m in range(1, 7):
        require(sphere_average_R2(nu, nu * m, 0) == 0,
                "L(P_nu^m) = 0 for all m >= 1")
        witness = sphere_average_R2(nu, nu * m, nu)
        require(witness == beta_value(nu * m) and witness != 0,
                "L(x^(2nu) P_nu^m) != 0 for all m >= 1")
    FAMILY_ROWS.append((nu, 8 * nu ** 2 + 4 * nu, 2 * nu, (4 * nu + 2) * nu))

# Direct polynomial confirmation of the nu=2 member, where the threshold bites.
R2_POLY = r_poly(2)
require(moment(R2_POLY) == Fraction(-8, 15) * double_factorial(21),
        "nu=2 sub-threshold value is nonzero at k=1")
require(moment(ppow(R2_POLY, 2)) == 0,
        "nu=2 at k=2 vanishes: P_2 = R_2^2 is the counterexample")
require(moment(pmul(ppow(X, 4), ppow(R2_POLY, 2))) != 0,
        "nu=2 witness Q = x^4 is detected")


# ------------------------------- 4.  the mechanism is the t-average, not termwise

def constant_term_in_z(nu, k, delta):
    """Coefficient list in t^(2j) of CT_z(x^(2delta) R_nu^k) on the sphere.

    On rho=1 one has A = 1+w (w = z^2) and C_nu = (1 - t^2 A^(2nu))/z, so
    CT_z = [w^(k-delta)] ( a^k (1 - t^2 a^(2nu))^(2k) ).
    """
    coefficients = []
    for j in range(2 * k + 1):
        term = ((-1) ** j * comb(2 * k, j)
                * generalized_binomial(k + 2 * nu * j, k - delta))
        coefficients.append(Fraction(term))
    return coefficients


MECHANISM_ROWS = []
for nu in (1, 2):
    for k in (nu, nu + 1, nu + 2):
        coefficients = constant_term_in_z(nu, k, 0)
        require(any(c != 0 for c in coefficients),
                "the z-constant term is a nonzero polynomial in t")
        integrated = sum(c / (2 * j + 1) for j, c in enumerate(coefficients))
        require(integrated == sphere_average_R2(nu, k, 0) == 0,
                "vanishing is created by the t-average, not termwise")
        MECHANISM_ROWS.append((nu, k, len([c for c in coefficients if c])))
require(MECHANISM_ROWS, "mechanism bank is non-vacuous")

# The nu=1 case is exactly the classical hypergeometric identity.
IDENTITY_MAX = 40
for m in range(1, IDENTITY_MAX + 1):
    total = sum(Fraction((-1) ** j * comb(2 * m, j) * comb(m + 2 * j, m),
                         2 * j + 1)
                for j in range(2 * m + 1))
    require(total == 0, "sum_j (-1)^j C(2m,j) C(m+2j,m)/(2j+1) = 0")
    witness = sum(Fraction((-1) ** j * comb(2 * m, j)
                           * comb(m + 2 * j, m - 1), 2 * j + 1)
                  for j in range(2 * m + 1))
    require(witness == beta_value(m),
            "the Q-identity evaluates to the Beta integral")

# Hostile control: the flatness order is sharp.  F_(2k) is flat to order 2k+1
# at a=1, and the prefactor a^(k-nu) has degree k-nu; raising the prefactor
# degree to k destroys the vanishing.
for k in (2, 3, 4):
    broken = generalized_binomial(k, k) * beta_value(k)
    require(broken != 0,
            "prefactor degree k (instead of < k) breaks the vanishing")


print("THM-3290 ARCHIMEDES FLATNESS GMC(3)/GVC(3) FAMILY EXACT CONTROL")
print("supplied_object=(rho=t^2+xy, A=rho+x^2, C=(rho^3-t^2A^2)/x, P=AC^2)")
print("supplied_P_degree_terms=(12,23)")
print("key_identity=x*C=rho^3-t^2*A^2")
print("supplied_verified_m=1.." + str(SUPPLIED_MAX))
print("supplied_closed_form=2^(8m+1)(6m+1)!(2m)!(12m+3)!!/(4m+1)!!  [m>=1]")
print("dictionary=Delta is the Laplacian of rho; L(rho^k)=(2k+1)!!")
print("master_formula=<x^(2delta) R_nu^k> = C(k-nu,k-delta) * 2^(2k)(2k)!/(4k+1)!!")
print("cross_route_checks=" + str(len(CROSS_ROWS)))
print("threshold=<R_nu^k> vanishes exactly for k>=nu; rows=" + repr(THRESHOLD_ROWS))
print("subthreshold_closed_form=(-1)^k C(nu-1,k) F_(2k)(1)")
print("family_rows=(nu, deg P_nu, deg Q_nu, operator power)=" + repr(FAMILY_ROWS))
print("mechanism=t-average, not termwise; rows=" + repr(MECHANISM_ROWS))
print("hypergeometric_identity_verified_m=1.." + str(IDENTITY_MAX))
print("scope=GMC(3)/GVC(3) only; no JC, JC(2), GMC(2) or NC2 consequence")
print("ALL EXACT CHECKS PASSED")
