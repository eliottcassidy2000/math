#!/usr/bin/env python3
"""Exact controls for THM-3300's factorial-Gaussian torus bridge and no-go.

Exact rational, cyclotomic Q(omega), and mod-p arithmetic only; no floating
point, no randomness, no imported executable.  Every gate is an explicit
``require`` so that ordinary and ``-O`` replay are byte-identical.
"""

from fractions import Fraction
from math import comb, factorial
from itertools import product


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


# ------------------------------------------------- real polynomials in R^N


def rmul(a, b, n):
    out = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            key = tuple(ka[i] + kb[i] for i in range(n))
            out[key] = out.get(key, 0) + va * vb
    return {k: v for k, v in out.items() if v}


def rpow(p, e, n):
    result = {tuple([0] * n): Fraction(1)}
    for _ in range(e):
        result = rmul(result, p, n)
    return result


def rlaplace(p, n):
    out = {}
    for key, value in p.items():
        for i in range(n):
            if key[i] >= 2:
                low = list(key)
                low[i] -= 2
                low = tuple(low)
                out[low] = out.get(low, 0) + key[i] * (key[i] - 1) * value
    return {k: v for k, v in out.items() if v}


def gaussian_moment(p, n, covariance):
    """E[f(G)] for G ~ N(0, covariance * I_n), i.e. (exp(cov*Delta/2) f)(0)."""
    total = Fraction(0)
    current = dict(p)
    order = 0
    while current:
        constant = current.get(tuple([0] * n), 0)
        if constant:
            total += (Fraction(constant) * Fraction(covariance) ** order
                      / Fraction(2 ** order * factorial(order)))
        current = rlaplace(current, n)
        order += 1
    return total


def modulus_squared_monomial(alpha, n):
    """prod_i (u_i^2 + v_i^2)^(alpha_i) as a polynomial in R^(2n)."""
    dim = 2 * n
    result = {tuple([0] * dim): Fraction(1)}
    for i, power in enumerate(alpha):
        block = {}
        first = [0] * dim
        first[2 * i] = 2
        block[tuple(first)] = Fraction(1)
        second = [0] * dim
        second[2 * i + 1] = 2
        block[tuple(second)] = block.get(tuple(second), 0) + Fraction(1)
        result = rmul(result, rpow(block, power, dim), dim)
    return result


# --------------------------- 1.  the bridge: L_fac is a Gaussian functional

BRIDGE_ROWS = 0
for n in (1, 2, 3):
    dim = 2 * n
    for alpha in product(range(4), repeat=n):
        if sum(alpha) > 5:
            continue
        factorial_value = Fraction(1)
        for power in alpha:
            factorial_value *= factorial(power)
        gaussian_value = gaussian_moment(modulus_squared_monomial(alpha, n),
                                         dim, Fraction(1, 2))
        require(factorial_value == gaussian_value,
                "L_fac(x^alpha) = alpha! = E[prod |z_i|^(2 alpha_i)]")
        BRIDGE_ROWS += 1
require(BRIDGE_ROWS > 50, "bridge bank is non-vacuous")

# Covariance 1/2 per real coordinate is forced: E|z|^2 = 1.
require(gaussian_moment(modulus_squared_monomial((1,), 1), 2,
                        Fraction(1, 2)) == 1,
        "normalization E|z|^2 = 1")
require(gaussian_moment(modulus_squared_monomial((1,), 1), 2,
                        Fraction(1)) != 1,
        "hostile: covariance 1 is the wrong normalization")


def simplex_monomial(alpha, n):
    """<x^alpha> over Delta_(n-1) with the uniform probability measure."""
    numerator = factorial(n - 1)
    for power in alpha:
        numerator *= factorial(power)
    return Fraction(numerator, factorial(n - 1 + sum(alpha)))


MOMENT_MAP_ROWS = 0
for n in (2, 3):
    dim = 2 * n
    radial = {}
    for i in range(dim):
        key = [0] * dim
        key[i] = 2
        radial[tuple(key)] = Fraction(1)
    for alpha in product(range(4), repeat=n):
        total = sum(alpha)
        if total == 0 or total > 4:
            continue
        numerator = gaussian_moment(modulus_squared_monomial(alpha, n), dim,
                                    Fraction(1))
        denominator = gaussian_moment(rpow(radial, total, dim), dim,
                                      Fraction(1))
        sphere = numerator / denominator
        require(sphere == simplex_monomial(alpha, n),
                "moment map pushes uniform sphere measure to uniform simplex")
        MOMENT_MAP_ROWS += 1
require(MOMENT_MAP_ROWS > 20, "moment-map bank is non-vacuous")


# ------------------------ 2.  exact Q(omega) arithmetic for the cyclic route


class Cyclo:
    """p + q*omega with omega^2 = -1 - omega."""

    __slots__ = ("p", "q")

    def __init__(self, p=0, q=0):
        self.p = Fraction(p)
        self.q = Fraction(q)

    def __add__(self, other):
        other = as_cyclo(other)
        return Cyclo(self.p + other.p, self.q + other.q)

    __radd__ = __add__

    def __neg__(self):
        return Cyclo(-self.p, -self.q)

    def __sub__(self, other):
        return self + (-as_cyclo(other))

    def __mul__(self, other):
        other = as_cyclo(other)
        return Cyclo(self.p * other.p - self.q * other.q,
                     self.p * other.q + self.q * other.p - self.q * other.q)

    __rmul__ = __mul__

    def inverse(self):
        norm = self.p * self.p - self.p * self.q + self.q * self.q
        require(norm != 0, "invertible cyclotomic")
        return Cyclo((self.p - self.q) / norm, -self.q / norm)

    def is_zero(self):
        return self.p == 0 and self.q == 0

    def __eq__(self, other):
        other = as_cyclo(other)
        return self.p == other.p and self.q == other.q


def as_cyclo(value):
    return value if isinstance(value, Cyclo) else Cyclo(value, 0)


OMEGA = Cyclo(0, 1)
require((OMEGA * OMEGA + OMEGA + Cyclo(1)).is_zero(), "omega^2+omega+1=0")


def cmul(a, b):
    out = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            key = tuple(ka[i] + kb[i] for i in range(len(ka)))
            out[key] = out.get(key, Cyclo()) + va * vb
    return {k: v for k, v in out.items() if not v.is_zero()}


def cpow(p, e):
    keylen = len(next(iter(p)))
    result = {tuple([0] * keylen): Cyclo(1)}
    for _ in range(e):
        result = cmul(result, p)
    return result


def cadd(*polys):
    out = {}
    for p in polys:
        for k, v in p.items():
            out[k] = out.get(k, Cyclo()) + v
    return {k: v for k, v in out.items() if not v.is_zero()}


def cscale(p, c):
    c = as_cyclo(c)
    return {k: v * c for k, v in p.items() if not (v * c).is_zero()}


def simplex_average(p, n):
    """<f> over Delta_(n-1), f homogeneous, coefficients in Q(omega)."""
    total = Cyclo()
    for key, value in p.items():
        numerator = factorial(n - 1)
        for power in key:
            numerator *= factorial(power)
        total = total + value * Cyclo(
            Fraction(numerator, factorial(n - 1 + sum(key))))
    return total


# The degree-one cyclic eigenvector has an exact closed form.
# g = s1 + omega s2 + omega^2 s3 on Delta_2 is the degree-one cyclic
# eigenvector; its powers have an exact closed form.
CLOSED_FORM_ROWS = []
CYCLIC_GENERATOR = {(1, 0, 0): Cyclo(1), (0, 1, 0): OMEGA,
                    (0, 0, 1): OMEGA * OMEGA}
for m in range(1, 13):
    value = simplex_average(cpow(CYCLIC_GENERATOR, m), 3)
    predicted = (Cyclo(Fraction(2, (m + 1) * (m + 2)))
                 if m % 3 == 0 else Cyclo())
    require(value == predicted,
            "<g^m> = (n-1)! m!/(m+n-1)! * [n divides m]")
    CLOSED_FORM_ROWS.append((m, m % 3 == 0))
require(len(CLOSED_FORM_ROWS) == 12, "closed-form bank complete")
require(sum(1 for _, hit in CLOSED_FORM_ROWS if hit) == 4,
        "exactly the multiples of 3 survive")


def cyclic_eigenbasis(degree):
    """Basis of the omega-eigenspace of degree-d forms in s1,s2,s3."""
    monomials = [(i, j, degree - i - j)
                 for i in range(degree + 1)
                 for j in range(degree + 1 - i)]
    projected = []
    for monomial in monomials:
        accumulator = {}
        current = {monomial: Cyclo(1)}
        weight = Cyclo(1)
        for _ in range(3):
            accumulator = cadd(accumulator, cscale(current, weight))
            current = {(k[2], k[0], k[1]): v for k, v in current.items()}
            weight = weight * OMEGA.inverse()
        accumulator = cscale(accumulator, Cyclo(Fraction(1, 3)))
        if accumulator:
            projected.append(accumulator)
    keys = sorted({k for p in projected for k in p})

    def rank(rows):
        matrix = [row[:] for row in rows]
        pivot_row = 0
        for column in range(len(keys)):
            pivot = next((r for r in range(pivot_row, len(matrix))
                          if not matrix[r][column].is_zero()), None)
            if pivot is None:
                continue
            matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
            scale = matrix[pivot_row][column].inverse()
            matrix[pivot_row] = [x * scale for x in matrix[pivot_row]]
            for r in range(len(matrix)):
                if r != pivot_row and not matrix[r][column].is_zero():
                    factor = matrix[r][column]
                    matrix[r] = [x - factor * y
                                 for x, y in zip(matrix[r], matrix[pivot_row])]
            pivot_row += 1
        return pivot_row

    rows = []
    basis = []
    for p in projected:
        row = [p.get(k, Cyclo()) for k in keys]
        if rank(rows + [row]) > len(rows):
            rows.append(row)
            basis.append(p)
    return basis


def univariate_gcd(a, b):
    def normalize(x):
        while x and x[-1].is_zero():
            x.pop()
        return x

    a = normalize(a[:])
    b = normalize(b[:])
    while b:
        scale = b[-1].inverse()
        while len(a) >= len(b) and a:
            factor = a[-1] * scale
            offset = len(a) - len(b)
            for i in range(len(b)):
                a[offset + i] = a[offset + i] - factor * b[i]
            a = normalize(a)
        a, b = b, a
    return a


EXCLUSION_ROWS = []
DEGREE_ONE = cyclic_eigenbasis(1)
require(len(DEGREE_ONE) == 1, "degree-1 eigenspace is a line")
require(not simplex_average(cpow(DEGREE_ONE[0], 3), 3).is_zero(),
        "degree 1: <g^3> != 0, so the line is excluded")
EXCLUSION_ROWS.append((1, 1, "explicit"))

DEGREE_TWO = cyclic_eigenbasis(2)
require(len(DEGREE_TWO) == 2, "degree-2 eigenspace is a plane")


def moment_in_parameter(basis, m):
    """Coefficients in the parameter a of <(B0 + a B1)^m>."""
    coefficients = [Cyclo() for _ in range(m + 1)]
    for j in range(m + 1):
        term = cmul(cpow(basis[1], j), cpow(basis[0], m - j))
        coefficients[j] = (coefficients[j]
                           + simplex_average(term, 3) * Cyclo(comb(m, j)))
    return coefficients


CUBIC = moment_in_parameter(DEGREE_TWO, 3)
SEXTIC = moment_in_parameter(DEGREE_TWO, 6)
require(any(not c.is_zero() for c in CUBIC), "the cubic is not identically 0")
GCD_TWO = univariate_gcd(CUBIC, SEXTIC)
require(len(GCD_TWO) <= 1,
        "degree 2: <g^3> and <g^6> share no root, family excluded")

# The affine parameter B0+a*B1 misses the projective point [0:1].  Audit it
# explicitly so the exclusion is genuinely projective rather than chart-local.
DEGREE_TWO_INFINITY_CUBIC = simplex_average(cpow(DEGREE_TWO[1], 3), 3)
DEGREE_TWO_INFINITY_SEXTIC = simplex_average(cpow(DEGREE_TWO[1], 6), 3)
require(DEGREE_TWO_INFINITY_CUBIC == Cyclo(Fraction(1, 11340)),
        "degree 2: exact cubic at the infinity point")
require(DEGREE_TWO_INFINITY_SEXTIC == Cyclo(Fraction(1, 43783740)),
        "degree 2: exact sextic at the infinity point")
EXCLUSION_ROWS.append((2, 2, "gcd"))

# Machinery control: a polynomial shares every root with itself.
require(len(univariate_gcd(CUBIC[:], CUBIC[:])) > 1,
        "positive control: gcd detects a shared root when one exists")


# ---------------------- 3.  mod-p certificate for the degree-three family

PRIME = 1000000009
require(PRIME % 3 == 1, "p = 1 mod 3 so omega exists in F_p")


def find_omega(prime):
    generator = 2
    while pow(generator, (prime - 1) // 3, prime) == 1:
        generator += 1
    return pow(generator, (prime - 1) // 3, prime)


OMEGA_P = find_omega(PRIME)
require((OMEGA_P * OMEGA_P + OMEGA_P + 1) % PRIME == 0, "omega mod p")


def minv(x):
    return pow(x, PRIME - 2, PRIME)


def fmul(a, b):
    out = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            key = (ka[0] + kb[0], ka[1] + kb[1], ka[2] + kb[2])
            out[key] = (out.get(key, 0) + va * vb) % PRIME
    return {k: v for k, v in out.items() if v}


def fpow(p, e):
    result = {(0, 0, 0): 1}
    for _ in range(e):
        result = fmul(result, p)
    return result


def fadd(*polys):
    out = {}
    for p in polys:
        for k, v in p.items():
            out[k] = (out.get(k, 0) + v) % PRIME
    return {k: v for k, v in out.items() if v}


def fscale(p, c):
    return {k: (v * c) % PRIME for k, v in p.items() if (v * c) % PRIME}


def faverage(p):
    total = 0
    for (i, j, k), v in p.items():
        total = (total + v * 2 * factorial(i) * factorial(j) * factorial(k)
                 % PRIME * minv(factorial(i + j + k + 2))) % PRIME
    return total % PRIME


def feigenbasis(degree):
    monomials = [(i, j, degree - i - j)
                 for i in range(degree + 1) for j in range(degree + 1 - i)]
    projected = []
    for monomial in monomials:
        accumulator = {}
        current = {monomial: 1}
        weight = 1
        for _ in range(3):
            accumulator = fadd(accumulator, fscale(current, weight))
            current = {(k[2], k[0], k[1]): v for k, v in current.items()}
            weight = weight * minv(OMEGA_P) % PRIME
        accumulator = fscale(accumulator, minv(3))
        if accumulator:
            projected.append(accumulator)
    keys = sorted({k for p in projected for k in p})
    rows = []
    basis = []

    def rank(candidate):
        matrix = [row[:] for row in candidate]
        pivot_row = 0
        for column in range(len(keys)):
            pivot = next((r for r in range(pivot_row, len(matrix))
                          if matrix[r][column]), None)
            if pivot is None:
                continue
            matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
            scale = minv(matrix[pivot_row][column])
            matrix[pivot_row] = [x * scale % PRIME
                                 for x in matrix[pivot_row]]
            for r in range(len(matrix)):
                if r != pivot_row and matrix[r][column]:
                    factor = matrix[r][column]
                    matrix[r] = [(x - factor * y) % PRIME
                                 for x, y in zip(matrix[r],
                                                 matrix[pivot_row])]
            pivot_row += 1
        return pivot_row

    for p in projected:
        row = [p.get(k, 0) for k in keys]
        if rank(rows + [row]) > len(rows):
            rows.append(row)
            basis.append(p)
    return basis


def resultant(a, b):
    da, db = len(a) - 1, len(b) - 1
    if da < 0 or db < 0:
        return 0
    size = da + db
    matrix = [[0] * size for _ in range(size)]
    for i in range(db):
        for j, c in enumerate(reversed(a)):
            matrix[i][i + j] = c
    for i in range(da):
        for j, c in enumerate(reversed(b)):
            matrix[db + i][i + j] = c
    determinant = 1
    for column in range(size):
        pivot = next((r for r in range(column, size)
                      if matrix[r][column]), None)
        if pivot is None:
            return 0
        if pivot != column:
            matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
            determinant = (-determinant) % PRIME
        determinant = determinant * matrix[column][column] % PRIME
        scale = minv(matrix[column][column])
        for r in range(column + 1, size):
            if matrix[r][column]:
                factor = matrix[r][column] * scale % PRIME
                for k in range(column, size):
                    matrix[r][k] = (matrix[r][k]
                                    - factor * matrix[column][k]) % PRIME
    return determinant % PRIME


def interpolate(values, nodes):
    size = len(nodes)
    coefficients = [0] * size
    for i in range(size):
        numerator = [1]
        denominator = 1
        for j in range(size):
            if i == j:
                continue
            numerator = ([(-nodes[j] * numerator[0]) % PRIME]
                         + [(numerator[k - 1] - nodes[j] * numerator[k])
                            % PRIME for k in range(1, len(numerator))]
                         + [numerator[-1]])
            denominator = denominator * (nodes[i] - nodes[j]) % PRIME
        scale = values[i] * minv(denominator) % PRIME
        for k in range(len(numerator)):
            coefficients[k] = (coefficients[k] + scale * numerator[k]) % PRIME
    while coefficients and coefficients[-1] == 0:
        coefficients.pop()
    return coefficients


def fgcd(a, b):
    def normalize(x):
        while x and x[-1] == 0:
            x.pop()
        return x

    a = normalize(a[:])
    b = normalize(b[:])
    while b:
        scale = minv(b[-1])
        while len(a) >= len(b) and a:
            factor = a[-1] * scale % PRIME
            offset = len(a) - len(b)
            for i in range(len(b)):
                a[offset + i] = (a[offset + i] - factor * b[i]) % PRIME
            a = normalize(a)
        a, b = b, a
    return a


BASIS_THREE = feigenbasis(3)
require(len(BASIS_THREE) == 3, "degree-3 eigenspace has dimension 3")


def moment_in_b(m, a_value):
    base = fadd(BASIS_THREE[0], fscale(BASIS_THREE[1], a_value))
    coefficients = [0] * (m + 1)
    for j in range(m + 1):
        term = fmul(fpow(BASIS_THREE[2], j), fpow(base, m - j))
        coefficients[j] = (coefficients[j]
                           + faverage(term) * comb(m, j)) % PRIME
    return coefficients


RESULTANTS = {}
for other, sample in ((6, 40), (9, 60)):
    nodes = list(range(1, sample + 1))
    values = [resultant(moment_in_b(3, x), moment_in_b(other, x))
              for x in nodes]
    RESULTANTS[other] = interpolate(values, nodes)
    require(len(RESULTANTS[other]) > 1,
            "the resultant is a nonconstant polynomial in a")
GCD_THREE = fgcd(RESULTANTS[6][:], RESULTANTS[9][:])
require(len(GCD_THREE) <= 1,
        "degree 3: the resultants share no root, family excluded")


def moment_on_infinity_line(exponent):
    """Coefficients in t of <(B1+t*B2)^exponent> on c0=0."""
    coefficients = []
    for j in range(exponent + 1):
        term = fmul(fpow(BASIS_THREE[2], j),
                    fpow(BASIS_THREE[1], exponent - j))
        coefficients.append(faverage(term) * comb(exponent, j) % PRIME)
    while coefficients and coefficients[-1] == 0:
        coefficients.pop()
    return coefficients


INFINITY_MOMENTS = {m: moment_on_infinity_line(m) for m in (3, 6, 9)}
require(tuple(len(INFINITY_MOMENTS[m]) - 1 for m in (3, 6, 9)) == (3, 6, 9),
        "degree 3: full moment degrees on the missing line c0=0")
require(all(INFINITY_MOMENTS[m][0] and INFINITY_MOMENTS[m][-1]
            for m in (3, 6)),
        "degree 3: both projective endpoints are nonzero")
INFINITY_GCD_36 = fgcd(INFINITY_MOMENTS[3][:], INFINITY_MOMENTS[6][:])
INFINITY_GCD_39 = fgcd(INFINITY_MOMENTS[3][:], INFINITY_MOMENTS[9][:])
require(len(INFINITY_GCD_36) == len(INFINITY_GCD_39) == 1,
        "degree 3: the missing projective line contains no common zero")
EXCLUSION_ROWS.append((3, 3, "resultant gcd mod p"))
require(len(fgcd(RESULTANTS[6][:], RESULTANTS[6][:])) > 1,
        "positive control: mod-p gcd detects a shared root")
require(len(fgcd(INFINITY_MOMENTS[3][:], INFINITY_MOMENTS[3][:])) > 1,
        "positive control: missing-line gcd detects a shared root")


print("THM-3300 FACTORIAL-GAUSSIAN TORUS BRIDGE EXACT CONTROL")
print("bridge=L_fac(x^alpha)=alpha!=E[prod |z_i|^(2 alpha_i)], cov 1/2")
print("bridge_rows=" + str(BRIDGE_ROWS))
print("moment_map=<f(|z|^2)>_{S^(2n-1)} = <f>_{Delta_(n-1)}")
print("moment_map_rows=" + str(MOMENT_MAP_ROWS))
print("cyclic_closed_form=<(sum omega^j s_j)^m> = (n-1)! m!/(m+n-1)! [n|m]")
print("closed_form_rows=" + str(len(CLOSED_FORM_ROWS)))
print("eigenspace_dims=(deg1,deg2,deg3)=("
      + str(len(DEGREE_ONE)) + "," + str(len(DEGREE_TWO)) + ","
      + str(len(BASIS_THREE)) + ")")
print("cyclic_exclusions=" + repr(EXCLUSION_ROWS))
print("resultant_degrees=(Res36,Res39)=("
      + str(len(RESULTANTS[6]) - 1) + "," + str(len(RESULTANTS[9]) - 1) + ")")
print("projective_controls=deg2_infinity=(1/11340,1/43783740);"
      " deg3_c0_zero_degrees=(3,6,9); gcds=(constant,constant);"
      " endpoints_nonzero")
print("scope=HFC route analysis only; no FC(n) proof or refutation")
print("ALL EXACT CHECKS PASSED")
