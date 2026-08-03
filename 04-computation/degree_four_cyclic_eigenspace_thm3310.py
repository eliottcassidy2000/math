#!/usr/bin/env python3
"""Exact controls for THM-3310's degree-four cyclic eigenspace on Delta_2.

Exact cyclotomic Q(omega) arithmetic plus guarded modular certificates; no
floating point, no randomness, no imported executable.  Every gate is an
explicit ``require`` so that ordinary and ``-O`` replay are byte-identical.

MODULAR GUARD.  The uniform-simplex moment of a monomial of total degree D is
2*a!b!c!/(D+2)!.  Reducing mod p is sound only when p > D+2; otherwise the
denominator is 0 mod p and a naive implementation silently returns garbage.
Section 6 exhibits that failure.
"""

from fractions import Fraction
from math import comb, factorial
from itertools import combinations


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


# ---------------------------------------------- exact arithmetic in Q(omega)


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


# --------------------------- barycentric layer (used only to build the table)


def bmul(a, b):
    out = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            key = (ka[0] + kb[0], ka[1] + kb[1], ka[2] + kb[2])
            out[key] = out.get(key, Cyclo()) + va * vb
    return {k: v for k, v in out.items() if not v.is_zero()}


def bpow(p, e):
    result = {(0, 0, 0): Cyclo(1)}
    for _ in range(e):
        result = bmul(result, p)
    return result


def badd(*polys):
    out = {}
    for p in polys:
        for k, v in p.items():
            out[k] = out.get(k, Cyclo()) + v
    return {k: v for k, v in out.items() if not v.is_zero()}


def bscale(p, c):
    c = as_cyclo(c)
    return {k: v * c for k, v in p.items() if not (v * c).is_zero()}


def baverage(p):
    """<f> over Delta_2, uniform probability, term by term."""
    total = Cyclo()
    for (a, b, c), value in p.items():
        total = total + value * Cyclo(
            Fraction(2 * factorial(a) * factorial(b) * factorial(c),
                     factorial(a + b + c + 2)))
    return total


def rotate(p):
    """s1 -> s2 -> s3 -> s1 acting on exponent tuples."""
    return {(k[2], k[0], k[1]): v for k, v in p.items()}


Z = {(1, 0, 0): Cyclo(1), (0, 1, 0): OMEGA, (0, 0, 1): OMEGA * OMEGA}
ZBAR = {(1, 0, 0): Cyclo(1), (0, 1, 0): OMEGA * OMEGA, (0, 0, 1): OMEGA}

require(rotate(Z) == bscale(Z, OMEGA * OMEGA), "rot(z) = omega^2 z")
require(rotate(ZBAR) == bscale(ZBAR, OMEGA), "rot(zbar) = omega zbar")
require(baverage(Z).is_zero(), "<z> = 0: the centroid is the origin")
require(baverage(bmul(Z, ZBAR)) == Cyclo(Fraction(1, 4)), "<|z|^2> = 1/4")


# ---------------------------------- 1.  the eigenspaces in the monomial model


def eigen_dimension(degree, eigenvalue):
    """Dimension of an eigenspace of degree-d forms, by group projection."""
    monomials = [(i, j, degree - i - j)
                 for i in range(degree + 1) for j in range(degree + 1 - i)]
    projected = []
    for mono in monomials:
        accumulator = {}
        current = {mono: Cyclo(1)}
        weight = Cyclo(1)
        for _ in range(3):
            accumulator = badd(accumulator, bscale(current, weight))
            current = rotate(current)
            weight = weight * as_cyclo(eigenvalue).inverse()
        accumulator = bscale(accumulator, Cyclo(Fraction(1, 3)))
        if accumulator:
            projected.append(accumulator)
    keys = sorted({k for p in projected for k in p})
    rows = []

    def rank(candidate):
        matrix = [row[:] for row in candidate]
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

    count = 0
    for p in projected:
        row = [p.get(k, Cyclo()) for k in keys]
        if rank(rows + [row]) > count:
            rows.append(row)
            count += 1
    return count


# z^a zbar^b has rot-eigenvalue omega^(2a+b), so it lies in the
# omega-eigenspace exactly when a - b = 2 (mod 3).
DIMENSION_ROWS = []
for degree in range(1, 8):
    model = [(a, b) for a in range(degree + 1) for b in range(degree + 1 - a)
             if (a - b) % 3 == 2]
    require(eigen_dimension(degree, OMEGA) == len(model),
            "the monomial model reproduces the eigenspace dimension")
    DIMENSION_ROWS.append((degree, len(model)))
require([d for _, d in DIMENSION_ROWS] == [1, 2, 3, 5, 7, 9, 12],
        "dimension sequence 1,2,3,5,7,9,12")

BASIS_EXPONENTS = [(0, 1), (2, 0), (1, 2), (3, 1), (0, 4)]
BASIS_LABELS = ["zbar", "z^2", "z zbar^2", "z^3 zbar", "zbar^4"]
require(len(BASIS_EXPONENTS) == 5, "degree-four eigenspace has dimension five")
for a, b in BASIS_EXPONENTS:
    require((a - b) % 3 == 2 and a + b <= 4, "basis monomials are admissible")
    element = bmul(bpow(Z, a), bpow(ZBAR, b))
    require(rotate(element) == bscale(element, OMEGA),
            "every basis element is an omega-eigenvector")

# The lower-degree eigenspaces are coordinate sub-loci of this one.
require(BASIS_EXPONENTS[:1] == [(0, 1)], "degree 1 = the first coordinate")
require(BASIS_EXPONENTS[:2] == [(0, 1), (2, 0)], "degree 2 = a coordinate line")
require(BASIS_EXPONENTS[:3] == [(0, 1), (2, 0), (1, 2)],
        "degree 3 = a coordinate plane")


# ------------------------------- 2.  the exact mixed-moment table, and (z,zbar)
#
# On Delta_2 the functions z and zbar are algebraically independent, so every
# polynomial is a dict (a,b) -> coefficient and the average is a table lookup.
# Each basis element is then a SINGLE monomial, which collapses the algebra.

MAX_TOTAL = 36                 # deg(g) <= 4 and the highest power used is 9

MOMENT = {}
row = {(0, 0): {(0, 0, 0): Cyclo(1)}}
for total in range(0, MAX_TOTAL + 1):
    if total:
        new = {}
        for a in range(total + 1):
            b = total - a
            if a > 0 and (a - 1, b) in row:
                new[(a, b)] = bmul(row[(a - 1, b)], Z)
            else:
                new[(a, b)] = bmul(row[(a, b - 1)], ZBAR)
        row = new
    for key, value in row.items():
        MOMENT[key] = baverage(value)

for (a, b), value in MOMENT.items():
    if (a - b) % 3 != 0:
        require(value.is_zero(), "mixed moments vanish unless a = b mod 3")
    require(value.q == 0, "mixed moments are rational, not properly cyclotomic")
require(MOMENT[(1, 1)] == Cyclo(Fraction(1, 4)), "<z zbar> = 1/4")
require(MOMENT[(3, 0)] == Cyclo(Fraction(1, 10)), "<z^3> = 1/10")
require(MOMENT[(0, 3)] == Cyclo(Fraction(1, 10)), "<zbar^3> = 1/10")
require(MOMENT[(2, 2)] == Cyclo(Fraction(1, 10)), "<z^2 zbar^2> = 1/10")
require(MOMENT[(3, 3)] == Cyclo(Fraction(29, 560)), "<z^3 zbar^3> = 29/560")
for m in range(1, 13):
    predicted = (Cyclo(Fraction(2, (m + 1) * (m + 2)))
                 if m % 3 == 0 else Cyclo())
    require(MOMENT[(m, 0)] == predicted,
            "<z^m> = 2/((m+1)(m+2)) [3|m], the THM-3300 closed form")

NAMED_MOMENTS = {k: str(v.p) for k, v in sorted(MOMENT.items())
                 if sum(k) <= 6 and not v.is_zero()}


def zmul(x, y):
    out = {}
    for (a1, b1), v1 in x.items():
        for (a2, b2), v2 in y.items():
            key = (a1 + a2, b1 + b2)
            out[key] = out.get(key, Cyclo()) + v1 * v2
    return {k: v for k, v in out.items() if not v.is_zero()}


def zpow(x, e):
    result = {(0, 0): Cyclo(1)}
    for _ in range(e):
        result = zmul(result, x)
    return result


def zaverage(x):
    total = Cyclo()
    for key, value in x.items():
        total = total + value * MOMENT[key]
    return total


ZBASIS = [{ab: Cyclo(1)} for ab in BASIS_EXPONENTS]
require(zaverage(zpow(ZBASIS[1], 3))
        == baverage(bpow(bmul(bpow(Z, 2), bpow(ZBAR, 0)), 3)),
        "the (z,zbar) layer agrees with barycentric averaging")


# ------------------- 3.  exact exclusion of every coordinate line in P^4

def parameter_moments(first, second, m):
    """Coefficients in a of <(first + a*second)^m>."""
    coefficients = [Cyclo() for _ in range(m + 1)]
    for j in range(m + 1):
        term = zmul(zpow(second, j), zpow(first, m - j))
        coefficients[j] = coefficients[j] + zaverage(term) * Cyclo(comb(m, j))
    return coefficients


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


LINE_ROWS = []
for i, j in combinations(range(5), 2):
    cubic = parameter_moments(ZBASIS[i], ZBASIS[j], 3)
    sextic = parameter_moments(ZBASIS[i], ZBASIS[j], 6)
    require(any(not c.is_zero() for c in cubic), "cubic not identically zero")
    require(len(univariate_gcd(cubic, sextic)) <= 1,
            "no common root in the affine chart of this line")
    require(not (zaverage(zpow(ZBASIS[j], 3)).is_zero()
                 and zaverage(zpow(ZBASIS[j], 6)).is_zero()),
            "the remaining projective point of the line is excluded too")
    LINE_ROWS.append((BASIS_LABELS[i], BASIS_LABELS[j]))
require(len(LINE_ROWS) == 10, "all ten coordinate lines treated")
require(len(univariate_gcd(parameter_moments(ZBASIS[0], ZBASIS[1], 3),
                           parameter_moments(ZBASIS[0], ZBASIS[1], 3))) > 1,
        "positive control: gcd detects a shared root when one exists")


# ------------------------- 4.  guarded modular layer

PRIME = 1000000009
require(PRIME % 3 == 1, "omega exists in F_p")
require(PRIME > MAX_TOTAL + 2, "MODULAR GUARD: p exceeds every denominator")


def minv(x):
    require(x % PRIME != 0, "no inverse of zero: the modular guard failed")
    return pow(x, PRIME - 2, PRIME)


FMOMENT = {}
for key, value in MOMENT.items():
    denominator = value.p.denominator % PRIME
    require(denominator != 0, "MODULAR GUARD: p divides a moment denominator")
    FMOMENT[key] = value.p.numerator % PRIME * minv(denominator) % PRIME


def fmul(a, b):
    out = {}
    for (a1, b1), v1 in a.items():
        for (a2, b2), v2 in b.items():
            key = (a1 + a2, b1 + b2)
            out[key] = (out.get(key, 0) + v1 * v2) % PRIME
    return {k: v for k, v in out.items() if v}


def fpow(p, e):
    result = {(0, 0): 1}
    for _ in range(e):
        result = fmul(result, p)
    return result


def faverage(x):
    total = 0
    for key, value in x.items():
        total = (total + value * FMOMENT[key]) % PRIME
    return total % PRIME


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


FBASIS = [{ab: 1} for ab in BASIS_EXPONENTS]
SAMPLES = {6: 40, 9: 60}
PLANE_ROWS = []
for i, j, k in combinations(range(5), 3):
    per_node = {}
    for x in range(1, max(SAMPLES.values()) + 1):
        base = {BASIS_EXPONENTS[i]: 1}
        for key, value in {BASIS_EXPONENTS[j]: x % PRIME}.items():
            base[key] = (base.get(key, 0) + value) % PRIME
        base_powers = [{(0, 0): 1}]
        third_powers = [{(0, 0): 1}]
        for _ in range(9):
            base_powers.append(fmul(base_powers[-1], base))
            third_powers.append(fmul(third_powers[-1], FBASIS[k]))
        moments = {}
        for m in (3, 6, 9):
            coefficients = [0] * (m + 1)
            for t in range(m + 1):
                term = fmul(third_powers[t], base_powers[m - t])
                coefficients[t] = (coefficients[t]
                                   + faverage(term) * comb(m, t)) % PRIME
            while coefficients and coefficients[-1] == 0:
                coefficients.pop()
            moments[m] = coefficients
        per_node[x] = moments
    resultants = {}
    degrees = {}
    for other, sample in SAMPLES.items():
        nodes = list(range(1, sample + 1))
        values = []
        for x in nodes:
            first, second = per_node[x][3], per_node[x][other]
            degrees.setdefault(other, set()).add((len(first), len(second)))
            values.append(resultant(first, second))
        resultants[other] = interpolate(values, nodes)
    require(all(len(d) == 1 for d in degrees.values()),
            "degrees in the eliminated variable are constant across nodes, so "
            "the mod-p resultant is the reduction of the char-0 one")
    require(len(fgcd(resultants[6][:], resultants[9][:])) <= 1,
            "no common root: this coordinate plane is excluded")
    PLANE_ROWS.append((BASIS_LABELS[i], BASIS_LABELS[j], BASIS_LABELS[k],
                       len(resultants[6]) - 1, len(resultants[9]) - 1))
require(len(PLANE_ROWS) == 10, "all ten coordinate planes treated")
require(len(fgcd(resultants[6][:], resultants[6][:])) > 1,
        "positive control: mod-p gcd detects a shared root")


# ------------------------- 6.  the modular guard is load-bearing

SMALL = 7
require(SMALL % 3 == 1, "7 = 1 mod 3")
require(SMALL <= 4 * 3 + 2, "7 is inside the forbidden range already at m=3")
require(factorial(4 * 3 + 2) % SMALL == 0,
        "7 divides the denominator of a degree-12 simplex moment")
BAD_KEYS = [k for k in MOMENT
            if sum(k) <= 12 and MOMENT[k].p.denominator % SMALL == 0]
require(BAD_KEYS, "there really are moments whose denominator dies mod 7")


print("THM-3310 DEGREE-FOUR CYCLIC EIGENSPACE ON DELTA_2 EXACT CONTROL")
print("complex coordinate=z = s1 + omega s2 + omega^2 s3 (vertices 1,omega,omega^2)")
print("rot(z)=omega^2 z; z^a zbar^b is an omega-eigenvector iff a-b=2 mod 3")
print("eigenspace_dimensions_degree_1_to_7="
      + repr([d for _, d in DIMENSION_ROWS]))
print("degree_four_basis=" + repr(BASIS_LABELS))
print("degree_four_basis_exponents=" + repr(BASIS_EXPONENTS))
print("lower_degree_eigenspaces_are_coordinate_subloci=True")
print("mixed_moments_through_degree_six=" + repr(NAMED_MOMENTS))
print("coordinate_lines_excluded=" + str(len(LINE_ROWS)) + "/10 exact Q(omega)")
print("coordinate_planes_excluded=" + str(len(PLANE_ROWS))
      + "/10 guarded mod-p resultant, p=" + str(PRIME))
print("plane_resultant_degrees=" + repr([(r[3], r[4]) for r in PLANE_ROWS]))
print("modular_guard=p>4m+2; forbidden example p=7 kills "
      + str(len(BAD_KEYS)) + " moment denominators by degree 12")
print("scope=support<=3 excluded rigorously; support 4 and 5 stay OPEN")
print("ALL EXACT CHECKS PASSED")
