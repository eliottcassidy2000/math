#!/usr/bin/env python3
"""Array-based independent audit of the affine exceptional retained packet."""

from fractions import Fraction
from hashlib import sha256
from math import comb

import sympy as sp


MOD = 137
CUT = 5
JET = 6
POINTS = (-1, 0, 1)
SPLIT_ROOTS = (44, 82, 92, 134)


def inv(value):
    return pow(value % MOD, -1, MOD)


def rat(a, b=1):
    return a % MOD * inv(b) % MOD


def zero():
    return [[0] * (JET + 1) for _ in range(JET + 1)]


def scalar(value):
    result = zero()
    result[0][0] = value % MOD
    return result


def add(*values):
    result = zero()
    for value in values:
        for i in range(JET + 1):
            for j in range(JET + 1 - i):
                result[i][j] = (result[i][j] + value[i][j]) % MOD
    return result


def scale(value, coefficient):
    return [[(coefficient * value[i][j]) % MOD if i + j <= JET else 0
             for j in range(JET + 1)] for i in range(JET + 1)]


def mul(left, right):
    result = zero()
    for i in range(JET + 1):
        for j in range(JET + 1 - i):
            if not left[i][j]:
                continue
            for u in range(JET + 1 - i):
                for v in range(JET + 1 - i - u - j):
                    if right[u][v]:
                        result[i + u][j + v] = (
                            result[i + u][j + v] + left[i][j] * right[u][v]
                        ) % MOD
    return result


def power(value, exponent):
    result = scalar(1)
    for _ in range(exponent):
        result = mul(result, value)
    return result


def diff(value, axis):
    result = zero()
    for i in range(JET + 1):
        for j in range(JET + 1 - i):
            exponent = i if axis == 0 else j
            if exponent:
                if axis == 0:
                    result[i - 1][j] = exponent * value[i][j] % MOD
                else:
                    result[i][j - 1] = exponent * value[i][j] % MOD
    return result


def inverse(value):
    unit = inv(value[0][0])
    remainder = add(scale(value, unit), scalar(-1))
    result = zero()
    term = scalar(1)
    for degree in range(JET + 1):
        result = add(result, scale(term, -1 if degree % 2 else 1))
        term = mul(term, remainder)
    return scale(result, unit)


HVAR = zero()
HVAR[1][0] = 1
TVAR = zero()
TVAR[0][1] = 1


def uadd(*polynomials):
    result = {}
    for polynomial in polynomials:
        for degree, coefficient in polynomial.items():
            result[degree] = (result.get(degree, 0) + coefficient) % MOD
    return {degree: coefficient for degree, coefficient in result.items() if coefficient}


def uscale(polynomial, coefficient):
    return {degree: coefficient * value % MOD for degree, value in polynomial.items()
            if coefficient * value % MOD}


def ushift(polynomial, amount):
    return {degree + amount: coefficient for degree, coefficient in polynomial.items()}


def q_polynomial(root):
    qbase = {0: rat(-3, 4), 1: 1, 2: rat(-27, 4), 3: -2 % MOD,
             4: rat(9, 2), 5: 1}
    pbase = {2: 1, 4: -2 % MOD, 6: 1}
    q6 = uadd(qbase, uscale(pbase, rat(-259, 36)))
    r1 = uadd(pbase, uscale(ushift(pbase, 2), -1))
    r2 = uadd(uscale(pbase, 4), uscale(ushift(pbase, 1), -9))
    parabola = (rat(520, 9) * root * root - rat(1688, 81) * root - rat(5717, 729)) % MOD
    return uadd(q6, uscale(r1, parabola), uscale(r2, root))


def around(polynomial, point):
    result = zero()
    for degree, coefficient in polynomial.items():
        for order in range(min(degree, JET) + 1):
            result[order][0] = (
                result[order][0]
                + coefficient * comb(degree, order) * point ** (degree - order)
            ) % MOD
    return result


ROWS = tuple((stable, point, source)
             for stable in range(CUT + 1)
             for source in range(CUT + 1 - stable)
             for point in POINTS)
CONSTANT = [1 if stable == 0 and source == 0 else 0
            for stable, _, source in ROWS]


def packets(root, mixed):
    Q = q_polynomial(root)
    result = {}
    for point in POINTS:
        qjet = add(around(Q, point), TVAR, power(TVAR, 2) if mixed else zero())
        xjet = add(scalar(point), HVAR)
        D = add(scalar(1), mul(power(xjet, 2), qjet))
        c = mul(mul(xjet, D), add(D, scalar(2)))
        y = scale(c, inv(3))
        z = add(mul(qjet, add(D, scalar(3))), scalar(3))
        Dinv = inverse(D)
        a = mul(qjet, power(Dinv, 2))
        A = add(a, scalar(rat(3, 4)))
        if y[0][0] or z[0][0] or A[0][0]:
            raise RuntimeError(("origin", root, point))
        result[point] = {
            "y": [power(y, n) for n in range(CUT + 2)],
            "z": [power(z, n) for n in range(CUT + 2)],
            "A": [power(A, n) for n in range(CUT + 2)],
            "c": [power(c, n) for n in range(CUT + 2)],
            "area": add(mul(diff(y, 0), diff(z, 1)), scale(mul(diff(y, 1), diff(z, 0)), -1)),
            "yx": diff(y, 0), "zx": diff(z, 0),
            "ax": diff(a, 0), "cx": diff(c, 0),
        }
    return result


def vector(values):
    return [values[point][source][stable] % MOD for stable, point, source in ROWS]


def rank(matrix):
    rows = [[entry % MOD for entry in row] for row in matrix]
    pivot_row = 0
    for column in range(len(rows[0])):
        pivot = next((r for r in range(pivot_row, len(rows)) if rows[r][column]), None)
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        coefficient = inv(rows[pivot_row][column])
        rows[pivot_row] = [coefficient * value % MOD for value in rows[pivot_row]]
        for r in range(len(rows)):
            if r != pivot_row and rows[r][column]:
                coefficient = rows[r][column]
                rows[r] = [(a - coefficient * b) % MOD
                           for a, b in zip(rows[r], rows[pivot_row])]
        pivot_row += 1
        if pivot_row == len(rows):
            break
    return pivot_row


def column_rank(columns):
    return rank([list(row) for row in zip(*columns)])


def contains(columns, target):
    rows = [list(row) for row in zip(*columns)]
    return rank(rows) == rank([row + [target[i]] for i, row in enumerate(rows)])


def build(root):
    mixed = packets(root, True)
    affine = packets(root, False)
    columns, labels = [], []
    for kind in range(3):
        for a in range(CUT + 1):
            for b in range(CUT + 1 - a):
                for cpower in range(CUT + 1 - a - b):
                    values = {}
                    for point in POINTS:
                        packet = mixed[point]
                        base = (packet["area"], packet["yx"], packet["zx"])[kind]
                        values[point] = mul(mul(mul(packet["y"][a], packet["z"][b]),
                                                      power(TVAR, cpower)), base)
                    columns.append(vector(values))
                    labels.append((kind, a, b, cpower))

    lookup = {label: columns[index] for index, label in enumerate(labels)}
    exact = []
    for kind in range(3):
        for a in range(CUT + 2):
            for b in range(CUT + 2 - a):
                for cpower in range(CUT + 2 - a - b):
                    terms = []
                    if kind == 0:
                        if b: terms.append((-b, (0, a, b - 1, cpower)))
                        if cpower: terms.append((-cpower, (1, a, b, cpower - 1)))
                    elif kind == 1:
                        if a: terms.append((a, (0, a - 1, b, cpower)))
                        if cpower: terms.append((-cpower, (2, a, b, cpower - 1)))
                    else:
                        if a: terms.append((a, (1, a - 1, b, cpower)))
                        if b: terms.append((b, (2, a, b - 1, cpower)))
                    if terms:
                        value = [0] * len(ROWS)
                        for coefficient, label in terms:
                            value = [(x + coefficient * y) % MOD
                                     for x, y in zip(value, lookup[label])]
                        if any(value): exact.append(value)

    fixed, linearized = [], []
    hprime = add(scalar(1), scale(TVAR, 2))
    for coordinate in ("fixed", "f", "g"):
        for a in range(CUT + 2):
            for b in range(CUT + 2 - a):
                for cpower in range(CUT + 2 - a - b):
                    if a + b + cpower == 0: continue
                    values = {}
                    for point in POINTS:
                        packet = mixed[point] if coordinate == "fixed" else affine[point]
                        density = zero()
                        if coordinate in ("fixed", "g"):
                            if b:
                                multiplier = hprime if coordinate == "fixed" else scalar(1)
                                density = add(density, scale(mul(mul(mul(multiplier, packet["A"][a]),
                                                                          packet["c"][b - 1]),
                                                                     power(TVAR, cpower)), -3 * b))
                            if cpower:
                                density = add(density, scale(mul(mul(mul(packet["ax"], packet["A"][a]),
                                                                          packet["c"][b]),
                                                                     power(TVAR, cpower - 1)), cpower))
                        else:
                            if a:
                                density = add(density, scale(mul(mul(packet["A"][a - 1], packet["c"][b]),
                                                                     power(TVAR, cpower)), 12 * a))
                            if cpower:
                                density = add(density, scale(mul(mul(mul(packet["cx"], packet["A"][a]),
                                                                          packet["c"][b]),
                                                                     power(TVAR, cpower - 1)), 4 * cpower))
                        values[point] = density
                    (fixed if coordinate == "fixed" else linearized).append(vector(values))
    deformation = [-24 % MOD if stable == 1 and source == 0 else 0
                   for stable, _, source in ROWS]
    return columns, exact, fixed, linearized, deformation


# Exact rational setup, with alpha-independent value and derivative packet.
qbase_Q = {0: Fraction(-3, 4), 1: Fraction(1), 2: Fraction(-27, 4),
           3: Fraction(-2), 4: Fraction(9, 2), 5: Fraction(1)}
pbase_Q = {2: Fraction(1), 4: Fraction(-2), 6: Fraction(1)}
q6_Q = {degree: qbase_Q.get(degree, 0) + Fraction(-259, 36) * pbase_Q.get(degree, 0)
        for degree in range(7)}


def evaluate(poly, point, derivative=False):
    if derivative:
        return sum(degree * coefficient * point ** (degree - 1)
                   for degree, coefficient in poly.items() if degree)
    return sum(coefficient * point**degree for degree, coefficient in poly.items())


qvalues = tuple(evaluate(q6_Q, point) for point in POINTS)
qderivatives = tuple(evaluate(q6_Q, point, True) for point in POINTS)
if qvalues != (Fraction(-3), Fraction(-3, 4), Fraction(-3)):
    raise RuntimeError(qvalues)
if qderivatives != (Fraction(-9, 2), Fraction(1), Fraction(9, 2)):
    raise RuntimeError(qderivatives)
if any((72783360*r**4-77822208*r**3-28419741*r**2+7849770*r-1276420) % MOD
       for r in SPLIT_ROOTS):
    raise RuntimeError("bad quartic roots")

xs, qs = sp.symbols("x q")
Ds = 1 + xs**2 * qs
as_ = qs / Ds**2
cs = xs * Ds * (Ds + 2)
bs = (Ds - 1) * (Ds + 2)**2
es = qs * (Ds + 3)
jac_ac = sp.factor(sp.diff(as_, xs) * sp.diff(cs, qs)
                   - sp.diff(as_, qs) * sp.diff(cs, xs))
if jac_ac != -3 or sp.factor(es / (bs + 4) - as_) != 0:
    raise RuntimeError(("local identities", jac_ac))

expected = (59, 59, 57, 59, True, True, False, True)
actual_by_root = []
for root in SPLIT_ROOTS:
    columns, exact, fixed, linearized, deformation = build(root)
    actual = (column_rank(columns), column_rank(exact), column_rank(fixed),
              column_rank(linearized), contains(columns, CONSTANT),
              contains(exact, CONSTANT), contains(fixed, CONSTANT),
              contains(linearized, deformation))
    if actual != expected:
        raise RuntimeError((root, actual, expected))
    actual_by_root.append((root, actual))

summary = "59,59,57,59;1,1,0,1"
print("AFFINE EXCEPTIONAL INDEPENDENT ARRAY AUDIT")
print(f"Qvalues={qvalues};Qderivatives={qderivatives};Qprime0=1=>non_even")
print("local_identity=Jac_xq(a,c)=-3, so (a,-4c) has Jac=12 for H=t")
print("cut5=rows63/all168/exact231/fixed83/linearized166")
print("roots44,82,92,134:ranks=59,59,57,59;spans(all,exact,fixed,mixed-tangent)=1,1,0,1")
print("fixed_a_obstruction=first appears at cut5 for H=t+t^2")
print("mixed_first_order=required -24t lies in affine Darboux tangent image")
print("scope=retained three-multigerm packet, not a global pair or all-order lift")
print("summary_sha256=" + sha256(summary.encode()).hexdigest())
print("RESULT=PASS")
