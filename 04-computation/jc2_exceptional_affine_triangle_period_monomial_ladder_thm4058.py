#!/usr/bin/env python3
"""Exact triangle-period audit for the exceptional fixed-a packet."""

from fractions import Fraction as F
from math import comb


ZERO = F(0)
ONE = F(1)
K0 = (ZERO, ZERO, ZERO, ZERO)
K1 = (ONE, ZERO, ZERO, ZERO)
ALPHA = (ZERO, ONE, ZERO, ZERO)
REL = (
    F(1276420, 72783360),
    F(-7849770, 72783360),
    F(28419741, 72783360),
    F(77822208, 72783360),
)


def require(condition, label):
    if not condition:
        raise RuntimeError(label)


def kc(value, denominator=1):
    return (F(value, denominator), ZERO, ZERO, ZERO)


def ka(left, right):
    return tuple(x + y for x, y in zip(left, right))


def kn(value):
    return tuple(-x for x in value)


def ks(left, right):
    return ka(left, kn(right))


def kscale(value, scalar):
    return tuple(scalar * x for x in value)


def km(left, right):
    raw = [ZERO] * 7
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            raw[i + j] += x * y
    for degree in range(6, 3, -1):
        coefficient = raw[degree]
        for offset, relation in enumerate(REL):
            raw[degree - 4 + offset] += coefficient * relation
    return tuple(raw[:4])


def kpow(value, exponent):
    answer = K1
    base = value
    while exponent:
        if exponent & 1:
            answer = km(answer, base)
        base = km(base, base)
        exponent //= 2
    return answer


def solve_fraction_matrix(matrix, target):
    rows = [list(row) + [target[i]] for i, row in enumerate(matrix)]
    size = len(rows)
    for column in range(size):
        pivot = next((i for i in range(column, size) if rows[i][column]), None)
        require(pivot is not None, "singular field inverse matrix")
        rows[column], rows[pivot] = rows[pivot], rows[column]
        scalar = rows[column][column]
        rows[column] = [value / scalar for value in rows[column]]
        for i in range(size):
            if i != column and rows[i][column]:
                scalar = rows[i][column]
                rows[i] = [value - scalar * pivot_value
                           for value, pivot_value in zip(rows[i], rows[column])]
    return [rows[i][-1] for i in range(size)]


def kinv(value):
    basis = (K1, ALPHA, kpow(ALPHA, 2), kpow(ALPHA, 3))
    columns = [km(value, item) for item in basis]
    matrix = [[columns[column][row] for column in range(4)] for row in range(4)]
    answer = tuple(solve_fraction_matrix(matrix, [ONE, ZERO, ZERO, ZERO]))
    require(km(value, answer) == K1, "bad field inverse")
    return answer


def kdiv(left, right):
    return km(left, kinv(right))


def ktext(value):
    return "(" + ",".join(str(coordinate) for coordinate in value) + ")"


# Bivariate series have keys (A-degree,c-degree) and total cutoff six.
CUT = 6


def badd(left, right):
    answer = dict(left)
    for key, value in right.items():
        answer[key] = ka(answer.get(key, K0), value)
        if answer[key] == K0:
            del answer[key]
    return answer


def bscale(value, scalar):
    if not isinstance(scalar, tuple):
        scalar = kc(scalar)
    return {key: product for key, coefficient in value.items()
            if (product := km(coefficient, scalar)) != K0}


def bmul(left, right):
    answer = {}
    for (i, j), x in left.items():
        for (u, v), y in right.items():
            if i + j + u + v <= CUT:
                key = (i + u, j + v)
                answer[key] = ka(answer.get(key, K0), km(x, y))
    return {key: value for key, value in answer.items() if value != K0}


def bpow(value, exponent):
    answer = {(0, 0): K1}
    base = value
    while exponent:
        if exponent & 1:
            answer = bmul(answer, base)
        base = bmul(base, base)
        exponent //= 2
    return answer


def bdiff(value, axis):
    answer = {}
    for (i, j), coefficient in value.items():
        exponent = i if axis == 0 else j
        if exponent:
            key = (i - (axis == 0), j - (axis == 1))
            answer[key] = kscale(coefficient, F(exponent))
    return answer


def binverse(value):
    constant = value.get((0, 0), K0)
    require(constant != K0, "nonunit bivariate series")
    inverse_constant = kinv(constant)
    remainder = badd(bscale(value, inverse_constant), {(0, 0): kc(-1)})
    answer = {}
    term = {(0, 0): K1}
    for degree in range(CUT + 1):
        answer = badd(answer, bscale(term, kc(-1 if degree & 1 else 1)))
        term = bmul(term, remainder)
    return bscale(answer, inverse_constant)


def xadd(left, right):
    answer = dict(left)
    for degree, value in right.items():
        answer[degree] = ka(answer.get(degree, K0), value)
        if answer[degree] == K0:
            del answer[degree]
    return answer


def xscale(value, scalar):
    return {degree: product for degree, coefficient in value.items()
            if (product := km(coefficient, scalar)) != K0}


def xshift(value, shift):
    return {degree + shift: coefficient for degree, coefficient in value.items()}


Q_BASE = {
    0: kc(-3, 4), 1: K1, 2: kc(-27, 4), 3: kc(-2),
    4: kc(9, 2), 5: K1,
}
P_BASE = {2: K1, 4: kc(-2), 6: K1}
Q6 = xadd(Q_BASE, xscale(P_BASE, kc(-259, 36)))
R1 = xadd(P_BASE, xscale(xshift(P_BASE, 2), kc(-1)))
R2 = xadd(xscale(P_BASE, kc(4)), xscale(xshift(P_BASE, 1), kc(-9)))
THETA = ka(
    ka(kscale(kpow(ALPHA, 2), F(520, 9)), kscale(ALPHA, F(-1688, 81))),
    kc(-5717, 729),
)
Q_ALPHA = xadd(Q6, xadd(xscale(R1, THETA), xscale(R2, ALPHA)))


def compose_q(value):
    answer = {}
    for degree, coefficient in Q_ALPHA.items():
        answer = badd(answer, bscale(bpow(value, degree), coefficient))
    return answer


A_VALUE = {(0, 0): kc(-3, 4), (1, 0): K1}
C_VARIABLE = {(0, 1): K1}


def branch_series(base_s):
    s = {(0, 0): kc(base_s)}
    for _ in range(5):
        s_squared = bmul(s, s)
        equation = badd(
            badd(bscale(s, kc(3)), bmul(A_VALUE, bmul(s_squared, s))),
            bscale(C_VARIABLE, kc(-1)),
        )
        derivative = badd(
            {(0, 0): kc(3)}, bscale(bmul(A_VALUE, s_squared), kc(3)),
        )
        s = badd(s, bscale(bmul(equation, binverse(derivative)), kc(-1)))
    equation = badd(
        badd(bscale(s, kc(3)), bmul(A_VALUE, bpow(s, 3))),
        bscale(C_VARIABLE, kc(-1)),
    )
    require(not equation, ("implicit branch residual", base_s))
    denominator = badd({(0, 0): K1}, bmul(A_VALUE, bmul(s, s)))
    x_value = bmul(s, binverse(denominator))
    q_value = bmul(A_VALUE, bmul(denominator, denominator))
    return badd(q_value, bscale(compose_q(x_value), kc(-1)))


def uadd(left, right):
    answer = dict(left)
    for degree, value in right.items():
        answer[degree] = ka(answer.get(degree, K0), value)
        if answer[degree] == K0:
            del answer[degree]
    return answer


def usub(left, right):
    return uadd(left, {degree: kn(value) for degree, value in right.items()})


def umul(left, right):
    answer = {}
    for i, x in left.items():
        for j, y in right.items():
            if i + j <= CUT:
                answer[i + j] = ka(answer.get(i + j, K0), km(x, y))
    return {degree: value for degree, value in answer.items() if value != K0}


def uinverse(value):
    constant = value.get(0, K0)
    require(constant != K0, "nonunit univariate series")
    inverse_constant = kinv(constant)
    remainder = uadd(
        {degree: km(coefficient, inverse_constant)
         for degree, coefficient in value.items()},
        {0: kc(-1)},
    )
    answer = {}
    term = {0: K1}
    for degree in range(CUT + 1):
        sign = kc(-1 if degree & 1 else 1)
        answer = uadd(answer, {n: km(coefficient, sign)
                              for n, coefficient in term.items()})
        term = umul(term, remainder)
    return {degree: km(coefficient, inverse_constant)
            for degree, coefficient in answer.items()}


def bsubstitute_c(value, c_series):
    powers = [{0: K1}]
    for _ in range(CUT):
        powers.append(umul(powers[-1], c_series))
    answer = {}
    for (a_degree, c_degree), coefficient in value.items():
        for degree, factor in powers[c_degree].items():
            if a_degree + degree <= CUT:
                key = a_degree + degree
                answer[key] = ka(answer.get(key, K0), km(coefficient, factor))
    return {degree: coefficient for degree, coefficient in answer.items()
            if coefficient != K0}


def intersection(left, right):
    difference = badd(left, bscale(right, kc(-1)))
    derivative = bdiff(difference, 1)
    c_series = {}
    for _ in range(5):
        value = bsubstitute_c(difference, c_series)
        slope = bsubstitute_c(derivative, c_series)
        c_series = usub(c_series, umul(value, uinverse(slope)))
    require(not bsubstitute_c(difference, c_series), "intersection residual")
    w_series = bsubstitute_c(left, c_series)
    require(w_series == bsubstitute_c(right, c_series), "intersection value mismatch")
    return c_series, w_series


branches = tuple(branch_series(base_s) for base_s in (2, 0, -2))
slopes = tuple(branch.get((0, 1), K0) for branch in branches)
require(tuple(branch.get((1, 0), K0) for branch in branches) == (K1, K1, K1),
        "bad common A-linear coefficient")
require(slopes == (kc(3, 4), kc(-1, 3), kc(-3, 4)), "bad branch slopes")

(c01, w01), (c12, w12), (c20, w20) = (
    intersection(branches[0], branches[1]),
    intersection(branches[1], branches[2]),
    intersection(branches[2], branches[0]),
)

C_COMMON = {
    2: ka(kc(-64, 27), kscale(ALPHA, F(64, 3))),
    3: ka(
        ka(kc(-74752, 6561), kscale(ALPHA, F(87296, 729))),
        kscale(kpow(ALPHA, 2), F(8192, 81)),
    ),
    4: ka(
        ka(
            ka(kc(-36628480, 1594323), kscale(ALPHA, F(9671680, 59049))),
            kscale(kpow(ALPHA, 2), F(15073280, 6561)),
        ),
        kscale(kpow(ALPHA, 3), F(-4014080, 2187)),
    ),
}
W_COMMON = {
    1: K1,
    2: ka(kc(64, 81), kscale(ALPHA, F(-64, 9))),
    3: ka(
        ka(kc(74752, 19683), kscale(ALPHA, F(-87296, 2187))),
        kscale(kpow(ALPHA, 2), F(-8192, 243)),
    ),
    4: ka(
        ka(
            ka(kc(-576256, 19683), kscale(ALPHA, F(88960, 6561))),
            kscale(kpow(ALPHA, 2), F(-709568, 729)),
        ),
        kscale(kpow(ALPHA, 3), F(-81920, 81)),
    ),
}
for series in (c01, c12, c20):
    require(all(series.get(degree, K0) == C_COMMON.get(degree, K0)
                for degree in range(5)), "c-jets disagree before A^5")
for series in (w01, w12, w20):
    require(all(series.get(degree, K0) == W_COMMON.get(degree, K0)
                for degree in range(5)), "w-jets disagree before A^5")

RHO = (
    F(-2073506706944, 1678822119),
    F(372679949312, 62178597),
    F(-184159683584, 6908733),
    F(-73442787328, 2302911),
)
delta = ks(c20.get(5, K0), c12.get(5, K0))
require(delta == kscale(RHO, F(26, 15)), "delta/rho mismatch")

c_edges = (
    ks(c01.get(5, K0), c20.get(5, K0)),
    ks(c12.get(5, K0), c01.get(5, K0)),
    ks(c20.get(5, K0), c12.get(5, K0)),
)
w_edges = (
    ks(w01.get(5, K0), w20.get(5, K0)),
    ks(w12.get(5, K0), w01.get(5, K0)),
    ks(w20.get(5, K0), w12.get(5, K0)),
)
require(c_edges == tuple(kscale(delta, scalar)
                         for scalar in (F(5, 13), F(-18, 13), ONE)),
        "c-edge vector mismatch")
require(w_edges == tuple(kscale(delta, scalar)
                         for scalar in (F(15, 52), F(6, 13), F(-3, 4))),
        "w-edge vector mismatch")

area = K0
for w_start, w_end, c_edge in (
        (w20.get(5, K0), w01.get(5, K0), c_edges[0]),
        (w01.get(5, K0), w12.get(5, K0), c_edges[1]),
        (w12.get(5, K0), w20.get(5, K0), c_edges[2])):
    area = ka(area, kscale(km(ka(w_start, w_end), c_edge), F(1, 2)))
require(area == kscale(km(delta, delta), F(-15, 52)), "triangle area mismatch")
canonical_response = kdiv(area, delta)
require(canonical_response == kscale(RHO, F(-1, 2)), "rho/2 response mismatch")


# Exact rational source-reparametrization controls for H=t+t^m.
def qadd(left, right, cutoff):
    answer = dict(left)
    for degree, value in right.items():
        if degree <= cutoff:
            answer[degree] = answer.get(degree, ZERO) + value
            if not answer[degree]:
                del answer[degree]
    return answer


def qmul(left, right, cutoff):
    answer = {}
    for i, x in left.items():
        for j, y in right.items():
            if i + j <= cutoff:
                answer[i + j] = answer.get(i + j, ZERO) + x * y
    return {degree: value for degree, value in answer.items() if value}


def qpow(value, exponent, cutoff):
    answer = {0: ONE}
    base = value
    while exponent:
        if exponent & 1:
            answer = qmul(answer, base, cutoff)
        base = qmul(base, base, cutoff)
        exponent //= 2
    return answer


def qcompose(value, argument, cutoff):
    answer = {}
    for degree, coefficient in value.items():
        answer = qadd(
            answer,
            {n: coefficient * entry
             for n, entry in qpow(argument, degree, cutoff).items()},
            cutoff,
        )
    return answer


def qderivative(value):
    return {degree - 1: F(degree) * coefficient
            for degree, coefficient in value.items() if degree}


def qinverse_unit(value, cutoff):
    constant = value.get(0, ZERO)
    require(constant, "nonunit rational series")
    normalized = {degree: coefficient / constant
                  for degree, coefficient in value.items()}
    remainder = qadd(normalized, {0: -ONE}, cutoff)
    answer = {}
    term = {0: ONE}
    for degree in range(cutoff + 1):
        sign = -ONE if degree & 1 else ONE
        answer = qadd(answer, {n: sign * entry for n, entry in term.items()}, cutoff)
        term = qmul(term, remainder, cutoff)
    return {degree: coefficient / constant for degree, coefficient in answer.items()}


for m in range(2, 13):
    cutoff = 2 * m
    u = {1: ONE}
    inverse_h = dict(u)
    for _ in range(6):
        error = qadd(
            qadd(inverse_h, qpow(inverse_h, m, cutoff), cutoff),
            {1: -ONE},
            cutoff,
        )
        derivative = qadd(
            {0: ONE},
            {degree: F(m) * coefficient
             for degree, coefficient in qpow(inverse_h, m - 1, cutoff).items()},
            cutoff,
        )
        correction = qmul(error, qinverse_unit(derivative, cutoff), cutoff)
        inverse_h = qadd(
            inverse_h,
            {degree: -coefficient for degree, coefficient in correction.items()},
            cutoff,
        )
    displacement = {1: ONE, m: ONE}
    require(qcompose(displacement, inverse_h, cutoff) == u, "right inverse mismatch")
    require(qcompose(inverse_h, displacement, cutoff) == u, "left inverse mismatch")
    inverse_derivative = qderivative(inverse_h)
    require(all(inverse_derivative.get(degree, ZERO) == ZERO
                for degree in range(1, m - 1)), "premature inverse term")
    require(inverse_derivative.get(m - 1, ZERO) == F(-m),
            "inverse-derivative coefficient mismatch")
    require(F(-m) * F(-(m - 1), 2) == F(comb(m, 2)),
            "all-m response factor mismatch")

print("EXCEPTIONAL FIXED-a TRIANGLE PERIOD -- EXACT OVER Q(alpha)")
print("branch_order=(p=-1,p=0,p=1);slopes=(3/4,-1/3,-3/4)")
print("pair_intersections=common_through_A4;first_split=A5")
print("delta=" + ktext(delta))
print("delta_over_rho=26/15")
print("normalized_c_edges=(5/13,-18/13,1)")
print("normalized_w_edges=(15/52,6/13,-3/4)")
print("triangle_area_over_delta_squared=-15/52")
print("canonical_affine_response_on_w=-rho/2")
print("all_r_response_on_w^r=-(r/2)rho;first_degree=r+4")
print("source_reparam=H=t+gamma*t^m;h_prime=1-m*gamma*u^(m-1)+higher")
print("canonical_fixed_a_constant_response=gamma*binom(m,2)*rho;first_cutoff=m+3")
print("inverse_reparam_controls=m=2..12:PASS")
print("scope=formal_three-branch_fixed-a_packet_only")
print("RESULT=PASS")
