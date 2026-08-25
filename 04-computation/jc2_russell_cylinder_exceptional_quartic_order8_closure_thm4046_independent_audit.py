"""Independent local-Taylor reconstruction of the THM-4046 order-eight window.

This audit intentionally does not import or execute the production companion.
It works directly in the exceptional quartic field and constructs local
bivariate Taylor series at x=-1,0,1 using sparse dictionary arithmetic.
"""

from __future__ import annotations

import hashlib
import itertools
import time
from math import comb

import sympy as sp
from sympy.polys.matrices import DomainMatrix


T, r, x = sp.symbols("T r x")
F6 = sp.Poly(
    72783360 * T**4
    - 77822208 * T**3
    - 28419741 * T**2
    + 7849770 * T
    - 1276420,
    T,
    domain=sp.QQ,
)
K = sp.QQ.alg_field_from_poly(F6)
alpha = K.convert(K.ext)
ZERO, ONE = K.zero, K.one
POINTS = (-1, 0, 1)
CUTOFF = 8
BUILD_CUTOFF = CUTOFF + 1


def kq(value):
    return K.convert(sp.Rational(value))


def eval_r(value):
    """Evaluate a QQ[r] polynomial at the field generator alpha."""
    poly = sp.Poly(value, r, domain=sp.QQ)
    out = ZERO
    for coefficient in poly.all_coeffs():
        out = out * alpha + K.convert(coefficient)
    return out


def anp_coordinates(value):
    """Return coefficients on (1,alpha,alpha^2,alpha^3)."""
    descending = list(value.to_list())
    descending = [K.dom.zero] * (4 - len(descending)) + descending
    return tuple(reversed(descending))


def sympy_q(value):
    return sp.Rational(value.numerator, value.denominator)


def coordinate_text(value):
    return ",".join(sp.sstr(sympy_q(v)) for v in anp_coordinates(value))


Q1 = (
    x**5
    + sp.Rational(9, 2) * x**4
    - 2 * x**3
    - sp.Rational(27, 4) * x**2
    + x
    - sp.Rational(3, 4)
)
P = x**2 * (x**2 - 1) ** 2
Q6 = sp.expand(Q1 - sp.Rational(259, 36) * P)
R1 = sp.expand(P * (1 - x**2))
R2 = sp.expand(P * (4 - 9 * x))
p_parabola = (
    sp.Rational(520, 9) * r**2
    - sp.Rational(1688, 81) * r
    - sp.Rational(5717, 729)
)
Q_expr = sp.Poly(sp.expand(Q6 + p_parabola * R1 + r * R2), x)
Q_coefficients = {
    degree: eval_r(Q_expr.nth(degree))
    for degree in range(Q_expr.degree() + 1)
}


def add(left, right, cutoff=BUILD_CUTOFF):
    out = dict(left)
    for key, value in right.items():
        if sum(key) <= cutoff:
            out[key] = out.get(key, ZERO) + value
            if not out[key]:
                del out[key]
    return out


def scale(value, scalar, cutoff=BUILD_CUTOFF):
    if not scalar:
        return {}
    return {
        key: coefficient * scalar
        for key, coefficient in value.items()
        if sum(key) <= cutoff and coefficient
    }


def multiply(left, right, cutoff=BUILD_CUTOFF):
    out = {}
    for (i, j), a in left.items():
        for (u, v), b in right.items():
            key = (i + u, j + v)
            if sum(key) > cutoff:
                continue
            out[key] = out.get(key, ZERO) + a * b
    return {key: value for key, value in out.items() if value}


def power(value, exponent, cutoff=CUTOFF):
    out = {(0, 0): ONE}
    for _ in range(exponent):
        out = multiply(out, value, cutoff=cutoff)
    return out


def derivative(value, axis, cutoff=CUTOFF):
    out = {}
    for (i, j), coefficient in value.items():
        exponents = (i, j)
        if not exponents[axis]:
            continue
        key = (i - (axis == 0), j - (axis == 1))
        if sum(key) <= cutoff:
            out[key] = coefficient * K.convert(exponents[axis])
    return out


def truncate(value, cutoff=CUTOFF):
    return {key: coefficient for key, coefficient in value.items() if sum(key) <= cutoff}


def local_q(point):
    out = {}
    for degree, coefficient in Q_coefficients.items():
        for local_degree in range(degree + 1):
            contribution = coefficient * K.convert(comb(degree, local_degree) * point ** (degree - local_degree))
            key = (local_degree, 0)
            out[key] = out.get(key, ZERO) + contribution
    out[(0, 2)] = out.get((0, 2), ZERO) + ONE
    return {key: value for key, value in out.items() if value}


def local_packet(point):
    Xlocal = {(0, 0): K.convert(point), (1, 0): ONE}
    qlocal = local_q(point)
    D = add({(0, 0): ONE}, multiply(power(Xlocal, 2, BUILD_CUTOFF), qlocal, BUILD_CUTOFF), BUILD_CUTOFF)
    Dplus2 = add(D, {(0, 0): K.convert(2)}, BUILD_CUTOFF)
    Dplus3 = add(D, {(0, 0): K.convert(3)}, BUILD_CUTOFF)
    y = scale(multiply(multiply(Xlocal, D, BUILD_CUTOFF), Dplus2, BUILD_CUTOFF), kq(sp.Rational(1, 3)), BUILD_CUTOFF)
    z = add(multiply(qlocal, Dplus3, BUILD_CUTOFF), {(0, 0): K.convert(3)}, BUILD_CUTOFF)

    yx = derivative(y, 0, CUTOFF)
    yt = derivative(y, 1, CUTOFF)
    zx = derivative(z, 0, CUTOFF)
    zt = derivative(z, 1, CUTOFF)
    area = add(
        multiply(yx, zt, CUTOFF),
        scale(multiply(yt, zx, CUTOFF), K.convert(-1), CUTOFF),
        CUTOFF,
    )
    y8, z8 = truncate(y), truncate(z)
    if y8.get((0, 0), ZERO) or z8.get((0, 0), ZERO):
        raise RuntimeError(f"retained origin failed at {point}")
    return (
        [power(y8, degree, CUTOFF) for degree in range(CUTOFF + 1)],
        [power(z8, degree, CUTOFF) for degree in range(CUTOFF + 1)],
        (area, yx, zx),
    )


started = time.time()
PACKETS = {point: local_packet(point) for point in POINTS}
print("packets=three_exact_local_branches", flush=True)


def monomial_pullback(kind, y_degree, z_degree, w_degree, point):
    yp, zp, bases = PACKETS[point]
    value = multiply(yp[y_degree], zp[z_degree], CUTOFF)
    value = multiply(value, {(0, w_degree): ONE}, CUTOFF)
    return multiply(value, bases[kind], CUTOFF)


ROWS = tuple(
    (stable_degree, point, source_degree)
    for stable_degree, max_source_degree in ((0, 4), (2, 3), (4, 2), (6, 1), (8, 0))
    for source_degree in range(max_source_degree + 1)
    for point in POINTS
)
LABELS = tuple(
    (kind, yd, zd, wd)
    for kind in range(3)
    for yd in range(CUTOFF + 1)
    for zd in range(CUTOFF + 1 - yd)
    for wd in range(CUTOFF + 1 - yd - zd)
)
COLUMNS = []
for number, label in enumerate(LABELS, 1):
    kind, yd, zd, wd = label
    branch = {
        point: monomial_pullback(kind, yd, zd, wd, point)
        for point in POINTS
    }
    COLUMNS.append(
        tuple(branch[point].get((source_degree, stable_degree), ZERO) for stable_degree, point, source_degree in ROWS)
    )
    if number % 100 == 0:
        print(f"columns_built={number}", flush=True)

if len(ROWS) != 45 or len(COLUMNS) != 495 or len(set(LABELS)) != 495:
    raise RuntimeError("universe count")

matrix = DomainMatrix([list(column) for column in COLUMNS], (len(COLUMNS), len(ROWS)), K, fmt="dense")
rank = matrix.rank()
kernel_dm = matrix.nullspace(divide_last=True)
kernel = kernel_dm.to_list()
print(f"matrix_rank={rank};nullity={len(kernel)}", flush=True)
if rank != 40 or len(kernel) != 5:
    raise RuntimeError("order-eight rank/nullity")

j8_indices = [i for i, row in enumerate(ROWS) if row[0] == 8]
projection = DomainMatrix(
    [[row[i] for i in j8_indices] for row in kernel],
    (len(kernel), len(j8_indices)),
    K,
    fmt="dense",
)
projection_rank = projection.rank()
if projection_rank != 1:
    raise RuntimeError("J8 projection rank")
active = next(row for row in kernel if any(row[i] for i in j8_indices))
relation = [value * kq(sp.Rational(5, 18)) / active[j8_indices[0]] for value in active]
lambda_expected = tuple(kq(v) for v in (sp.Rational(5, 18), -1, sp.Rational(13, 18)))
if tuple(relation[i] for i in j8_indices) != lambda_expected:
    raise RuntimeError("normalized Lambda block")
for column in COLUMNS:
    if sum((relation[i] * column[i] for i in range(len(ROWS))), ZERO):
        raise RuntimeError("active relation residual")

nonzero = [(ROWS[i], relation[i]) for i in range(len(ROWS)) if relation[i]]
if len(nonzero) != 35:
    raise RuntimeError("active support size")

positive_value_sums = {}
for stable_degree in (2, 4, 6, 8):
    indices = [
        i
        for i, row in enumerate(ROWS)
        if row[0] == stable_degree and row[2] == 0
    ]
    positive_value_sums[stable_degree] = sum((relation[i] for i in indices), ZERO)
    if positive_value_sums[stable_degree]:
        raise RuntimeError(f"J{stable_degree} value block does not kill constants")

constant_indices = [i for i, row in enumerate(ROWS) if row[0] == 0 and row[2] == 0]
kappa = sum((relation[i] for i in constant_indices), ZERO)
expected_kappa = (
    -kq(sp.Rational(5183766767360, 3**19))
    + kq(sp.Rational(931699873280, 3**16)) * alpha
    - kq(sp.Rational(460399208960, 3**14)) * alpha**2
    - kq(sp.Rational(183606968320, 3**13)) * alpha**3
)
if kappa != expected_kappa:
    raise RuntimeError(f"kappa mismatch {coordinate_text(kappa)}")
if not kappa:
    raise RuntimeError("zero kappa")

# Determinant of multiplication by kappa on the power basis is the field norm.
basis = [ONE, alpha, alpha**2, alpha**3]
mult_matrix = sp.Matrix(
    [
        [sympy_q(anp_coordinates(kappa * basis[column])[row]) for column in range(4)]
        for row in range(4)
    ]
)
norm = sp.factor(mult_matrix.det())
if norm == 0:
    raise RuntimeError("zero norm")

# Reconstruct the complete inherited order-six cokernel inside this order-eight window.
ROWS6 = tuple(
    (stable_degree, point, source_degree)
    for stable_degree, max_source_degree in ((0, 3), (2, 2), (4, 1), (6, 0))
    for source_degree in range(max_source_degree + 1)
    for point in POINTS
)
row6_indices = [ROWS.index(row) for row in ROWS6]
column6_indices = [i for i, label in enumerate(LABELS) if sum(label[1:]) <= 6]
matrix6 = DomainMatrix(
    [[COLUMNS[column][row] for row in row6_indices] for column in column6_indices],
    (len(column6_indices), len(ROWS6)),
    K,
    fmt="dense",
)
if matrix6.rank() != 26:
    raise RuntimeError("inherited order-six rank")
kernel6 = matrix6.nullspace(divide_last=True).to_list()
if len(kernel6) != 4:
    raise RuntimeError("inherited order-six nullity")
embedded = []
for row6 in kernel6:
    row8 = [ZERO] * len(ROWS)
    for local, global_index in enumerate(row6_indices):
        row8[global_index] = row6[local]
    for column in COLUMNS:
        if sum((row8[i] * column[i] for i in range(len(ROWS))), ZERO):
            raise RuntimeError("inherited relation does not extend")
    embedded.append(row8)
span = DomainMatrix([*embedded, relation], (5, len(ROWS)), K, fmt="dense")
if span.rank() != 5:
    raise RuntimeError("active quotient independence")

payload = ";".join(
    f"{stable},{point},{source}:{coordinate_text(value)}"
    for (stable, point, source), value in nonzero
)
relation_hash = hashlib.sha256(payload.encode("ascii")).hexdigest()
primary_gauge_payload = ";".join(
    f"J{stable}@{point}:d{source}=({coordinate_text(value)})"
    for (stable, point, source), value in nonzero
)
primary_gauge_hash = hashlib.sha256(primary_gauge_payload.encode("ascii")).hexdigest()

print("universe=three_two_form_kinds_times_target_monomials_total_degree_le_8")
print("rows=J0_xle4,J2_xle3,J4_xle2,J6_xle1,J8_values")
print(f"order8=rows45_columns495_rank{rank}_nullity{len(kernel)}_J8projection{projection_rank}_active_nonzero{len(nonzero)}")
print("J8_block=5/18,-1,13/18")
print(f"kappa_coordinates={coordinate_text(kappa)}")
print(f"kappa_norm={sp.sstr(norm)}")
print(f"kappa_norm_factor={sp.factorint(abs(sp.numer(norm)))} / {sp.factorint(sp.denom(norm))}")
print("positive_value_block_sums=J2:0,J4:0,J6:0,J8:0")
print("order6_embedded=rank26_nullity4;active_quotient_dimension1")
print(f"active_relation_sha256={relation_hash}")
print(f"primary_serialization_sha256={primary_gauge_hash}")
print("RESULT=PASS")
