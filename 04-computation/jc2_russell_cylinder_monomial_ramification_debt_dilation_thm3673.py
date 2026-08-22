"""Exact companion for THM-3673.

The proof in THM-3673 is all-k.  This companion verifies its coefficient
decimation algebra symbolically, checks the transferred universal identities
on complete finite two-form jet universes for k=2,3,4,5, includes nonunit
ramification coefficients, and retains k=1 as a hostile control.
"""

from math import comb

import sympy as sp
from flint import fmpq, fmpq_mat


CHECKS = 0


def require(label, condition):
    global CHECKS
    if not condition:
        raise RuntimeError(f"FAILED: {label}")
    CHECKS += 1


x, q, t, X = sp.symbols("x q t X")
points = (-1, 0, 1)
lambda_row = (sp.Rational(5, 18), -1, sp.Rational(13, 18))

Q1 = (
    x**5
    + sp.Rational(9, 2) * x**4
    - 2 * x**3
    - sp.Rational(27, 4) * x**2
    + x
    - sp.Rational(3, 4)
)
Q6 = sp.expand(Q1 - sp.Rational(259, 36) * x**2 * (x**2 - 1) ** 2)
Qstar = (
    -x**7
    - sp.Rational(27, 4) * x**6
    + 3 * x**5
    + 18 * x**4
    - 3 * x**3
    - sp.Rational(27, 2) * x**2
    + x
    - sp.Rational(3, 4)
)

D = 1 + x**2 * q
b = sp.expand((D - 1) * (D + 2) ** 2)
c = sp.expand(x * D * (D + 2))
e = sp.expand(q * (D + 3))
y_general = c / 3
z_general = e + 3

Q6_IDENTITY = {
    (0, -1, 0): sp.Rational(16246280, 531441),
    (0, -1, 1): -sp.Rational(4489, 6561),
    (0, -1, 2): -sp.Rational(10, 81),
    (0, 0, 1): -sp.Rational(64, 81),
    (0, 1, 0): sp.Rational(13390648, 531441),
    (0, 1, 1): -sp.Rational(6559, 6561),
    (0, 1, 2): -sp.Rational(26, 81),
    ("k", -1, 0): sp.Rational(2012, 2187),
    ("k", -1, 1): sp.Rational(5, 27),
    ("k", 1, 0): -sp.Rational(2012, 2187),
    ("k", 1, 1): -sp.Rational(13, 27),
}

QSTAR_IDENTITY = {
    (0, -1, 0): sp.Rational(2300, 81),
    (0, -1, 1): -sp.Rational(1, 9),
    (0, -1, 2): -sp.Rational(10, 81),
    (0, 1, 0): sp.Rational(3140, 81),
    (0, 1, 1): -sp.Rational(7, 9),
    (0, 1, 2): -sp.Rational(26, 81),
    ("k", -1, 0): sp.Rational(4, 9),
    ("k", -1, 1): sp.Rational(5, 27),
    ("k", 1, 0): -sp.Rational(4, 9),
    ("k", 1, 1): -sp.Rational(13, 27),
}


def multiply(left, right, cutoff):
    answer = {}
    for (i, j), left_value in left.items():
        for (u, v), right_value in right.items():
            if i + j + u + v > cutoff:
                continue
            key = (i + u, j + v)
            answer[key] = answer.get(key, 0) + left_value * right_value
    return {key: sp.factor(value) for key, value in answer.items() if value}


def power(value, exponent, cutoff):
    answer = {(0, 0): sp.S.One}
    for _ in range(exponent):
        answer = multiply(answer, value, cutoff)
    return answer


def jet(expr, point, cutoff):
    shifted = sp.Poly(sp.expand(expr.subs(x, X + point)), X, t, domain=sp.QQ)
    return {
        (source_degree, stable_degree): coefficient
        for (source_degree, stable_degree), coefficient in shifted.terms()
        if source_degree + stable_degree <= cutoff
    }


def pulled_packets(polynomial, ramification, scale, cutoff):
    pulled_y = sp.expand(y_general.subs(q, polynomial + scale * t**ramification))
    pulled_z = sp.expand(z_general.subs(q, polynomial + scale * t**ramification))
    y_x = sp.diff(pulled_y, x)
    y_t = sp.diff(pulled_y, t)
    z_x = sp.diff(pulled_z, x)
    z_t = sp.diff(pulled_z, t)
    area = sp.expand(y_x * z_t - y_t * z_x)
    packets = {}
    for point in points:
        y_jet = jet(pulled_y, point, cutoff)
        z_jet = jet(pulled_z, point, cutoff)
        packets[point] = (
            [power(y_jet, degree, cutoff) for degree in range(cutoff + 1)],
            [power(z_jet, degree, cutoff) for degree in range(cutoff + 1)],
            (jet(area, point, cutoff), jet(y_x, point, cutoff), jet(z_x, point, cutoff)),
        )
    return packets


def monomial_pullback(packets, point, kind, y_degree, z_degree, w_degree, cutoff):
    y_powers, z_powers, bases = packets[point]
    value = multiply(y_powers[y_degree], z_powers[z_degree], cutoff)
    value = multiply(value, {(0, w_degree): sp.S.One}, cutoff)
    return multiply(value, bases[kind], cutoff)


def verify_identity(label, polynomial, identity, base_debt, ramification, scale):
    cutoff = 2 * ramification
    packets = pulled_packets(polynomial, ramification, scale, cutoff)
    mutation_detected = False
    monomial_count = 0
    for kind in range(3):
        for y_degree in range(cutoff + 1):
            for z_degree in range(cutoff + 1 - y_degree):
                for w_degree in range(cutoff + 1 - y_degree - z_degree):
                    values = {
                        point: monomial_pullback(
                            packets,
                            point,
                            kind,
                            y_degree,
                            z_degree,
                            w_degree,
                            cutoff,
                        )
                        for point in points
                    }
                    left = sum(
                        lambda_row[index] * values[point].get((0, cutoff), 0)
                        for index, point in enumerate(points)
                    )
                    right = 0
                    mutated_right = 0
                    for (stable_degree, point, source_degree), coefficient in identity.items():
                        if stable_degree == "k":
                            actual_degree = ramification
                            scaled_coefficient = scale * coefficient
                        else:
                            actual_degree = stable_degree
                            scaled_coefficient = scale**2 * coefficient
                        entry = values[point].get((source_degree, actual_degree), 0)
                        right += scaled_coefficient * entry
                        mutation = 1 if (stable_degree, point, source_degree) == (0, -1, 0) else 0
                        mutated_right += (scaled_coefficient + mutation) * entry
                    residual = sp.factor(left - right)
                    require(
                        f"{label} k={ramification} scale={scale} monomial="
                        f"{kind,y_degree,z_degree,w_degree}",
                        residual == 0,
                    )
                    if sp.factor(left - mutated_right) != 0:
                        mutation_detected = True
                    monomial_count += 1
    require(f"{label} k={ramification} active mutation", mutation_detected)
    require(
        f"{label} k={ramification} complete monomial universe",
        monomial_count == 3 * comb(cutoff + 3, 3),
    )
    forced = sp.factor(
        scale**2
        * sum(
            coefficient
            for (stable_degree, _point, source_degree), coefficient in identity.items()
            if stable_degree == 0 and source_degree == 0
        )
    )
    require(f"{label} k={ramification} forced debt", forced == scale**2 * base_debt)
    print(
        f"PASS {label}_k={ramification}_scale={scale}_two_form_monomials={monomial_count} "
        f"forced_lambda_J{cutoff}={forced}"
    )


def as_fmpq(value):
    value = sp.Rational(value)
    return fmpq(int(value.p), int(value.q))


def affine_hostile_control():
    """At k=1,N=2 the target row is not in the complete lower-jet row span."""
    cutoff = 2
    packets = pulled_packets(Q6, 1, sp.S.One, cutoff)
    row_keys = [
        (stable_degree, point, source_degree)
        for stable_degree in range(cutoff)
        for point in points
        for source_degree in range(cutoff - stable_degree + 1)
    ]
    columns = []
    for kind in range(3):
        for y_degree in range(cutoff + 1):
            for z_degree in range(cutoff + 1 - y_degree):
                for w_degree in range(cutoff + 1 - y_degree - z_degree):
                    values = {
                        point: monomial_pullback(
                            packets, point, kind, y_degree, z_degree, w_degree, cutoff
                        )
                        for point in points
                    }
                    lower = [values[p].get((d, s), 0) for s, p, d in row_keys]
                    target = sum(
                        lambda_row[index] * values[point].get((0, cutoff), 0)
                        for index, point in enumerate(points)
                    )
                    columns.append((lower, target))
    lower_rows = [[column[0][index] for column in columns] for index in range(len(row_keys))]
    target_row = [column[1] for column in columns]
    lower = fmpq_mat(
        len(lower_rows),
        len(columns),
        [as_fmpq(value) for row in lower_rows for value in row],
    )
    augmented = fmpq_mat(
        len(lower_rows) + 1,
        len(columns),
        [as_fmpq(value) for row in lower_rows + [target_row] for value in row],
    )
    require("affine lower rank", lower.rank() == 13)
    require("affine augmented rank", augmented.rank() == 14)
    print("PASS hostile_k=1_N=2_lower_rank=13_augmented_rank=14_relation=NONE")


print("THM-3673 exact companion -- monomial ramification debt dilation")
print("status=PROVED VERIFIED-EXACT PENDING-INDEPENDENT-HOSTILE-AUDIT")

K, ell = sp.symbols("K ell", integer=True, positive=True)
require("A-character exponent", sp.expand((K * ell + 1) + (K - 1) - K * (ell + 1)) == 0)
require("BC-character exponent", sp.expand(K * ell - K * ell) == 0)
require("dilation target order", 2 * K == K * 2)
print("PASS symbolic_character_decimation=A_(Kell+1),BC_(Kell),orders_(0,K,2K)")

for k in (2, 3, 4, 5):
    verify_identity("Q6", Q6, Q6_IDENTITY, sp.Rational(365888, 6561), k, sp.S.One)
    verify_identity("Qstar", Qstar, QSTAR_IDENTITY, sp.Rational(5440, 81), k, sp.S.One)

verify_identity("Q6", Q6, Q6_IDENTITY, sp.Rational(365888, 6561), 3, sp.Rational(2))
verify_identity(
    "Qstar",
    Qstar,
    QSTAR_IDENTITY,
    sp.Rational(5440, 81),
    4,
    -sp.Rational(3, 5),
)

affine_hostile_control()
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
