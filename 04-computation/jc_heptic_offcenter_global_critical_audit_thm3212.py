#!/usr/bin/env python3
"""Exact quotient-algebra audit of the THM-3212 off-center control.

No repository executable is imported.  Arithmetic is implemented directly in
Q[u]/(q(u)) and then in that algebra's polynomial ring in x.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import hashlib

import sympy as sp


Q = Fraction


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


@dataclass(frozen=True)
class Alg:
    c: tuple[Q, Q, Q]

    # Filled per passport before constructing elements.
    relation: tuple[Q, Q, Q] = (Q(0), Q(0), Q(0))

    @staticmethod
    def make(a=0, b=0, c=0) -> "Alg":
        return Alg((Q(a), Q(b), Q(c)))

    def __add__(self, other) -> "Alg":
        other = as_alg(other)
        return Alg(tuple(self.c[i] + other.c[i] for i in range(3)))

    __radd__ = __add__

    def __neg__(self) -> "Alg":
        return Alg(tuple(-value for value in self.c))

    def __sub__(self, other) -> "Alg":
        return self + (-as_alg(other))

    def __rsub__(self, other) -> "Alg":
        return as_alg(other) - self

    def __mul__(self, other) -> "Alg":
        other = as_alg(other)
        raw = [Q(0)] * 5
        for i, left in enumerate(self.c):
            for j, right in enumerate(other.c):
                raw[i + j] += left * right
        # u^3 = relation[0] + relation[1] u + relation[2] u^2.
        r0, r1, r2 = Alg.relation
        for degree in (4, 3):
            value = raw[degree]
            raw[degree] = Q(0)
            raw[degree - 3] += value * r0
            raw[degree - 2] += value * r1
            raw[degree - 1] += value * r2
        return Alg(tuple(raw[:3]))

    __rmul__ = __mul__

    def __pow__(self, exponent: int) -> "Alg":
        require(exponent >= 0, "negative algebra exponent")
        result = ONE_A
        base = self
        power = exponent
        while power:
            if power & 1:
                result = result * base
            base = base * base
            power //= 2
        return result

    def inverse(self) -> "Alg":
        require(self != ZERO_A, "division by zero in accessory algebra")
        basis = (ONE_A, U_A, U_A * U_A)
        columns = [(self * item).c for item in basis]
        matrix = [
            [columns[column][row] for column in range(3)]
            + [Q(1 if row == 0 else 0)]
            for row in range(3)
        ]
        for pivot in range(3):
            row = next((r for r in range(pivot, 3) if matrix[r][pivot]), None)
            require(row is not None, "nonunit accessory element")
            matrix[pivot], matrix[row] = matrix[row], matrix[pivot]
            scale = matrix[pivot][pivot]
            matrix[pivot] = [entry / scale for entry in matrix[pivot]]
            for r in range(3):
                if r == pivot:
                    continue
                scale = matrix[r][pivot]
                if scale:
                    matrix[r] = [
                        matrix[r][j] - scale * matrix[pivot][j]
                        for j in range(4)
                    ]
        return Alg(tuple(matrix[row][3] for row in range(3)))

    def __truediv__(self, other) -> "Alg":
        return self * as_alg(other).inverse()

    def __rtruediv__(self, other) -> "Alg":
        return as_alg(other) / self

    def is_zero(self) -> bool:
        return self == ZERO_A


def as_alg(value) -> Alg:
    if isinstance(value, Alg):
        return value
    return Alg.make(value)


ZERO_A = Alg.make(0)
ONE_A = Alg.make(1)
U_A = Alg.make(0, 1)


@dataclass(frozen=True)
class Poly:
    c: tuple[Alg, ...]

    @staticmethod
    def make(coefficients) -> "Poly":
        values = [as_alg(value) for value in coefficients]
        while values and values[-1].is_zero():
            values.pop()
        return Poly(tuple(values))

    @property
    def degree(self) -> int:
        return len(self.c) - 1

    @property
    def lc(self) -> Alg:
        require(self.degree >= 0, "zero polynomial has no leading coefficient")
        return self.c[-1]

    def is_zero(self) -> bool:
        return self.degree < 0

    def __add__(self, other) -> "Poly":
        other = as_poly(other)
        size = max(len(self.c), len(other.c))
        return Poly.make(
            (self.c[i] if i < len(self.c) else ZERO_A)
            + (other.c[i] if i < len(other.c) else ZERO_A)
            for i in range(size)
        )

    __radd__ = __add__

    def __neg__(self) -> "Poly":
        return Poly.make(-value for value in self.c)

    def __sub__(self, other) -> "Poly":
        return self + (-as_poly(other))

    def __rsub__(self, other) -> "Poly":
        return as_poly(other) - self

    def __mul__(self, other) -> "Poly":
        other = as_poly(other)
        if self.is_zero() or other.is_zero():
            return ZERO_P
        raw = [ZERO_A] * (len(self.c) + len(other.c) - 1)
        for i, left in enumerate(self.c):
            for j, right in enumerate(other.c):
                raw[i + j] = raw[i + j] + left * right
        return Poly.make(raw)

    __rmul__ = __mul__

    def __pow__(self, exponent: int) -> "Poly":
        require(exponent >= 0, "negative polynomial exponent")
        result = ONE_P
        base = self
        power = exponent
        while power:
            if power & 1:
                result = result * base
            base = base * base
            power //= 2
        return result

    def scalar(self, value) -> "Poly":
        factor = as_alg(value)
        return Poly.make(coefficient * factor for coefficient in self.c)

    def derivative(self) -> "Poly":
        return Poly.make(self.c[i] * i for i in range(1, len(self.c)))

    def divmod(self, divisor: "Poly") -> tuple["Poly", "Poly"]:
        require(not divisor.is_zero(), "polynomial division by zero")
        remainder = list(self.c)
        quotient = [ZERO_A] * max(0, self.degree - divisor.degree + 1)
        inverse_lc = divisor.lc.inverse()
        while remainder and len(remainder) - 1 >= divisor.degree:
            shift = len(remainder) - 1 - divisor.degree
            factor = remainder[-1] * inverse_lc
            quotient[shift] = factor
            for index, coefficient in enumerate(divisor.c):
                remainder[index + shift] = remainder[index + shift] - factor * coefficient
            while remainder and remainder[-1].is_zero():
                remainder.pop()
        return Poly.make(quotient), Poly.make(remainder)

    def exact_quotient(self, divisor: "Poly") -> "Poly":
        quotient, remainder = self.divmod(divisor)
        require(remainder.is_zero(), "expected exact polynomial quotient")
        return quotient

    def monic(self) -> "Poly":
        if self.is_zero():
            return self
        return self.scalar(self.lc.inverse())


def as_poly(value) -> Poly:
    if isinstance(value, Poly):
        return value
    return Poly.make([value])


ZERO_P = Poly.make([])
ONE_P = Poly.make([1])
X_P = Poly.make([0, 1])


def gcd(left: Poly, right: Poly) -> Poly:
    a, b = left, right
    while not b.is_zero():
        _, remainder = a.divmod(b)
        a, b = b, remainder
    return a.monic()


def serialize(poly: Poly) -> bytes:
    rows = []
    for coefficient in poly.c:
        rows.append(
            ",".join(f"{entry.numerator}/{entry.denominator}" for entry in coefficient.c)
        )
    return (";".join(rows) + "\n").encode()


def configure_relation(c3: int, c2: int, c1: int, c0: int) -> None:
    # q=c3*u^3+c2*u^2+c1*u+c0.
    Alg.relation = (Q(-c0, c3), Q(-c1, c3), Q(-c2, c3))


def symbolic_resultant_identity_audit() -> None:
    """Independently derive the universal constant-lane resultant in y."""
    y, V, W, dV, b, c, k = sp.symbols("y V W dV b c k")
    L = y**2 + b * y + c * V
    R1 = 2 * L * (2 * y + b) + W
    R2 = k * V**3 + V**2 * y - b * dV * y * L
    derived = sp.expand(sp.resultant(R1, R2, y))
    expected = sp.expand(
        64 * V**9 * k**3
        - 96 * V**8 * b * k**2
        + 64 * V**8 * c * k
        + 32 * V**7 * b**2 * c * dV * k**2
        + 32 * V**7 * b**2 * k
        - 32 * V**7 * b * c
        + 48 * V**6 * W * b * dV * k**2
        - 16 * V**6 * W
        - 8 * V**6 * b**4 * dV * k**2
        - 32 * V**6 * b**3 * c * dV * k
        + 32 * V**6 * b**2 * c**2 * dV
        - 24 * V**5 * W * b**2 * dV * k
        + 16 * V**5 * W * b * c * dV
        + 8 * V**5 * b**5 * dV * k
        - 8 * V**5 * b**4 * c * dV
        + 16 * V**4 * W * b**3 * c * dV**2 * k
        - 4 * V**4 * W * b**3 * dV
        + 12 * V**3 * W**2 * b**2 * dV**2 * k
        - 4 * V**3 * W * b**5 * dV**2 * k
        + 2 * V * W**2 * b**4 * c * dV**3
        + W**3 * b**3 * dV**3
    )
    require(derived == expected, "universal symbolic resultant identity differs")
    quotient = sp.cancel(derived.subs(W, V * sp.Symbol("A")) / V**3)
    require(quotient.is_polynomial(), "saturated universal resultant is not polynomial")
    print("universal_symbolic_resultant_identity=PASS")


def modular_squarefree_audit(
    name: str,
    a: int,
    b: int,
    prime: int,
    root_u: int,
    case_index: int,
    constant_c: int = 0,
) -> tuple[int, int, tuple[int, ...]]:
    x = sp.symbols("x")
    p = prime

    def inverse(value: int) -> int:
        return pow(value % p, -1, p)

    if case_index == 1:
        v = (8 * root_u**2 + 9 * root_u + 8) * inverse(7) % p
        s0 = 5 * (root_u + 1) * inverse(7) % p
        A0 = 80 * v**2 * (root_u + 1) * inverse(343) % p
    else:
        v = (24 * root_u**2 - 16 * root_u - 16) * inverse(21) % p
        s0 = (5 * root_u - 4) * inverse(7) % p
        A0 = 9 * v**2 * (5 * root_u - 4) * inverse(343) % p
    C = -7 * A0 % p
    require(C != 0, f"bad modular C for {name}")

    def polynomial(expression) -> sp.Poly:
        return sp.Poly(expression, x, modulus=p)

    q2 = polynomial(x**2 - root_u * x + v)
    D = polynomial(x**a * (x - 1) ** b) * q2
    T = polynomial(x * (x - 1)) * q2
    E = polynomial(
        (
            a * (x - 1) * (x**2 - root_u * x + v)
            + b * x * (x**2 - root_u * x + v)
            + x * (x - 1) * (2 * x - root_u)
        )
        * inverse(7)
    )
    S = polynomial(x + s0)
    V = (S * D * T**2).mul_ground(4 * inverse(C**2) % p)
    A = (S * E * T).mul_ground(2 * inverse(C) % p)
    dV = V.diff()
    g = S * T
    c = constant_c
    K = (
        64 * V**6
        - 96 * V**5
        + 64 * V**5 * c
        + 32 * V**4 * c * dV
        + 32 * V**4
        - 32 * V**4 * c
        + 48 * V**4 * A * dV
        - 16 * V**4 * A
        - 8 * V**3 * dV
        - 32 * V**3 * c * dV
        + 32 * V**3 * c**2 * dV
        - 24 * V**3 * A * dV
        + 16 * V**3 * A * c * dV
        + 8 * V**2 * dV
        - 8 * V**2 * c * dV
        + 16 * V**2 * A * c * dV**2
        - 4 * V**2 * A * dV
        + 12 * V**2 * A**2 * dV**2
        - 4 * V * A * dV**2
        + 2 * A**2 * c * dV**3
        + A**3 * dV**3
    )
    residual = K
    removed = []
    while True:
        common = sp.gcd(residual, g)
        if common.degree() <= 0:
            break
        removed.append(common.degree())
        residual = sp.exquo(residual, common)
    derivative_gcd = sp.gcd(residual, residual.diff())
    require(residual.degree() == 52, f"bad modular residual degree for {name}")
    require(derivative_gcd.degree() == 0, f"modular residual not squarefree for {name}")
    require(sp.gcd(residual, V).degree() == 0, f"modular residual meets V for {name}")
    return residual.degree(), derivative_gcd.degree(), tuple(removed)


def constant_lane_resultant(V: Poly, A: Poly, dV: Poly, b, c, k) -> Poly:
    """Res_y(R1,R2)/V^3 for constant B=b,C0=c,E0'=k."""
    b, c, k = as_alg(b), as_alg(c), as_alg(k)
    return (
        (V**6).scalar(64 * k**3)
        + (V**5).scalar(-96 * b * k**2 + 64 * c * k)
        + (V**4 * dV).scalar(32 * b**2 * c * k**2)
        + (V**4).scalar(32 * b**2 * k - 32 * b * c)
        + (V**4 * A * dV).scalar(48 * b * k**2)
        - 16 * V**4 * A
        + (V**3 * dV).scalar(
            -8 * b**4 * k**2 - 32 * b**3 * c * k + 32 * b**2 * c**2
        )
        + (V**3 * A * dV).scalar(-24 * b**2 * k + 16 * b * c)
        + (V**2 * dV).scalar(8 * b**5 * k - 8 * b**4 * c)
        + (V**2 * A * dV**2).scalar(16 * b**3 * c * k)
        + (V**2 * A * dV).scalar(-4 * b**3)
        + (V**2 * A**2 * dV**2).scalar(12 * b**2 * k)
        + (V * A * dV**2).scalar(-4 * b**5 * k)
        + (A**2 * dV**3).scalar(2 * b**4 * c)
        + (A**3 * dV**3).scalar(b**3)
    )


def saturate_by(poly: Poly, boundary: Poly) -> tuple[Poly, tuple[int, ...]]:
    residual = poly
    removed = []
    while True:
        common = gcd(residual, boundary)
        if common.degree <= 0:
            break
        removed.append(common.degree)
        residual = residual.exact_quotient(common)
    return residual, tuple(removed)


def audit_case(
    name: str,
    a: int,
    b: int,
    q_coeffs,
    v: Alg,
    s0: Alg,
    A0: Alg,
    extra_x: int,
    extra_x_minus_one: int,
    modular_prime: int,
    modular_root: int,
    case_index: int,
) -> None:
    print(f"begin_passport={name}", flush=True)
    configure_relation(*q_coeffs)
    # Recreate inputs after relation configuration where multiplication occurs.
    v = eval_alg(v.c)
    s0 = eval_alg(s0.c)
    A0 = eval_alg(A0.c)
    C = -7 * A0
    Q2 = X_P**2 - X_P.scalar(U_A) + v
    D = X_P**a * (X_P - 1) ** b * Q2
    T = X_P * (X_P - 1) * Q2
    E = (
        (X_P - 1) * Q2 * a
        + X_P * Q2 * b
        + X_P * (X_P - 1) * (2 * X_P - U_A)
    ).scalar(Q(1, 7))
    S = X_P + s0
    V = (S * D * T**2).scalar(4 / C**2)
    A = (S * E * T).scalar(2 / C)
    dV = V.derivative()
    g = S * T
    print(f"built_passport={name} degree_V={V.degree} degree_A={A.degree}", flush=True)

    # Universal saturated resultant after y=Vz, divided by its forced V^3.
    K = constant_lane_resultant(V, A, dV, 1, 0, 1)
    require(K.degree == 96, "unexpected critical resultant degree")
    print(f"built_resultant={name} degree={K.degree}", flush=True)
    residual, removed = saturate_by(K, g)
    running_degree = K.degree
    for removed_degree in removed:
        running_degree -= removed_degree
        print(f"removed_g_factor={name}:{removed_degree}->{running_degree}", flush=True)
    forced_boundary = S**3 * T**8 * X_P**extra_x * (X_P - 1) ** extra_x_minus_one
    require(
        forced_boundary * residual == K,
        "explicit boundary-factor formula differs",
    )
    # V and g have the same radical, so the terminal gcd(residual,g)=1
    # already proves that the residual lies wholly in V!=0.  A separate
    # finite-field squarefreeness audit is cheaper than a coefficient-heavy
    # characteristic-zero gcd with the degree-51/.. derivative.
    modular_degree, repeated_degree, modular_removed = modular_squarefree_audit(
        name, a, b, modular_prime, modular_root, case_index
    )
    require(tuple(removed) == modular_removed, "modular saturation pattern differs")
    print(f"finished_saturation={name}", flush=True)

    print(f"passport={name}")
    print(f"degree_V={V.degree} degree_A={A.degree} degree_g={g.degree}")
    print(f"raw_saturated_resultant_degree={K.degree}")
    print(f"successive_g_gcd_degrees={tuple(removed)}")
    print(
        f"forced_boundary_factor=S^3*T^8*x^{extra_x}*(x-1)^{extra_x_minus_one}"
    )
    print(f"away_from_g_resultant_degree={residual.degree}")
    print(f"away_from_g_resultant_derivative_gcd_degree={repeated_degree}")
    print(
        f"squarefree_good_reduction=(p={modular_prime},u={modular_root},degree={modular_degree})"
    )
    print(f"away_from_g_resultant_sha256={hashlib.sha256(serialize(residual)).hexdigest()}")
    print(f"away_from_g_resultant_lc={residual.lc.c}")
    print("gcd_away_resultant_with_V=1")

    # Cheapest coefficient perturbation C0=1, with B=1 and E0=x retained.
    # Its local clutch is Delta=1-A', and its global resultant is the c=1
    # specialization of the same universal 21-term formula.
    delta_c1 = ONE_P - A.derivative()
    require(gcd(delta_c1, g).degree == 0, "C0=1 fails the local clutch")
    K_c1 = constant_lane_resultant(V, A, dV, 1, 1, 1)
    residual_c1, removed_c1 = saturate_by(K_c1, g)
    require(forced_boundary * residual_c1 == K_c1, "C0=1 boundary factor differs")
    require(residual_c1.degree == 52, "C0=1 residual degree differs")
    mod_degree_c1, mod_gcd_c1, mod_removed_c1 = modular_squarefree_audit(
        name, a, b, modular_prime, modular_root, case_index, constant_c=1
    )
    require(removed_c1 == mod_removed_c1, "C0=1 modular saturation differs")
    require(mod_degree_c1 == 52 and mod_gcd_c1 == 0, "C0=1 modular gate failed")
    print("nearby_control=(B=1,C0=1,E0=x)")
    print("nearby_local_Delta=1-A_src_prime_is_unit_mod_g")
    print("nearby_away_resultant_degree=52")
    print(
        f"nearby_away_resultant_sha256={hashlib.sha256(serialize(residual_c1)).hexdigest()}"
    )
    print("nearby_away_resultant_squarefree_good_reduction=PASS")


def eval_alg(coefficients) -> Alg:
    # Re-evaluate a raw degree <=2 tuple after relation configuration.
    return Alg(tuple(Q(value) for value in coefficients))


def main() -> None:
    symbolic_resultant_identity_audit()
    # Values requiring multiplication (A0) are constructed after configuring
    # the relation, so each block is explicit rather than sharing stale state.
    configure_relation(100, 244, 237, 44)
    v1 = (8 * U_A**2 + 9 * U_A + 8) / 7
    s1 = 5 * (U_A + 1) / 7
    A01 = 80 * v1**2 * (U_A + 1) / 343
    audit_case(
        "(4,1,1,1)",
        4,
        1,
        (100, 244, 237, 44),
        v1,
        s1,
        A01,
        9,
        0,
        113,
        85,
        1,
    )

    configure_relation(75, -89, -31, 61)
    v2 = (24 * U_A**2 - 16 * U_A - 16) / 21
    s2 = (5 * U_A - 4) / 7
    A02 = 9 * v2**2 * (5 * U_A - 4) / 343
    audit_case(
        "(3,2,1,1)",
        3,
        2,
        (75, -89, -31, 61),
        v2,
        s2,
        A02,
        6,
        3,
        101,
        64,
        2,
    )
    print("FAILED_CHECKS=NONE")


if __name__ == "__main__":
    main()
