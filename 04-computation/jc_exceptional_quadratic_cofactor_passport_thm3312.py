#!/usr/bin/env python3
"""Independent C2 trace/norm replay over the THM-3309 exceptional quadratic.

This is a scratch hostile continuation of the fixed C=c+x, d=k=1 slice.
It never chooses a root of F0 and keeps the degree-36 linear pivot distinct
from the degree-32 quadratic pivot P2.
"""

from __future__ import annotations

import importlib.util
from hashlib import sha256
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "04-computation" / "jc_affine_c_exceptional_quadratic_blowup_scout_20260803.py"
SPEC = importlib.util.spec_from_file_location("exceptional_prior", SOURCE)
prior = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(prior)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def digest(element):
    text = "|".join(
        f"{monomial}:{coefficient}"
        for monomial, coefficient in sorted(element.to_dict().items())
    )
    return sha256(text.encode("utf-8")).hexdigest()


def packet_digest(elements):
    return sha256("\n".join(digest(item) for item in elements).encode("ascii")).hexdigest()


def reconstruct(case):
    """Return the exact A=K[x]/(linear_a) objects used by the C2 passport."""

    name = case["name"]
    field = case["field"]
    ring = case["ring"]
    X = case["X"]
    V = case["V"]
    Aresponse = case["A"]
    boundary = case["boundary"]
    Vp = V.diff(X)
    Vpp = Vp.diff(X)

    def specialize_linear(parameter):
        # Direct closed PRS row from the independently proved universal
        # formulas.  The tracked predecessor obtains the same row by first
        # building the generic SymPy PRS; avoiding that generic expansion
        # makes this next-operation audit substantially cheaper.
        C = X + field.convert(parameter)
        P = 2 * Vp + 8 * V**2
        Q = 2 * Vp + 12 * V**2
        R = V * (4 * V**2 + 8 * C * V**2 + Vp * (2 * C + Aresponse))
        q = 2 + 4 * C * V
        r = V * (2 * C + Aresponse)
        D = 6 * P - 4 * Q
        ell1 = P * (P * q - 4 * R) - D * Q
        ell0 = P**2 * r - D * R
        return ell1.exquo(boundary), ell0.exquo(boundary)

    raw_a0, raw_b0 = specialize_linear(0)
    raw_a1, raw_b1 = specialize_linear(1)
    require(raw_a0 == raw_a1, f"{name}: linear coefficient depends on c")
    common_unit = raw_a0.LC
    linear_a = raw_a0 / common_unit
    b0 = raw_b0 / common_unit
    b1 = (raw_b1 - raw_b0) / common_unit
    require(linear_a.degree() == 36, f"{name}: degree-36 linear pivot")
    require(prior.factor_profile(linear_a) == ((36, 1),), f"{name}: A is a field")
    require(linear_a.gcd(b1).degree() == 0, f"{name}: transverse b1")

    inv_b1 = prior.inverse_mod(b1, linear_a)
    cstar = (-b0 * inv_b1) % linear_a
    Cstar = X + cstar

    P = (2 * Vp + 8 * V**2) % linear_a
    Q = (2 * Vp + 12 * V**2) % linear_a
    R = (
        V
        * (4 * V**2 + 8 * Cstar * V**2 + Vp * (2 * Cstar + Aresponse))
    ) % linear_a
    require((2 * Vp + 8 * V**2).degree() == 32, f"{name}: degree-32 P2")
    require(
        all(linear_a.gcd(item).degree() == 0 for item in (P, Q, R, V, Vp)),
        f"{name}: exceptional coefficients/V/Vprime are units",
    )
    delta = (Q**2 - 4 * P * R) % linear_a
    require(delta != ring.zero, f"{name}: separable quadratic")

    a_prime = linear_a.diff(X) % linear_a
    b_x = (b0.diff(X) + cstar * b1.diff(X)) % linear_a
    inv_a_prime = prior.inverse_mod(a_prime, linear_a)
    inv_P = prior.inverse_mod(P, linear_a)
    inv_delta = prior.inverse_mod(delta, linear_a)

    P_x = (2 * Vpp + 16 * V * Vp) % linear_a
    Q_x = (2 * Vpp + 24 * V * Vp) % linear_a
    bracket = 4 * V**2 + 8 * Cstar * V**2 + Vp * (2 * Cstar + Aresponse)
    R_x = (
        Vp * bracket
        + V
        * (
            8 * V * Vp
            + 8 * V**2
            + 16 * Cstar * V * Vp
            + Vpp * (2 * Cstar + Aresponse)
            + Vp * (2 + Aresponse.diff(X))
        )
    ) % linear_a
    R_c = (V * (8 * V**2 + 2 * Vp)) % linear_a
    F1_2 = (P_x * inv_a_prime) % linear_a
    F1_1 = (-Q_x * inv_a_prime + R_c * inv_b1) % linear_a
    F1_0 = (R_x * inv_a_prime - R_c * b_x * inv_a_prime * inv_b1) % linear_a
    leading_ratio = (F1_2 * inv_P) % linear_a
    rem1 = (F1_1 + leading_ratio * Q) % linear_a
    rem0 = (F1_0 - leading_ratio * R) % linear_a
    velocity_1 = (-(rem1 * Q + 2 * P * rem0) * inv_delta) % linear_a
    velocity_0 = ((2 * rem1 * R + rem0 * Q) * inv_delta) % linear_a
    require(velocity_1 != ring.zero, f"{name}: velocity separates branches")

    return {
        "name": name,
        "field": field,
        "ring": ring,
        "linear_a": linear_a,
        "P": P,
        "Q": Q,
        "R": R,
        "delta": delta,
        "V": V % linear_a,
        "Vp": Vp % linear_a,
        "Aresponse": Aresponse % linear_a,
        "Cstar": Cstar % linear_a,
        "velocity_0": velocity_0,
        "velocity_1": velocity_1,
    }


class QuadraticAlgebra:
    """B=A[t]/(P*t^2-Q*t+R), represented in the basis (1,t)."""

    def __init__(self, data):
        self.data = data
        self.ring = data["ring"]
        self.modulus = data["linear_a"]
        self.P = data["P"]
        self.Q = data["Q"]
        self.R = data["R"]
        self.delta = data["delta"]
        self.invP = prior.inverse_mod(self.P, self.modulus)
        self.sum_roots = self.red(self.Q * self.invP)
        self.product_roots = self.red(self.R * self.invP)

    def red(self, value):
        return value % self.modulus

    def add(self, left, right):
        return (self.red(left[0] + right[0]), self.red(left[1] + right[1]))

    def neg(self, value):
        return (self.red(-value[0]), self.red(-value[1]))

    def sub(self, left, right):
        return self.add(left, self.neg(right))

    def scale(self, scalar, value):
        return (self.red(scalar * value[0]), self.red(scalar * value[1]))

    def mul(self, left, right):
        a, b = left
        c, d = right
        return (
            self.red(a * c - b * d * self.product_roots),
            self.red(a * d + b * c + b * d * self.sum_roots),
        )

    def conjugate(self, value):
        a, b = value
        return (self.red(a + b * self.sum_roots), self.red(-b))

    def trace(self, value):
        return self.red(2 * value[0] + value[1] * self.sum_roots)

    def norm(self, value):
        product = self.mul(value, self.conjugate(value))
        require(product[1] == self.ring.zero, "norm is not invariant")
        return product[0]

    def inverse(self, value):
        norm = self.norm(value)
        require(norm != self.ring.zero, "division by zero in B")
        inv_norm = prior.inverse_mod(norm, self.modulus)
        return self.scale(inv_norm, self.conjugate(value))

    def divide(self, numerator, denominator):
        return self.mul(numerator, self.inverse(denominator))

    def square(self, value):
        return self.mul(value, value)

    def diff_square(self, value):
        difference = self.sub(value, self.conjugate(value))
        square = self.square(difference)
        require(square[1] == self.ring.zero, "conjugate difference square")
        expected = self.red(
            value[1] ** 2 * self.delta * prior.inverse_mod(self.P**2, self.modulus)
        )
        require(square[0] == expected, "quadratic discriminant scaling")
        return square[0]


def report_A(label, value, data):
    require(value != data["ring"].zero, f"{data['name']}: {label} vanished")
    print(
        f"case={data['name']};A_invariant={label};degree={value.degree()};"
        f"digest={digest(value)};nonzero_unit=1;in_accessory_field={int(value.degree() == 0)}"
    )


def audit_passport(data):
    name = data["name"]
    ring = data["ring"]
    P, Q, R = data["P"], data["Q"], data["R"]
    V, Vp = data["V"], data["Vp"]
    algebra = QuadraticAlgebra(data)
    zero, one = ring.zero, ring.one
    t = (zero, one)
    y = (zero, -one)

    # F0(t)=0 and the two original critical equations vanish without choosing
    # either geometric root.
    F0 = algebra.add(
        algebra.add(algebra.scale(P, algebra.square(t)), algebra.scale(-Q, t)),
        (R, zero),
    )
    require(F0 == (zero, zero), f"{name}: F0 relation")
    L = algebra.add(algebra.add(algebra.square(y), y), (data["Cstar"] * V, zero))
    R1 = algebra.add(
        algebra.scale(2, algebra.mul(L, algebra.add(algebra.scale(2, y), (one, zero)))),
        (V * data["Aresponse"], zero),
    )
    R2 = algebra.add(
        algebra.add((V**3, zero), algebra.scale(V**2, y)),
        algebra.mul(L, algebra.add(algebra.scale(-Vp, y), (2 * V**2, zero))),
    )
    require(R1 == (zero, zero) and R2 == (zero, zero), f"{name}: critical pair")

    D = algebra.red(6 * P - 4 * Q)
    require(D == algebra.red(4 * Vp), f"{name}: D=4Vprime")
    quotient = (D, algebra.red(-4 * P))       # 4*P*y+D at y=-t
    W = algebra.scale(-4, quotient)
    U = algebra.add((P**2, zero), algebra.scale(-Vp, quotient))
    identity = algebra.sub(U, algebra.scale(Vp / 4, W))
    require(identity == (algebra.red(P**2), zero), f"{name}: cofactor identity")

    velocity = (data["velocity_0"], data["velocity_1"])
    ratio = algebra.divide(U, W)

    objects = {
        "selected_y": y,
        "normal_velocity": velocity,
        "elimination_U": U,
        "elimination_W": W,
        "projective_ratio_U_over_W": ratio,
    }
    for label, value in objects.items():
        trace = algebra.trace(value)
        norm = algebra.norm(value)
        difference_square = algebra.diff_square(value)
        report_A(label + "_trace", trace, data)
        report_A(label + "_norm", norm, data)
        report_A(label + "_conjugate_difference_square", difference_square, data)
        print(
            f"case={name};B_object={label};pair_digest={packet_digest(value)};"
            "relative_trace_norm_descend_to_A=1;branches_distinct=1"
        )

    # The two projective cofactor ratios are genuinely different.  Their
    # alternating cross determinant has a normalization-free square.
    cross = algebra.sub(
        algebra.mul(U, algebra.conjugate(W)),
        algebra.mul(algebra.conjugate(U), W),
    )
    cross_square = algebra.square(cross)
    require(cross_square[1] == zero, f"{name}: cross square invariant")
    expected_cross_square = algebra.red(256 * P**4 * data["delta"])
    require(cross_square[0] == expected_cross_square, f"{name}: cofactor cross formula")
    report_A("projective_cofactor_cross_square_256P4delta", cross_square[0], data)

    # Exact closed relative passports for the critical y-root and quotient pair.
    require(algebra.trace(y) == algebra.red(-Q * algebra.invP), f"{name}: trace y")
    require(algebra.norm(y) == algebra.red(R * algebra.invP), f"{name}: norm y")
    require(algebra.trace(quotient) == algebra.red(-48 * V**2), f"{name}: trace quotient")
    require(algebra.trace(W) == algebra.red(192 * V**2), f"{name}: trace W")

    # A Keller/mate Bezout row cannot even enter this fibre: both gradient
    # components vanish, while Tr_{B/A}(1)=2 and Nm_{B/A}(1)=1.
    require(algebra.trace((one, zero)) == 2 * one, f"{name}: trace one")
    require(algebra.norm((one, zero)) == one, f"{name}: norm one")
    print(
        f"case={name};mate_gate=BLOCKED_BEFORE_mu;R1=R2=0_in_B;"
        "no_gradient_Bezout_row=1;trace_obstruction=0_ne_2;norm_obstruction=0_ne_1"
    )
    print(
        f"case={name};passport_summary="
        "critical_y_root_and_velocity_and_elimination_cofactor_ratio_are_C2_conjugate;"
        "unordered_passport_trace_norm_and_difference_square_descends_to_A;"
        "projective_ratio_generates_B_over_A;no_branch_label_descends_to_A"
    )


def main():
    print("THM-3309 quadratic branch trace/norm independent replay")
    print("linear_row_route=direct_closed_PRS_formula;no_generic_root_choice")
    for name in ("4111", "3211"):
        audit_passport(reconstruct(prior.response_case(name)))
    print("scope=fixed_C=c+x_d=k=1_two_accessory_fields;no_root_choice;no_Keller_mate;no_JC2")
    print("ALL EXACT PASSPORT CHECKS PASSED")


if __name__ == "__main__":
    main()
