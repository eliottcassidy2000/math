#!/usr/bin/env python3
"""Exact companion for THM-3784's different/codifferent trace duality."""

from __future__ import annotations

import ast
import hashlib
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


T, U, P = sp.symbols("T U P")
m = sp.symbols("m", integer=True, positive=True)

f_symbolic = T ** (m + 1) - m * P * T + m * U
g_symbolic = P - (m + 1) * T**m / m
gate(sp.factor(g_symbolic + sp.diff(f_symbolic, T) / m) == 0,
     "symbolic different identity")

degree_controls = []
norm_controls = []
trace_controls = []
determinant_controls = []
packet_rows = []

for degree_parameter in range(1, 11):
    n = degree_parameter + 1
    f = T**n - degree_parameter * P * T + degree_parameter * U
    g = P - sp.Rational(n, degree_parameter) * T**degree_parameter
    branch = (
        n**n * U**degree_parameter
        - degree_parameter**n * P**n
    )
    expected_discriminant = sp.factor(
        (-1) ** (degree_parameter * n // 2)
        * degree_parameter**degree_parameter
        * branch
    )
    discriminant = sp.factor(sp.discriminant(f, T))
    gate(sp.expand(discriminant - expected_discriminant) == 0,
         f"discriminant m={degree_parameter}")

    expected_norm = sp.factor(
        sp.Rational((-1) ** n, degree_parameter) * branch
    )
    resultant_norm = sp.factor(sp.resultant(f, g, T))
    gate(sp.expand(resultant_norm - expected_norm) == 0,
         f"resultant norm m={degree_parameter}")

    # Multiplication by t in the power basis 1,t,...,t^m.
    multiplication_t = sp.zeros(n)
    for column in range(n - 1):
        multiplication_t[column + 1, column] = 1
    multiplication_t[0, n - 1] = -degree_parameter * U
    multiplication_t[1, n - 1] = degree_parameter * P

    multiplication_g = (
        P * sp.eye(n)
        - sp.Rational(n, degree_parameter)
        * multiplication_t**degree_parameter
    )
    determinant_norm = sp.factor(multiplication_g.det())
    gate(sp.expand(determinant_norm - expected_norm) == 0,
         f"matrix norm m={degree_parameter}")
    gate(sp.expand(determinant_norm - resultant_norm) == 0,
         f"independent norm paths m={degree_parameter}")

    multiplication_x = multiplication_g.inv()
    traces = [
        sp.factor(sp.trace(multiplication_t**power * multiplication_x))
        for power in range(3 * n)
    ]
    for power in range(degree_parameter):
        gate(traces[power] == 0,
             f"trace-zero codifferent rung m={degree_parameter},k={power}")
    gate(traces[degree_parameter] == -degree_parameter,
         f"top trace coefficient m={degree_parameter}")

    coefficients = sp.symbols(f"a0:{n}")
    coefficient_functional = sp.expand(
        sum(coefficients[power] * traces[power] for power in range(n))
    )
    gate(
        sp.expand(
            coefficient_functional
            + degree_parameter * coefficients[degree_parameter]
        )
        == 0,
        f"coefficient-of-T^m trace functional m={degree_parameter}",
    )
    gate(
        traces[0] != -degree_parameter,
        f"constant-polynomial hostile to leading-coefficient wording m={degree_parameter}",
    )

    for power in range(n, 3 * n):
        recurrence = sp.factor(
            traces[power]
            - degree_parameter * P * traces[power - degree_parameter]
            + degree_parameter * U * traces[power - degree_parameter - 1]
        )
        gate(recurrence == 0,
             f"trace recurrence m={degree_parameter},k={power}")

    trace_pairing = sp.Matrix(
        n,
        n,
        lambda row, column: traces[row + column],
    )
    expected_pairing_determinant = (
        (-1) ** (n * (n + 1) // 2) * degree_parameter**n
    )
    pairing_determinant = sp.factor(trace_pairing.det())
    gate(pairing_determinant == expected_pairing_determinant,
         f"codifferent basis determinant m={degree_parameter}")

    norm_x = sp.factor(multiplication_x.det())
    expected_norm_x = sp.factor(
        sp.Rational((-1) ** n * degree_parameter, 1) / branch
    )
    gate(sp.factor(norm_x - expected_norm_x) == 0,
         f"inverse different norm m={degree_parameter}")
    gate(sp.factor(norm_x * determinant_norm - 1) == 0,
         f"reciprocal norm m={degree_parameter}")

    degree_controls.append(str(degree_parameter))
    norm_controls.append(f"m{degree_parameter}:{sp.sstr(expected_norm)}")
    trace_controls.append(
        f"m{degree_parameter}:"
        + ",".join(sp.sstr(value) for value in traces[:n])
    )
    determinant_controls.append(
        f"m{degree_parameter}:{sp.sstr(pairing_determinant)}"
    )
    packet_rows.extend(
        (
            f"m={degree_parameter};f={sp.sstr(f)}",
            f"m={degree_parameter};g={sp.sstr(g)}",
            f"m={degree_parameter};Norm_g={sp.sstr(expected_norm)}",
            f"m={degree_parameter};trace_packet={trace_controls[-1]}",
        )
    )

# Direct source checks guard the typing of the inverse different as an
# honest polynomial on the opposite affine-plane filling.
x, y = sp.symbols("x y")
A_source = 1 + x * y
source_controls = []
for degree_parameter in range(1, 9):
    tau = sp.expand(x**2 * A_source)
    B_source = sp.expand(
        1 + x ** (2 * degree_parameter + 1) * A_source**degree_parameter
    )
    U_source = sp.expand(x * A_source * B_source)
    P_source = sp.cancel(
        ((degree_parameter + 1) * B_source - 1) / (degree_parameter * x)
    )
    fprime_source = sp.expand(
        (degree_parameter + 1) * tau**degree_parameter
        - degree_parameter * P_source
    )
    g_source = sp.cancel(
        P_source
        - sp.Rational(degree_parameter + 1, degree_parameter)
        * tau**degree_parameter
    )
    gate(sp.cancel(g_source - 1 / x) == 0,
         f"source different unit m={degree_parameter}")
    gate(sp.cancel(x + degree_parameter / fprime_source) == 0,
         f"source inverse different m={degree_parameter}")
    gate(
        sp.cancel(
            tau ** (degree_parameter + 1)
            - degree_parameter * P_source * tau
            + degree_parameter * U_source
        )
        == 0,
        f"source minimal relation m={degree_parameter}",
    )
    for power in range(degree_parameter + 1):
        gate(
            sp.Poly(sp.expand(x * tau**power), x, y).as_expr()
            == sp.expand(x * tau**power),
            f"source codifferent ladder polynomial m={degree_parameter},k={power}",
        )
    source_controls.append(str(degree_parameter))

semantic_rows = (
    "cover:f_m=t^(m+1)-mPt+mU;different_unit:g=-f_m_prime(t)/m",
    "source_inverse_different:x=1/g=-m/f_m_prime(t)",
    "norm_g=(-1)^(m+1)/m*((m+1)^(m+1)U^m-m^(m+1)P^(m+1))",
    "norm_x=reciprocal_branch_factor",
    "trace:Tr(t^k*x)=0_for_0<=k<m;Tr(t^m*x)=-m",
    "coefficient_functional:Tr(x*q(t))=-m*coefficient_of_T^m_for_deg_q<=m",
    "codifferent:x*S=f_prime(t)^-1*S;different=g*S",
    "pairing_determinant=(-1)^((m+1)(m+2)/2)*m^(m+1)",
    "repair_lane:trace_zero_codifferent_rungs_open_no_counterexample_claim",
)
semantic = hashlib.sha256("\n".join(semantic_rows).encode()).hexdigest()
packet = hashlib.sha256("\n".join(packet_rows).encode()).hexdigest()

source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "inactive Python assert",
)

print("theorem=THM-3784-rational-keller-tower-different-codifferent-trace-duality")
print("field=algebraically_closed_characteristic_zero;m_integer_at_least_1")
print("minimal=f_m(t)=t^(m+1)-m*P*t+m*U")
print("different_unit=g=P-(m+1)*t^m/m=-f_m_prime(t)/m")
print("inverse_different=x=1/g=-m/f_m_prime(t)")
print("norm_g=(-1)^(m+1)*((m+1)^(m+1)*U^m-m^(m+1)*P^(m+1))/m")
print("norm_x=(-1)^(m+1)*m/((m+1)^(m+1)*U^m-m^(m+1)*P^(m+1))")
print("trace_packet=Tr(t^k*x)=0_for_0<=k<m;Tr(t^m*x)=-m")
print("coefficient_extraction=Tr(x*q(t))=-m*[T^m]q(T)_for_deg(q)<=m")
print("codifferent=x*k[U,P,t];different=g*k[U,P,t]")
print("pairing_determinant=(-1)^((m+1)*(m+2)/2)*m^(m+1)")
print("geometric_identity=axis_pole_inverse_different;branch_cusp_norm_of_different")
print("degree_controls_m=" + ",".join(degree_controls))
print("source_controls_m=" + ",".join(source_controls))
print("norm_controls=" + ",".join(norm_controls))
print("trace_controls=" + ",".join(trace_controls))
print("pairing_determinant_controls=" + ",".join(determinant_controls))
print(f"packet_sha256={packet}")
print(f"semantic_sha256={semantic}")
print(f"CHECKS={CHECKS}")
print("OPEN_REPAIR=trace_zero_codifferent_combinations_across_sheets")
print("NO_CLAIM=polynomial_Keller_pair_or_planar_Jacobian_counterexample")
print("RESULT=PASS")
