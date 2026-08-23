#!/usr/bin/env python3
"""Exact companion for THM-3780's all-affine-plane filling obstruction."""

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


def jac(first: sp.Expr, second: sp.Expr, left: sp.Symbol, right: sp.Symbol) -> sp.Expr:
    return sp.expand(
        sp.diff(first, left) * sp.diff(second, right)
        - sp.diff(first, right) * sp.diff(second, left)
    )


x, y = sp.symbols("x y", nonzero=True)
t, p, g, L = sp.symbols("t p g L")
m = sp.symbols("m", integer=True, positive=True)

A_source = 1 + x * y
t_source = x**2 * A_source
c_m = (m + 1) / m
P_source = 1 / x + c_m * t_source**m
g_source = sp.factor(P_source - c_m * t_source**m)
U_source = sp.expand(t_source * P_source - t_source ** (m + 1) / m)
B_source = 1 + x * t_source**m
U_tower = sp.expand(x * A_source * B_source)

gate(sp.factor(g_source - 1 / x) == 0, "opposite-boundary unit")
gate(sp.expand(U_source - U_tower) == 0, "source tower identity")
gate(sp.factor(jac(t_source, P_source, x, y) - x) == 0,
     "birational chart Jacobian")
gate(sp.factor(jac(U_source, P_source, x, y) - 1) == 0,
     "rational Keller response")
gate(
    sp.factor(
        t_source ** (m + 1) - m * P_source * t_source + m * U_source
    )
    == 0,
    "monic integral equation",
)

# The inverse formulas prove equality of the two localized cylinder rings.
x_inverse = 1 / g
y_inverse = t * g**3 - g
A_inverse = sp.expand(1 + x_inverse * y_inverse)
t_roundtrip = sp.factor(x_inverse**2 * A_inverse)
P_roundtrip = sp.factor(g + c_m * t**m)
gate(sp.factor(A_inverse - t * g**2) == 0, "inverse A formula")
gate(sp.factor(t_roundtrip - t) == 0, "inverse t formula")
gate(sp.factor(P_roundtrip - g - c_m * t**m) == 0, "inverse P formula")
gate(
    sp.factor((t_source * g_source**3 - g_source) - y) == 0,
    "source y recovery",
)

# The alternative affine-plane completion is the normalization plane
# k[t,p].  It makes both target functions polynomial but exposes their
# ramification divisor g=0.
U_normal = t * p - t ** (m + 1) / m
g_normal = p - c_m * t**m
gate(sp.factor(sp.diff(U_normal, t) - g_normal) == 0,
     "normalization Jacobian coefficient")
gate(sp.diff(U_normal, p) == t, "normalization transverse coefficient")
gate(sp.factor(jac(U_normal, p, t, p) - g_normal) == 0,
     "normalization ramification divisor")
gate(
    sp.factor(t ** (m + 1) - m * p * t + m * U_normal) == 0,
    "normalization monic relation",
)

# The critical graph normalizes the discriminant cusp.
critical_substitution = {p: c_m * t**m}
branch_U = sp.factor(U_normal.subs(critical_substitution))
gate(sp.factor(branch_U - t ** (m + 1)) == 0, "critical divisor target U")
branch_binomial = (
    (m + 1) ** (m + 1) * L**m
    - m ** (m + 1) * p ** (m + 1)
)
gate(
    sp.factor(
        sp.expand_power_base(
            sp.expand_power_exp(
                branch_binomial.subs({L: t ** (m + 1), p: c_m * t**m})
            ),
            force=True,
        )
    )
    == 0,
    "branch cusp normalization",
)

# The differential identity used in the abstract normal-model proof is the
# literal coefficient identity dU wedge dp = g dt wedge dp.
gate(sp.factor(sp.diff(U_normal, t) - g_normal) == 0,
     "cotangent forced-unit identity")

direct_controls = []
branch_controls = []
packet_rows = []
for degree_parameter in range(1, 9):
    coefficient = sp.Rational(degree_parameter + 1, degree_parameter)
    Axy = 1 + x * y
    tau = sp.expand(x**2 * Axy)
    Bxy = sp.expand(1 + x ** (2 * degree_parameter + 1) * Axy**degree_parameter)
    Uxy = sp.expand(x * Axy * Bxy)
    Pxy = sp.cancel(1 / x + coefficient * tau**degree_parameter)
    Gxy = sp.cancel(Pxy - coefficient * tau**degree_parameter)

    gate(Gxy == 1 / x, f"source unit m={degree_parameter}")
    gate(sp.cancel(jac(Uxy, Pxy, x, y) - 1) == 0,
         f"source response m={degree_parameter}")
    gate(
        sp.cancel(
            tau ** (degree_parameter + 1)
            - degree_parameter * Pxy * tau
            + degree_parameter * Uxy
        )
        == 0,
        f"source integral equation m={degree_parameter}",
    )
    gate(sp.cancel(1 / Gxy - x) == 0,
         f"opposite completion inverse m={degree_parameter}")

    U_alt = t * p - sp.Rational(1, degree_parameter) * t ** (degree_parameter + 1)
    g_alt = p - coefficient * t**degree_parameter
    gate(jac(U_alt, p, t, p) == g_alt,
         f"alternative completion Jacobian m={degree_parameter}")
    gate(
        sp.expand(U_alt.subs(p, coefficient * t**degree_parameter)
                  - t ** (degree_parameter + 1))
        == 0,
        f"critical graph image m={degree_parameter}",
    )
    discriminant_factor = (
        (degree_parameter + 1) ** (degree_parameter + 1) * L**degree_parameter
        - degree_parameter ** (degree_parameter + 1) * p ** (degree_parameter + 1)
    )
    gate(
        sp.factor(
            discriminant_factor.subs(
                {
                    L: t ** (degree_parameter + 1),
                    p: coefficient * t**degree_parameter,
                }
            )
        )
        == 0,
        f"branch parametrization m={degree_parameter}",
    )

    direct_controls.append(str(degree_parameter))
    branch_controls.append(
        f"m{degree_parameter}:g=p-{sp.sstr(coefficient)}*t^{degree_parameter}"
    )
    packet_rows.extend(
        (
            f"m={degree_parameter};U={sp.sstr(U_alt)}",
            f"m={degree_parameter};g={sp.sstr(g_alt)}",
            f"m={degree_parameter};source_t={sp.sstr(tau)}",
            f"m={degree_parameter};source_P={sp.sstr(Pxy)}",
        )
    )

semantic_rows = (
    "tower:t=x^2(1+xy);P=1/x+(m+1)t^m/m;U=tP-t^(m+1)/m",
    "unit:g=P-(m+1)t^m/m=1/x;inverse:x=1/g;y=tg^3-g",
    "open:k[x,x^-1,y]=k[t,g,g^-1]=k[t,P,g^-1]",
    "normalization:k[t,P];monic:t^(m+1)-mPt+mU;ramification:g=0",
    "differential:dU_wedge_dP=g*dt_wedge_dP",
    "normal_model:integrality_puts_t_in_R;etale_forces_g_unit",
    "affine_plane:no_nonconstant_units;therefore_no_same-pair_Keller_filling",
    "duality:x_completion_keeps_volume_loses_P;g_completion_keeps_P_loses_etale",
)
semantic = hashlib.sha256("\n".join(semantic_rows).encode()).hexdigest()
packet = hashlib.sha256("\n".join(packet_rows).encode()).hexdigest()

source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "inactive Python assert",
)

print("theorem=THM-3780-rational-keller-tower-all-affine-plane-fillings-obstruction")
print("field=algebraically_closed_characteristic_zero;m_integer_at_least_1")
print("tower=t=x^2*(1+x*y);P=1/x+(m+1)*t^m/m;U=t*P-t^(m+1)/m")
print("opposite_boundary_unit=g=P-(m+1)*t^m/m=1/x")
print("birational_inverse=x=1/g;y=t*g^3-g")
print("shared_open=k[x,x^-1,y]=k[t,g,g^-1]=k[t,P,g^-1]=Gm_times_A1")
print("normalization_ring=k[t,P];minimal=t^(m+1)-m*P*t+m*U")
print("normalization_map=(t,P)->(U=t*P-t^(m+1)/m,P);Jacobian=g")
print("ramification=g=0;branch_param=P=(m+1)*t^m/m,U=t^(m+1)")
print("differential_identity=dU_wedge_dP=g*dt_wedge_dP")
print("normal_model_gate=U,P_in_R_and_etale_implies_t_in_R_and_g_in_R_units")
print("affine_plane_exit=units_are_constants_but_g_is_nonconstant")
print("completion_duality=x=0_keeps_volume_loses_P;g=0_keeps_P_loses_etale")
print("sharp_model=Spec(k[t,P,g^-1])_is_etale_and_has_nonconstant_unit_g")
print("direct_controls_m=" + ",".join(direct_controls))
print("branch_controls=" + ",".join(branch_controls))
print(f"packet_sha256={packet}")
print(f"semantic_sha256={semantic}")
print(f"CHECKS={CHECKS}")
print("NO_CLAIM=obstruction_for_changed_pair_or_planar_Jacobian_conjecture")
print("RESULT=PASS")
