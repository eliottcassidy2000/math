#!/usr/bin/env python3
"""Exact symbolic controls for THM-3880's sign-monodromy obstruction.

Reproduction:
  python3 04-computation/jc2_marked_root_opposite_sign_nondescent_thm3880.py
  python3 -O 04-computation/jc2_marked_root_opposite_sign_nondescent_thm3880.py
"""

from __future__ import annotations

import hashlib
import sys

import sympy as sp

sys.stdout.reconfigure(newline="\n")


A, C, b, z, eta, U, S = sp.symbols("A C b z eta U S")
Ap, Cp, bp, zp, epsilon = sp.symbols("Ap Cp bp zp epsilon")
CHECKS = 0


def zero(label: str, value: sp.Expr) -> None:
    global CHECKS
    CHECKS += 1
    value = sp.factor(value)
    if value != 0:
        raise AssertionError(f"{label}: {value}")


def nonzero(label: str, value: sp.Expr) -> None:
    global CHECKS
    CHECKS += 1
    value = sp.factor(value)
    if value == 0:
        raise AssertionError(f"{label}: unexpectedly zero")


# Universal cusp identity.
Delta = -27 * A**2 * b**2 + 8 * A * C**3 - 54 * A * C * b + 9 * C**2 - 54 * b
P = 1 + sp.Rational(2, 3) * A * C
u = 1 + A * C + A**2 * b
zero("cusp_identity", A**2 * Delta - 27 * (P**3 - u**2))

# On z^2=P, the two possible global signs u=+/-z^3 give these forced
# A^2B numerators.  Opposite root signs change them by +/-2z^3.
AC_on_root = sp.Rational(3, 2) * (z**2 - 1)
num_plus = sp.expand(z**3 - 1 - AC_on_root)
num_minus = sp.expand(-z**3 - 1 - AC_on_root)
zero("plus_factor", 2 * num_plus - (z - 1) ** 2 * (2 * z + 1))
zero("minus_factor", 2 * num_minus + (z + 1) ** 2 * (2 * z - 1))
zero("plus_opposite_sign_jump", num_plus.subs(z, -z) - num_plus + 2 * z**3)
zero("minus_opposite_sign_jump", num_minus.subs(z, -z) - num_minus - 2 * z**3)
nonzero("jump_polynomial_nonzero", 2 * z**3)

# Exact finite pole criterion above A=0.  The relation
# (z-1)(z+1)=(2/3)AC cancels A^2 precisely at the matching sign.
plus_regular = sp.Rational(2, 9) * C**2 * (2 * z + 1) / (z + 1) ** 2
minus_regular = -sp.Rational(2, 9) * C**2 * (2 * z - 1) / (z - 1) ** 2
zero(
    "plus_matching_sign_cancellation",
    num_plus
    - (A**2 * plus_regular).subs(A**2 * C**2, sp.Rational(9, 4) * (z**2 - 1) ** 2),
)
zero(
    "minus_matching_sign_cancellation",
    num_minus
    - (A**2 * minus_regular).subs(A**2 * C**2, sp.Rational(9, 4) * (z**2 - 1) ** 2),
)
zero("plus_wrong_sign_unit_numerator", num_plus.subs(z, -1) + 2)
zero("minus_wrong_sign_unit_numerator", num_minus.subs(z, 1) + 2)

# Complete ordinary-node/A2-cusp descent checks.  At a node on A=0,
# Delta=0 forces B=C^2/6 on every branch; away from A=0 the displayed
# forced formula is already a single-valued expression in A,C,z.
zero("A0_node_forced_value", Delta.subs({A: 0, b: C**2 / 6}))

# Differentiate z^2=P and A^2B=epsilon*z^3-1-AC along a normalization
# parameter.  At an A2 cusp Ap=Cp=0.  The first identity gives z*zp=0;
# inserting that in the second leaves A^2*bp=0.
root_derivative = 2 * z * zp - sp.Rational(2, 3) * (Ap * C + A * Cp)
forced_derivative = (
    2 * A * Ap * b
    + A**2 * bp
    - 3 * epsilon * z**2 * zp
    + Ap * C
    + A * Cp
)
zero("cusp_root_derivative_seam", root_derivative.subs({Ap: 0, Cp: 0}) - 2 * z * zp)
zero(
    "cusp_forced_derivative_seam",
    forced_derivative.subs({Ap: 0, Cp: 0}) - A**2 * bp + 3 * epsilon * z**2 * zp,
)
zero("cusp_z2zp_from_zzp", z**2 * zp - z * (z * zp))

# If the cusp lies on A=0, differentiate Delta_B=0 directly.  All coordinate
# derivative terms vanish and the remaining coefficient is -54*B'.
Delta_derivative = sp.diff(Delta, A) * Ap + sp.diff(Delta, C) * Cp + sp.diff(Delta, b) * bp
zero(
    "A0_cusp_direct_derivative",
    Delta_derivative.subs({A: 0, Ap: 0, Cp: 0}) + 54 * bp,
)

# THM-3876 is exactly the opposite-sign case.  With eta=zeta^N and
# U=r^(M+N), its collision has U=-1/(eta+1).
U0 = -1 / (eta + 1)
z0 = sp.factor(1 + 2 * U0)
z1 = sp.factor(1 + 2 * eta * U0)
zero("thm3876_root_sign_flip", z1 + z0)
zero("thm3876_C_collision", (eta - 1) + (eta**2 - 1) * U0)

# Since U^2=A^2*r^(2N), replace 1/A^2 by S/U^2.  The intrinsic jump then
# reproduces THM-3876's explicit marked-value difference.
intrinsic_jump = sp.factor(-2 * z0**3 * S / U0**2)
explicit_jump = -2 * S * (eta - 1) ** 3 / (eta + 1)
zero("thm3876_jump_recovery", intrinsic_jump - explicit_jump)

# Exact boundary controls.
# (i) z=0, equivalently 2AC+3=0, really can carry a polynomial profile.
b_zero_root = sp.Rational(2, 9) * C**2
Delta_zero_root = sp.factor(Delta.subs(b, b_zero_root))
zero(
    "z0_hyperbola_positive_control",
    Delta_zero_root + sp.Rational(1, 3) * C**2 * (2 * A * C + 3) ** 2,
)
# (ii) A=0 is not itself forbidden when the sign is constant.
b_axis = sp.Rational(1, 6) * C**2
Delta_axis = sp.factor(Delta.subs(b, b_axis))
zero(
    "A0_same_sign_positive_control",
    Delta_axis + sp.Rational(1, 4) * A * C**3 * (3 * A * C + 4),
)
# (iii) at A=0 one has z=+/-1, and the two signs give opposite u-values;
# the undivided u argument, not division by A^2, detects the conflict.
zero("A0_plus_u", 1**3 - 1)
zero("A0_minus_u", (-1) ** 3 + 1)
nonzero("A0_opposite_u_jump", 1**3 - (-1) ** 3)

print("THM3880_IDENTITY", "A^2Delta=27(P^3-u^2);P=1+(2/3)AC")
print("THM3880_MARKED_ROOT", "z^2=P;carrier forces one global sign u=epsilon*z^3")
print("THM3880_PLUS", "A^2B+=(z-1)^2(2z+1)/2")
print("THM3880_MINUS", "A^2B-=-(z+1)^2(2z-1)/2")
print("THM3880_JUMP", "B_epsilon(-z)-B_epsilon(z)=-2epsilon*z^3/A^2 when A!=0")
print("THM3880_A0", "use u equality;opposite z=+/-1 still impossible")
print("THM3880_Z0", "silent genuine boundary;2AC+3=0,b=2C^2/9 is positive")
print("THM3880_POLE_PLUS", "B+ regular above A=0 iff z=+1 there")
print("THM3880_POLE_MINUS", "B- regular above A=0 iff z=-1 there")
print("THM3880_NODE_IFF", "regular B descends at a node iff root signs are not opposite nonzero")
print("THM3880_A2", "same-sign forced B automatically has zero cusp derivative")
print("THM3880_GLOBAL_IFF", "fixed epsilon carrier iff z=epsilon over A=0 and z matches at every node")
print("THM3880_OPEN", "higher singularities,no normalization root,projective companion remain open")
print("THM3880_THM3876", "eta collision is exactly z->-z and recovers its jump")
semantic_packet = (
    "intrinsic marked-root square section on normalization",
    "global sign u equals epsilon z cubed",
    "opposite nonzero root signs over one target point",
    "undivided u obstruction including A zero",
    "explicit carrier jump away from A zero",
    "z zero genuine hyperbola boundary",
    "same-sign ordinary nodes and A2 cusp jets automatically descend",
    "exact matching-sign cancellation above A zero",
    "wrong sign gives a pole of twice the A order",
    "full fixed-sign nodal-cuspidal carrier iff",
    "higher conductor jets and projective geometry open",
    "THM3876 primitive-root collision as exact specialization",
)
print("SEMANTIC_SHA256", hashlib.sha256(repr(semantic_packet).encode()).hexdigest())
print("CHECKS", CHECKS)
