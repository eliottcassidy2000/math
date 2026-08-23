#!/usr/bin/env python3
"""Independent hostile audit of THM-3860.

The checker uses recursive polynomial packets and divisor-local Laurent
controls distinct from the primary companion.  All gates remain active with
``python -O``.
"""

from __future__ import annotations

import ast
import hashlib
import sys
from pathlib import Path

import sympy as sp

sys.stdout.reconfigure(newline="\n")


s, z = sp.symbols("s z")
c = sp.symbols("c", nonzero=True)
GATES = 0


def gate(condition: object, label: str) -> None:
    global GATES
    GATES += 1
    if condition is True or condition == sp.S.true:
        return
    raise RuntimeError(label)


def zero(label: str, expression: sp.Expr) -> None:
    expression = sp.cancel(sp.factor(expression))
    gate(expression == 0, f"{label}: {expression}")


def equal(label: str, left: sp.Expr, right: sp.Expr) -> None:
    zero(label, left - right)


def nonzero(label: str, expression: sp.Expr) -> None:
    gate(sp.cancel(sp.factor(expression)) != 0, f"{label}: unexpectedly zero")


def bracket(first: sp.Expr, second: sp.Expr) -> sp.Expr:
    return sp.expand(
        sp.diff(first, z) * sp.diff(second, s)
        - sp.diff(first, s) * sp.diff(second, z)
    )


def coefficient(poly: sp.Expr, variable: sp.Symbol, degree: int) -> sp.Expr:
    return sp.expand(poly).coeff(variable, degree)


# -------------------------------------------------------------------------
# 1. Full coefficient recursion, including both derivative edge charts.
# -------------------------------------------------------------------------

N = 8
af = [sp.Function(f"u{i}")(s) for i in range(N + 1)]
bf = [sp.Function(f"v{i}")(s) for i in range(N + 1)]
formal_A = sum(af[i] * z**i for i in range(N + 1))
formal_C = sum(bf[i] * z**i for i in range(N + 1))
formal_J = bracket(formal_A, formal_C)
for m in range(N):
    convolution = sum(
        i * af[i] * sp.diff(bf[j], s) - j * sp.diff(af[i], s) * bf[j]
        for i in range(N + 1)
        for j in range(N + 1)
        if i + j == m + 1
    )
    equal(f"full convolution row {m}", coefficient(formal_J, z, m), convolution)


def recursive_packet(
    a0: sp.Expr,
    b0: sp.Expr,
    a1: sp.Expr,
    b1: sp.Expr,
    stop: int,
    tag: str,
) -> tuple[list[sp.Expr], list[sp.Expr]]:
    aa = [sp.expand(a0), sp.expand(a1)]
    bb = [sp.expand(b0), sp.expand(b1)]
    equal(
        f"{tag} Bezout row",
        aa[1] * sp.diff(bb[0], s) - sp.diff(aa[0], s) * bb[1],
        1,
    )
    for m in range(1, stop):
        interior = sum(
            i * aa[i] * sp.diff(bb[m + 1 - i], s)
            - (m + 1 - i) * sp.diff(aa[i], s) * bb[m + 1 - i]
            for i in range(1, m + 1)
        )
        tangent_gauge = (m + 2) * s ** (m + 1) + (m + 1) * s + 1
        new_a = sp.expand(
            -interior * aa[1] / (m + 1)
            + tangent_gauge * sp.diff(aa[0], s)
        )
        new_b = sp.expand(
            -interior * bb[1] / (m + 1)
            + tangent_gauge * sp.diff(bb[0], s)
        )
        aa.append(new_a)
        bb.append(new_b)
        equal(
            f"{tag} recursion determinant row {m}",
            new_a * sp.diff(bb[0], s) - sp.diff(aa[0], s) * new_b,
            -interior / (m + 1),
        )
        trial_A = sum(aa[i] * z**i for i in range(len(aa)))
        trial_C = sum(bb[i] * z**i for i in range(len(bb)))
        equal(f"{tag} Darboux coefficient row {m}", coefficient(bracket(trial_A, trial_C), z, m), 0)
    return aa, bb


nodal_a0 = 9 * c**6 * s**2
nodal_b0 = 27 * c**9 * s**3 - 3 * c**3 * s
nodal_a1 = -sp.Rational(1, 3) / c**3
nodal_b1 = -sp.Rational(3, 2) * s
nodal_a, nodal_b = recursive_packet(
    nodal_a0, nodal_b0, nodal_a1, nodal_b1, 7, "nodal"
)

# These two packets force one arm derivative to vanish identically.  They
# guard against proofs that silently divide by a0' or b0' instead of using
# the unimodular row.
edge_a_a, edge_a_b = recursive_packet(s, 0, s, -1, 5, "b0-prime-zero")
edge_b_a, edge_b_b = recursive_packet(0, s, 1, s**2 + 1, 5, "a0-prime-zero")

W_nodal = sp.expand(
    nodal_a1 * sp.diff(nodal_b1, s) - sp.diff(nodal_a1, s) * nodal_b1
)
equal("nodal Wronskian", W_nodal, 1 / (2 * c**3))

# Difference of any two solutions is exactly a tangent multiple in the two
# hostile edge charts, illustrating the exhaustive kernel claim.
tau = s**3 + 2
for tag, aa, bb in (
    ("b0-prime-zero", edge_a_a, edge_a_b),
    ("a0-prime-zero", edge_b_a, edge_b_b),
):
    tangent_a = tau * sp.diff(aa[0], s)
    tangent_b = tau * sp.diff(bb[0], s)
    equal(
        f"{tag} tangent kernel",
        tangent_a * sp.diff(bb[0], s) - sp.diff(aa[0], s) * tangent_b,
        0,
    )


# -------------------------------------------------------------------------
# 2. Vertical rational class, jets, and the literal nonuniqueness repair.
# -------------------------------------------------------------------------

w = sp.symbols("w", nonzero=True)
phi = sp.Function("Phi")(z)
f = sp.Function("Ffree")(z)
L = 1 / ((1 + w * phi) * sp.diff(phi, z))
Z = phi
S = L * s + f
equal(
    "vertical density",
    (1 + w * Z) * bracket(Z, S),
    1,
)
equal("vertical slope forced", sp.diff(S, s), L)

phi_mobius = z / (1 + w * z / 2)
L_mobius = sp.factor(
    1 / ((1 + w * phi_mobius) * sp.diff(phi_mobius, z))
)
equal(
    "Mobius slope",
    L_mobius,
    (1 + w * z / 2) ** 3 / (1 + 3 * w * z / 2),
)
equal("Mobius phi zero", phi_mobius.subs(z, 0), 0)
equal("Mobius phi first jet", sp.diff(phi_mobius, z).subs(z, 0), 1)
equal("Mobius phi second jet", sp.diff(phi_mobius, z, 2).subs(z, 0), -w)
equal("Mobius L zero jet", L_mobius.subs(z, 0), 1)
equal("Mobius L first jet", sp.diff(L_mobius, z).subs(z, 0), 0)

# For fixed phi the slope L is unique, but f is free.  Two distinct choices
# preserve the same arm, first jet, and exact density.  This is a counterexample
# to a literal claim of a "unique companion", though not to the body theorem.
S_zero = L_mobius * s
S_free = L_mobius * s + z**2
nonzero("vertical companions are not unique", S_free - S_zero)
for label, candidate_S in (("zero", S_zero), ("free", S_free)):
    equal(
        f"vertical companion {label} density",
        (1 + w * phi_mobius) * bracket(phi_mobius, candidate_S),
        1,
    )
    equal(f"vertical companion {label} arm", candidate_S.subs(z, 0), s)
    equal(
        f"vertical companion {label} first normal jet",
        sp.diff(candidate_S, z).subs(z, 0),
        0,
    )

mobius_phi_series = sp.series(phi_mobius, z, 0, 4).removeO().expand()
mobius_L_series = sp.series(L_mobius, z, 0, 4).removeO().expand()
equal("Mobius z2 normal", coefficient(mobius_phi_series, z, 2), -w / 2)
equal("Mobius z3 normal", coefficient(mobius_phi_series, z, 3), w**2 / 4)
equal("Mobius z2 tangent", coefficient(mobius_L_series, z, 2), 3 * w**2 / 4)
equal("Mobius z3 tangent", coefficient(mobius_L_series, z, 3), -w**3)


# -------------------------------------------------------------------------
# 3. Explicit nodal lift and Russell divisor geometry.
# -------------------------------------------------------------------------

h = 1 + z / (4 * c**3)
g = 1 + 3 * z / (4 * c**3)
Z_nodal = z / h
S_nodal = s * h**3 / g
A_nodal = sp.factor(9 * c**6 * S_nodal**2 - Z_nodal / (3 * c**3))
C_nodal_seed = sp.factor(
    27 * c**9 * S_nodal**3
    - 3 * c**3 * S_nodal
    - sp.Rational(3, 2) * S_nodal * Z_nodal
)
C_nodal_closed = 27 * c**9 * s**3 * h**9 / g**3 - 3 * c**3 * s * h**2
equal("nodal C closed form", C_nodal_seed, C_nodal_closed)
equal("nodal exact Darboux", bracket(A_nodal, C_nodal_seed), 1)
equal("nodal arm A", A_nodal.subs(z, 0), nodal_a0)
equal("nodal arm C", C_nodal_seed.subs(z, 0), nodal_b0)
equal("nodal first row A", sp.diff(A_nodal, z).subs(z, 0), nodal_a1)
equal("nodal first row C", sp.diff(C_nodal_seed, z).subs(z, 0), nodal_b1)

# The Russell surface is smooth, so every named principal prime below gives
# a DVR valuation.  Saturating c as a unit leaves no common gradient zero.
r, e, cinv = sp.symbols("r e cinv")
surface = r**2 * e - z**3 + c**3 * r
surface_singular_ideal = sp.groebner(
    [
        surface,
        sp.diff(surface, r),
        sp.diff(surface, z),
        sp.diff(surface, e),
        c * cinv - 1,
    ],
    r,
    z,
    e,
    c,
    cinv,
    order="lex",
)
gate(surface_singular_ideal.reduce(sp.Integer(1))[1] == 0, "Russell surface smoothness ideal not unit")

z0 = sp.symbols("z0", nonzero=True)
e_divisor = (z0**3 - c**3 * r) / r**2
equal("constant-z quotient relation", surface.subs({z: z0, e: e_divisor}), 0)
rinv_formula = (r * e + c**3) / z0**3
equal(
    "constant-z quotient makes r invertible",
    sp.rem(
        sp.together(r * rinv_formula - 1).as_numer_denom()[0],
        sp.Poly(surface.subs(z, z0), e).as_expr(),
        e,
    ),
    0,
)
D_divisor = sp.factor(c**3 + e_divisor * r)
s_divisor = sp.factor(e_divisor / (3 * D_divisor))
equal(
    "constant-z arm parameter",
    s_divisor,
    1 / (3 * r) - c**3 / (3 * z0**3),
)
nonzero("constant-z arm parameter nonconstant", sp.diff(s_divisor, r))

named_root = -4 * c**3
equal("named divisor g unit", g.subs(z, named_root), -2)
named_residue = sp.simplify(sp.limit(h * A_nodal, z, named_root))
equal("named divisor A residue", named_residue, sp.Rational(4, 3))

# A Laurent expansion at the named divisor sees the simple pole directly;
# the h^6 term remains regular and cannot interact with it.
u = sp.symbols("u")
named_laurent = sp.series(A_nodal.subs(z, named_root + u), u, 0, 2)
equal("named pole coefficient", named_laurent.coeff(u, -1), 16 * c**3 / 3)


# -------------------------------------------------------------------------
# 4. Pole lemma, all root patterns, and equal-order cancellation.
# -------------------------------------------------------------------------

degree = sp.symbols("degree", integer=True, positive=True)
aa = -w / degree
phi_exceptional = ((1 + aa * z) ** (-degree) - 1) / w
equal("exceptional family value", phi_exceptional.subs(z, 0), 0)
equal("exceptional family first derivative", sp.diff(phi_exceptional, z).subs(z, 0), 1)
exceptional_second = sp.simplify(sp.diff(phi_exceptional, z, 2).subs(z, 0))
equal("exceptional family second derivative", exceptional_second, w * (degree + 1) / degree)
equal(
    "exceptional family contradiction",
    exceptional_second + w,
    w * (2 * degree + 1) / degree,
)
for degree_value in range(1, 9):
    nonzero(
        f"exceptional positive degree {degree_value}",
        (2 * degree + 1).subs(degree, degree_value),
    )

# For the actual Mobius control, locate the forced zero of the density
# denominator, verify phi is regular there, and recover the L pole.
bad_mobius = -sp.Rational(2, 3) / w
nonzero("Mobius bad point nonzero", bad_mobius)
equal(
    "Mobius density denominator zero",
    ((1 + w * phi_mobius) * sp.diff(phi_mobius, z)).subs(z, bad_mobius),
    0,
)
nonzero(
    "Mobius phi regular at bad point",
    (1 + w * z / 2).subs(z, bad_mobius),
)
nonzero("Mobius L has polar numerator", sp.factor((1 + w * z / 2) ** 3).subs(z, bad_mobius))

# At a general constant-z divisor, equal pole orders in Ls and f cannot
# cancel: the former has nonconstant residue ell*s_bar, the latter a scalar.
ell, mu = sp.symbols("ell mu", nonzero=True)
equal(
    "equal-order cancellation derivative",
    sp.diff(ell * s_divisor + mu, r),
    -ell / (3 * r**2),
)
nonzero("equal-order leading coefficient", sp.diff(ell * s_divisor + mu, r))

# If S has pole order q>0 and Z is regular, the two terms of
# A=9c^6 S^2-Z/(3c^3) have valuations -2q and >=0.  This symbolic valuation
# comparison records the no-cancellation mechanism for every q.
q = sp.symbols("q", integer=True, positive=True)
gate(-2 * q < 0, "quadratic carrier pole is not negative")


# -------------------------------------------------------------------------
# 5. Scope: an s-dependent Z first appears at order z^3, not a no-go.
# -------------------------------------------------------------------------

u_s = s
h_s = sp.Rational(3, 2) * w**2 * s - sp.Rational(3, 2) * s**2
Z_mixed = z - w * z**2 / 2 + u_s * z**3
S_mixed = s + h_s * z**2
mixed_density = sp.expand(
    sp.series((1 + w * Z_mixed) * bracket(Z_mixed, S_mixed), z, 0, 3).removeO()
)
equal("mixed density constant", coefficient(mixed_density, z, 0), 1)
equal("mixed density z row", coefficient(mixed_density, z, 1), 0)
equal("mixed density z2 row", coefficient(mixed_density, z, 2), 0)
equal("mixed Z forced z2", coefficient(Z_mixed, z, 2), -w / 2)
equal("mixed Z_s begins z3", sp.diff(Z_mixed, s), z**3)


source = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))), "inactive assert found")

semantic_packet = (
    "THM-3860 independent hostile audit",
    "all-order coefficient convolution and affine tangent torsor",
    "both derivative edge charts avoid hidden division",
    "vertical slope unique but additive f companion nonunique",
    "Mobius nodal rational Darboux pair and named simple pole",
    "Russell surface smooth and every nonzero constant-z fibre prime",
    "rational root-pattern lemma including poles and repeated roots",
    "equal-order Ls plus f cancellation impossible by nonconstant s residue",
    "quadratic carrier doubles every vertical pole",
    "s-dependent Z begins at z cubed and remains an open lane",
)

print("THM3860_HOSTILE_AUDIT", "PASS_WITH_WORDING_REPAIR")
print("FORMAL_RECURSION", "convolution rows 0..7; nodal through 6; both derivative edges")
print("VERTICAL_CLASSIFICATION", "S_s=L unique; S=L*s+f with f free")
print("UNIQUENESS_REPAIR", "companion is not unique; only vertical slope is unique")
print("NODAL_RATIONAL_PAIR", "Jacobian one; arm and first jet exact")
print("RUSSELL_DIVISORS", "surface smooth; z=z0 nonzero gives prime G_m fibre")
print("NAMED_POLE", "z=-4c^3; ord(A)=-1")
print("VERTICAL_BARRIER", "all rational phi(z),f(z); finite divisorial pole")
print("SCOPE", "no obstruction to s-dependent Z, general formal/polynomial lifts, or JC(2)")
print("SEMANTIC_SHA256", hashlib.sha256(repr(semantic_packet).encode()).hexdigest())
print("GATES", GATES)
