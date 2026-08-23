#!/usr/bin/env python3
"""Independent exact audit for THM-3815's quadratic-r/constant-z^2 cell.

The script independently derives the slope equations, exact resultant,
repeated-root-safe boundary divisibility obstruction, five-term Groebner
certificate, and denominator-free reconstruction.  It contains no Python
assert, so normal and optimized modes execute the same gates.
"""

from __future__ import annotations

import ast
import hashlib
import itertools
import json
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: object, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(label)


def same(left: object, right: object, label: str) -> None:
    difference = left - right  # type: ignore[operator]
    if isinstance(difference, sp.MatrixBase):
        gate(difference.applyfunc(lambda q: sp.factor(sp.cancel(q)))
             == sp.zeros(*difference.shape), label)
        return
    gate(sp.factor(sp.cancel(difference)) == 0, label)


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[0]
THM3810 = REPO / "01-canon/theorems/THM-3810-affine-r-profile-constant-z2-nodal-carriers-have-critical-points.md"
THM3805 = REPO / "01-canon/theorems/THM-3805-quadratic-r-repairs-of-nodal-carriers-have-critical-points.md"
THM3810_SHA256 = "b9c4acfa84df98e050623e15642b8a38d0d703ef84bf9da0de0d0546a9ee52a3"
THM3805_SHA256 = "d2868ba6cfcae9d6114c9329cb822a77a2875fd9dc86f32b2baedeaf76234eca"


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


thm3810_bytes = THM3810.read_bytes()
thm3805_bytes = THM3805.read_bytes()
gate(sha(thm3810_bytes) == THM3810_SHA256, "current THM-3810 hash")
gate(sha(thm3805_bytes) == THM3805_SHA256, "current THM-3805 hash")
gate("leading coefficient is the moving factor `1+a_2r`" in
     thm3810_bytes.decode("utf-8"), "THM-3810 sharpened boundary missing")


# ---------------------------------------------------------------------------
# 1. Hamiltonian packet on the c=1 cubic pseudoplane.
# ---------------------------------------------------------------------------

r, z, e, u = sp.symbols("r z e u")
a0, a1, a2, h = sp.symbols("a0 a1 a2 h")
g = a0 + a1*e + a2*e**2
gp = sp.diff(g, e)
carrier = e**2 - z/sp.Integer(3) + r*g + h*z**2
surface = r**2*e - z**3 + r
L_rz = 1 - 6*h*z
K_re = 1 + 2*r*e


def bracket(left: sp.Expr, right: sp.Expr) -> sp.Expr:
    # Coordinates are ordered (r,z,e), with the canonical brackets used in
    # THM-3810.
    left_r, left_z, left_e = (sp.diff(left, variable) for variable in (r,z,e))
    right_r, right_z, right_e = (sp.diff(right, variable) for variable in (r,z,e))
    return sp.expand(
        3*r**2*(left_r*right_z-left_z*right_r)
        + 9*z**2*(left_r*right_e-left_e*right_r)
        + 3*K_re*(left_z*right_e-left_e*right_z)
    )


Cr = sp.factor(bracket(carrier, r))
Cz = sp.factor(bracket(carrier, z))
Ce = sp.factor(bracket(carrier, e))
Cr_expected = L_rz*r**2 - 9*z**2*(2*e+r*gp)
Cz_expected = 3*g*r**2 - 3*K_re*(2*e+r*gp)
Ce_expected = 9*g*z**2 - L_rz*K_re
same(Cr, Cr_expected, "Hamiltonian r component")
same(Cz, Cz_expected, "Hamiltonian z component")
same(Ce, Ce_expected, "Hamiltonian e component")
same(K_re*Cr-3*z**2*Cz+r**2*Ce, 0, "Casimir syzygy")

# A critical point cannot have z=0: Cr=r^2 first forces r=0, while Ce=-1.
same(Cr.subs(z, 0), r**2, "z=0 first component")
same(Ce.subs({z: 0, r: 0}), -1, "z=0 terminal contradiction")


# ---------------------------------------------------------------------------
# 2. Slope compression with the moving coefficient W=1+a2*r.
# ---------------------------------------------------------------------------

L = 1 - 6*h*z
W = 1 + a2*u*z
M = L*u**2 - 9*a1*u*z
e_slope = M/(18*W)
r_slope = u*z

# Cr=z^2(M-18We), so away from W=0 it gives e=M/(18W).
same(Cr.subs(r, r_slope), z**2*(M-18*W*e),
     "moving-coefficient slope form")

# Clear W from the surface equation and 36W^2 from Ce.  The raw Ce equation
# is quartic; subtracting 3*a2*L*z times p1 cancels its u^4 term.
p1 = sp.expand(18*W*z**2-u**2*z*M-18*W*u)
p2_raw = sp.expand(
    324*a0*z**2*W**2
    + 18*a1*M*z**2*W
    + a2*M**2*z**2
    - 4*L*W*(9*W+u*z*M)
)
p2 = sp.expand(p2_raw-3*a2*L*z*p1)

same(surface.subs({r: r_slope, e: e_slope}), -z*p1/(18*W),
     "surface reconstruction identity")
same(36*W**2*Ce.subs({r: r_slope, e: e_slope}), p2_raw,
     "cleared Ce identity")
same(p2_raw, p2+3*a2*L*z*p1, "quartic-to-cubic row operation")

gate(sp.degree(p1, u) == 4, "p1 generic u degree")
gate(sp.degree(p2_raw, u) == 4, "raw Ce generic u degree")
gate(sp.degree(p2, u) == 3, "reduced Ce generic u degree")

N = 4*L-9*a1*a2*z**2
same(sp.LC(sp.Poly(p1, u)), -z*L, "p1 projective leading coefficient")
same(sp.LC(sp.Poly(p2, u)), -z*L*N, "p2 projective leading coefficient")

# N=0 is only a degree drop of p2, not a forbidden divisor: p1 retains
# degree four whenever zL is nonzero and therefore has no root at infinity.
gate(not sp.LC(sp.Poly(p1, u)).has(a1, a2),
     "p1 leading coefficient is independent of the p2 degree-drop divisor")

# At W=0 one has u=-1/(a2*z).  A common root of p1,p2 then forces the
# projected moving-denominator boundary B=0.
B = L + 9*a1*a2*z**2
u_W = -1/(a2*z)
same(p1.subs(u, u_W), -B/(a2**4*z**3), "W boundary projection through p1")
same(M.subs(u, u_W), B/(a2**2*z**2), "W boundary projection through M")


# ---------------------------------------------------------------------------
# 3. Exact resultant and degree-sixteen residual.
# ---------------------------------------------------------------------------

resultant = sp.factor(sp.resultant(p1, p2, u))
residual = sp.cancel(resultant/(11664*z**3*L*B**2))
gate(sp.denom(residual) == 1, "resultant residual denominator")
residual = sp.expand(residual)
residual_poly = sp.Poly(residual, z)
gate(residual_poly.degree() == 16, "residual degree sixteen")
same(residual_poly.LC(), 314928*a2**4*h**4,
     "residual fixed leading coefficient")
same(residual_poly.eval(0), 144, "residual nonzero constant")
same(resultant, 11664*z**3*L*B**2*residual,
     "complete resultant factorization")

# The a2=0 compression is exactly twice the THM-3810 p2 convention, while
# h=0 is the already audited THM-3805 family.  The proof below inverts
# neither a0 nor a1, so their zero strata remain inside the genuine cell.
M_affine = L*u**2-9*a1*u*z
p1_affine = 18*z**2-u**2*z*M_affine-18*u
p2_affine = 162*a0*z**2+9*a1*M_affine*z**2-18*L-2*L*u*z*M_affine
same(p1.subs(a2, 0), p1_affine, "a2=0 p1 inheritance")
same(p2.subs(a2, 0), 2*p2_affine, "a2=0 p2 inheritance")


# ---------------------------------------------------------------------------
# 4. Repeated-root-safe boundary exclusion and exact Groebner certificate.
# ---------------------------------------------------------------------------

# If every root of residual were supported on V(L*B), factorization over an
# algebraically closed field would imply residual | L*B*residual'.  Divide
# symbolically over Q(a0,a1,a2,h).  In the genuine cell, 8*a2*h^3 clears
# every denominator in the remainder.
parameter_field = sp.QQ.frac_field(a0, a1, a2, h)
dividend = sp.Poly(sp.expand(L*B*sp.diff(residual, z)), z,
                   domain=parameter_field)
divisor = sp.Poly(residual, z, domain=parameter_field)
quotient, remainder = sp.div(dividend, divisor)
rho = sp.Poly(sp.cancel(8*a2*h**3*remainder.as_expr()), z)
gate(all(sp.denom(coefficient) == 1 for coefficient in rho.all_coeffs()),
     "denominator-free logarithmic remainder")
gate(rho.degree() == 15, "logarithmic remainder degree")

# Only the first five low coefficients are needed.  Exact Buchberger
# reduction gives a2^2*h^5 in their ideal.  Thus they cannot all vanish when
# a2*h is nonzero.
E = []
for exponent in range(5):
    coefficient = rho.coeff_monomial(z**exponent)
    primitive = sp.Poly(coefficient, a0, a1, a2, h).primitive()[1]
    E.append(sp.expand(primitive.as_expr()))
groebner = sp.groebner(E, a0, a1, a2, h, order="grevlex")
membership_remainder = groebner.reduce(a2**2*h**5)[1]
same(membership_remainder, 0, "five-coefficient Groebner monomial certificate")
gate(any(sp.Poly(poly, a0, a1, a2, h).monoms() == [(0,0,2,5)]
         for poly in groebner.polys),
     "explicit a2^2*h^5 basis element")

residual_sha = sha(str(residual_poly.as_expr()).encode())
five_coefficient_sha = sha("\n".join(str(item) for item in E).encode())
groebner_sha = sha("\n".join(str(poly.as_expr()) for poly in groebner.polys).encode())


# ---------------------------------------------------------------------------
# 5. Finite-field hostile scans (controls only, not proof dependencies).
# ---------------------------------------------------------------------------

residual_coefficients = residual_poly.all_coeffs()
finite_field_cases = 0
for prime in (5, 7):
    survivors = []
    for values in itertools.product(range(prime), range(prime),
                                    range(1, prime), range(1, prime)):
        finite_field_cases += 1
        substitution = dict(zip((a0, a1, a2, h), values))
        specialized = sp.Poly.from_list(
            [int(coefficient.subs(substitution)) % prime
             for coefficient in residual_coefficients],
            gens=z, modulus=prime,
        )
        L_specialized = sp.Poly(1-6*values[3]*z, z, modulus=prime)
        B_specialized = sp.Poly(
            1-6*values[3]*z+9*values[1]*values[2]*z**2,
            z, modulus=prime,
        )
        boundary_remainder = (L_specialized*B_specialized
                              * specialized.diff()).rem(specialized)
        if boundary_remainder.is_zero:
            survivors.append(values)
    gate(survivors == [], f"finite-field boundary survivors p={prime}")


# ---------------------------------------------------------------------------
# 6. Off-boundary reconstruction and logical conclusion.
# ---------------------------------------------------------------------------

# At a residual root z0 off z*L*B, the resultant supplies a finite common
# u-root: p1 has nonzero quartic leading coefficient even on N=0.  B!=0
# prevents W=0.  These symbolic identities then reconstruct the source and
# all Hamiltonian equations without dividing by N, u, or K.
same(Cr.subs({r: r_slope, e: e_slope}), 0,
     "reconstructed Cr identity")
same(Ce.subs({r: r_slope, e: e_slope}),
     (p2+3*a2*L*z*p1)/(36*W**2),
     "reconstructed Ce identity")
same(surface.subs({r: r_slope, e: e_slope}), -z*p1/(18*W),
     "reconstructed surface identity")
same(p1.subs(u, 0), 18*z**2, "common slope root is nonzero")

audit_tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(audit_tree)),
     "probe contains optimization-erased assert")

semantic = {
    "carrier": "e^2-z/3+r(a0+a1e+a2e^2)+hz^2",
    "compression": "W=1+a2uz; M=Lu^2-9a1uz; quartic p1 plus cubic p2",
    "resultant": "11664*z^3*L*(L+9a1a2z^2)^2*H16",
    "residual": "LC=314928*a2^4*h^4; H(0)=144",
    "boundary": "H|LBH' would force five coefficients whose ideal contains a2^2*h^5",
    "reconstruction": "off zLB, p1 keeps degree4, W is nonzero, and source/Cr/Ce/Cz reconstruct",
    "status": "independent hostile audit; a2=0 THM3810 and h=0 THM3805",
}
semantic_sha = sha(json.dumps(semantic, sort_keys=True,
                              separators=(",", ":")).encode())

print("AUDIT=THM-3815-quadratic_r_profile_plus_constant_z2")
print("CELL=genuine_a2_times_h_nonzero;axes_inherit_THM3805_THM3810")
print("SLOPE=W=1+a2*u*z;M=L*u^2-9*a1*u*z;p1_degree4;p2_degree3")
print("FORBIDDEN=z=0;L=0;W=0_projects_to_B=L+9*a1*a2*z^2=0")
print("HARMLESS=p2_leading_N=4L-9a1a2z^2;u=0;K=0")
print("RESULTANT=11664*z^3*L*B^2*H16")
print("H16=degree16;LC=314928*a2^4*h^4;constant=144")
print("BOUNDARY_CERT=first_five_log_remainders_generate_a2^2*h^5")
print(f"RESIDUAL_SHA256={residual_sha}")
print(f"FIVE_COEFFICIENT_SHA256={five_coefficient_sha}")
print(f"GROEBNER_SHA256={groebner_sha}")
print(f"FINITE_FIELD_HOSTILES={finite_field_cases};survivors=0")
print(f"CHECKS={CHECKS}")
print(f"SEMANTIC_SHA256={semantic_sha}")
print("STATUS=INDEPENDENT_HOSTILE_AUDIT_PASS")
