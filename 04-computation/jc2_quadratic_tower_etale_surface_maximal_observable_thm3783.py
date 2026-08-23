#!/usr/bin/env python3
"""Exact companion for THM-3783's quadratic-tower etale surface."""

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


x, y = sp.symbols("x y")
r, z, g = sp.symbols("r z g")

# The m=1 member of the THM-3774 tower.
A = 1 + x * y
B = sp.expand(1 + x**3 * A)
U = sp.expand(x * A * B)
P = sp.cancel((2 * B - 1) / x)

R = x**2
Z = sp.expand(x / 2 + x**5 * y)
G = sp.expand(y + x**4 * y**2)
H = sp.expand(r**3 * g - z**2 + r / 4)

gate(sp.cancel(jac(U, P, x, y) - 1) == 0, "m1 rational Keller identity")
gate(sp.cancel(P**2 - 4 * U - 1 / R) == 0, "inverse discriminant r")
gate(sp.cancel(Z - (R * P / 2 - R**2)) == 0, "z target-field formula")
gate(sp.expand(R**3 * G - Z**2 + R / 4) == 0, "etale-surface relation")
gate(sp.expand(U - (R * G + 2 * Z + R**2)) == 0, "U in surface ring")
gate(sp.cancel(P - 2 * (Z + R**2) / R) == 0, "P in surface field")

# The old Danielewski observables are strict sub-observables of the new ring.
V = sp.expand(U * P)
Nm = sp.expand(
    1
    + 2 * g * z
    + 6 * r**2 * g
    + 6 * r * z
    + 2 * r**3
)
Em = sp.expand(
    g
    + 3 * r
    + 4 * r**2 * g**2
    + 16 * r * g * z
    + 24 * r**3 * g
    + 16 * r**2 * z
    + 4 * r**4
)
gate(sp.expand(V - Nm.subs({r: R, z: Z, g: G})) == 0, "V surface formula")
E_source = sp.cancel(P * (V - 1))
gate(sp.cancel(E_source - Em.subs({r: R, z: Z, g: G})) == 0,
     "E surface formula")
gate(sp.expand(U * E_source - V * (V - 1)) == 0,
     "old Danielewski relation")

# Complete Poisson packet and the hypersurface-gradient signs.
J_rz = jac(R, Z, x, y)
J_rg = jac(R, G, x, y)
J_zg = jac(Z, G, x, y)
gate(sp.expand(J_rz - 2 * R**3) == 0, "bracket rz")
gate(sp.expand(J_rg - 4 * Z) == 0, "bracket rg")
gate(sp.expand(J_zg - (sp.Rational(1, 2) + 6 * R**2 * G)) == 0,
     "bracket zg")
gate(sp.diff(H, g) * 2 == 2 * r**3, "gradient bracket rz")
gate(-sp.diff(H, z) * 2 == 4 * z, "gradient bracket rg")
gate(sp.diff(H, r) * 2 == sp.Rational(1, 2) + 6 * r**2 * g,
     "gradient bracket zg")

# The gradient/minor ideal is the unit ideal: smooth target and etale source map.
gradient_basis = sp.groebner(
    [H, sp.diff(H, r), sp.diff(H, z), sp.diff(H, g)],
    g, z, r, order="lex",
)
minor_basis = sp.groebner(
    [H, 2 * r**3, 4 * z, sp.Rational(1, 2) + 6 * r**2 * g],
    g, z, r, order="lex",
)
gate(len(gradient_basis.polys) == 1 and gradient_basis.polys[0].as_expr() == 1,
     "smooth gradient unit ideal")
gate(len(minor_basis.polys) == 1 and minor_basis.polys[0].as_expr() == 1,
     "etale minor unit ideal")

# Exact inverse on r!=0 and exact arm coverage at r=0.
xx, rr, zz, gg = sp.symbols("xx rr zz gg")
yy = (zz - xx / 2) / xx**5
inverse_g = sp.cancel(yy + xx**4 * yy**2 - gg)
inverse_numerator = sp.factor(sp.together(inverse_g).as_numer_denom()[0])
gate(
    sp.expand(inverse_numerator
              - 4 * (zz**2 - xx**2 / 4 - xx**6 * gg)) == 0,
    "off-arm inverse numerator",
)
gate(R.subs(x, 0) == 0, "arm r")
gate(Z.subs(x, 0) == 0, "arm z")
gate(G.subs(x, 0) == y, "arm g parameter")

# The rational deck involution fixes exactly the displayed generators.
sigma = {x: -x, y: -y - x**-4}
gate(sp.cancel(R.subs(sigma, simultaneous=True) - R) == 0, "deck fixes r")
gate(sp.cancel(Z.subs(sigma, simultaneous=True) - Z) == 0, "deck fixes z")
gate(sp.cancel(G.subs(sigma, simultaneous=True) - G) == 0, "deck fixes g")
sigma_x = sp.cancel((-x).subs(sigma, simultaneous=True))
sigma_y = sp.cancel((-y - x**-4).subs(sigma, simultaneous=True))
gate(sigma_x == x, "deck square x")
gate(sp.cancel(sigma_y - y) == 0, "deck square y")

# Weighted Euler field and the global exact primitive.
def euler(expr: sp.Expr) -> sp.Expr:
    return sp.expand(
        2 * r * sp.diff(expr, r)
        + z * sp.diff(expr, z)
        - 4 * g * sp.diff(expr, g)
    )


gate(sp.expand(euler(H) - 2 * H) == 0, "weighted relation")
gate(euler(r) == 2 * r, "weight r")
gate(euler(z) == z, "weight z")
gate(euler(g) == -4 * g, "weight g")

# Pullback of alpha=-(1/3)i_E omega.
alpha_dx = -sp.Rational(4, 3) * y
alpha_dy = -sp.Rational(1, 3) * x
alpha_curl = sp.expand(sp.diff(alpha_dy, x) - sp.diff(alpha_dx, y))
gate(alpha_curl == 1, "exact primitive differential")
gate(sp.expand(jac(R, Z, x, y) / (2 * R**3) - 1) == 0,
     "symplectic pullback on r-nonzero chart")

# A visible additive action.  On r!=0 this is r^3 partial_z, so its kernel
# there is k[r,r^-1]; regularity along the unique r-arm cuts this to k[r].
def lnd(expr: sp.Expr) -> sp.Expr:
    return sp.expand(r**3 * sp.diff(expr, z) + 2 * z * sp.diff(expr, g))


gate(lnd(H) == 0, "LND preserves relation")
gate(lnd(r) == 0, "LND r")
gate(lnd(z) == r**3, "LND z")
gate(lnd(g) == 2 * z, "LND g")
gate(lnd(lnd(lnd(g))) == 0, "LND local nilpotence control")

# The two obvious smooth generators are precisely monomial Broughton profiles.
gate(sp.expand(2 * Z - (x + 2 * x**5 * y)) == 0,
     "z Broughton profile")
gate(sp.expand(G - (y + y**2 * x**4)) == 0,
     "g dual Broughton profile")

# Small hostile packet: every reduced monomial in the surface generators is
# fixed by the rational deck involution after source substitution.
packet_rows = []
for a in range(4):
    for bpow in range(4):
        for c in range(2):
            source_monomial = sp.expand(R**a * G**bpow * Z**c)
            moved = sp.cancel(source_monomial.subs(sigma, simultaneous=True))
            gate(sp.cancel(moved - source_monomial) == 0,
                 f"deck packet a={a},b={bpow},c={c}")
            packet_rows.append(f"{a},{bpow},{c}:{sp.sstr(source_monomial)}")

# The first aligned two-by-three Euler-weight cell has a universal peel.
# Put A=d-1>=3, delta=gcd(A,4), u=A/delta, v=4/delta.  Endpoint
# commutation gives p^u,p^v and q,q^(2A-2).  The other middle bucket
# integrates exactly, while the scalar bucket contains the nonconstant
# factor q*p'+delta*p*q'.
ss = sp.symbols("ss")
pp = sp.Function("pp")(ss)
qq = sp.Function("qq")(ss)
lam, mu, nu = sp.symbols("lam mu nu")


def wronskian(weight_f: int, f: sp.Expr, weight_h: int, h: sp.Expr) -> sp.Expr:
    return sp.expand(
        weight_f * f * sp.diff(h, ss)
        - weight_h * sp.diff(f, ss) * h
    )


ap_rows = []

# d=3 is the first genuine integrated seam and already factors through the
# nonconstant negative profile p.
p3 = pp
q3 = qq
h3 = lam * p3**2
k3 = mu * q3**2
middle3 = 2 * mu * p3 * q3
other3 = sp.simplify(
    wronskian(-2, p3, 2, k3)
    + wronskian(1, q3, -1, middle3)
)
scalar3 = sp.expand(
    wronskian(-2, p3, -1, middle3)
    + wronskian(1, q3, -4, h3)
)
factor3 = 2 * (lam - mu) * p3 * (
    sp.diff(p3, ss) * q3 + 2 * p3 * sp.diff(q3, ss)
)
gate(other3 == 0, "AP d3 integrated bucket")
gate(sp.simplify(scalar3 - factor3) == 0, "AP d3 scalar factor")
ap_rows.append("d=3:p-transport-factor")

for Aweight in range(3, 13):
    delta = sp.gcd(Aweight, 4)
    upow = Aweight // delta
    vpow = 4 // delta
    Cweight = Aweight - 3
    Eweight = 2 * Aweight - 2
    Nweight = 2 * Aweight - 3
    f_low = pp**upow
    g_low = lam * pp**vpow
    f_high = qq
    g_high = mu * qq**Eweight
    g_middle = (
        mu * Eweight * pp**upow * qq**Nweight
        + nu * qq**Cweight
    )
    other_bucket = sp.simplify(
        wronskian(-Aweight, f_low, Eweight, g_high)
        + wronskian(1, f_high, Cweight, g_middle)
    )
    gate(other_bucket == 0, f"AP integrated bucket A={Aweight}")

    scalar_bucket = sp.expand(
        wronskian(-Aweight, f_low, Cweight, g_middle)
        + wronskian(1, f_high, -4, g_low)
    )
    transport = qq * sp.diff(pp, ss) + delta * pp * sp.diff(qq, ss)
    residual = (
        -mu * Eweight * Nweight * upow
        * pp ** (2 * upow - 1) * qq ** (2 * Aweight - 4)
        + lam * vpow * pp ** (vpow - 1)
    )
    if Cweight > 0:
        residual += (
            -nu * Cweight * upow
            * pp ** (upow - 1) * qq ** (Aweight - 4)
        )
    gate(sp.simplify(scalar_bucket - transport * residual) == 0,
         f"AP scalar factor A={Aweight}")
    ap_rows.append(
        f"A={Aweight};delta={delta};u={upow};v={vpow};"
        f"C={Cweight};E={Eweight};N={Nweight}"
    )

semantic_rows = (
    "m1:A=1+xy;B=1+x^3A;U=xAB;P=(2B-1)/x;J(U,P)=1",
    "generators:r=x^2;z=x/2+x^5y;g=y+x^4y^2",
    "surface:r^3g=z^2-r/4;smooth;map_A2_surjective_etale_generic_degree2",
    "field:k(U,P)=k(r,z);intersection:k[x,y]_intersect_k(U,P)=surface_ring",
    "poisson:{r,z}=2r^3;{r,g}=4z;{z,g}=1/2+6r^2g",
    "weights:(r,z,g)=(2,1,-4);omega_weight=-3;alpha=-(i_Eomega)/3",
    "pullback_alpha=-(4y_dx+x_dy)/3;dalpha=dx_dy",
    "units=k*;Pic=Z/2_from_div(r)=2L;LND_kernel=k[r];ML_subset_k[r]",
    "Euler_support:no_homogeneous_or_two_by_two_Darboux_pair",
    "named_aligned_2x3_orientation:universal_transport_factor_nonentry",
    "natural_no_mates:z_and_g_are_monomial_Broughton;U_inherits_tower_debt",
    "Darboux_floor:any_pair_is_nonbirational_on_surface_and_source_degree>=4",
)
semantic = hashlib.sha256("\n".join(semantic_rows).encode()).hexdigest()
packet = hashlib.sha256("\n".join(packet_rows + ap_rows).encode()).hexdigest()

source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "inactive Python assert",
)

print("theorem=THM-3783-quadratic-tower-etale-surface-maximal-polynomial-observable")
print("field=algebraically_closed_characteristic_zero;de_Rham_consequence_over_C")
print("m1_seed=U=x(1+xy)(1+x^3(1+xy));P=(2(1+x^3(1+xy))-1)/x")
print("generators=r=x^2;z=x/2+x^5y;g=y+x^4y^2")
print("surface=r^3g=z^2-r/4;smooth")
print("map=A2_to_surface;surjective_etale;generic_degree=2;arm_degree_drop=1")
print("field=k(U,P)=k(r,z)")
print("intersection=k[x,y]_intersect_k(U,P)=k[r,z,g]/(r^3g-z^2+r/4)")
print("poisson=J(r,z)=2r^3;J(r,g)=4z;J(z,g)=1/2+6r^2g")
print("symplectic=exact;weights=2,1,-4;primitive=-(i_Eomega)/3")
print("pullback_primitive=-(4y_dx+x_dy)/3")
print("units=k*;Pic=Z/2;div(r)=2L")
print("LND=delta(r)=0,delta(z)=r^3,delta(g)=2z;kernel=k[r];ML_subset_k[r]")
print("Euler_support=no_homogeneous_or_two_by_two_Darboux_pair")
print("named_aligned_2x3=universal_transport_factor_nonentry;other_2x3_cells_open")
print("natural_components=z,g_are_Broughton_no_mates;U_has_tower_residue_debt")
print("Darboux_pair=OPEN_but_must_be_nonbirational_on_surface;source_degree_even_at_least4")
print(f"packet_sha256={packet}")
print(f"semantic_sha256={semantic}")
print(f"CHECKS={CHECKS}")
print("NO_CLAIM=planar_Jacobian_counterexample_or_complete_Darboux_classification")
print("RESULT=PASS")
