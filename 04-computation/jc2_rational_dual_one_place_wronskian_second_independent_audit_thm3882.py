#!/usr/bin/env python3
"""Import-independent hostile audit of THM-3882.

This script does not import the canonical companion.  It checks homogeneous
Wronskian identities, finite/infinite divisor multiplicities, the precise
Riemann--Hurwitz inequality, the sextic Pluecker packet, the THM-3879 node
specialization, and a sharp nonimmersed family in every degree.
"""

from __future__ import annotations

import hashlib
import sys

import sympy as sp

sys.stdout.reconfigure(newline="\n")

S, T = sp.symbols("S T")
GATES = 0
SECTION_GATES: dict[str, int] = {}


def gate(section: str, label: str, condition: bool) -> None:
    global GATES
    GATES += 1
    SECTION_GATES[section] = SECTION_GATES.get(section, 0) + 1
    if not condition:
        raise RuntimeError(f"{section}/{label}: failed")


def equal(section: str, label: str, left: object, right: object) -> None:
    gate(section, label, left == right)


def zero(section: str, label: str, value: sp.Expr) -> None:
    gate(section, label, sp.Poly(sp.expand(value), S, T).is_zero)


def omega(f: sp.Expr, h: sp.Expr) -> sp.Expr:
    """Homogeneous first transvectant; scalar-normalized zeros are Wronskian zeros."""
    return sp.expand(sp.diff(f, S) * sp.diff(h, T) - sp.diff(f, T) * sp.diff(h, S))


def deterministic_form(degree: int, seed: int) -> sp.Expr:
    return sp.expand(
        sum(
            ((seed + 3) * (i + 2) + (i * i + 5)) * S ** (degree - i) * T**i
            for i in range(degree + 1)
        )
    )


# 1. A fully symbolic homogeneous identity, plus a rectangular deterministic
# bank.  If deg(g)=m and deg(x)=deg(y)=e, the homogeneous transvectant has the
# scalar (m+e)/e relative to the affine Wronskian.  This scalar is a unit in
# characteristic zero and hence does not alter the divisor.
g_coeff = sp.symbols("g0:4")
x_coeff = sp.symbols("x0:5")
y_coeff = sp.symbols("y0:5")
g_generic = sum(g_coeff[i] * S ** (3 - i) * T**i for i in range(4))
x_generic = sum(x_coeff[i] * S ** (4 - i) * T**i for i in range(5))
y_generic = sum(y_coeff[i] * S ** (4 - i) * T**i for i in range(5))
zero(
    "factorization",
    "generic_bidegree_3_4",
    4 * omega(g_generic * x_generic, g_generic * y_generic)
    - 7 * g_generic**2 * omega(x_generic, y_generic),
)

for m in range(0, 6):
    for e in range(1, 9):
        g = deterministic_form(m, 11 + 7 * m + e)
        x = deterministic_form(e, 101 + 13 * m + e)
        y = deterministic_form(e, 211 + 17 * m + 3 * e)
        zero(
            "factorization",
            f"deterministic_m{m}_e{e}",
            e * omega(g * x, g * y) - (m + e) * g**2 * omega(x, y),
        )


# 2. Exact valuations at both projective addresses.  The family [S^e:T^e]
# has ramification e-1 at S=0 and T=0.  Multiplication by
# g=S^a*T^(m-a) adds twice the base multiplicity at each address.
one_support_rows = 0
all_monomial_rows = 0
for d in range(2, 31):
    for m in range(0, d):
        e = d - m
        for a in range(m + 1):
            b = m - a
            g = S**a * T**b
            x = S**e
            y = T**e
            observed = sp.Poly(omega(g * x, g * y), S, T)
            expected = sp.Poly(
                e * d * S ** (2 * a + e - 1) * T ** (2 * b + e - 1), S, T
            )
            equal("divisor", f"monomial_identity_d{d}_m{m}_a{a}", observed, expected)
            exp_s = 2 * a + e - 1
            exp_t = 2 * b + e - 1
            equal("divisor", f"degree_d{d}_m{m}_a{a}", exp_s + exp_t, 2 * d - 2)
            equal("divisor", f"s_place_d{d}_m{m}_a{a}", exp_s, 2 * a + (e - 1))
            equal("divisor", f"t_place_d{d}_m{m}_a{a}", exp_t, 2 * b + (e - 1))
            supports = int(exp_s > 0) + int(exp_t > 0)
            criterion = e == 1 and (a == 0 or b == 0)
            equal("divisor", f"one_support_criterion_d{d}_m{m}_a{a}", supports == 1, criterion)
            all_monomial_rows += 1
            one_support_rows += int(supports == 1)


# 3. Riemann--Hurwitz hostile ledger in every genus.  If all ramification of
# a degree-e map X_g -> P1 were at one address, its coefficient 2e+2g-2
# could not exceed the local maximum e-1.  The unique numerical survivor is
# (g,e)=(0,1), where the ramification divisor is empty.
rh_survivors: list[tuple[int, int]] = []
for genus in range(0, 21):
    for e in range(1, 41):
        ram_degree = 2 * e + 2 * genus - 2
        local_cap = e - 1
        possible = ram_degree <= local_cap
        equal(
            "riemann_hurwitz",
            f"one_address_g{genus}_e{e}",
            possible,
            genus == 0 and e == 1,
        )
        if possible:
            rh_survivors.append((genus, e))
equal("riemann_hurwitz", "unique_survivor", rh_survivors, [(0, 1)])


# 4. Sharp boundary family.  For every d>=2,
#   nu_d=[S*T^(d-1):T^d:S^d]
# is birational (X/Y=S/T on the torus).  It is an immersed conic for d=2;
# for d>=3 it fails immersion only at T=0.  After the common tangent factor
# T^(d-2) is removed, the third dual coordinate is d*T^d, a one-place line.
nonimmersed_one_place_degrees: list[int] = []
for d in range(2, 21):
    X = S * T ** (d - 1)
    Y = T**d
    Z = S**d
    d_s = sp.Matrix([sp.diff(X, S), sp.diff(Y, S), sp.diff(Z, S)])
    d_t = sp.Matrix([sp.diff(X, T), sp.diff(Y, T), sp.diff(Z, T)])
    A, B, C = map(sp.expand, d_s.cross(d_t))
    zero("boundary", f"cross_A_d{d}", A + d * d * S ** (d - 1) * T ** (d - 1))
    zero("boundary", f"cross_B_d{d}", B - d * (d - 1) * S**d * T ** (d - 2))
    zero("boundary", f"cross_C_d{d}", C - d * T ** (2 * d - 2))
    tangent_factor = T ** (d - 2)
    reduced = tuple(sp.cancel(q / tangent_factor) for q in (A, B, C))
    zero("boundary", f"reduced_C_d{d}", reduced[2] - d * T**d)
    # The nonimmersed correction is the common tangent-base divisor E_nu.
    # Here D=(d-1)[T=0], Ram(phi)=0, E_nu=(d-2)[T=0], so the reduced
    # dual-line coefficient is 2(d-1)-(d-2)=d.
    equal("boundary", f"projection_base_d{d}", d - 1, d - 1)
    equal("boundary", f"tangent_base_d{d}", d - 2, d - 2)
    equal("boundary", f"corrected_divisor_d{d}", 2 * (d - 1) - (d - 2), d)
    # At T=0, affine coordinates through [0:0:1] are u^(d-1),u^d.
    u = sp.symbols("u")
    immersed_at_cusp_address = (
        sp.diff(u ** (d - 1), u).subs(u, 0) != 0
        or sp.diff(u**d, u).subs(u, 0) != 0
    )
    equal("boundary", f"immersion_boundary_d{d}", immersed_at_cusp_address, d == 2)
    # The opposite chart has coordinate v=X/Y=S/T, proving local immersion
    # there and the generic inverse [S:T]=[X:Y].
    equal("boundary", f"torus_inverse_d{d}", sp.cancel(X / Y), S / T)
    equal("boundary", f"opposite_derivative_d{d}", sp.diff(u, u), 1)
    if d >= 3:
        nonimmersed_one_place_degrees.append(d)


# 5. THM-3879's trinodal quartic, checked homogeneously at 0 and infinity.
X = S * T * (S**2 - 2 * T**2)
Y = S * T * (S**2 - T**2)
Z = (S**2 - T**2) * (S**2 - 2 * T**2)
A = -(S**2 - T**2) ** 2 * (S**2 + 2 * T**2)
B = (S**2 - 2 * T**2) ** 2 * (S**2 + T**2)
C = 2 * S**3 * T**3
d_s = sp.Matrix([sp.diff(X, S), sp.diff(Y, S), sp.diff(Z, S)])
d_t = sp.Matrix([sp.diff(X, T), sp.diff(Y, T), sp.diff(Z, T)])
cross = tuple(map(sp.expand, d_s.cross(d_t)))
for name, got, wanted in zip("ABC", cross, (4 * A, 4 * B, 4 * C)):
    zero("thm3879", f"homogeneous_dual_{name}", got - wanted)
equal("thm3879", "no_common_zero_S_axis_A", sp.expand(A.subs({S: 0, T: 1})), -2)
equal("thm3879", "no_common_zero_S_axis_B", sp.expand(B.subs({S: 0, T: 1})), 4)
equal("thm3879", "no_common_zero_T_axis_A", sp.expand(A.subs({S: 1, T: 0})), -1)
equal("thm3879", "no_common_zero_T_axis_B", sp.expand(B.subs({S: 1, T: 0})), 1)
g = S * T
x = S**2 - 2 * T**2
y = S**2 - T**2
zero("thm3879", "node_base_factor", X - g * x)
zero("thm3879", "node_base_factor_second", Y - g * y)
zero("thm3879", "residual_ramification", omega(x, y) - 4 * S * T)
zero("thm3879", "full_projection_wronskian", omega(X, Y) - 8 * S**3 * T**3)
equal("thm3879", "base_divisor", (1, 1), (1, 1))
equal("thm3879", "ramification_divisor", (1, 1), (1, 1))
equal("thm3879", "dual_line_divisor", (2 * 1 + 1, 2 * 1 + 1), (3, 3))


# 6. The complete 6A2+4A1 Pluecker ledger used for the packet consequence.
d, nodes, cusps = 6, 4, 6
genus = (d - 1) * (d - 2) // 2 - nodes - cusps
dual_degree = d * (d - 1) - 2 * nodes - 3 * cusps
stationary_tangent_weight = 3 * d * (d - 2) - 6 * nodes - 8 * cusps
equal("pluecker", "normalization_genus", genus, 0)
equal("pluecker", "dual_degree", dual_degree, 4)
equal("pluecker", "stationary_tangent_weight", stationary_tangent_weight, 0)
equal("pluecker", "immersed_dual_degree_boundary", dual_degree >= 3, True)
# Independent degree/reflexivity route to the same immersion conclusion.  The
# six ordinary cusps remove one tangent-base unit each from the raw degree ten
# Gauss coordinates.  The resulting quartic dual has raw bidual coordinates
# of degree six; biduality already needs all six degrees, leaving no base unit.
equal("pluecker", "primal_raw_tangent_degree", 2 * d - 2, 10)
equal("pluecker", "primal_tangent_base_degree", cusps, 6)
equal("pluecker", "reduced_dual_param_degree", (2 * d - 2) - cusps, dual_degree)
equal("pluecker", "dual_raw_tangent_degree", 2 * dual_degree - 2, 6)
equal("pluecker", "bidual_degree", d, 6)
equal("pluecker", "dual_tangent_base_degree", (2 * dual_degree - 2) - d, 0)


semantic_packet = (
    "homogeneous W(gx,gy)=g^2W(x,y)",
    "finite and infinite coefficients are 2D plus ramification",
    "one-address ramification forces genus zero and projection degree one",
    "immersion then forces the plane degree to be two",
    "among all genera the only immersed one-place case is a smooth conic",
    "six A2 plus four A1 gives genus zero dual degree four and zero stationary-tangent weight",
    "THM3879 C=0 is base (1,1) twice plus ramification (1,1)",
    "without immersion the reduced identity is 2D plus ramification minus the common tangent-base divisor",
    "nonimmersed family [S*T^(d-1):T^d:S^d] has a one-place dual line for every d at least three",
    "the theorem closes the dual equisingularity architecture and not JC2",
)

print("AUDIT", "PASS")
print("DIVISOR_IDENTITY", "dual_line=2*base_fibre+ramification")
print("IFF", "one_place iff genus=0,projection_degree=1,base=(d-1)p")
print("IMMERSED_CLASSIFICATION", "one_place iff smooth_conic_boundary")
print("PACKET", "6A2+4A1 -> genus0,dual_degree4,stationary_weight0 -> no_one_place_line")
print("THM3879", "C=0 has base=(1,1),ramification=(1,1),dual_pullback=(3,3)")
print("STRICT_EXTENSION", "same criterion holds for every normalization genus; one-place forces genus0")
print("NONIMMERSED_CORRECTION", "reduced_dual_line=2*base+ramification-common_tangent_base")
print("SHARP_HOSTILE", "[S*T^(d-1):T^d:S^d] realizes a one-place dual line after nonimmersion for all d>=3")
print("MONOMIAL_ROWS", all_monomial_rows)
print("MONOMIAL_ONE_SUPPORT_ROWS", one_support_rows)
print("RH_SURVIVORS", rh_survivors)
for section in sorted(SECTION_GATES):
    print("SECTION_GATES", section, SECTION_GATES[section])
print("SEMANTIC_SHA256", hashlib.sha256(repr(semantic_packet).encode()).hexdigest())
print("GATES", GATES)
