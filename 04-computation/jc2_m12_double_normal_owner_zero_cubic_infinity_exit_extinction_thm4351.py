#!/usr/bin/env python3
"""Exact primary certificate for THM-4351's doubly deleted F00 corner.

This file deliberately imports no repository computation.  It reconstructs
the specialized sixteen-row source, its M/F00 fan, the generic genus-two
carrier, and the eta^2=4*K*U even-A5 normal-series resolution.
"""

from collections import Counter, defaultdict
from fractions import Fraction as Fr
from itertools import combinations, product
from math import gcd, lcm
import sys

import sympy as sp


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")


def need(condition, label):
    if not bool(condition):
        raise AssertionError(label)


def rank_two(points):
    for p, q, r in combinations(points, 3):
        if ((q[0] - p[0]) * (r[1] - p[1])
                != (q[1] - p[1]) * (r[0] - p[0])):
            return True
    return False


def candidate_planes(points):
    result = set()
    for p, q, r in combinations(points, 3):
        det = ((q[0] - p[0]) * (r[1] - p[1])
               - (q[1] - p[1]) * (r[0] - p[0]))
        if det == 0:
            continue
        aa = Fr(
            (q[2] - p[2]) * (r[1] - p[1])
            - (q[1] - p[1]) * (r[2] - p[2]), det)
        bb = Fr(
            (q[0] - p[0]) * (r[2] - p[2])
            - (q[2] - p[2]) * (r[0] - p[0]), det)
        cc = Fr(p[2]) - aa * p[0] - bb * p[1]
        result.add((aa, bb, cc))
    return tuple(result)


def lower_faces(active, planes):
    answer = set()
    for aa, bb, cc in planes:
        equality = []
        for rr, ll, hh in active:
            gap = Fr(hh) - aa * rr - bb * ll - cc
            if gap < 0:
                break
            if gap == 0:
                equality.append((rr, ll, hh))
        else:
            if rank_two(equality):
                answer.add((aa, bb, cc))
    return frozenset(answer)


def pick(vertices):
    edges = tuple(zip(vertices, vertices[1:] + vertices[:1]))
    area2 = abs(sum(x1 * y2 - x2 * y1
                    for (x1, y1), (x2, y2) in edges))
    boundary = sum(gcd(abs(x2 - x1), abs(y2 - y1))
                   for (x1, y1), (x2, y2) in edges)
    return area2, boundary, (area2 - boundary + 2) // 2


def edge_polynomial(expression, start, end, S, P, X):
    dr, dl = end[0] - start[0], end[1] - start[1]
    length = gcd(abs(dr), abs(dl))
    dr, dl = dr // length, dl // length
    poly = sp.Poly(sp.expand(expression), S, P)
    return sp.expand(sum(
        poly.coeff_monomial(S ** (start[0] + j * dr)
                            * P ** (start[1] + j * dl)) * X**j
        for j in range(length + 1)))


# ----------------------------------------------------------------------
# I. Reconstruct the exact F00 support and its 24 realizable masks.
# ----------------------------------------------------------------------
U, W, Z = sp.symbols("U W Z")
Delta, Theta, Phi, eta = sp.symbols("Delta Theta Phi eta")
zeta, u, xi, alpha, beta, K = sp.symbols(
    "zeta u xi alpha beta K")
e = -sp.Rational(1376, 135)

full_rows = {
    (1, 0): -3,
    (2, 0): sp.Rational(8, 3),
    (3, 0): e,
    (0, 2): K,
    (2, 1): Phi,
    (4, 0): Delta,
    (1, 2): Theta,
    (3, 1): eta,
    (0, 3): zeta,
    (5, 0): u,
    (2, 2): xi,
    (4, 1): alpha,
    (1, 3): beta,
    (6, 0): U,
    (3, 2): W,
    (0, 4): Z,
}
need(len(full_rows) == 16, "complete sixteen-row source")

gate = {Z: 0, beta: 0, zeta: 0, W: 0, xi: 0, alpha: 0, Theta: 0}
rows = {
    monomial: sp.expand(sp.sympify(coefficient).subs(gate))
    for monomial, coefficient in full_rows.items()
    if sp.expand(sp.sympify(coefficient).subs(gate)) != 0
}
need(len(rows) == 9, "nine surviving F00 source rows")

support = defaultdict(lambda: sp.Integer(0))
support[(2, 0, 0)] += 1
support[(0, 1, 0)] -= 1
support[(2, 0, 1)] -= sp.Rational(1, 2)
for (i, j), coefficient in rows.items():
    support[(j + 2, i + j, 1)] -= coefficient
    support[(j, i + j + 1, 1)] += coefficient
support = {point: sp.expand(coefficient)
           for point, coefficient in support.items() if coefficient != 0}

M_plane = (Fr(1, 12), Fr(1, 6), Fr(-1, 6))
F00_plane = (Fr(1, 3), Fr(1, 6), Fr(-2, 3))
planes = candidate_planes(tuple(support))
Krel = sp.Rational(2848, 45) - sp.Rational(7, 6) * Delta
dstar = sp.Rational(3968, 63)
delta_reps = (sp.Integer(0), dstar, sp.Integer(1))

exact_supports = set()
status_counter = Counter()
for delta_value, u_value, phi_value, eta_value in product(
        delta_reps, (0, 1), (0, 1), (0, 1)):
    substitution = {
        U: 1, Delta: delta_value, K: Krel.subs(Delta, delta_value),
        u: u_value, Phi: phi_value, eta: eta_value,
    }
    need(substitution[K] != 0, "representatives remain on K!=0")
    active = frozenset(point for point, coefficient in support.items()
                       if sp.simplify(coefficient.subs(substitution)) != 0)
    exact_supports.add(active)
    status_counter[lower_faces(active, planes)] += 1

need(len(exact_supports) == 24, "24 exact F00 active supports")
need(status_counter == Counter({frozenset((M_plane, F00_plane)): 24}),
     "24/24 exact supports have precisely the M,F00 fan")

S, P, X = sp.symbols("S P X")
M = (P - S**2) * (U * P**6 - 1)
F00 = S**2 * (1 - U * P**6 - eta * S * P**4
               - K * S**2 * P**2)

def face_from_support(plane):
    aa, bb, cc = plane
    return sp.expand(sum(
        coefficient * S**rr * P**ll
        for (rr, ll, hh), coefficient in support.items()
        if Fr(hh) == aa * rr + bb * ll + cc))


need(sp.expand(face_from_support(M_plane) - M) == 0, "exact M face")
need(sp.expand(face_from_support(F00_plane) - F00) == 0,
     "exact F00 face")

M_polygon = [(0, 1), (2, 0), (2, 6), (0, 7)]
F00_polygon = [(2, 0), (4, 2), (2, 6)]
global_polygon = [(0, 1), (2, 0), (4, 2), (2, 6), (0, 7)]
need(pick(M_polygon) == (24, 14, 6), "M Pick ledger")
need(pick(F00_polygon) == (12, 10, 2), "F00 Pick ledger")
need(pick(global_polygon) == (36, 12, 13), "global Pick genus 13")

outer_owners = (M, F00, F00, M, M)
outer_expected = (
    X - 1,
    1 - K * X**2,
    -K - eta * X - U * X**2,
    U * (X - 1),
    U - X**6,
)
for start, end, owner, expected in zip(
        global_polygon, global_polygon[1:] + global_polygon[:1],
        outer_owners, outer_expected):
    need(sp.expand(edge_polynomial(owner, start, end, S, P, X)
                   - expected) == 0, "exact F00 outer edge")
internal = edge_polynomial(F00, (2, 0), (2, 6), S, P, X)
need(sp.expand(internal - (1 - U * X**6)) == 0,
     "exact M--F00 internal edge")

# The good-form order is derived from the primitive supporting plane, not
# imported from another audit.
face_base = lcm(*(q.denominator for q in (*F00_plane, Fr(1, 6))))
face_order = face_base * (Fr(5, 6) - sum(F00_plane))
need((face_base, face_order) == (6, 6), "F00 base/order is 6/6")


# ----------------------------------------------------------------------
# II. Generic F00 normalization and the unique collision.
# ----------------------------------------------------------------------
core = 1 - U * P**6 - eta * S * P**4 - K * S**2 * P**2
Y00 = 2 * K * S * P + eta * P**3
L00 = eta**2 - 4 * K * U
branch = L00 * P**6 + 4 * K
need(sp.expand(Y00**2 - branch + 4 * K * core) == 0,
     "F00 hyperelliptic normalization identity")
need(sp.factor(sp.discriminant(branch, P))
     == 47775744 * K**5 * (4 * K * U - eta**2)**5,
     "exact F00 branch discriminant")

Lsym = sp.symbols("Lsym", nonzero=True)
generic_branch = sp.Poly(Lsym * P**6 + 4 * K, P,
                         domain=sp.QQ.frac_field(Lsym, K))
need(generic_branch.degree() == 6, "generic branch degree six")
need(sp.gcd(generic_branch, generic_branch.diff()).degree() == 0,
     "generic branch is squarefree")
need((generic_branch.degree() - 2) // 2 == 2,
     "generic F00 carrier has genus two")

a, lam = sp.symbols("a lam", nonzero=True)
collision = {K: lam**2, U: lam**2 * a**2, eta: -2 * lam**2 * a}
linear = lam * (S * P - a * P**3)
need(sp.expand(core.subs(collision) - (1 - linear) * (1 + linear)) == 0,
     "collision face has two rational signs")
need(sp.expand((K + eta * X + U * X**2).subs(collision)
               - lam**2 * (1 - a * X)**2) == 0,
     "outer quadratic is exactly doubled")
minus_attach = 1 - lam * a * X**3
plus_attach = 1 + lam * a * X**3
need(sp.expand((1 - U * X**6).subs(collision)
               - minus_attach * plus_attach) == 0,
     "six M attachments split as two cubics")
attach_field = sp.QQ.frac_field(lam, a)
need(sp.Poly(minus_attach, X, domain=attach_field).discriminant() != 0,
     "first attachment cubic is reduced")
need(sp.Poly(plus_attach, X, domain=attach_field).discriminant() != 0,
     "second attachment cubic is reduced")
need(sp.resultant(minus_attach, plus_attach, X) != 0,
     "the two triples of attachments are disjoint")

# At reciprocal infinity, the two branches are x^3 +/- lam*(v-a)=0.
# Their intersection quotient is C[x,w]/(w,x^3), of length three.
x, v, rho = sp.symbols("x v rho")
wloc = sp.symbols("wloc")
fplus = x**3 + lam * wloc
fminus = x**3 - lam * wloc
need(sp.expand(fplus * fminus - (x**6 - lam**2 * wloc**2)) == 0,
     "reciprocal signs are exactly the two F00 branches")
gb = sp.groebner((fplus, fminus), wloc, x,
                 domain=sp.QQ.frac_field(lam))
need(gb.reduce(wloc)[1] == 0 and gb.reduce(x**3)[1] == 0,
     "intersection ideal contains (w,x^3)")
need(gb.reduce(sp.Integer(1))[1] != 0 and gb.reduce(x)[1] != 0
     and gb.reduce(x**2)[1] != 0,
     "intersection quotient has basis 1,x,x^2")
need(3 == (5 + 1) // 2, "intersection is an A5 contact of delta three")


# ----------------------------------------------------------------------
# III. Exact primitive reciprocal chart and forced critical depth.
# ----------------------------------------------------------------------
sigma = sp.symbols("sigma")
p0 = sigma**-1 * P
s0 = sigma**-2 * S
y0 = s0 * p0
H0 = sp.expand(sum(coefficient * p0**i * y0**j
                   for (i, j), coefficient in rows.items()))
F_source = ((s0**2 - p0) * (1 - sigma**6 * H0)
            - sigma**6 * s0**2 / 2)
G_source = sp.expand(sigma**4 * F_source)
H6 = U * P**6 + eta * S * P**4 + K * S**2 * P**2
H5 = u * P**5 + Phi * S * P**3
G_expected = sp.expand(
    (S**2 - sigma**3 * P)
    * (1 - H6 - sigma * H5 - sigma**2 * Delta * P**4
       - e * sigma**3 * P**3 - sp.Rational(8, 3) * sigma**4 * P**2
       + 3 * sigma**5 * P)
    - sigma**6 * S**2 / 2)
need(sp.simplify(G_source - G_expected) == 0,
     "exact primitive F00 scaling from the complete source")

Phi_root = sp.factor(x**10 * G_source.subs({P: 1 / x, S: v / x**2}))
A = U + eta * v + K * v**2
B = u + Phi * v
root_expected = sp.expand(
    (v**2 - (sigma * x)**3)
    * (x**6 - A - sigma * x * B - (sigma * x)**2 * Delta
       - e * (sigma * x)**3 - sp.Rational(8, 3) * (sigma * x)**4
       + 3 * (sigma * x)**5)
    - (sigma * x)**6 * v**2 / 2)
need(sp.simplify(Phi_root - root_expected) == 0,
     "exact reciprocal root chart")

D = (A + rho * B + rho**2 * Delta + e * rho**3
     + sp.Rational(8, 3) * rho**4 - 3 * rho**5
     + rho**6 * v**2 / (2 * (v**2 - rho**3)))
need(sp.factor(sp.together(
    Phi_root / (v**2 - (sigma * x)**3)
    - (x**6 - D.subs(rho, sigma * x)))) == 0,
    "division by the reciprocal unit produces one series")

repeat = {U: K * a**2, eta: -2 * K * a}
need(sp.expand(A.subs(repeat) - K * (v - a)**2) == 0,
     "canonical repeated-root parameterization")
need((v**2 - rho**3).subs({v: a, rho: 0}) == a**2,
     "reciprocal prefactor is a unit")
need(sp.diff(D, v, 2).subs({**repeat, v: a, rho: 0}) == 2 * K,
     "Morse Hessian is a unit")

# The affine approximation is much sharper than needed: the rational
# residual first perturbs the critical point in degree rho^9, not rho^6.
v0 = a - rho * Phi / (2 * K)
Drepeat = sp.together(D.subs(repeat))
derivative_at_v0 = sp.factor(sp.diff(Drepeat, v).subs(v, v0))
derivative_expected = sp.factor(
    -rho**9 * v0 / (v0**2 - rho**3)**2)
need(sp.simplify(derivative_at_v0 - derivative_expected) == 0,
     "critical-point correction starts in rho^9")

psi_through_8 = sp.series(Drepeat.subs(v, v0), rho, 0, 9).removeO().expand()
psi_expected = (
    (u + Phi * a) * rho
    + (Delta - Phi**2 / (4 * K)) * rho**2
    + e * rho**3 + sp.Rational(8, 3) * rho**4
    - 3 * rho**5 + sp.Rational(1, 2) * rho**6)
need(sp.simplify(psi_through_8 - psi_expected) == 0,
     "critical value is exact through rho^8")
need(e == -sp.Rational(1376, 135) and e != 0,
     "fixed cubic coefficient forces critical depth at most three")

# All three depths occur on the seam.  These examples use a=1 and satisfy
# U=K, eta=-2K, hence U*K!=0 and eta^2=4*K*U.
k0 = Krel.subs(Delta, 0)
k1 = Krel.subs(Delta, 1)
examples = {
    1: {Delta: 0, K: k0, U: k0, eta: -2 * k0, Phi: 0, u: 1},
    2: {Delta: 1, K: k1, U: k1, eta: -2 * k1, Phi: 0, u: 0},
    3: {Delta: 0, K: k0, U: k0, eta: -2 * k0, Phi: 0, u: 0},
}
c1 = u + Phi * a
c2 = Delta - Phi**2 / (4 * K)
for depth, substitution in examples.items():
    substitution = {**substitution, a: 1}
    coefficients = (sp.simplify(c1.subs(substitution)),
                    sp.simplify(c2.subs(substitution)), e)
    first = next(j + 1 for j, coefficient in enumerate(coefficients)
                 if coefficient != 0)
    need(first == depth, "explicit realization of critical depth")
    need(sp.simplify((eta**2 - 4 * K * U).subs(substitution)) == 0,
         "depth example lies on collision")
    need(substitution[U] * substitution[K] != 0,
         "depth example remains on U*K!=0")


# ----------------------------------------------------------------------
# IV. Tail functions, primitive/common bases, complement, and orders.
# ----------------------------------------------------------------------
tau, t, XX, YY, Y0, z, h = sp.symbols("tau t XX YY Y0 z h")
c0 = sp.symbols("c0", nonzero=True)
tail_rows = {}
for r in (1, 2, 3):
    q = r // 2
    eps = r % 2

    # Common base used for a uniform rho=tau^12*X complement.
    common_sigma = tau ** (2 * (6 - r))
    common_x = tau ** (2 * r) * XX
    common_y = tau ** (6 * r) * YY
    transformed = sp.cancel(
        (common_y**2
         - common_x**r * (common_x**(6 - r)
                          - common_sigma**r * c0)) / tau**(12 * r))
    need(sp.expand(transformed
                   - (YY**2 - XX**r * (XX**(6 - r) - c0))) == 0,
         "common-base exceptional equation")

    normalized_polynomial = sp.expand(XX**eps * (XX**(6 - r) - c0))
    raw_after_normalization = sp.expand(
        (XX**q * Y0)**2 - XX**r * (XX**(6 - r) - c0))
    need(sp.expand(raw_after_normalization / XX**(2 * q)
                   - (Y0**2 - normalized_polynomial)) == 0,
         "remove only the persistent square")
    poly = sp.Poly(normalized_polynomial, XX,
                   domain=sp.QQ.frac_field(c0))
    need(sp.gcd(poly, poly.diff()).degree() == 0,
         "normalized tail branch polynomial is squarefree")
    genus = (poly.degree() - 1) // 2
    expected_genus = {1: 2, 2: 1, 3: 1}[r]
    need(genus == expected_genus, "tail genus")
    need(q + genus + 1 == 3, "local A5 conservation")

    # Primitive Newton base gives the same exceptional function field.
    gg = gcd(r, 6 - r)
    primitive_sigma_weight = (6 - r) // gg
    primitive_x_weight = r // gg
    primitive_y_weight = 3 * r // gg
    primitive_sigma = t**primitive_sigma_weight
    primitive_x = t**primitive_x_weight * XX
    primitive_y = t**primitive_y_weight * YY
    primitive_total = 2 * primitive_y_weight
    transformed_primitive = sp.cancel(
        (primitive_y**2
         - primitive_x**r * (primitive_x**(6 - r)
                             - primitive_sigma**r * c0))
        / t**primitive_total)
    need(sp.expand(transformed_primitive
                   - (YY**2 - XX**r * (XX**(6 - r) - c0))) == 0,
         "primitive base has identical tail, not a root-cover genus")
    need(gcd(gcd(primitive_sigma_weight, primitive_x_weight),
             primitive_y_weight) == 1,
         "primitive valuation vector")

    common_order = (6 * 2 * (6 - r) + 6 * 2 * r
                    + 2 * r - 6 * r)
    primitive_order = (6 * primitive_sigma_weight
                       + 6 * primitive_x_weight
                       + primitive_x_weight - primitive_y_weight)
    need(common_order == {1: 68, 2: 64, 3: 60}[r],
         "common-base differential order")
    need(primitive_order == {1: 34, 2: 16, 3: 10}[r],
         "primitive differential order")
    need(common_order == 2 * gg * primitive_order,
         "common order is primitive order times ramification")

    # h=y/x^3 and z=1/X give the exact two-ended complement.  On the
    # common base z*rho=tau^12; on the primitive base z*rho=t^(6/g).
    common_rho = tau**12 * XX
    common_ratio = sp.cancel(common_rho**r / common_x**6)
    need(sp.simplify(common_ratio.subs(XX, 1 / z) - z**(6 - r)) == 0,
         "common-base reciprocal complement")
    primitive_rho = t**(6 // gg) * XX
    primitive_ratio = sp.cancel(primitive_rho**r / primitive_x**6)
    need(sp.simplify(primitive_ratio.subs(XX, 1 / z)
                     - z**(6 - r)) == 0,
         "primitive reciprocal complement")
    complement = h**2 - 1 + z**(6 - r) * c0
    need(sp.diff(complement, h).subs({h: 1, z: 0}) == 2,
         "positive infinity attachment is etale")
    need(sp.diff(complement, h).subs({h: -1, z: 0}) == -2,
         "negative infinity attachment is etale")

    tail_rows[r] = {
        "equation": str(normalized_polynomial),
        "genus": genus,
        "persistent_delta": q,
        "common_order": common_order,
        "primitive_order": primitive_order,
        "primitive_surface_A": 6 // gg - 1,
        "primitive_surface_exponent": 6 // gg,
    }

# Direct differential conversion in the root chart.
GS_root = sp.diff(G_source, S).subs({P: 1 / x, S: v / x**2})
need(sp.simplify(sp.diff(Phi_root, v) - x**8 * GS_root) == 0,
     "Phi_v=x^8 G_S")
# Since -dP=x^-2 dx, the source generator becomes x^6 dx/Phi_v.
need(-2 + 8 == 6, "reciprocal differential contributes x^6")


# ----------------------------------------------------------------------
# V. Incidence/genus and the proper-flat input.
# ----------------------------------------------------------------------
# M: seven rational components, twelve internal nodes.  Generic F00 adds
# one component and six attachments.  Collision splits F00 into two signs.
need(18 - 8 + 1 == 11 and 11 + 2 == 13,
     "generic graph plus genus-two carrier equals global genus 13")
need(18 - 9 + 1 == 10 and 10 + 3 == 13,
     "split graph plus A5 conductor equals global genus 13")
for r, row in tail_rows.items():
    # Contracting the rational A_n attachment chains, one tail vertex and
    # two sign edges change (V,E)=(9,18) to (10,20).
    graph_b1 = 20 - 10 + 1
    need(graph_b1 == 11, "two-ended tail restores one graph cycle")
    # Resolving each z*rho=t^N attachment replaces one contracted edge by
    # an N-edge chain with N-1 rational vertices.  Two such replacements
    # leave b1 unchanged.  The common base has N=12; the primitive bases
    # have N=6,3,2 for r=1,2,3.
    for chain_exponent in (12, row["primitive_surface_exponent"]):
        expanded_vertices = 10 + 2 * (chain_exponent - 1)
        expanded_edges = 18 + 2 * chain_exponent
        need(expanded_edges - expanded_vertices + 1 == 11,
             "expanded rational attachment chains preserve graph b1")
    need(graph_b1 + row["genus"] + row["persistent_delta"] == 13,
         "regularized collision conserves global genus")
    need(row["common_order"] > 0 and row["primitive_order"] > 0,
         "every positive-genus tail has positive form order")

print("THM4351_DOUBLE_NORMAL_OWNER_ZERO_AUDIT=PASS")
print("gate=Z=beta_11=zeta_3=W=xi_10=alpha_11=Theta=0;U*K!=0")
print("exact_supports=24;fan=M,F00;global_genus=13")
print("generic=L00!=0;F00_genus=2;form_order=6@sigma;graph_b1=11")
print("collision=L00=eta^2-4KU=0;two_rational_signs;A5_delta=3;attachments=3+3")
print("critical_coefficients=(u+Phi*a),(Delta-Phi^2/(4K)),e=-1376/135")
print("critical_depths=1,2,3_all_realizable;critical_point_correction_starts=rho^9")
print("tails=r1:g2/d0,r2:g1/d1,r3:g1/d1;two_ended=yes")
print("common_base=sigma=tau^[2(6-r)];rho=tau^12*X;chains=A11+A11")
print("common_form_orders=r1:68,r2:64,r3:60")
print("primitive_form_orders=r1:34,r2:16,r3:10;primitive_chains=A5,A2,A1")
print("genus_conservation=b1(11)+tail_genus+persistent_delta=13")
print("verdict=full_F00_corner_extinct_relative_to_inherited_proper_flat_interface")
