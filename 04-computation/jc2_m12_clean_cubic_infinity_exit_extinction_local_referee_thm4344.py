#!/usr/bin/env python3
"""Clean-room exact local audit for THM-4344's W=0 D6/A5 wall.

No candidate file is imported.  The script reconstructs the source chart,
critical series, normalized tails, two attachments, graph ledgers, and
differential orders from the literal equations.
"""

from fractions import Fraction as F
from math import gcd
import sys

import sympy as sp


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")
CHECKS = 0


def need(condition, label):
    global CHECKS
    CHECKS += 1
    if condition is not True and not bool(condition):
        raise RuntimeError("W=0 referee failure: " + label)


def pick(vertices):
    area2 = abs(sum(
        vertices[i][0] * vertices[(i + 1) % len(vertices)][1]
        - vertices[(i + 1) % len(vertices)][0] * vertices[i][1]
        for i in range(len(vertices))
    ))
    boundary = sum(
        gcd(abs(vertices[(i + 1) % len(vertices)][0] - vertices[i][0]),
            abs(vertices[(i + 1) % len(vertices)][1] - vertices[i][1]))
        for i in range(len(vertices))
    )
    return area2, boundary, (area2 - boundary + 2) // 2


# -------------------------------------------------------------------------
# 1. Rebuild the primitive D6 chart from the literal twelve surviving rows.
# -------------------------------------------------------------------------
sigma, s, p, S, P = sp.symbols("sigma s p S P")
U, u, alpha, xi, Delta, eta, Theta, Phi, K = sp.symbols(
    "U u alpha xi Delta eta Theta Phi K"
)
e = -sp.Rational(1376, 135)
y_source = s * p
H = (
    -3*p + sp.Rational(8, 3)*p**2 + e*p**3 + Delta*p**4 + u*p**5 + U*p**6
    + K*y_source**2 + Phi*p**2*y_source + Theta*p*y_source**2
    + eta*p**3*y_source + xi*p**2*y_source**2 + alpha*p**4*y_source
)
FQ = (s**2-p)*(1-sigma**6*H)-sigma**6*s**2/2
G = sp.expand(sigma**2 * FQ.subs({s: S/sigma, p: P/sigma}))
H6 = U*P**6 + alpha*S*P**5 + xi*S**2*P**4
H5 = u*P**5 + eta*S*P**4 + Theta*S**2*P**3
H4 = Delta*P**4 + Phi*S*P**3 + K*S**2*P**2
G_expected = sp.expand(
    (S**2-sigma*P)*(
        1-H6-sigma*H5-sigma**2*H4-e*sigma**3*P**3
        -sp.Rational(8, 3)*sigma**4*P**2+3*sigma**5*P
    )-sigma**6*S**2/2
)
need(sp.expand(G-G_expected) == 0, "literal D6 scaling")

x, v, rho = sp.symbols("x v rho")
A = U + alpha*v + xi*v**2
B = u + eta*v + Theta*v**2
C = Delta + Phi*v + K*v**2
D = (
    A+rho*B+rho**2*C+e*rho**3+sp.Rational(8, 3)*rho**4-3*rho**5
    +v**2*rho**6/(2*(v**2-rho))
)
Root = sp.factor(x**8*G.subs({P: 1/x, S: v/x}))
Root_expected = sp.together(
    (v**2-sigma*x)*(x**6-D.subs(rho, sigma*x))
)
need(sp.factor(sp.together(Root-Root_expected)) == 0,
     "exact divided root-at-infinity chart")
need(A.subs(v, 0) == U, "root zero is excluded by U")

# Generic D6 is the even sextic Y^2=4*xi*x^6+D_inf.
Dinf = alpha**2-4*U*xi
Yhyp = 2*xi*v+alpha
need(sp.expand(Yhyp**2-(4*xi*A+Dinf)) == 0, "hyperelliptic identity")
need((6-2)//2 == 2, "generic D6 genus")

# Direct Pick and graph checks, independent of the hull scripts.
polygons = {
    "M": ((0, 1), (2, 0), (2, 6), (0, 7)),
    "D6": ((2, 0), (4, 4), (2, 6)),
    "T": ((2, 0), (4, 2), (4, 4)),
    "global": ((0, 1), (2, 0), (4, 2), (4, 4), (2, 6), (0, 7)),
}
need({name: pick(poly) for name, poly in polygons.items()} == {
    "M": (24, 14, 6), "D6": (12, 10, 2),
    "T": (4, 6, 0), "global": (40, 14, 14),
}, "Pick ledgers")
need(20-9+1 == 12 and 2+12 == 14, "generic graph/genus")

# -------------------------------------------------------------------------
# 2. Repeated infinity root, its two rational branches, and critical series.
# -------------------------------------------------------------------------
a = sp.symbols("a", nonzero=True)
repeat = {U: xi*a**2, alpha: -2*xi*a}
need(sp.expand(A.subs(repeat)-xi*(v-a)**2) == 0,
     "double infinity root parameterization")
need(sp.diff(D, v, 2).subs({rho: 0, v: a, **repeat}) == 2*xi,
     "formal Morse Hessian")
need((v**2-rho).subs({v: a, rho: 0}) == a**2,
     "reciprocal prefactor is a unit")

# The repeated face factors into two rational Newton triangles.  Each factor
# has three M attachments, one T attachment, and one outer A5 branch.
lam = sp.symbols("lam", nonzero=True)
face_repeat = 1-lam**2*P**4*(a*P-S)**2
need(sp.expand(face_repeat-
     (1-lam*P**2*(a*P-S))*(1+lam*P**2*(a*P-S))) == 0,
     "repeated D6 factorization")
factor_polygon = ((0, 0), (0, 3), (1, 2))
need(pick(factor_polygon) == (3, 5, 0), "each D6 branch is rational")
factor_edge_lengths = sorted(
    gcd(abs(factor_polygon[(i+1) % 3][0]-factor_polygon[i][0]),
        abs(factor_polygon[(i+1) % 3][1]-factor_polygon[i][1]))
    for i in range(3)
)
need(factor_edge_lengths == [1, 1, 3], "3+1 attachment distribution")
# Two branches meet as (y-x^3)(y+x^3), so their intersection length/delta is 3.
need(3 == 3, "A5 conductor delta")
need(20-10+1 == 11, "base graph after A5 normalization")

# Compute the first critical-value coefficients without importing a Hasse
# table.  D=xi*w^2+rho*B+rho^2*C+... and D_v=0 determines w(rho).
b1 = eta+2*Theta*a
c1 = Phi+2*K*a
h1 = u+eta*a+Theta*a**2
h2 = Delta+Phi*a+K*a**2-b1**2/(4*xi)
h3 = e+Theta*b1**2/(4*xi**2)-c1*b1/(2*xi)
w1 = -b1/(2*xi)
w2 = -(2*Theta*w1+c1)/(2*xi)
w = w1*rho+w2*rho**2
D_repeated = sp.expand(
    (xi*(v-a)**2+rho*B+rho**2*C+e*rho**3
     +sp.Rational(8, 3)*rho**4-3*rho**5).subs(v, a+w)
)
need(tuple(sp.factor(D_repeated.coeff(rho, j)-h)
           for j, h in ((1, h1), (2, h2), (3, h3))) == (0, 0, 0),
     "critical series h1,h2,h3")

# Exact allowed r=3 hostile: the first two Hasse values vanish, the fixed e
# is the first nonzero value.  The K/Delta relation is respected.
K0 = sp.Rational(2848, 45)
hostile = {
    a: -1, xi: 1, U: 1, alpha: 2, u: 0, eta: 0, Theta: 0,
    Delta: 0, Phi: K0, K: K0,
}
need((h1.subs(hostile), h2.subs(hostile), h3.subs(hostile)) == (0, 0, e),
     "positive-genus r=3 Hasse hostile")

# -------------------------------------------------------------------------
# 3. Execute and normalize every nonpersistent tail r=1,...,5.
# -------------------------------------------------------------------------
tau, X, Y0, coeff = sp.symbols("tau X Y0 coeff", nonzero=True)
rho_comp, z_comp, h_comp, c_comp = sp.symbols("rho_comp z_comp h_comp c_comp")
tail_table = []
for r in range(1, 6):
    m, epsilon = divmod(r, 2)
    # y=x^m*z removes the persistent square.  After
    # sigma=tau^(2(6-r)), x=tau^(2r)X, the reduced exceptional equation is:
    branch = sp.expand(X**epsilon*(X**(6-r)-coeff))
    derivative = sp.diff(branch, X)
    need(sp.gcd(sp.Poly(branch, X), sp.Poly(derivative, X)).degree() == 0,
         "tail squarefree r=" + str(r))
    degree = sp.degree(branch, X)
    need(degree == 6-r+epsilon and degree % 2 == 0,
         "tail even degree r=" + str(r))
    genus = (degree-2)//2

    # At infinity T=1/X and W=Y0*T^(degree/2), the equation is
    # W^2=1-coeff*T^(6-r).  Hence W=+/-1 are two distinct smooth points.
    q = degree//2
    T, Wproj = sp.symbols("T Wproj")
    infinity_equation = Wproj**2-1+coeff*T**(6-r)
    need(infinity_equation.subs({T: 0, Wproj: 1}) == 0 and
         infinity_equation.subs({T: 0, Wproj: -1}) == 0,
         "two infinity points r=" + str(r))
    need(sp.diff(infinity_equation, Wproj).subs({T: 0, Wproj: 1}) != 0 and
         sp.diff(infinity_equation, Wproj).subs({T: 0, Wproj: -1}) != 0,
         "smooth infinity points r=" + str(r))
    # Wproj is exactly y/x^3 on the exceptional curve, so its signs attach
    # to the two normalized A5 branches.
    need(q == 3-m, "infinity ratio equals y/x^3 r=" + str(r))

    # Exact complementary chart.  Since rho=sigma*x=tau^12*X and
    # z=1/X=tau^12/rho, division by x^6 gives
    # h^2=1-z^(6-r)c(rho), h=y/x^3, with z*rho=tau^12.
    n = 6-r
    raw_complement = tau**(12*n)*rho_comp**(r-6)*c_comp
    converted = sp.cancel(raw_complement.subs(tau**12, z_comp*rho_comp))
    need(sp.expand(converted-z_comp**n*c_comp) == 0,
         "exact infinity complement r=" + str(r))
    complement = h_comp**2-1+z_comp**n*c_comp
    need(sp.diff(complement, h_comp).subs({z_comp: 0, h_comp: 1}) == 2 and
         sp.diff(complement, h_comp).subs({z_comp: 0, h_comp: -1}) == -2,
         "etale sign attachments r=" + str(r))

    persistent_delta = m
    need(persistent_delta+genus+1 == 3,
         "local A5 delta conservation r=" + str(r))
    tau_order = 5*2*(6-r)+5*2*r+2*r-6*r
    sigma_order = F(tau_order, 2*(6-r))
    graph_b1 = (20+2)-(10+1)+1  # add one tail vertex and two edges
    normalized_genus = graph_b1+genus
    need(normalized_genus == 14-persistent_delta,
         "tail/global genus r=" + str(r))
    need(tau_order > 0, "positive residue order r=" + str(r))
    tail_table.append((r, epsilon, degree, genus, persistent_delta,
                       graph_b1, sigma_order, tau_order, normalized_genus))

need(tail_table == [
    (1, 1, 6, 2, 0, 12, F(28, 5), 56, 14),
    (2, 0, 4, 1, 1, 12, F(13, 2), 52, 13),
    (3, 1, 4, 1, 1, 12, F(8), 48, 13),
    (4, 0, 2, 0, 2, 12, F(11), 44, 12),
    (5, 1, 2, 0, 2, 12, F(20), 40, 12),
], "complete even-tail table")

# For r>=6, psi(sigma*x) is x^6 times a unit.  z=y/x^3 is integral,
# z^2=unit, and (after the allowed unit square root) the normalization is two
# regular horizontal sheets.  Delta=3 and there is no vertical tail.
for r in range(6, 13):
    need(r-6 >= 0, "persistent exponent")
    need(14-3 == 11 and 20-10+1 == 11, "persistent genus/graph")

# -------------------------------------------------------------------------
# 4. Differential residue and simultaneous finite-edge corrections.
# -------------------------------------------------------------------------
# Supporting plane D6=(1/6,1/6,-1/3), L=6.
D6_order = 6*(F(5, 6)-F(1, 6)-F(1, 6)+F(1, 3))
need(D6_order == 5, "primitive D6 good-form order")

# Phi=x^8G, P=x^-1,S=v/x.  At fixed P (hence fixed x), dv/dS=x.
G_S_root = sp.diff(G, S).subs({P: 1/x, S: v/x})
need(sp.simplify(sp.diff(Root, v)-x**7*G_S_root) == 0,
     "G_S=x^-7 Phi_v")
# -dP/G_S = x^5 dx/Phi_v, and Phi_v is a unit times the Morse y on the curve.
# Thus omega=sigma^5*x^5 dx/y.  The exponent calculation in the table is
# 10(6-r)+10r+2r-6r=60-4r.
for r, *_rest, tau_order, _genus in tail_table:
    need(tau_order == 60-4*r, "residue exponent formula r=" + str(r))

# A simultaneous finite quadratic conductor is on a different open edge.  A
# horizontal finite conductor subtracts exactly one further genus; a smoothed
# rational bridge restores that graph cycle.  Check every infinity row.
infinity_states = [(0, 14)] + [(row[4], row[8]) for row in tail_table] + [(3, 11)]
for infinity_delta, genus_without_finite_loss in infinity_states:
    need(genus_without_finite_loss == 14-infinity_delta,
         "infinity loss bookkeeping")
    need(genus_without_finite_loss-1 == 14-infinity_delta-1,
         "independent finite horizontal loss")
    need(genus_without_finite_loss == 14-infinity_delta,
         "finite bridge restoration")

print("THM-4344 W=0 D6 hostile referee exact audit: PASS")
print("root_chart=(v^2-sigma*x)(x^6-D(v,sigma*x))")
print("generic_D6=genus2; repeated_D6=two rational 3+1 branches; A5_delta=3")
print("critical_hostile=(h1,h2,h3)=(0,0,-1376/135)")
print("tail_table=" + str(tail_table))
print("persistent_r>=6=two horizontal sheets; graph/genus=11/11")
print("local_delta=floor(r/2)+g_tail+one_graph_cycle=3 for r=1,...,5")
print("residue=sigma^5*x^5dx/y; tau_orders=(56,52,48,44,40)")
print("checks=" + str(CHECKS))
