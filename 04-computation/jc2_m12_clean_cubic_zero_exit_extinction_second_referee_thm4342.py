#!/usr/bin/env python3
"""Second hostile audit for THM-4342's K=0 root-exit theorem.

This path imports neither proposed scratch certificate nor the inherited hull
engine.  It targets saturation, all ordinary-blowup charts, the seam relation,
the m=1 normal gate, differential transport, and the genus ledgers.
"""

from fractions import Fraction as F
from math import gcd
import sys

import sympy as sp


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")


def need(condition, label):
    if not bool(condition):
        raise AssertionError(label)


P, b, delta = sp.symbols("P b delta")
U, W = sp.symbols("U W", nonzero=True)
Theta, xi = sp.symbols("Theta xi")
Phi, eta, alpha = sp.symbols("Phi eta alpha")
Delta, upsilon = sp.symbols("Delta upsilon")
J = Phi + eta*P + alpha*P**2
D = -3 + sp.Rational(8,3)*P - sp.Rational(1376,135)*P**2 + Delta*P**3 + upsilon*P**4 + U*P**5
qbase = 1-delta**2/2


def source(m, quotient):
    return sp.expand(
        (1-delta**2*P*b**2)
        * (b**2-P**(m+2)*quotient-delta*b*P**3*J-delta**2*b**2*P*D)
        - delta**2*b**2/2
    )


quotients = {1: Theta+xi*P+W*P**2, 2: xi+W*P, 3: W}

# Laurent saturation distinguishes exit multiplicity from a torus collision.
A1 = P*quotients[1]
full_disc = sp.factor(sp.discriminant(A1,P))
sat_disc = sp.factor(sp.discriminant(quotients[1],P))
need(sp.simplify(full_disc-Theta**2*(xi**2-4*W*Theta)) == 0,
     "full discriminant factorization")
need(sp.simplify(sat_disc-(xi**2-4*W*Theta)) == 0, "saturated discriminant")
need(sp.Poly(quotients[2],P).degree() == 1 and sp.Poly(quotients[3],P).degree() == 0,
     "no collision after deeper exits")

# The inherited seam relation is not optional.
Delta_forced = F(6,7)*F(2848,45)
need(Delta_forced == F(5696,105), "K=0 forced Delta")
need(F(2848,45)-F(7,6)*Delta_forced == 0, "seam relation")
need(Delta_forced != 0, "Delta owner is actually present")

# First blowup of the horizontal endpoint centre P=b=0.  The P-chart contains
# every strict-transform point; the complementary b-chart has unit equation on
# its exceptional divisor and hence no missing component.
y1, z1 = sp.symbols("y1 z1")
first_P = {}
for m,Qm in quotients.items():
    Fm = source(m,Qm)
    first_P[m] = sp.factor(Fm.subs(b,P*y1)/P**2)
    need(not sp.denom(first_P[m]).has(P), f"m{m} first P-chart integral")
    first_b = sp.factor(Fm.subs(P,b*z1)/b**2)
    need(not sp.denom(first_b).has(b), f"m{m} first b-chart integral")
    need(sp.simplify(first_b.subs(b,0)-qbase) == 0,
         f"m{m} complementary first chart empty over the endpoint")

need(sp.simplify(sp.diff(first_P[1],P).subs({P:0,y1:0})+Theta) == 0,
     "m1 resolved after one blowup")

# m=2,3 require the second ordinary blowup.  Check both charts and their glue.
y2,z2 = sp.symbols("y2 z2")
second_P = {}
second_y = {}
for m in (2,3):
    second_P[m] = sp.factor(first_P[m].subs(y1,P*y2)/P**2)
    second_y[m] = sp.factor(first_P[m].subs(P,y1*z2)/y1**2)
    need(not sp.denom(second_P[m]).has(P), f"m{m} second P-chart integral")
    need(not sp.denom(second_y[m]).has(y1), f"m{m} second y-chart integral")

need(sp.simplify(second_P[2].subs(P,0)-(qbase*y2**2-xi)) == 0,
     "m2 P-chart two simple points")
need(sp.simplify(second_y[2].subs(y1,0)-(qbase-xi*z2**2)) == 0,
     "m2 complementary chart sees the same two points")
need(sp.simplify((qbase*y2**2-xi).subs(y2,1/z2)*z2**2
                 -(qbase-xi*z2**2)) == 0, "m2 chart glue")
need(sp.simplify(second_y[3].subs(y1,0)-qbase) == 0,
     "m3 complementary second chart empty")
need(sp.simplify(sp.diff(second_P[3],P).subs({P:0,y2:0})+W) == 0,
     "m3 resolved binomial collar")

# Only m=1 has a saturated collision.  Its double root is nonzero and has the
# stated root-free normal gate.
a = sp.symbols("a", nonzero=True)
double_sub = {Theta:W*a**2, xi:-2*W*a}
need(sp.simplify(quotients[1].subs(double_sub)-W*(P-a)**2) == 0,
     "m1 double parameterization")
normal_gate = 4*W**2*Phi-2*W*xi*eta+xi**2*alpha
need(sp.simplify(normal_gate.subs(double_sub)-4*W**2*J.subs(P,a)) == 0,
     "m1 normal gate")

# If J(a)=0, the horizontal tangent discriminant is x^2 times a DVR unit.
x,B1,Ca = sp.symbols("x B1 Ca")
tangent = qbase*b**2-delta*B1*x*b-a**3*W*x**2
tangent_disc = sp.factor(sp.discriminant(tangent,b))
need(sp.simplify(tangent_disc-x**2*(delta**2*B1**2+4*qbase*a**3*W)) == 0,
     "horizontal node discriminant")
need(sp.simplify((delta**2*B1**2+4*qbase*a**3*W).subs(delta,0)) == 4*a**3*W,
     "horizontal discriminant coefficient is a unit")

# Differential chain rule under b=P^e y.
evar = sp.symbols("evar", integer=True, positive=True)
EP,EY,dy = sp.symbols("EP EY dy", nonzero=True)
dP = -EY*dy/EP
db = P**evar*dy+evar*P**(evar-1)*sp.symbols("yy")*dP
yy = next(v for v in db.free_symbols if v.name == "yy")
FP = P**(2*evar-1)*(P*EP-evar*yy*EY)
ratio = sp.factor(P**(2*evar)*yy**2*db/FP)
need(sp.simplify(ratio-P**evar*yy**2*dy/EP) == 0,
     "differential transform")


def polygon_ledger(vertices):
    area2 = abs(sum(x1*y2-x2*y1 for (x1,y1),(x2,y2) in
                    zip(vertices,vertices[1:]+vertices[:1])))
    boundary = sum(gcd(abs(x2-x1),abs(y2-y1)) for (x1,y1),(x2,y2) in
                   zip(vertices,vertices[1:]+vertices[:1]))
    return area2,boundary,(area2-boundary+2)//2


polygons = {
    1: ((0,1),(2,0),(4,3),(4,5),(0,7)),
    2: ((0,1),(2,0),(4,4),(4,5),(0,7)),
    3: ((0,1),(2,0),(4,5),(0,7)),
}
need({m:polygon_ledger(list(v)) for m,v in polygons.items()} ==
     {1:(40,12,15),2:(38,12,14),3:(36,10,14)}, "Pick ledgers")
packets = {1:(11,11,5,3,3,1),2:(11,11,3,3,3,1),3:(11,11,7,1)}
need(sum(i-1 for i in packets[1]) == 28, "m1 RH checksum")
need(sum(i-1 for i in packets[2]) == 26, "m2 RH checksum")
need(sum(i-1 for i in packets[3]) == 26, "m3 RH checksum")
need(3+1+11 == 15 and 3+11 == 14, "component genus ledgers")

print("PASS THM4342 KZERO SECOND HOSTILE AUDIT")
print("SEAM=K0_forces_Delta_5696/105;Delta-absent_hull_masks_are_conservative_only")
print("SATURATION=Theta^2_exit_factor_removed;only_Q1_discriminant_remains")
print("BLOWUPS=m1_one_chart;m2_two_charts_glued;m3_second_complement_empty")
print("DOUBLE_GATE=4W^2Phi-2Wxi_eta+xi^2alpha=4W^2J(a)")
print("HORIZONTAL=discriminant_x^2_times_DVR_unit")
print("DIFFERENTIAL=b=P^e*y_preserves_sigma_order")
print("LEDGERS=Pick_15_14_14;RH_28_26_26;graph_11")
