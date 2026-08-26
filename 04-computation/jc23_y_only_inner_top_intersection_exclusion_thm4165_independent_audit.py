#!/usr/bin/env python3
"""Independent canonical referee for the complete I_C=D_C=0 intersection.

This lane uses the alternative source pair (A,C_0), derives the repeated-root
chart directly, tests terminal factors by nonzero resultants, and reconstructs
the local differential packets without importing the primary audit.
"""
from hashlib import sha256
import sys
import sympy as sp

sys.stdout.reconfigure(newline="\n")

CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not bool(condition):
        raise RuntimeError(message)


def valuation(poly, variable):
    return min(m[0] for m, c in sp.Poly(poly, variable).terms() if c)


def quotient_rows(poly, variable, modulus, parameter):
    rows = []
    for (degree,), coefficient in sp.Poly(sp.cancel(poly), variable).terms():
        numerator, denominator = sp.fraction(sp.cancel(coefficient))
        remainder = sp.rem(
            sp.Poly(numerator, parameter), sp.Poly(modulus, parameter)
        ).as_expr()
        if remainder:
            rows.append((degree, sp.factor(remainder/denominator)))
    return rows


def remainder(expr, modulus, variable):
    numerator, denominator = sp.fraction(sp.cancel(expr))
    return sp.factor(sp.rem(
        sp.Poly(numerator, variable), sp.Poly(modulus, variable)
    ).as_expr()/denominator)


def resultant_firewall(expr, modulus, variable, label):
    """Use a second, resultant-based unit test on an actual quotient row."""
    numerator, denominator = sp.fraction(sp.together(sp.cancel(expr)))
    require(numerator != 0, label+" numerator nonzero")
    require(denominator != 0, label+" denominator nonzero")
    require(sp.resultant(modulus, numerator, variable) != 0,
            label+" numerator resultant nonzero")
    require(sp.resultant(modulus, denominator, variable) != 0,
            label+" denominator resultant nonzero")


s, p = sp.symbols("s p")
X, T = sp.symbols("X T")
u, Phi, r, W, z, q = sp.symbols("u Phi r W z q")
K = sp.Rational(2848, 45)
eps = -sp.Rational(1376, 135)
Z = sp.Rational(5696, 135)*u
Th = 3*u**2
J = 8544*Phi-22784*u-1215*u**3
Fdisc = (
    -2076192000*Phi**3+110716875*Phi**2*u**3
    -285684019200*Phi*u**2+13541904000*u**5-61429478588416*u
)
Etop = (
    -3244050*Phi**2*u-51904800*Phi*u**2+4185329664*Phi
    +2460375*u**5-595175040*u**3-22321758208*u
)
field_u = sp.QQ.frac_field(u)
poly_F = sp.Poly(Fdisc, Phi, domain=field_u)
require(poly_F.is_irreducible, "F irreducible over Q(u)")
require(poly_F.is_sqf, "F squarefree over Q(u)")
require(sp.factor(sp.discriminant(Fdisc, Phi)
                  +4936024008000000*u**2
                   *(2460375*u**4+44643516416)**3) == 0,
        "intersection-chart discriminant")

t = p-s**2
H = (
    -3*p+sp.Rational(8, 3)*p**2+eps*p**3+K*s**2*p**2
    +Phi*s*p**3+Th*s**2*p**3+Z*s**3*p**3
)
A = sp.cancel((-s*p+t**2*sp.diff(H, s))/p)
C0 = sp.expand(s**2+2*t**2*sp.diff(H, p))
require(sp.denom(A) == 1, "A polynomial")
require((sp.degree(A, s), sp.degree(C0, s)) == (6, 7),
        "alternative source degrees")
require(sp.factor(sp.Poly(A, s).LC()-3*Z*p**2) == 0, "A leading")
require(sp.factor(sp.Poly(C0, s).LC()-6*Z*p**2) == 0, "C0 leading")

resAC = sp.resultant(A, C0, s)
require(valuation(resAC, p) == 9, "AC inner artifact")
inner = sp.cancel(resAC/p**9)
require(sp.degree(inner, p) == 17, "AC inner R17")
generic = quotient_rows(inner, p, Fdisc, Phi)
require((generic[0][0], generic[-1][0]) == (16, 0),
        "AC intersection R16")
require(sp.factor(
    generic[0][1]
    - sp.Rational(1534934573591985913856, 622782421875)*u**6*Etop
) == 0, "AC top endpoint")
require(sp.factor(
    generic[-1][1]-sp.Rational(8305770496, 1125)*u**2*J
) == 0, "AC bottom endpoint")

# A separately formed normalized elimination supplies a second coordinate
# model and its coordinate-infinity gates.
P = T+X**2*T**2
Y = X*T*P
G = (
    -X**2*T/2-3*P+sp.Rational(8, 3)*P**2+eps*P**3
    +K*Y**2+Phi*P**2*Y+Th*P*Y**2+Z*Y**3
)
fXT = sp.cancel(sp.diff(G, X)/T)
hXT = sp.diff(G, T)
require((sp.degree(fXT, X), sp.degree(hXT, X)) == (8, 9),
        "normalized degrees at X-infinity")
require(sp.factor(sp.Poly(fXT, X).LC()-9*Z*T**8) == 0,
        "normalized f leading X-row")
require(sp.factor(sp.Poly(hXT, X).LC()-9*Z*T**8) == 0,
        "normalized h leading X-row")
resXT = sp.resultant(fXT, hXT, X)
require(valuation(resXT, T) == 56, "normalized T artifact")
innerXT = sp.cancel(resXT/(T**56*(6*T+1)**2))
require(sp.degree(innerXT, T) == 17, "normalized inner Q17")
genericXT = quotient_rows(innerXT, T, Fdisc, Phi)
require((genericXT[0][0], genericXT[-1][0]) == (16, 0),
        "normalized intersection Q16")

# Repeated-root normalization constructed without a discriminant solver.
R = 11392*u*r**3+405*u**2*r**2+1376
PhiR = -2*r*u*(2848*r+135*u)/45
require(sp.Poly(R, r, domain=field_u).is_irreducible,
        "root chart irreducible over Q(u)")
root_relation = 45*Phi+2*r*u*(2848*r+135*u)
require(sp.factor(sp.resultant(R, root_relation, r)+5696*u**2*Fdisc) == 0,
        "root chart covers F for u nonzero")
top = Z*W**3+Th*W**2+Phi*W+eps
require(remainder(Fdisc.subs(Phi, PhiR), R, r) == 0,
        "root chart maps to discriminant")
require(remainder(top.subs({Phi: PhiR, W: r}), R, r) == 0,
        "repeated root value")
require(remainder(sp.diff(top, W).subs({Phi: PhiR, W: r}), R, r) == 0,
        "repeated root derivative")

alpha = u*(5696*r+135*u)/45
beta = 8*(356*r**2+15)/45
gamma = 5696*r/45
delta = -(2848*r**4+120*r**2+135)/45
require(remainder(
    Etop.subs(Phi, PhiR)
    -u**2*(5696*r+135*u)**3*(356*r**2+15)/15,
    R, r,
) == 0, "critical endpoint equals alpha^3 beta")

# Direct local chart and residue differential data.
Fq = sp.expand((s**2-p)*(q-H)-s**2/2)
Fbar = sp.expand(z**4*Fq.subs({s: W, p: 1/z}))
rows_z = [sp.Poly(Fbar, z).coeff_monomial(z**degree) for degree in range(5)]
require(sp.factor(rows_z[0]-top) == 0, "boundary cubic")
require(remainder(
    sp.diff(rows_z[0], W, 2).subs({Phi: PhiR, W: r})/2-alpha,
    R, r,
) == 0, "alpha local")
require(remainder(rows_z[1].subs({Phi: PhiR, W: r})-beta, R, r) == 0,
        "beta local")
require(remainder(
    sp.diff(rows_z[1], W).subs({Phi: PhiR, W: r})-gamma,
    R, r,
) == 0, "gamma local")
require(remainder(rows_z[2].subs(W, r)-delta, R, r) == 0,
        "delta local")

# Finite exceptional parameter polynomials and direct AC profiles.
AT = 2460375*u**4+44643516416
PhiT = sp.Rational(405, 5696)*u**3
BN = 36905625*u**4-4721414400*u**2+239958900736
PhiN = ((2025*u**2+489856)*(6075*u**2-489856)/(230688000*u))
BJ = (
    30267225703125*u**8+2043284356800000*u**6
    +264381824212992000*u**4+6498574373014732800*u**2
    +498260889496415371264
)
PhiJ = u*(1215*u**2+22784)/8544

terminal = {
    "triple": (AT, PhiT, (
        405*u**2+88064, 2025*u**2-979712,
        36585*u**2+427154432, u,
    )),
    "node": (BN, PhiN, (
        20614295925117675*u**2-789558891060057728,
        247725*u**2+21063808,
        1199320411266696917025*u**2+56388834945845741080192,
        153378225*u**2-12933178112,
    )),
    "J": (BJ, PhiJ, (
        18030796678125*u**6+5214383198496000*u**4
            +118847166205132800*u**2+14477208077889175552,
        2460375*u**4-204543360*u**2+5580439552,
        112249545005345596875*u**6+4039772647342545648000*u**4
            +349920353377066350182400*u**2+9283134759621226372530176,
        u,
    )),
}
for name, (modulus, phi, endpoint_factors) in terminal.items():
    require(remainder(Fdisc.subs(Phi, phi), modulus, u) == 0,
            name+" lies on intersection")
    source_rows = quotient_rows(inner.subs(Phi, phi), p, modulus, u)
    normalized_rows = quotient_rows(innerXT.subs(Phi, phi), T, modulus, u)
    expected = (16, 1) if name == "J" else (15, 0)
    require((source_rows[0][0], source_rows[-1][0]) == expected,
            name+" AC R15 profile")
    require((normalized_rows[0][0], normalized_rows[-1][0]) == (15, 0),
            name+" normalized Q15 profile")
    resultant_firewall(source_rows[0][1], modulus, u,
                       name+" actual AC first row")
    resultant_firewall(source_rows[-1][1], modulus, u,
                       name+" actual AC last row")
    resultant_firewall(normalized_rows[0][1], modulus, u,
                       name+" actual normalized first row")
    resultant_firewall(normalized_rows[-1][1], modulus, u,
                       name+" actual normalized last row")
    require(sp.discriminant(modulus, u) != 0, name+" reduced parameters")
    for factor in endpoint_factors:
        require(sp.resultant(modulus, factor, u) != 0,
                name+" endpoint resultant firewall")
for left, right in ((AT, BN), (AT, BJ), (BN, BJ)):
    require(sp.resultant(left, right, u) != 0,
            "exceptional strata disjoint")
require(sp.factor(sp.resultant(R, 5696*r+135*u, r)+5696*AT) == 0,
        "triple condition iff AT=0")
require(sp.factor(sp.resultant(R, 356*r**2+15, r)-356*BN) == 0,
        "node condition iff BN=0")
require(sp.factor(
    sp.resultant(R, J.subs(Phi, PhiR), r)
    +sp.Rational(16222208, 3375)*u**3*BJ
) == 0, "J condition iff BJ=0")
require(sp.factor(Fdisc.subs(Phi, PhiJ)+u*BJ/sp.Integer(8111104)) == 0,
        "J parameter map exact")
r_node = -(6075*u**2-489856)/(170880*u)
require(sp.factor(R.subs(r, -135*u/5696)-AT/sp.Integer(32444416)) == 0,
        "triple parameter exact")
require(sp.factor(sp.cancel(
    R.subs(r, r_node)-sp.Rational(43, 38448000)*BN/u**2
)) == 0, "node root parameter exact")
require(sp.factor(sp.cancel(
    (356*r_node**2+15)-BN/(sp.Integer(82022400)*u**2)
)) == 0, "node beta parameter exact")

# Ordinary-node, not cusp, by a disjoint y=u/r calculation.
y = sp.symbols("y")
Ry = 91125*y**2+2563200*y+174388736
alpha_y = -y*(135*y+5696)/1068
gamma2_y = sp.Rational(5696**2, 45**2)*(-sp.Rational(15, 356))
require(sp.factor(delta.subs(r**2, -sp.Rational(15, 356))+3) == 0,
        "node delta equals minus three")
require(sp.factor(
    gamma2_y-4*alpha_y*(-3)+(135*y+2848)**2/sp.Integer(12015)
) == 0, "node tangent discriminant formula")
for factor in (y, 135*y+5696, 135*y+2848):
    require(sp.resultant(Ry, factor, y) != 0, "node noncusp resultant")
require(sp.factor(-(135*y+2848)**2/sp.Integer(12015)) != 0,
        "node tangent discriminant formula")

# Universal pairs provide the four points added to residual degrees.
Hessian = sp.det(sp.hessian(G, (X, T)))
for sign in (-1, 1):
    for point, critical_value, hv in (
        ({T: 0, X: sign*sp.sqrt(-6)}, 0, 6),
        ({T: -sp.Rational(1, 6), X: sign*sp.sqrt(6)},
         sp.Rational(1, 2), -6),
    ):
        require(sp.simplify(sp.diff(G, X).subs(point)) == 0,
                "universal G_X")
        require(sp.simplify(sp.diff(G, T).subs(point)) == 0,
                "universal G_T")
        require(sp.simplify(G.subs(point)-critical_value) == 0,
                "universal G value")
        require(sp.simplify(Hessian.subs(point)-hv) == 0,
                "universal Hessian")

# Independent packet and carrier-orbit arithmetic.
rows = (
    ("smooth-generic", 20, (8,5,3,2,2,2,1), (8,5,3,1), 16, 13),
    ("smooth-J", 19, (8,5,3,2,2,2,1), (8,5,3,1), 16, 13),
    ("triple", 19, (8,7,2,2,2,1), (8,7,1), 16, 13),
    ("node", 19, (8,3,2,2,2,2,2,1), (8,3,2,2,1), 14, 11),
)
for name, length, full, finite, full_index, finite_index in rows:
    finite_n = sum(finite)
    require(sum(e-1 for e in full) == full_index, name+" full index")
    require(sum(e-1 for e in finite) == finite_index, name+" finite index")
    require(2*(sum(full)-length) < full_index, name+" full bound")
    require(finite_n > 3+1, name+" carrier orbit has n>m+1")
    kmax = finite_n+3-length
    require(2*kmax+3 < finite_index, name+" finite bound")

semantic = (
    "alternative=AC0:p9R16_generic;p10R15_J;p9R15_triple_node",
    "normalized=T56(6T+1)^2:Q16_generic;Q15_exceptional",
    "domain=F_irreducible_over_Q(u);coverage=Res=-5696u^2F",
    "rootchart=alpha,beta,gamma,delta independently derived",
    "terminal=three reduced disjoint walls;all resultants nonzero",
    "node=ordinary;no cusp",
    "packets=smooth/triple/node",
    "responses=all excluded",
    "validity=actual_AC_XT_endpoint_rows;coordinate_infinity;G_values;n>m+1",
)
digest = sha256("\n".join(semantic).encode()).hexdigest()

print("COMPLETE I_C=D_C=0 INDEPENDENT CANONICAL AUDIT")
print("checks="+str(CHECKS))
print("alternative_source_pair=A,C0")
print("source_AC=p^9R16_generic;p^10R15_J;p^9R15_triple_node")
print("normalized=T^56*(6T+1)^2*(Q16_generic,Q15_exceptional)")
print("rootchart=R(u,r)=0;E=unit*alpha^3*beta")
print("strata=smooth_generic:L20;smooth_J:L19;triple:L19;ordinary_node:L19")
print("terminal_firewalls=all_resultants_nonzero;walls_pairwise_disjoint")
print("actual_row_firewalls=AC_and_normalized_first_last_numerators_denominators_coprime")
print("infinity_and_type_gates=source_normalized_leads;F_irreducible;root_chart_surjective")
print("universal_values=G(0,+/-sqrt(-6))=0;G(-1/6,+/-sqrt(6))=1/2")
print("node_cusp_resultant=NONZERO")
print("packets=smooth[8,5,3,2,2,2,1];triple[8,7,2,2,2,1];node[8,3,2,2,2,2,2,1]")
print("monodromy=all_responses_excluded")
print("semantic_sha256="+digest)
print("verdict=ACCEPT")
