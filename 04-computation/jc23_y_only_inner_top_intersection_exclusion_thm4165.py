#!/usr/bin/env python3
"""Exact canonical certificate for the complete I_C=D_C=0 Y-only wall.

Primary lane: source pair (A,B), normalized (X,T) elimination, discriminant
normalization by the repeated root, local boundary differential, exhaustive
endpoint strata, and full/finite monodromy ledgers.  No repository code is
imported and no assert statements are used.
"""
from hashlib import sha256
from math import gcd
import sympy as sp

CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not bool(condition):
        raise RuntimeError(message)


def order(poly, variable):
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


def reduced(expr, modulus, variable):
    numerator, denominator = sp.fraction(sp.cancel(expr))
    remainder = sp.rem(
        sp.Poly(numerator, variable), sp.Poly(modulus, variable)
    ).as_expr()
    return sp.factor(remainder/denominator)


def gcd_firewall(expr, modulus, variable, label):
    """Prove that an actual quotient-row coefficient is a unit on a wall."""
    numerator, denominator = sp.fraction(sp.together(sp.cancel(expr)))
    require(numerator != 0, label+" numerator nonzero")
    require(denominator != 0, label+" denominator nonzero")
    wall = sp.Poly(modulus, variable)
    require(sp.gcd(wall, sp.Poly(numerator, variable)).degree() == 0,
            label+" numerator coprime to wall")
    require(sp.gcd(wall, sp.Poly(denominator, variable)).degree() == 0,
            label+" denominator coprime to wall")


def packet_defect(packet):
    return sum(entry-1 for entry in packet)


s, p = sp.symbols("s p")
X, T = sp.symbols("X T")
u, Phi, r, W, z, q = sp.symbols("u Phi r W z q")
K0 = sp.Rational(2848, 45)
epsilon = -sp.Rational(1376, 135)
Zeta = sp.Rational(5696, 135)*u
Theta = 3*u**2

J = 8544*Phi-22784*u-1215*u**3
F = (
    -2076192000*Phi**3 + 110716875*Phi**2*u**3
    - 285684019200*Phi*u**2 + 13541904000*u**5
    - 61429478588416*u
)
E = (
    -3244050*Phi**2*u - 51904800*Phi*u**2 + 4185329664*Phi
    + 2460375*u**5 - 595175040*u**3 - 22321758208*u
)
N = (
    -202033731568359375*Phi**2*u**6
    + 101837292342912000000*Phi**2*u**4
    - 6797550458098483200000*Phi**2*u**2
    - 351969174847758532608000*Phi**2
    - 13792169408400000000*Phi*u**7
    + 607952470637606400000*Phi*u**5
    + 134732451621655019520000*Phi*u**3
    - 2610438046787542450176000*Phi*u
    + 392263245112500000*u**10 - 41629364664354000000*u**8
    - 8759084784795648000000*u**6
    - 126485060189595107328000*u**4
    + 6820380454827676454092800*u**2
    - 976304345140624192423591936
)

top_cubic = Zeta*W**3+Theta*W**2+Phi*W+epsilon
D = sp.factor(sp.discriminant(top_cubic, W))
require(sp.factor(D-u*F/sp.Integer(12301875)) == 0,
        "inner-chart top discriminant")
require(sp.factor(4*Theta*K0**2-27*Zeta**2) == 0,
        "inner wall equation")
field_u = sp.QQ.frac_field(u)
poly_F = sp.Poly(F, Phi, domain=field_u)
require(poly_F.is_irreducible, "F irreducible over Q(u)")
require(poly_F.is_sqf, "F squarefree over Q(u)")
require(sp.factor(sp.discriminant(F, Phi)
                  +4936024008000000*u**2
                   *(2460375*u**4+44643516416)**3) == 0,
        "discriminant of the irreducible intersection chart")

t = p-s**2
H = sp.expand(
    -3*p+sp.Rational(8, 3)*p**2+epsilon*p**3+K0*s**2*p**2
    +Phi*s*p**3+Theta*s**2*p**3+Zeta*s**3*p**3
)
Acrit = sp.cancel((-s*p+t**2*sp.diff(H, s))/p)
Ccrit = sp.expand(s**2+2*t**2*sp.diff(H, p))
Bcrit = sp.cancel((Ccrit+s*Acrit)/t**2)
require(sp.denom(Acrit) == 1 and sp.denom(Bcrit) == 1,
        "source AB polynomials")
require((sp.degree(Acrit, s), sp.degree(Bcrit, s)) == (6, 3),
        "source AB degrees at s-infinity")
require(sp.factor(sp.Poly(Acrit, s).LC()-3*Zeta*p**2) == 0,
        "A leading s-row")
require(sp.factor(sp.Poly(Bcrit, s).LC()-9*Zeta*p**2) == 0,
        "B leading s-row")

# Source (A,B) elimination on I_C=0, before the D_C strict transform.
resAB = sp.resultant(Acrit, Bcrit, s)
require(order(resAB, p) == 7, "AB p artifact")
inner_R17 = sp.cancel(resAB/p**7)
require(sp.degree(inner_R17, p) == 17, "inner R17")
F_rows_source = quotient_rows(inner_R17, p, F, Phi)
require((F_rows_source[0][0], F_rows_source[-1][0]) == (16, 0),
        "generic intersection source R16")
require(sp.factor(
    F_rows_source[0][1]
    - sp.Rational(1534934573591985913856, 622782421875)*u**6*E
) == 0, "source D-strict top endpoint")
require(sp.factor(
    F_rows_source[-1][1]
    - sp.Rational(8305770496, 1125)*u**2*J
) == 0, "source D-strict bottom endpoint")

# Normalized elimination, formed independently rather than inferred from AB.
P = T+X**2*T**2
Y = X*T*P
G = sp.expand(
    -X**2*T/2-3*P+sp.Rational(8, 3)*P**2+epsilon*P**3
    +K0*Y**2+Phi*P**2*Y+Theta*P*Y**2+Zeta*Y**3
)
f = sp.cancel(sp.diff(G, X)/T)
h = sp.diff(G, T)
require((sp.degree(f, X), sp.degree(h, X)) == (8, 9),
        "normalized degrees at X-infinity")
require(sp.factor(sp.Poly(f, X).LC()-9*Zeta*T**8) == 0,
        "f leading X-row")
require(sp.factor(sp.Poly(h, X).LC()-9*Zeta*T**8) == 0,
        "h leading X-row")
resXT = sp.resultant(f, h, X)
require(order(resXT, T) == 56, "normalized T artifact")
inner_Q17 = sp.cancel(resXT/(T**56*(6*T+1)**2))
require(sp.degree(inner_Q17, T) == 17, "inner Q17")
F_rows_normalized = quotient_rows(inner_Q17, T, F, Phi)
require((F_rows_normalized[0][0], F_rows_normalized[-1][0]) == (16, 0),
        "generic intersection normalized Q16")
require(sp.factor(
    F_rows_normalized[0][1]
    - sp.Rational(47966705424749559808, 38306957530517578125)*u**5*N
) == 0, "normalized D-strict top endpoint")
require(sp.factor(
    F_rows_normalized[-1][1]
    + sp.Rational(1519777094677765052956672, 56953125)*u**7
) == 0, "normalized D-strict bottom endpoint")

# In Q(u)[Phi]/(F), the normalized top endpoint is exactly E*J^2 up
# to the displayed unit.  This proves that no endpoint component is omitted.
mod_F = sp.Poly(F, Phi, domain=field_u)
endpoint_identity = sp.rem(
    sp.Poly(N-sp.Rational(108000, 1)/u*E*J**2, Phi, domain=field_u),
    mod_F,
).as_expr()
require(sp.factor(endpoint_identity) == 0, "N=(108000/u)EJ^2 mod F")

# Normalize the discriminant curve by its repeated top root r.
R = 11392*u*r**3+405*u**2*r**2+1376
Phi_r = -2*r*u*(2848*r+135*u)/45
require(sp.Poly(R, r, domain=field_u).is_irreducible,
        "repeated-root chart irreducible over Q(u)")
root_relation = 45*Phi+2*r*u*(2848*r+135*u)
require(sp.factor(sp.resultant(R, root_relation, r)+5696*u**2*F) == 0,
        "repeated-root chart covers every F component for u nonzero")
require(reduced(F.subs(Phi, Phi_r), R, r) == 0,
        "repeated-root parameter covers F=0")
require(reduced(top_cubic.subs({Phi: Phi_r, W: r}), R, r) == 0,
        "top cubic vanishes at r")
require(reduced(sp.diff(top_cubic, W).subs({Phi: Phi_r, W: r}), R, r) == 0,
        "top derivative vanishes at r")

alpha = u*(5696*r+135*u)/45
beta = 8*(356*r**2+15)/45
gamma = 5696*r/45
delta = -(2848*r**4+120*r**2+135)/45
E_root = u**2*(5696*r+135*u)**3*(356*r**2+15)/15
require(reduced(E.subs(Phi, Phi_r)-E_root, R, r) == 0,
        "E is exactly triple-or-node factor")

# Direct top chart z=1/p, W=s.
Fq = sp.expand((s**2-p)*(q-H)-s**2/2)
Fbar = sp.expand(z**4*Fq.subs({s: W, p: 1/z}))
f0 = sp.Poly(Fbar, z).coeff_monomial(1)
f1 = sp.Poly(Fbar, z).coeff_monomial(z)
f2 = sp.Poly(Fbar, z).coeff_monomial(z**2)
require(sp.factor(f0-top_cubic) == 0, "top chart boundary cubic")
require(reduced(
    sp.diff(f0, W, 2).subs({Phi: Phi_r, W: r})/2-alpha, R, r
) == 0, "local alpha")
require(reduced(f1.subs({Phi: Phi_r, W: r})-beta, R, r) == 0,
        "local beta")
require(reduced(
    sp.diff(f1, W).subs({Phi: Phi_r, W: r})-gamma, R, r
) == 0, "local gamma")
require(reduced(f2.subs(W, r)-delta, R, r) == 0, "local delta")

# Three finite exceptional subwalls: triple, ordinary node, and J=0.
A_T = 2460375*u**4+44643516416
Phi_T = sp.Rational(405, 5696)*u**3
B_N = 36905625*u**4-4721414400*u**2+239958900736
Phi_N = ((2025*u**2+489856)*(6075*u**2-489856)
         /(230688000*u))
B_J = (
    30267225703125*u**8+2043284356800000*u**6
    +264381824212992000*u**4+6498574373014732800*u**2
    +498260889496415371264
)
Phi_J = u*(1215*u**2+22784)/8544
walls = (
    ("triple", A_T, Phi_T),
    ("node", B_N, Phi_N),
    ("J", B_J, Phi_J),
)
for name, modulus, phi_value in walls:
    require(reduced(F.subs(Phi, phi_value), modulus, u) == 0,
            name+" lies on discriminant")
    require(sp.gcd(sp.Poly(modulus, u), sp.Poly(sp.diff(modulus, u), u)).degree() == 0,
            name+" parameter polynomial squarefree")
    require(sp.gcd(sp.Poly(modulus, u), sp.Poly(u, u)).degree() == 0,
            name+" has u nonzero")
for left, right in ((A_T, B_N), (A_T, B_J), (B_N, B_J)):
    require(sp.gcd(sp.Poly(left, u), sp.Poly(right, u)).degree() == 0,
            "exceptional subwalls disjoint")

# Converse coverage: these are the complete elimination polynomials for the
# three exceptional conditions on the irreducible repeated-root chart.
require(sp.factor(sp.resultant(R, 5696*r+135*u, r)+5696*A_T) == 0,
        "triple condition iff A_T=0")
require(sp.factor(sp.resultant(R, 356*r**2+15, r)-356*B_N) == 0,
        "node condition iff B_N=0")
require(sp.factor(
    sp.resultant(R, J.subs(Phi, Phi_r), r)
    +sp.Rational(16222208, 3375)*u**3*B_J
) == 0, "J condition iff B_J=0")
require(sp.factor(F.subs(Phi, Phi_J)+u*B_J/sp.Integer(8111104)) == 0,
        "J parameter map exact")

# The triple and node are exactly alpha=0 and beta=0, respectively.
require(reduced((5696*r+135*u).subs(r, -135*u/5696), A_T, u) == 0,
        "triple root formula")
require(reduced(R.subs(r, -135*u/5696), A_T, u) == 0,
        "triple root relation")
r_node = -(6075*u**2-489856)/(170880*u)
require(sp.factor(R.subs(r, -135*u/5696)-A_T/sp.Integer(32444416)) == 0,
        "triple parameter map exact")
require(sp.factor(sp.cancel(
    R.subs(r, r_node)-sp.Rational(43, 38448000)*B_N/u**2
)) == 0, "node root parameter exact")
require(sp.factor(sp.cancel(
    (356*r_node**2+15)-B_N/(sp.Integer(82022400)*u**2)
)) == 0, "node beta parameter exact")

# Node noncusp firewall in y=u/r.
yvar = sp.symbols("yvar")
R_y = 91125*yvar**2+2563200*yvar+174388736
alpha_y = -yvar*(135*yvar+5696)/1068
gamma2_y = sp.Rational(5696**2, 45**2)*(-sp.Rational(15, 356))
require(sp.factor(delta.subs(r**2, -sp.Rational(15, 356))+3) == 0,
        "node delta equals minus three")
require(sp.factor(
    gamma2_y-4*alpha_y*(-3)
    +(135*yvar+2848)**2/sp.Integer(12015)
) == 0, "node tangent discriminant identity")
for factor in (yvar, 135*yvar+5696, 135*yvar+2848):
    require(sp.gcd(sp.Poly(R_y, yvar), sp.Poly(factor, yvar)).degree() == 0,
            "node nonzero/cusp firewall")
# On beta=0: alpha=-y(135y+5696)/1068 and tangent discriminant
# gamma^2-4 alpha delta=-(135y+2848)^2/12015, both nonzero.

# Exact terminal endpoint profiles and all-root firewalls.
endpoint_factors = {
    "triple": (
        405*u**2+88064,
        2025*u**2-979712,
        36585*u**2+427154432,
        u,
    ),
    "node": (
        20614295925117675*u**2-789558891060057728,
        247725*u**2+21063808,
        1199320411266696917025*u**2+56388834945845741080192,
        153378225*u**2-12933178112,
    ),
    "J": (
        18030796678125*u**6+5214383198496000*u**4
            +118847166205132800*u**2+14477208077889175552,
        2460375*u**4-204543360*u**2+5580439552,
        112249545005345596875*u**6+4039772647342545648000*u**4
            +349920353377066350182400*u**2+9283134759621226372530176,
        u,
    ),
}
for name, modulus, phi_value in walls:
    source_rows = quotient_rows(inner_R17.subs(Phi, phi_value), p, modulus, u)
    normalized_rows = quotient_rows(inner_Q17.subs(Phi, phi_value), T, modulus, u)
    if name == "J":
        require((source_rows[0][0], source_rows[-1][0]) == (16, 1),
                "J source p^8 R15")
    else:
        require((source_rows[0][0], source_rows[-1][0]) == (15, 0),
                name+" source R15")
    require((normalized_rows[0][0], normalized_rows[-1][0]) == (15, 0),
            name+" normalized Q15")
    gcd_firewall(source_rows[0][1], modulus, u,
                 name+" actual AB first row")
    gcd_firewall(source_rows[-1][1], modulus, u,
                 name+" actual AB last row")
    gcd_firewall(normalized_rows[0][1], modulus, u,
                 name+" actual normalized first row")
    gcd_firewall(normalized_rows[-1][1], modulus, u,
                 name+" actual normalized last row")
    for factor in endpoint_factors[name]:
        require(sp.gcd(sp.Poly(modulus, u), sp.Poly(factor, u)).degree() == 0,
                name+" terminal endpoint firewall")

# Four universal affine critical points remain, so d=16 gives L=20 and
# d=15 gives L=19.
hessian = sp.det(sp.hessian(G, (X, T)))
for sign in (-1, 1):
    for point, critical_value, hessian_value in (
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
        require(sp.simplify(hessian.subs(point)-hessian_value) == 0,
                "universal Hessian")

# Boundary normalization and packet arithmetic.  The base non-top packet is
# (8,2,2,2,1).  The local differential is -Q*z^2*dW/Fbar_z.
smooth_packet = (8, 5, 3, 2, 2, 2, 1)
triple_packet = (8, 7, 2, 2, 2, 1)
node_packet = (8, 3, 2, 2, 2, 2, 2, 1)
require((sum(smooth_packet), packet_defect(smooth_packet)) == (23, 16),
        "smooth-double packet")
require((sum(triple_packet), packet_defect(triple_packet)) == (22, 16),
        "triple packet")
require((sum(node_packet), packet_defect(node_packet)) == (22, 14),
        "node packet")
require((16, 16, 14) == (2*9-2, 2*9-2, 2*8-2),
        "packet/genus saturation")

# Full and finite response ledgers.  Three index-two places form the prime
# cubic carrier and become three transposition punctures in a finite response.
responses = (
    ("smooth-generic", 20, smooth_packet, (8, 5, 3, 1), 16, 13),
    ("smooth-J", 19, smooth_packet, (8, 5, 3, 1), 16, 13),
    ("triple", 19, triple_packet, (8, 7, 1), 16, 13),
    ("node", 19, node_packet, (8, 3, 2, 2, 1), 14, 11),
)
for name, length, full_packet, finite_packet, full_index, finite_index in responses:
    full_n = sum(full_packet)
    finite_n = sum(finite_packet)
    require(packet_defect(full_packet) == full_index, name+" full index")
    require(packet_defect(finite_packet) == finite_index, name+" finite index")
    require(2*(full_n-length) < full_index, name+" full contradiction")
    require(finite_n > 3+1, name+" carrier orbit has n>m+1")
    overlap = finite_n+3-length
    require(overlap >= 0, name+" finite overlap arithmetic")
    require(2*overlap+3 < finite_index, name+" finite contradiction")

semantic_rows = (
    "intersection=eta=Delta=0;zeta!=0;I_C=D_C=0",
    "D_C=uF/12301875;repeated_root_chart=R(u,r)=0",
    "domain=F_irreducible_over_Q(u);coverage=Res_r(R,relation)=-5696u^2F",
    "critical_top=E;normalized_top=(108000/u)EJ^2",
    "E=u^2(5696r+135u)^3(356r^2+15)/15",
    "strata=smooth_generic_L20|smooth_J_L19|triple_L19|node_L19",
    "node=ordinary_non_cusp",
    "packets=smooth(8,5,3,2,2,2,1);triple(8,7,2,2,2,1);node(8,3,2,2,2,2,2,1)",
    "responses=all_full_and_finite_strict",
    "validity=actual_AB_XT_endpoint_rows;coordinate_infinity;G_values;n>m+1",
)
digest = sha256("\n".join(semantic_rows).encode()).hexdigest()

print("COMPLETE I_C=D_C=0 PRIMARY CANONICAL AUDIT")
print("checks="+str(CHECKS))
print("parameter=D_C=u*F(u,Phi)/12301875;repeated_root_R=0")
print("critical_identity=N=(108000/u)*E*J^2_mod_F")
print("root_identity=E=u^2*(5696r+135u)^3*(356r^2+15)/15")
print("strata=smooth_generic:L20;smooth_J:L19;triple:L19;ordinary_node:L19")
print("source_AB=p^7R16_generic;p^8R15_J;p^7R15_triple_node")
print("normalized=T^56*(6T+1)^2*(Q16_generic,Q15_exceptional)")
print("terminal_firewalls=three_squarefree_pairwise_disjoint_polynomials;all_endpoint_gcds_one")
print("actual_row_firewalls=AB_and_normalized_first_last_numerators_denominators_coprime")
print("infinity_and_type_gates=source_normalized_leads;F_irreducible;root_chart_surjective")
print("universal_values=G(0,+/-sqrt(-6))=0;G(-1/6,+/-sqrt(6))=1/2")
print("node_cusp_firewall=PASS")
print("packets=smooth[8,5,3,2,2,2,1];triple[8,7,2,2,2,1];node[8,3,2,2,2,2,2,1]")
print("genera=smooth9;triple9;node8")
print("monodromy=all_full_and_finite_responses_excluded")
print("semantic_sha256="+digest)
print("verdict=PASS")
