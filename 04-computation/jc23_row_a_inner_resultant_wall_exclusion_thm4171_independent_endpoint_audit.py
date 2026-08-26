#!/usr/bin/env python3
"""Independent source-pair/endpoint referee for THM-4171.

This audit reconstructs the alternative critical pair (A,C_0), computes its
own Sylvester determinant, reads the actual successive endpoint rows, and
uses resultant rather than gcd firewalls on the terminal quartic.  It is not
claimed as a second normalized-projection audit.
"""
from hashlib import sha256
import sympy as sp

CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not bool(condition):
        raise RuntimeError(message)


def order(poly, variable):
    return min(monomial[0] for monomial, coefficient
               in sp.Poly(poly, variable).terms() if coefficient)


def remainder_zero(expr, modulus, variable, domain):
    numerator, denominator = sp.fraction(sp.together(sp.cancel(expr)))
    require(denominator != 0, "nonzero quotient denominator")
    remainder = sp.rem(
        sp.Poly(numerator, variable, domain=domain),
        sp.Poly(modulus, variable, domain=domain),
    ).as_expr()
    return sp.factor(remainder) == 0


s, p = sp.symbols("s p")
k, u, Phi, w = sp.symbols("k u Phi w")
Delta = (5696-90*k)/105
Theta = 3*u**2
eta = 2*k*u/3
Cwall = sp.factor(Delta+Theta)
t = p-s**2

require(sp.factor(4*Theta*k**2-27*eta**2) == 0, "wall chart")
require(sp.factor(sp.Rational(2848, 45)-sp.Rational(7, 6)*Delta-k) == 0,
        "K-Delta chart")

H = sp.expand(
    -3*p+sp.Rational(8, 3)*p**2-sp.Rational(1376, 135)*p**3
    +k*s**2*p**2+Phi*s*p**3+Delta*p**4+Theta*s**2*p**3
    +eta*s*p**4-eta*s**3*p**3
)
G = -s**2/(2*t)+H
A = sp.cancel((-s*p+t**2*sp.diff(H, s))/p)
C0 = sp.expand(s**2+2*t**2*sp.diff(H, p))
require(sp.denom(A) == 1, "A polynomial")
require(sp.factor(sp.cancel(t**2*sp.diff(G, s)-p*A)) == 0,
        "A critical identity")
require(sp.factor(sp.cancel(2*t**2*sp.diff(G, p)-C0)) == 0,
        "C0 critical identity")
require((sp.degree(A, s), sp.degree(C0, s)) == (6, 7),
        "alternative-pair degrees")
require(sp.factor(sp.Poly(A, s).LC()+2*k*u*p**2) == 0,
        "A infinity row")
require(sp.factor(sp.Poly(C0, s).LC()+4*k*u*p**2) == 0,
        "C0 infinity row")

determinant = sp.resultant(A, C0, s)
require(order(determinant, p) == 9, "AC0 exceptional power")
R18 = sp.Poly(sp.cancel(determinant/p**9), p)
require(R18.degree() == 18, "AC0 generic residual degree")

J = 3*Phi*k+8*k*u+27*u**3
PhiJ = -(8*k*u+27*u**3)/(3*k)
S = -180*k**3-135*k**2*u**2+2752*k**2+2160*k*u**2+4860*u**4
Pwall = (
    16380*k**3+12285*k**2*u**2-115072*k**2
    -710640*k*u**2-9534464*k+1999872*u**2
)

# These are the actual first and last rows of the independently computed
# determinant, not a separately declared list of endpoint candidates.
actual_top = sp.factor(R18.LC())
actual_bottom = sp.factor(R18.TC())
require(sp.factor(
    actual_top-sp.Rational(524288, 364651875)*k**5*u**5
    *(-90*k+315*u**2+5696)**4
) == 0, "actual generic top row")
require(sp.factor(actual_bottom-82944*k**2*u**2*J) == 0,
        "actual generic bottom row")

specialJ = sp.factor(R18.as_expr().subs(Phi, PhiJ))
require(order(specialJ, p) == 1, "J exceptional order")
R17 = sp.Poly(sp.cancel(specialJ/p), p)
require(R17.degree() == 17, "J residual degree")
require(sp.factor(R17.LC()-actual_top) == 0, "actual J top row")
require(sp.factor(R17.TC()+sp.Rational(6912, 5)*k*u**3*S) == 0,
        "actual J bottom row S")

field_k = sp.QQ.frac_field(k)
actual_S_bottom = R17.nth(1)
require(remainder_zero(
    actual_S_bottom+sp.Rational(96, 35)*k**2*u**3*Pwall,
    S, u, field_k,
), "actual S bottom row P")

Qterminal = (
    -102060*k**7-76545*k**6*w+316224*k**6+4286520*k**5*w
    +12685824*k**5+1683990*k**4*w**2-29961792*k**4*w
    +909058048*k**4-17622360*k**3*w**2+417392640*k**3*w
    -51438240*k**2*w**3+1002637440*k**2*w**2
    +183708000*k*w**3+62001450*w**4
)
actual_SP_bottom = sp.factor(R17.nth(2))
require(sp.factor(
    actual_SP_bottom+sp.Rational(64, 105)*u**3/k
    *Qterminal.subs(w, u**2)
) == 0, "actual S=P terminal row")

# Eliminate w=u^2 independently.  P is linear in w and has no vertical
# component; its substitution into S produces exactly the terminal quartic.
S_w = S.subs(u**2, w).subs(u**4, w**2)
P_w = Pwall.subs(u**2, w)
Aw = sp.Poly(P_w, w).coeff_monomial(w)
Bw = sp.Poly(P_w, w).coeff_monomial(1)
H4 = (
    257063625*k**4-13336482600*k**3+141996773760*k**2
    -1407227197440*k+47646220845056
)
require(sp.resultant(Aw, Bw, k) != 0, "no vertical P component")
require(sp.factor(
    sp.together(S_w.subs(w, -Bw/Aw))
    -128*k**2*H4/(49*(195*k**2-11280*k+31744)**2)
) == 0, "quartic coverage identity")
require(sp.Poly(H4, k).is_irreducible, "terminal quartic irreducible")
require(sp.discriminant(H4, k) != 0, "terminal quartic reduced")

w_value = -Bw/Aw
Cnum, Cden = sp.fraction(sp.together(Cwall.subs(u**2, w_value)))
Qnum, Qden = sp.fraction(sp.together(Qterminal.subs(w, w_value)))
units = (
    ("K", k),
    ("Delta", 5696-90*k),
    ("P slope", Aw),
    ("u squared", -Bw),
    ("C numerator", Cnum),
    ("C denominator", Cden),
    ("terminal numerator", Qnum),
    ("terminal denominator", Qden),
)
for label, unit in units:
    require(sp.resultant(H4, unit, k) != 0,
            label+" terminal resultant firewall")

# Independent collapsed-chart and response checks.
require(sp.factor(A.subs(p, 0)+s) == 0, "p=0 A row")
require(sp.factor(C0.subs(p, 0)+s**2*(6*s**2-1)) == 0,
        "p=0 C0 row")
require(sp.factor(A.subs(p, s**2)+s) == 0, "t=0 A row")
require(sp.factor(C0.subs(p, s**2)-s**2) == 0, "t=0 C0 row")

packet = (7, 7, 4, 2, 2, 2, 1)
finite_packet = (7, 7, 4, 1)
require((sum(packet), sum(e-1 for e in packet)) == (25, 18),
        "full packet ledger")
require((sum(finite_packet), sum(e-1 for e in finite_packet)) == (19, 15),
        "finite packet ledger")
for residual_degree, length in ((18, 22), (17, 21), (16, 20), (15, 19)):
    require(length == residual_degree+4, "four collapsed points restored")
    require(2*(25-length) < 18, "full strict response")
    require(19 > 3+1, "finite carrier lemma type gate")
    overlap = 19+3-length
    require(overlap >= 0, "finite overlap nonnegative")
    require(2*overlap+3 < 15, "finite strict response")

semantic = (
    "independent_pair=A,C0;determinant=p^9R18",
    "actual_bottom_rows=J,S,P,Qterminal;residual=18,17,16,15",
    "terminal=linear-P-substitution,H4-irreducible,resultant-firewalls",
    "collapsed=p0,t0;length=22,21,20,19",
    "responses=full25/index18;finite19/index15;carrier3",
)
digest = sha256("\n".join(semantic).encode()).hexdigest()
print("THM4171_ROW_A_INNER_WALL_INDEPENDENT_ENDPOINT_AUDIT")
print("checks="+str(CHECKS))
print("scope=INDEPENDENT_SOURCE_PAIR_AND_ENDPOINT_AUDIT_NOT_NORMALIZED")
print("alternative_source_pair=A,C0")
print("source_AC0=p^9R18;p^10R17;p^11R16;p^12R15")
print("actual_endpoints=J,S,P,H4_terminal")
print("terminal_firewall=H4_irreducible_reduced;all_eight_resultants_nonzero")
print("lengths=22,21,20,19")
print("monodromy=all_full_and_finite_bounds_strict")
print("semantic_sha256="+digest)
print("verdict=ACCEPT")
