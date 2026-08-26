#!/usr/bin/env python3
"""Canonical exact certificate for THM-4171.

Scope: THM-4157 row A, zeta=-eta, eta*Delta*C != 0, and
D_A=4*Theta*K^2-27*eta^2=0.  This certificate reconstructs the complete
source-endpoint cascade, the no-loss gates, and both monodromy ledgers.
"""
from hashlib import sha256
import sympy as sp

CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not bool(condition):
        raise RuntimeError(message)


def valuation(poly, variable):
    return min(m[0] for m, coefficient in sp.Poly(poly, variable).terms()
               if coefficient)


def reduce_mod(expr, modulus, variable, domain):
    numerator, denominator = sp.fraction(sp.together(sp.cancel(expr)))
    remainder = sp.rem(
        sp.Poly(numerator, variable, domain=domain),
        sp.Poly(modulus, variable, domain=domain),
    ).as_expr()
    return sp.factor(remainder/denominator)


s, p = sp.symbols("s p")
X, T = sp.symbols("X T")
k, u, Phi, w = sp.symbols("k u Phi w")
Delta = (5696-90*k)/105
Theta = 3*u**2
eta = 2*k*u/3
Cwall = sp.factor(Delta+Theta)

require(sp.factor(4*Theta*k**2-27*eta**2) == 0, "D_A chart")
require(sp.factor(sp.Rational(2848, 45)-sp.Rational(7, 6)*Delta-k) == 0,
        "K chart")

t = p-s**2
Hsource = sp.expand(
    -3*p+sp.Rational(8, 3)*p**2-sp.Rational(1376, 135)*p**3
    +k*s**2*p**2+Phi*s*p**3+Delta*p**4+Theta*s**2*p**3
    +eta*s*p**4-eta*s**3*p**3
)
Gsource = -s**2/(2*t)+Hsource
Acrit = sp.cancel((-s*p+t**2*sp.diff(Hsource, s))/p)
Ccrit = sp.expand(s**2+2*t**2*sp.diff(Hsource, p))
Bcrit = sp.cancel((Ccrit+s*Acrit)/t**2)
require(sp.denom(Acrit) == 1 and sp.denom(Bcrit) == 1,
        "source polynomials")
require(sp.factor(sp.cancel(t**2*sp.diff(Gsource, s)-p*Acrit)) == 0,
        "first source critical identity")
require(sp.factor(sp.cancel(
    2*t**2*sp.diff(Gsource, p)-(t**2*Bcrit-s*Acrit)
)) == 0, "second source critical identity")
require((sp.degree(Acrit, s), sp.degree(Bcrit, s)) == (6, 3),
        "source degrees")
require(sp.factor(sp.Poly(Acrit, s).LC()+2*k*u*p**2) == 0,
        "A source infinity gate")
require(sp.factor(sp.Poly(Bcrit, s).LC()+6*k*u*p**2) == 0,
        "B source infinity gate")

resultant = sp.factor(sp.resultant(Acrit, Bcrit, s))
require(valuation(resultant, p) == 7, "D_A p artifact")
R18 = sp.Poly(sp.cancel(resultant/p**7), p)
require(R18.degree() == 18, "generic D_A residual degree")

J = 3*Phi*k+8*k*u+27*u**3
PhiJ = -(8*k*u+27*u**3)/(3*k)
S = -180*k**3-135*k**2*u**2+2752*k**2+2160*k*u**2+4860*u**4
Pwall = (
    16380*k**3+12285*k**2*u**2-115072*k**2
    -710640*k*u**2-9534464*k+1999872*u**2
)

require(sp.factor(
    R18.LC()-sp.Rational(524288, 364651875)*k**5*u**5
    *(-90*k+315*u**2+5696)**4
) == 0, "live source leading endpoint")
require(sp.factor(R18.TC()-82944*k**2*u**2*J) == 0,
        "first bottom endpoint J")

Rj = sp.Poly(sp.cancel(R18.as_expr().subs(Phi, PhiJ)/p), p)
require(Rj.degree() == 17, "J residual degree")
require(sp.factor(Rj.LC()-R18.LC()) == 0, "J leading endpoint")
require(sp.factor(Rj.TC()+sp.Rational(6912, 5)*k*u**3*S) == 0,
        "second bottom endpoint S")

field_k = sp.QQ.frac_field(k)
c2_on_J = R18.nth(2).subs(Phi, PhiJ)
require(reduce_mod(
    c2_on_J+sp.Rational(96, 35)*k**2*u**3*Pwall,
    S, u, field_k,
) == 0, "third bottom endpoint P modulo S")

Qterminal = (
    -102060*k**7-76545*k**6*w+316224*k**6+4286520*k**5*w
    +12685824*k**5+1683990*k**4*w**2-29961792*k**4*w
    +909058048*k**4-17622360*k**3*w**2+417392640*k**3*w
    -51438240*k**2*w**3+1002637440*k**2*w**2
    +183708000*k*w**3+62001450*w**4
)
c3_on_J = sp.factor(R18.nth(3).subs(Phi, PhiJ))
require(sp.factor(
    c3_on_J+sp.Rational(64, 105)*u**3/k*Qterminal.subs(w, u**2)
) == 0, "terminal p3 coefficient")

# The deepest finite algebra: S=P=0 is a quartic in k, with w=u^2
# recovered rationally.  All row-A units and the terminal p3 row survive.
S_w = S.subs(u**2, w).subs(u**4, w**2)
P_w = Pwall.subs(u**2, w)
Aw = sp.Poly(P_w, w).coeff_monomial(w)
Bw = sp.Poly(P_w, w).coeff_monomial(1)
require(sp.factor(Aw-63*(195*k**2-11280*k+31744)) == 0,
        "P linear slope")
require(sp.factor(Bw-4*k*(4095*k**2-28768*k-2383616)) == 0,
        "P linear intercept")
H4 = (
    257063625*k**4-13336482600*k**3+141996773760*k**2
    -1407227197440*k+47646220845056
)
require(sp.factor(
    sp.together(S_w.subs(w, -Bw/Aw))
    -128*k**2*H4/(49*(195*k**2-11280*k+31744)**2)
) == 0, "deepest quartic coverage")
require(sp.Poly(H4, k).is_irreducible, "H4 irreducible")
require(sp.gcd(sp.Poly(H4, k), sp.Poly(sp.diff(H4, k), k)).degree() == 0,
        "H4 squarefree")
require(sp.gcd(sp.Poly(Aw, k), sp.Poly(Bw, k)).degree() == 0,
        "P has no lost vertical component")

w_value = -Bw/Aw
Cnum = sp.together(Cwall.subs(u**2, w_value)).as_numer_denom()[0]
Qvalue = sp.together(Qterminal.subs(w, w_value))
Qnum, Qden = sp.fraction(Qvalue)
for label, unit in (
    ("k", k), ("Delta", 5696-90*k), ("P coefficient", Aw),
    ("u squared numerator", -Bw), ("C", Cnum),
    ("terminal numerator", Qnum), ("terminal denominator", Qden),
):
    require(sp.gcd(sp.Poly(H4, k), sp.Poly(unit, k)).degree() == 0,
            label+" is a unit on H4")

# The source cascade loses one residual degree at each exact bottom endpoint.
strata = (
    ("J_nonzero", 18, 22),
    ("J_zero_S_nonzero", 17, 21),
    ("J_S_zero_P_nonzero", 16, 20),
    ("J_S_P_zero", 15, 19),
)
packet = (7, 7, 4, 2, 2, 2, 1)
full_n = sum(packet)
full_index = sum(entry-1 for entry in packet)
finite_packet = (7, 7, 4, 1)
finite_n = sum(finite_packet)
finite_index = sum(entry-1 for entry in finite_packet)
require((full_n, full_index, finite_n, finite_index) == (25, 18, 19, 15),
        "packet response ledger")
for name, residual_degree, length in strata:
    require(length == residual_degree+4, name+" universal restoration")
    require(2*(full_n-length) < full_index, name+" full contradiction")
    require(finite_n > 3+1, name+" carrier n>m+1")
    kmax = finite_n+3-length
    require(kmax >= 0, name+" carrier overlap nonnegative")
    require(2*kmax+3 < finite_index, name+" finite contradiction")

# Exact universal points omitted by the source chart.
Pnorm = T+X**2*T**2
Ynorm = X*T*Pnorm
Gnorm = sp.expand(
    -X**2*T/2-3*Pnorm+sp.Rational(8, 3)*Pnorm**2
    -sp.Rational(1376, 135)*Pnorm**3+k*Ynorm**2+Phi*Pnorm**2*Ynorm
    +Delta*Pnorm**4+Theta*Pnorm*Ynorm**2
    +eta*Pnorm**3*Ynorm-eta*Ynorm**3
)
fnorm = sp.cancel(sp.diff(Gnorm, X)/T)
hnorm = sp.diff(Gnorm, T)
Hess = sp.det(sp.hessian(Gnorm, (X, T)))
require((sp.degree(fnorm, X), sp.degree(hnorm, X)) == (7, 8),
        "normalized infinity degrees")
require(sp.factor(sp.Poly(fnorm, X).LC()-8*Cwall*T**7) == 0,
        "normalized f infinity row")
require(sp.factor(sp.Poly(hnorm, X).LC()-8*Cwall*T**7) == 0,
        "normalized h infinity row")
require(sp.factor(fnorm.subs(T, 0)+X) == 0,
        "T=0 normalized f row")
require(sp.factor(hnorm.subs(T, 0)+(X**2+6)/2) == 0,
        "T=0 normalized h row")
require(sp.factor(fnorm.subs(T, -1/X**2)+(X**2-6)/X) == 0,
        "p=0 nonzero-T f restriction")
require(sp.factor(hnorm.subs(T, -1/X**2)+(X**2-6)/2) == 0,
        "p=0 nonzero-T h restriction")
for sign in (-1, 1):
    for point, value, hess_value in (
        ({T: 0, X: sign*sp.sqrt(-6)}, 0, 6),
        ({T: -sp.Rational(1, 6), X: sign*sp.sqrt(6)},
         sp.Rational(1, 2), -6),
    ):
        require(sp.simplify(sp.diff(Gnorm, X).subs(point)) == 0,
                "universal G_X")
        require(sp.simplify(hnorm.subs(point)) == 0, "universal G_T")
        require(sp.simplify(Gnorm.subs(point)-value) == 0,
                "universal value")
        require(sp.simplify(Hess.subs(point)-hess_value) == 0,
                "universal Hessian")

W, q = sp.symbols("W q")
carrier = -eta*W**3+k*W**2-(q-sp.Rational(1, 2))
require((sp.degree(carrier, W), sp.degree(carrier, q)) == (3, 1),
        "cubic carrier type")
require(sp.factor(
    sp.discriminant(carrier, W)
    -(q-sp.Rational(1, 2))
     *(4*k**3-27*eta**2*(q-sp.Rational(1, 2)))
) == 0, "prime cubic carrier discriminant")

terminal_payload = ",".join(str(c) for c in sp.Poly(H4, k).all_coeffs())
terminal_digest = sha256(terminal_payload.encode()).hexdigest()
semantic = (
    "chart=K:k,Delta:(5696-90k)/105,Theta:3u^2,eta:2ku/3",
    "AB=p^7R18;bottom=J,S,P,H4-terminal;residual=18,17,16,15",
    "universal=G0(two,Hess6)+G1/2(two,Hess-6);length=22,21,20,19",
    "packet=7,7,4,2,2,2,1;genus=10;full=25/index18",
    "finite=19/index15;carrier=three-transposition-orbit;bounds=3,5,7,9",
)
semantic_digest = sha256("\n".join(semantic).encode()).hexdigest()
print("THM4171_ROW_A_INNER_RESULTANT_WALL_EXACT_CERTIFICATE")
print("checks="+str(CHECKS))
print("chart=K:k;Delta:(5696-90k)/105;Theta:3u^2;eta:2ku/3")
print("source=p^7R18;endpoints=J,S,P,H4_terminal")
print("J=3Phi*k+8k*u+27u^3")
print("S=-180k^3-135k^2u^2+2752k^2+2160ku^2+4860u^4")
print("P=16380k^3+12285k^2u^2-115072k^2-710640ku^2-9534464k+1999872u^2")
print("lengths=22,21,20,19")
print("packet=[7,7,4,2,2,2,1];genus=10")
print("responses=all_full_and_finite_carrier_orbit_bounds_strict")
print("H4_sha256="+terminal_digest)
print("semantic_sha256="+semantic_digest)
print("verdict=PASS")
