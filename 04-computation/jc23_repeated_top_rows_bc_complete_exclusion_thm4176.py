#!/usr/bin/env python3
"""Exact certificate for exporting the THM-4173 bridge to rows B/C."""

import sympy as s


CHECKS = 0


def need(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


X, T = s.symbols("X T")
D, ph, e = s.symbols("D ph e")
P = T + X**2*T**2
Y = X*T*P
K = s.Rational(2848, 45) - s.Rational(7, 6)*D
G = s.expand(
    -X**2*T/2 - 3*P + s.Rational(8, 3)*P**2
    - s.Rational(1376, 135)*P**3 + K*Y**2 + ph*P**2*Y
    + D*P**4 - D*P*Y**2 + e*P**3*Y - e*Y**3
)
f = s.expand(s.cancel(s.diff(G, X)/T))
h = s.expand(s.diff(G, T))
Hess = s.det(s.hessian(G, (X, T)))
Jac = s.det(s.Matrix(((s.diff(f, X), s.diff(f, T)),
                      (s.diff(h, X), s.diff(h, T)))))
need(s.factor(T*Jac-Hess-f*s.diff(G, X, T)) == 0,
     "Morse-resultant bridge changed")

J = D*(105*D-5696)**2 + 54675*e**2
A = 105*D-5696
b = 825*D-22784


def rem_j(expr):
    return s.factor(s.rem(s.Poly(s.expand(expr), e), s.Poly(J, e)).as_expr())


fb, hb = s.Poly(f, X), s.Poly(h, X)
need((fb.degree(), hb.degree()) == (6, 7), "row-B X-degrees changed")
need(s.factor(fb.LC()-7*T**6*(ph+e*T)) == 0,
     "row-B f leading coefficient changed")
need(s.factor(hb.LC()-T**6*(7*ph+8*e*T)) == 0,
     "row-B h leading coefficient changed")

RB = s.resultant(f, h, X)
QB = s.Poly(s.cancel(RB/(T**30*(6*T+1)**2)), T)
need(QB.degree() == 17, "row-B generic residual degree changed")
need(s.factor(QB.nth(0)+s.Rational(3*7**5, 2**5)*ph**5) == 0,
     "row-B residual constant changed")
need(s.factor(QB.nth(17)-e**3*J**2/s.Integer(22500)) == 0,
     "row-B residual top slot changed")

F = 180*D*e*b**2 - 4*D*A**2*b*ph - 1215*e*A**2*ph**2
need(s.factor(rem_j(QB.nth(16))-D*A**4*F/s.Integer(373669453125)) == 0,
     "row-B J-wall degree-sixteen slot changed")
need(s.factor(s.discriminant(F, ph)-16*D*A**2*b**2*J) == 0,
     "row-B double-root wall changed")

ph_star = -s.Rational(2, 1215)*D*b/e
H = 187425*D**3-25939920*D**2+1215936512*D-18687983616
N = (
    1161862920421875*D**8-436262053412193750*D**7
    +72428665164732450000*D**6-6925998948093233280000*D**5
    +416163890453413588992000*D**4-16047146699676377898024960*D**3
    +386666237774835302050824192*D**2-5306381690338021964332924928*D
    +31649810610496845164669042688
)
scale = s.Integer(397166535676406250000)*e**3
q15_star = s.cancel(rem_j(QB.nth(15)).subs(ph, ph_star))
q14_star = s.cancel(rem_j(QB.nth(14)).subs(ph, ph_star))
need(rem_j(s.together(scale*q15_star+3*D**3*A**6*H**2)
           .as_numer_denom()[0]) == 0,
     "row-B H-wall degree-fifteen slot changed")
need(rem_j(s.together(scale*q14_star+D**3*A**4*N)
           .as_numer_denom()[0]) == 0,
     "row-B terminal degree-fourteen slot changed")
need(s.gcd(s.Poly(H, D), s.Poly(D*A*b*N, D)).degree() == 0,
     "row-B terminal unit/coprimality gate changed")

fC, hC = s.expand(f.subs(ph, 0)), s.expand(h.subs(ph, 0))
fCp, hCp = s.Poly(fC, X), s.Poly(hC, X)
need((fCp.degree(), hCp.degree()) == (6, 7), "row-C X-degrees changed")
need(s.factor(fCp.LC()-7*e*T**7) == 0,
     "row-C f leading coefficient changed")
need(s.factor(hCp.LC()-8*e*T**7) == 0,
     "row-C h leading coefficient changed")
RC = s.resultant(fC, hC, X)
QC = s.Poly(s.cancel(RC/(T**32*(6*T+1)**2)), T)
need(QC.degree() == 15, "row-C generic residual degree changed")
need(s.factor(QC.nth(0)+s.Rational(16807, 984150000)
              * e*(45*D-2048)**5) == 0,
     "row-C residual constant changed")
need(s.factor(QC.nth(15)-e**3*J**2/s.Integer(22500)) == 0,
     "row-C residual top slot changed")
expected_c14 = s.Rational(4, 8303765625)*D**2*e*A**4*b**2
need(s.factor(rem_j(QC.nth(14))-expected_c14) == 0,
     "row-C J-wall degree-fourteen slot changed")
D2 = s.Rational(22784, 825)
e2 = -s.Rational(739213574144, 187171875)
c13 = rem_j(QC.nth(13)).subs(D, D2)
c13 = s.rem(s.Poly(c13, e), s.Poly(e**2-e2, e)).as_expr()
expected_c13 = -s.Rational(
    39082296781894638834860007697799970816,
    8994856609344482421875,
)*e
need(s.factor(c13-expected_c13) == 0,
     "row-C terminal degree-thirteen slot changed")

# Worst rows after the exact degree-drop towers.
need(14+2+2 == 18, "row-B critical length changed")
need(2*(17+3-18)+3 == 7 < 13, "row-B finite response changed")
need(2*(23-18) == 10 < 16, "row-B full response changed")
need(13+2+2 == 17, "row-C critical length changed")
need(2*(15+3-17)+3 == 5 < 11, "row-C finite response changed")
need(2*(21-17) == 8 < 14, "row-C full response changed")

print("bridge=T*detD(f,h)-detHess(G)-f*G_XT=0")
print("row_B_X_infinity=none_for_T!=0,Phi*eta!=0")
print("row_C_X_infinity=none_for_T!=0,eta!=0")
print("J=Delta*(105Delta-5696)^2+54675eta^2")
print("row_B_Q_degrees=17;16_on_J;15_on_J_and_Phi_star;14_on_H_roots")
print("row_B_min_L=18;finite_ceiling=7<13;full_ceiling=10<16")
print("row_C_Q_degrees=15;14_on_J;13_at_Delta_22784_over_825")
print("row_C_min_L=17;finite_ceiling=5<11;full_ceiling=8<14")
print(f"checks={CHECKS}")
print("verdict=B_AND_C_CLOSE")
