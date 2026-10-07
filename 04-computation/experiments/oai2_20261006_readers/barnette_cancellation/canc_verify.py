#!/usr/bin/env python3
"""Exact checks of the explicit identities in 'An explicit failure of complex affine-space
cancellation' (openai/math, 2026-09-23).  We verify the stabilization A[w] = C^[5] side
(sections 2, 4 identities) symbolically; the non-polynomiality proof (sections 3-6) is read,
not machine-verified."""
from sympy import symbols, expand, simplify, Poly, diff, cancel, factor, together, S, Matrix

p, s, u, F, J, w = symbols('p s u F J w')
x = s**2 + u**3 + p**2*F
y = s + x*(x - u**3)
z = s*x + p**2*J
H = x**2*F - (1 + 2*s*x)*J - p**2*J**2 - p*u

ok = lambda e: expand(e) == 0
print("(con:identity) xy - z(z+1) = p^2 (H + p u):", ok(x*y - z*(z+1) - p**2*(H + p*u)))

# Delta = -p^2 d/du in coordinates (p,x,y,z,u); compute its values on s,u,F,J via the inverse formulas
X, Y, Z, U = symbols('X Y Z U')
s_of = Y - X*(X - U**3)
F_of = (X - s_of**2 - U**3)/p**2
J_of = (Z - s_of*X)/p**2
Dl = lambda e: -p**2*diff(e, U)
sub = {X: x, Y: y, Z: z, U: u}
Ds = expand(Dl(s_of).subs(sub)); DF = expand(Dl(F_of).subs(sub)); DJ = expand(Dl(J_of).subs(sub))
print("Delta s = -3p^2 x u^2:", ok(Ds - (-3*p**2*x*u**2)))
print("Delta F = (6sx+3)u^2:", ok(DF - (6*s*x + 3)*u**2))
print("Delta J = 3x^2u^2:", ok(DJ - 3*x**2*u**2))
# Delta as a derivation of P = C[p,s,u,F,J]
def Delta(e):
    return expand(diff(e, s)*(-3*p**2*x*u**2) + diff(e, u)*(-p**2) + diff(e, F)*((6*s*x + 3)*u**2)
                  + diff(e, J)*(3*x**2*u**2))
print("Delta H = p^3:", ok(Delta(H) - p**3))
print("Delta^2 H = 0:", ok(Delta(Delta(H))))
# local nilpotence on generators (finite iterates)
for name, g in [("s", s), ("u", u), ("F", F), ("J", J)]:
    e, k = g, 0
    while expand(e) != 0 and k < 30:
        e = Delta(e); k += 1
    print(f"  Delta^{k} {name} = 0" if expand(e) == 0 else f"  {name}: not nilpotent within 30")

# linear change L, M (determinant 1) and expansion
x0 = s**2 + u**3
L_, M_ = symbols('L M')
L = x0**2*F - (1 + 2*s*x0)*J
M = (1 - 2*s*x0)*F + 4*s**2*J
print("det of (F,J)->(L,M) = 1:", ok(Matrix([[x0**2, -(1+2*s*x0)], [1-2*s*x0, 4*s**2]]).det() - 1))
FL = 4*s**2*L_ + (1 + 2*s*x0)*M_
JL = -(1 - 2*s*x0)*L_ + x0**2*M_
print("inverse formulas:", ok(L.subs({F: FL, J: JL}) - L_) and ok(M.subs({F: FL, J: JL}) - M_))
HL = expand(H.subs({F: FL, J: JL}, simultaneous=True))
xL = x0 + p**2*FL
Q = 2*x0*FL**2 - 2*s*FL*JL - JL**2
print("(con:h-expansion) H(L) = L - p u + p^2 Q(L) + p^4 F(L)^3:", ok(HL - (L_ - p*u + p**2*Q + p**4*FL**3)))
C = expand(Q.subs(L_, 0))
print("C = Q(0) = M^2(2x0+6sx0^2+4s^2x0^3-x0^4):", ok(C - M_**2*(2*x0 + 6*s*x0**2 + 4*s**2*x0**3 - x0**4)))
Q1 = cancel((Q - C)/L_)
print("Q1 polynomial:", Poly(expand(Q1), L_, M_, s, u, p).is_polynomial if hasattr(Poly, 'is_polynomial') else expand(Q1*L_ - (Q - C)) == 0)
Ls = p*u - p**2*C
HLs = expand(HL.subs(L_, Ls))
r = cancel(HLs/p**3)
print("(con:root-divisibility) H(L*) in p^3 B:", ok(HLs - p**3*expand(r)) and Poly(expand(r), p, s, u, M_).as_expr() == expand(r))
h = u - p*Q - p**3*FL**3 - p**2*w
e0 = -w - p*FL**3 - Q1*h
print("H(L) + p^3 w = L - p h:", ok(HL + p**3*w - (L_ - p*h)))
print("(con:e-certificate) p^3 e0 - (L - L*) = (p^2 Q1 - 1)(H(L) + p^3 w):",
      ok(p**3*e0 - (L_ - Ls) - (p**2*Q1 - 1)*(HL + p**3*w)))
Zs = symbols('Zs')
Wpoly = cancel(-HL.subs(L_, Ls + p**3*Zs)/p**3)
print("W(Z) = -H(L*+p^3 Z)/p^3 is a polynomial:", expand(Wpoly*p**3 + HL.subs(L_, Ls + p**3*Zs)) == 0 and
      all(d >= 0 for d in Poly(expand(Wpoly), p, Zs, s, u, M_).monoms()[0]))

# Section 4: principalization in R~ = C[a,d,b,c,u]/(ac-bd-1)
a, d, b, c = symbols('a d b c')
xr, yr, zr = a*b, d*c, d*b
sR = yr - xr*(xr - u**3)
f = xr - sR**2 - u**3
g = zr - sR*xr
v = a**3*b - a**2*u**3 - d**2
hh = xr - u**3
C1 = c**2 - b**2*hh
al = a**2*(1 - 2*b*d)
be = d**2*(3 + 2*b*d) + al*hh
rel = a*c - b*d - 1
def zero_mod_rel(e):
    # substitute c = (1+bd)/a and clear a-denominators
    return simplify(together(expand(e).subs(c, (1 + b*d)/a))) == 0
print("xy = z(z+1) in R~:", zero_mod_rel(xr*yr - zr*(zr + 1)))
print("f = C1 v:", zero_mod_rel(f - C1*v))
print("g = b^2 v:", zero_mod_rel(g - b**2*v))
print("v = alpha f + beta g:", zero_mod_rel(v - al*f - be*g))
print("alpha C1 + beta b^2 = 1:", zero_mod_rel(al*C1 + be*b**2 - 1))
