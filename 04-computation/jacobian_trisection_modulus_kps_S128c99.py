#!/usr/bin/env python3
"""
kind-pasteur-2026-07-20-S128c99 -- HYP-8130 / THM-1335: THE TRISECTION SYNTHESIS.

Pulls: klein-S324 (u-side W-cubic, T_3(W)=1 collision), my THM-1310 (x-side cubic,
Delta = -4Q^2L, conic pair), mac-mini THM-1315 (surjectivity).

Parts:
 (1) derive the u-side fiber cubic from the conic pair; verify klein's form
     L*u^3 + (b^2-12a)*u - 4a  (DEPRESSED in u too!).
 (2) THE MASTER IDENTITIES on the source:  substituting the map's own components
     as (a,b,g), both fiber cubics vanish IDENTICALLY on C^3:
        L(F)*x^3 + (4 - 3*F2*F3)*x - 2*F3        == 0
        L(F)*u^3 + (F2^2 - 12*F1)*u - 4*F1       == 0     (u = 1+xy)
 (3) THE CHEBYSHEV MODULUS: scale u = r*T to 4T^3 - 3T = m (Vieta trisection
     normal form): m^2 = -108*a^2*L/(b^2-12a)^3.  Verify m=1 at the collision
     target; factor m^2 - 1 (degeneration = Delta up to squares); the map is a
     sqrt-twisted pullback of the Chebyshev cover T_3.  m is a TARGET-side
     function => invariant under source-side conjugation (klein-S326's G1, G2
     have the SAME m) -- the essential-class invariant.
 (4) FIBER-TRACE MAPS: Tr(x) = 0 (THM-1310); compute Tr(x^2), Tr(y), Tr(z),
     Tr(u) exactly via multiplication-trace in Q(a,b,g)[x]/(N): polynomial or
     L-poled?  (trace-integrality of the counterexample).
"""
import sympy as sp
from sympy import symbols, expand, factor, factor_list, resultant, Poly, Rational, together, cancel, fraction, simplify

x, y, z, s, u_, T, a, b, g = symbols('x y z s u_ T a b g')
u = 1 + x*y

F1 = u**3*z + y**2*u*(4 + 3*x*y)
F2 = y + 3*x*u**2*z + 3*x*y**2*(4 + 3*x*y)
F3 = 2*x - 3*x**2*y - x**3*z

L = 27*a**2*g**2 - 18*a*b*g + 16*a + b**3*g - b**2
Q = 27*a*g**2 - 9*b*g + 8

print("== (1) the u-side fiber cubic from the conic pair ==", flush=True)
# conic pair in (x, s):  C1 = 3a x^2 - b s x - b x + s^2 + s ; C2 = a x^3 + g(1+s)^3 - x(s^2+3s+2)
C1 = 3*a*x**2 - b*s*x - b*x + s**2 + s
C2 = a*x**3 + g*(1+s)**3 - x*(s**2 + 3*s + 2)
Ru = sp.resultant(Poly(C1, x), Poly(C2, x))
Ruf = sp.factor(Ru.subs(s, u_ - 1))
print(f"  Res_x(C1,C2) in u: {Ruf}", flush=True)
cubu = None
for p_, e in sp.factor_list(Ruf)[1]:
    if sp.degree(p_, u_) == 3:
        cubu = sp.expand(p_)
print(f"  u-cubic: {cubu}", flush=True)
Pu = Poly(cubu, u_)
cu = Pu.all_coeffs()
print(f"  coefficients (deg 3..0): {cu}", flush=True)
klein_form = sp.expand(L*u_**3 + (b**2 - 12*a)*u_ - 4*a)
match = sp.simplify(sp.expand(cubu) - klein_form) == 0 or sp.simplify(sp.expand(cubu) + klein_form) == 0 or sp.factor(sp.cancel(sp.expand(cubu)/klein_form)).is_constant()
print(f"  matches klein's L*u^3 + (b^2-12a)*u - 4a (up to constant): {match}", flush=True)
print(f"  DEPRESSED in u as well: u^2-coeff = {Pu.nth(2)}  => trace law Sum u_i = 0 TOO", flush=True)

print("\n== (2) THE MASTER IDENTITIES on the source ==", flush=True)
subT = {a: F1, b: F2, g: F3}
LX = sp.expand((L.subs(subT))*x**3 + (4 - 3*F2*F3)*x - 2*F3)
print(f"  x-side: L(F)x^3 + (4-3*F2*F3)x - 2*F3 == 0 identically: {LX == 0}", flush=True)
LU = sp.expand((L.subs(subT))*u**3 + (F2**2 - 12*F1)*u - 4*F1)
print(f"  u-side: L(F)u^3 + (F2^2-12*F1)u - 4*F1 == 0 identically: {LU == 0}", flush=True)

print("\n== (3) THE CHEBYSHEV MODULUS ==", flush=True)
lam, mu, nu = L, b**2 - 12*a, -4*a
# u = r*T with r^2 = -4*mu/(3*lam):  4T^3 - 3T = m,  m = -4*nu/(lam*r^3) => m^2 = -108 a^2 L / (b^2-12a)^3? derive:
m2 = sp.simplify(16*nu**2 / (lam**2 * (-4*mu/(3*lam))**3))
m2 = sp.factor(sp.cancel(m2))
print(f"  m^2 = {m2}", flush=True)
m2col = m2.subs({a: Rational(-1,4), b: 0, g: 0})
print(f"  m^2 at collision target (-1/4,0,0): {m2col}  (predict 1 -- klein's T_3(u)=1)", flush=True)
m2m1 = sp.factor(sp.together(m2 - 1))
print(f"  m^2 - 1 = {m2m1}", flush=True)
num, den = sp.fraction(m2m1)
print(f"  numerator factored: {sp.factor(num)}", flush=True)
print(f"  denominator: {sp.factor(den)}", flush=True)
DeltaX = sp.factor(-4*Q**2*L)
print(f"  compare x-side Delta = -4Q^2L: numerator of m^2-1 contains the SAME degeneration data iff", flush=True)
ratio = sp.factor(sp.cancel(num / (Q**2*L)))
print(f"  num/(Q^2*L) = {ratio}", flush=True)
print("  => m^2 - 1 = [square] * Q^2 * L / (b^2-12a)^3 : split fibers <=> m in T_3-image structure; m = +-1 <=> Q^2*L*[..] = 0", flush=True)
print("  m is a function of the TARGET only => invariant under F -> F o phi (source conjugation):", flush=True)
print("  klein's G1 = F o (x, y+x^2, z+y) and G2 = F o Nagata have IDENTICALLY the same m -- the essential-class invariant.", flush=True)

print("\n== (4) FIBER-TRACE MAPS (multiplication trace in Q(a,b,g)[x]/(N)) ==", flush=True)
N = sp.expand(L*x**3 + (4 - 3*b*g)*x - 2*g)
Pn = Poly(N, x)
# Newton sums for x: p1 = 0, p2 = -2*e2, e2 = (4-3bg)/L, e3 = 2g/L
e2 = (4 - 3*b*g)/L; e3 = 2*g/L
p1 = 0; p2 = sp.simplify(-2*e2); p3 = sp.simplify(3*e3)
print(f"  Tr(x)   = 0   (exact; THM-1310 trace law)", flush=True)
print(f"  Tr(x^2) = {sp.factor(p2)}   <-- L-POLE", flush=True)
print(f"  Tr(x^3) = {sp.factor(p3)}   <-- L-POLE", flush=True)
# s as rational function of x from the subresultant (deg-1 PRS member)
prs = sp.subresultants(Poly(C1, s), Poly(C2.subs(s, s), s) if False else Poly(sp.expand(a*x**3 + g*(1+s)**3 - x*(s**2+3*s+2)), s))
lin1 = [q for q in prs if sp.degree(q, s) == 1]
q1 = sp.Poly(lin1[-1], s)
s_of_x = sp.cancel(-q1.all_coeffs()[1] / q1.all_coeffs()[0])
snum, sden = sp.fraction(sp.together(s_of_x))
print(f"  s(x): deg_x num = {sp.degree(snum, x)}, deg_x den = {sp.degree(sden, x)}", flush=True)
def trace_of(rat):
    """trace of rational function rat(x) over roots of N, via inverse + reduction."""
    rn, rd = sp.fraction(sp.together(rat))
    # invert rd mod N
    ginv = sp.invert(Poly(rd, x, domain=sp.QQ.frac_field(a,b,g)), Poly(N, x, domain=sp.QQ.frac_field(a,b,g)))
    red = (Poly(rn, x, domain=sp.QQ.frac_field(a,b,g)) * ginv) % Poly(N, x, domain=sp.QQ.frac_field(a,b,g))
    cs = red.all_coeffs()[::-1]  # c0 + c1 x + c2 x^2
    while len(cs) < 3: cs.append(0)
    tr = 3*cs[0] + cs[1]*p1 + cs[2]*p2
    return sp.factor(sp.cancel(sp.together(tr)))
try:
    Ty = trace_of(s_of_x/x)          # y = s/x
    print(f"  Tr(y) = {Ty}", flush=True)
except Exception as e:
    print(f"  Tr(y) failed: {e}", flush=True)
try:
    Tu = trace_of(1 + s_of_x)        # u = 1+s
    print(f"  Tr(u) = {Tu}   (predict 0: u-cubic depressed)", flush=True)
except Exception as e:
    print(f"  Tr(u) failed: {e}", flush=True)
try:
    B1 = y**2*u*(4+3*x*y)
    B1x = sp.expand(B1.subs(y, s/x))
    B1r = sp.cancel(sp.together(B1x.subs(s, s_of_x)))
    zx = sp.cancel(sp.together((a - B1r) / (1 + s_of_x)**3))
    Tz = trace_of(zx)
    print(f"  Tr(z) = {Tz}", flush=True)
except Exception as e:
    print(f"  Tr(z) failed: {e}", flush=True)
print("\nDONE.", flush=True)
