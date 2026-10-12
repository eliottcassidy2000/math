"""Gao's degree-four Keller map G (arXiv:2608.00222, section 3.5) on its coordinate-pair fibrations
(HYP-9255 test): is any generic coordinate-pair fibre disconnected (a union of lines)?
G = (gamma x, (p(w)+2 gamma)/(gamma x), (q(w)+gamma w)/(gamma x)^2), p = w^3-6w^2+6w, q = 3/8 w^4 - 2w^3 + 3/2 w^2,
gamma = 2 - 4xy - x^2 z, u = 1+xy, w = gamma u.  Session mac-mini-2026-10-12-logic4."""
import sympy as sp
x,y,z,c,X,Y,w = sp.symbols('x y z c X Y w')
p = lambda t: t**3 - 6*t**2 + 6*t
q = lambda t: sp.Rational(3,8)*t**4 - 2*t**3 + sp.Rational(3,2)*t**2
gam = 2 - 4*x*y - x**2*z; u = 1 + x*y; W = gam*u
G1 = sp.expand(gam*x); G2 = sp.cancel((p(W)+2*gam)/(gam*x)); G3 = sp.cancel((q(W)+gam*W)/(gam*x)**2)
G = [G1,G2,G3]
print("polynomial components:", [sp.denom(sp.together(g))==1 for g in G], " degrees:", [sp.Poly(g,x,y,z).total_degree() for g in G])
J = sp.Matrix([[sp.diff(g,v) for v in (x,y,z)] for g in G]); print("det JG =", sp.expand(J.det()))
print("q' = (w/2) p' :", sp.expand(sp.diff(q(w),w) - w/2*sp.diff(p(w),w))==0)
# (G1,G2)-fibre over (c,X): gamma = c/x, z = (2x-4x^2y-c)/x^3 ; w = c(1+xy)/x ; p(w)+2c/x = cX => x = 2c/(cX - p(w))
xw = 2*c/(c*X - p(w)); yw = sp.cancel(w/c - 1/xw); zw = sp.cancel((2*xw - 4*xw**2*yw - c)/xw**3)
sub = {x:xw, y:yw, z:zw}
print("\n(G1,G2)-fibre: G1,G2 =", sp.cancel(G1.subs(sub)), sp.cancel(G2.subs(sub)))
G3w = sp.cancel(G3.subs(sub)); print("   G3|fibre =", sp.factor(G3w), " deg_w =", sp.Poly(sp.numer(sp.together(G3w)),w).degree(), " denominator:", sp.denom(sp.together(G3w)))
print("   dG3/dw =", sp.factor(sp.diff(G3w,w)), "  (= (cX - p(w))/(2c^2)?", sp.cancel(sp.diff(G3w,w) - (c*X-p(w))/(2*c**2))==0, ")")
print("   punctures: w = infinity and the 3 roots of cX - p(w); connected (rational, parametrized by w); b_1 = 3")
# (G1,G3)-fibre over (c,Y): q(w) + c w/x = c^2 Y => x = c w/(c^2 Y - q(w))
xw2 = c*w/(c**2*Y - q(w)); yw2 = sp.cancel(w/c - 1/xw2); zw2 = sp.cancel((2*xw2 - 4*xw2**2*yw2 - c)/xw2**3)
sub2 = {x:xw2, y:yw2, z:zw2}
print("\n(G1,G3)-fibre: G1,G3 =", sp.cancel(G1.subs(sub2)), sp.cancel(G3.subs(sub2)))
G2w = sp.cancel(G2.subs(sub2)); print("   G2|fibre =", sp.factor(G2w))
print("   dG2/dw =", sp.factor(sp.diff(G2w,w)))
print("   punctures: w = 0 (x=0), w = infinity, the 4 roots of c^2 Y - q(w); connected; b_1 = 5")
# (G2,G3)-fibre over (X,Y): eliminate c: c = p(w) x/(X x - 2); then q(w)(Xx-2)^2 + w p(w)(Xx-2) = Y p(w)^2 x^2
A = sp.expand(q(w)*X**2 - Y*p(w)**2); B = sp.expand(-4*q(w)*X + w*p(w)*X); C = sp.expand(4*q(w) - 2*w*p(w))
quad = sp.expand(q(w)*(X*x-2)**2 + w*p(w)*(X*x-2) - Y*p(w)**2*x**2)
print("\n(G2,G3)-fibre: quadratic in x: A x^2 + B x + C with A,B,C as computed:", sp.expand(quad - (A*x**2+B*x+C))==0)
disc = sp.expand(B**2 - 4*A*C); fl = sp.factor_list(disc)
print("   disc degree in w:", sp.Poly(disc,w).degree(), " factorization:", fl)
sqf = sp.sqf_part(sp.Poly(disc, w)); print("   squarefree part degree:", sqf.degree(), " => genus of the double cover =", (sqf.degree()-1)//2 if sqf.degree()%2==1 else (sqf.degree()-2)//2)
print("   disc a perfect square in Q(X,Y)[w]?", all(e%2==0 for _,e in fl[1]), " (if yes the generic fibre would be disconnected)")
# where is the fibre punctured? x = infinity at roots of A (deg in w), x = 0 at roots of C, Xx - 2 = 0 at ...
print("   A =", sp.factor(A)); print("   C =", sp.factor(C))
