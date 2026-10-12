"""General tangent-sweep family, THM-1300 component order F = (F1, F2, F3) = ((q+gamma w)/(gamma x)^2, (p+2gamma)/(gamma x), gamma x).
(F1,F2)-fibre over (a,b): with c = F3 free, x = 2c/(cb - p(w)), c w/x = c^2 a - q(w) => c = p x/(b x - 2) and
   q (b x - 2)^2 + w p (b x - 2) - a p^2 x^2 = 0, a quadratic in x whose discriminant is p(w)^2 (b^2 w^2 - 8 a (w p - 2 q)).
Checks for generic p of degree d = 2..6: the discriminant identity; w p - 2 q = P(w) := int_0^w p; squarefreeness of
b^2 - 8 a P(w)/w^2 generically (so the double cover t^2 = b^2 - 8aP/w^2 is irreducible of genus floor((d-2)/2));
and for Gao's degree-four map: the critical values of F1 on the (F2,F3)-fibre are pairwise distinct => monodromy S_4."""
import sympy as sp
w,s,x,a,b,c = sp.symbols('w s x a b c')
for d in range(2,7):
    ps = sp.symbols('p1:%d'%(d+1)); p = sum(ps[i-1]*w**i for i in range(1,d+1))
    q = sp.integrate((s/2)*sp.diff(p,w).subs(w,s),(s,0,w)); P = sp.integrate(p.subs(w,s),(s,0,w))
    quad = sp.expand(q*(b*x-2)**2 + w*p*(b*x-2) - a*p**2*x**2)
    A,B,C = [sp.expand(quad.coeff(x,k)) for k in (2,1,0)]
    disc = sp.expand(B**2 - 4*A*C)
    ok1 = sp.expand(disc - p**2*(b**2*w**2 - 8*a*(w*p-2*q)))==0
    ok2 = sp.expand(w*p - 2*q - P)==0
    h = sp.expand(sp.cancel((b**2 - 8*a*P/w**2)))
    D = sp.discriminant(h, w)
    print("d=%d: disc = p^2 (b^2 w^2 - 8a(wp-2q)): %s ; wp - 2q = int_0^w p: %s ; deg_w(b^2 - 8aP/w^2) = %d ; its w-discriminant is a nonzero polynomial in (a,b,p_i): %s => generic genus %d"%(
        d, ok1, ok2, sp.Poly(h,w).degree(), sp.expand(D)!=0, (sp.Poly(h,w).degree()-1)//2))
# Gao's G (d=3): F1|(F2,F3)-fibre = (q(w) + w(cb - p(w))/2)/c^2 with p = w^3-6w^2+6w ; critical points = roots of cb - p
p3 = w**3 - 6*w**2 + 6*w; q3 = sp.integrate((s/2)*sp.diff(p3,w).subs(w,s),(s,0,w))
F1w = (q3 + w*(c*b - p3)/2)/c**2
r = sp.symbols('r')
# critical values: F1w at roots r of c b - p(r) = 0; two critical values coincide iff resultant condition; compute the polynomial whose roots are the critical values
# use b = p(r)/c to express: critical value as function of r: V(r) = (q3(r) + r(p3(r) - p3(r))/2)/c^2 = q3(r)/c^2 ; so distinct critical points r1 != r2 with p(r1) = p(r2) (= cb) give equal critical values iff q3(r1) = q3(r2)
# => check resultant of (p(r1)-p(r2))/(r1-r2) and (q3(r1)-q3(r2))/(r1-r2) in r2 is a nonzero polynomial in r1 (then for generic b no coincidence)
r1,r2 = sp.symbols('r1 r2')
e1 = sp.cancel((p3.subs(w,r1)-p3.subs(w,r2))/(r1-r2)); e2 = sp.cancel((q3.subs(w,r1)-q3.subs(w,r2))/(r1-r2))
res = sp.resultant(e1, e2, r2)
print("Gao G: resultant of (p(r1)-p(r2))/(r1-r2) and (q(r1)-q(r2))/(r1-r2) in r2 =", sp.factor(res), "=> critical values generically distinct:", sp.expand(res)!=0)
print("   (critical value at a critical point r is q(r)/c^2; two coincide iff p(r1)=p(r2) and q(r1)=q(r2))")
