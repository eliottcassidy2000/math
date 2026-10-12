"""General tangent-sweep family (Gao-Gallagher): F = (gamma x, (p(w)+2gamma)/(gamma x), (q(w)+gamma w)/(gamma x)^2),
w = gamma u, u = 1+xy, q' = (w/2) p', p(0) = 0.  On the (F1,F2)-fibre over (c,X): gamma = c/x and x = 2c/(cX - p(w));
on the (F1,F3)-fibre over (c,Y): x = c w/(c^2 Y - q(w)).  Verify, for generic p of degree d = 2..6:
  F3|fibre = (q(w) + w(cX - p(w))/2)/c^2 has dF3/dw = (cX - p(w))/(2c^2)   (degree d+1 polynomial, critical points = punctures)
  F2|fibre = (p(w) + 2(c^2 Y - q(w))/w)/c  has dF2/dw = -2(c^2 Y - q(w))/(c w^2) (degree d+1 rational function, critical points = punctures)
and that the leading coefficient of F2|fibre at infinity is lead(p)/((d+1)c) (nonzero)."""
import sympy as sp
w,c,X,Y,s = sp.symbols('w c X Y s')
for d in range(2,7):
    ps = sp.symbols('p1:%d'%(d+1))
    p = sum(ps[i-1]*w**i for i in range(1,d+1))
    q = sp.integrate((s/2)*sp.diff(p,w).subs(w,s),(s,0,w))
    F3 = (q + w*(c*X - p)/2)/c**2
    F2 = (p + 2*(c**2*Y - q)/w)/c
    ok3 = sp.expand(sp.diff(F3,w) - (c*X-p)/(2*c**2))==0
    ok2 = sp.cancel(sp.diff(F2,w) + 2*(c**2*Y-q)/(c*w**2))==0
    deg3 = sp.Poly(sp.expand(F3*c**2),w).degree()
    lead2 = sp.Poly(sp.expand((F2*c*w).subs(Y,0)),w).LC()  # leading coeff of w*F2*c
    print("d=%d: dF3/dw identity %s (deg_w F3 = %d); dF2/dw identity %s; leading coefficient of F2 at infinity = %s"%(d,ok3,deg3,ok2,sp.factor(lead2/c)))
