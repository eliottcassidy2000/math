#!/usr/bin/env python3
"""
paley_starstar_bivariate_eqn_monad.py
monad-explorer-2026-06-07 (deep-research, 11th session)

Functional equation for G(x,y) = 1 + U(x,y),  U = sum_{k,m>=1} t(k,m) x^k y^m.
G(x,-1) = F(x) = (sqrt(1+4x)-1)/(2x),  x F^2 + F - 1 = 0  (the (**) loop eqn).
Both THM-438 handoffs are Taylor coeffs of V(s,y)=U(x,y) (s=x/(1-x)) at s=-1:
   #1 V(-1,y)=-y ;  #2 V_s(-1,y)=y/((1-y)(1-2y)).
Columns known EXACTLY (all k) for m=1..4; partial m=5,6.
"""
import sympy as sp

x, y, s = sp.symbols('x y s')
c2, c3, c4 = sp.symbols('c2 c3 c4')

T = {}
T[1] = x/(1-x)
T[2] = 3*x**2/(1-x)**3
T[3] = (13+7*x)*x**3/(1-x)**5
T[4] = (69+97*x+15*x**2)*x**4/(1-x)**7
P5 = 421+1056*x+c2*x**2+c3*x**3+c4*x**4
T[5] = P5*x**5/(1-x)**9
T[6] = sp.Integer(2867)*x**6/(1-x)**11   # only constant of P6 known

U = sum(T[m]*y**m for m in range(1,7))
G = 1 + U

print("="*72); print("PART 1: confirm G(x,-1)=F, loop eqn"); print("="*72)
F = (sp.sqrt(1+4*x)-1)/(2*x)
chk = sp.series(G.subs(y,-1) - F, x, 0, 7).removeO()
print("G(x,-1)-F up to x^6:", sp.simplify(chk))

def yser_coeffs(expr, ny, nx):
    """return list of x-series for [y^0..y^ny]."""
    e = sp.expand(expr)
    out = []
    for j in range(ny+1):
        cj = e.coeff(y, j)
        out.append(sp.series(cj, x, 0, nx).removeO())
    return out

print("="*72); print("PART 2: E := x G^2 + G - 1 ; divide by (1+y); is quotient nice?"); print("="*72)
E = x*G**2 + G - 1
Ej = yser_coeffs(E, 6, 7)
for j,c in enumerate(Ej):
    print(f"  [y^{j}] E:", sp.factor(c))
# quotient R = E/(1+y) as power series in y: R_j = sum_{i<=j} (-1)^(j-i) E_i
print("  --- quotient R=E/(1+y) coeffs ---")
R = []
for j in range(len(Ej)):
    rj = sum((-1)**(j-i)*Ej[i] for i in range(j+1))
    rj = sp.series(rj, x, 0, 7).removeO()
    R.append(rj)
    print(f"  [y^{j}] R:", sp.factor(rj))

print("="*72); print("PART 3: try x*y*G^2 + G - 1 and variants"); print("="*72)
for name, E2 in [("xyG^2+G-1", x*y*G**2+G-1),
                 ("xG^2+G-1-x", x*G**2+G-1-x),
                 ("xG^2 + G(1) ...", None)]:
    if E2 is None: continue
    Ej2 = yser_coeffs(E2, 6, 6)
    print(f"  {name}: [y^j] =", [sp.factor(c) for c in Ej2])

print("="*72); print("PART 4: Q_m(s) handoffs (fixed cancellation)"); print("="*72)
def Qm(m, Pm):
    u = s/(1+s)
    return sp.cancel(s**m*(1+s)**(m-1)*Pm.subs(x,u))
Ps = {1:sp.Integer(1),2:sp.Integer(3),3:13+7*x,4:69+97*x+15*x**2,5:P5}
for m in range(1,6):
    Q = sp.expand(Qm(m, Ps[m]))
    print(f"  m={m}: Q={Q}")
    print(f"        Q(-1)={sp.expand(Q.subs(s,-1))}  Q'(-1)={sp.expand(sp.diff(Q,s).subs(s,-1))}")
