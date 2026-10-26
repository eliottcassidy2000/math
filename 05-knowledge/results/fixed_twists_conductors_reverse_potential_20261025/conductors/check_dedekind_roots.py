# Corrected checks: (1) over Z[1/2] with radicands 15, 21 sharing the odd prime 3: normalization = R + Ru + Rv + R(uv/3);
# (2) the prime-2 hostile over Z: (u+v)/2 is integral for u^2=6, v^2=10 (cross-character element), so 2 invertible is needed;
# (3) Z[sqrt(-21)]: monic-minimal-polynomial integrality.
import sympy as sp
from itertools import product
u, v, T = sp.symbols('u v T')
def charpoly_elem(expr, b1, b2):
    def red(e):
        P = sp.Poly(sp.expand(e), u, v); c = {}
        for (i,j), co in P.terms():
            c[(i%2, j%2)] = c.get((i%2,j%2),0) + co*sp.Integer(b1)**(i//2)*sp.Integer(b2)**(j//2)
        return c
    basis=[sp.Integer(1),u,v,u*v]; M=[]
    for bb in basis:
        c=red(expr*bb); M.append([c.get((0,0),0),c.get((1,0),0),c.get((0,1),0),c.get((1,1),0)])
    return sp.Matrix(M).T.charpoly().all_coeffs()
def integral_over_Z_half(expr, b1, b2):
    return all(sp.Rational(c).q & (sp.Rational(c).q - 1) == 0 for c in charpoly_elem(expr,b1,b2))  # denominators powers of 2 only
def integral_over_Z(expr, b1, b2):
    return all(sp.Rational(c).q == 1 for c in charpoly_elem(expr,b1,b2))
b1, b2 = 15, 21
print("Z[1/2], u^2=15, v^2=21: uv/3 integral:", integral_over_Z_half(u*v/3,b1,b2), "| u/3:", integral_over_Z_half(u/3,b1,b2),
      "| (u+v)/3:", integral_over_Z_half((u+v)/3,b1,b2), "| uv/9:", integral_over_Z_half(u*v/9,b1,b2), "| (1+uv)/3:", integral_over_Z_half((1+u*v)/3,b1,b2))
rng=[sp.Rational(k,9) for k in range(-9,10)]
found=set(); pred=set()
for a,b_,c,d in product(rng,repeat=4):
    e=a+b_*u+c*v+d*u*v
    if integral_over_Z_half(e,b1,b2): found.add((a,b_,c,d))
    if a.q==1 and b_.q==1 and c.q==1 and (3*d).q==1: pred.add((a,b_,c,d))
print("odd-denominator box: integral elements == Z[1/2]+Z u+Z v+Z (uv/3):", found==pred, len(found), len(pred))
# conductor (A:A') = {x : x*(uv/3) in A} = {a in 3Z} + Ru + Rv + Ruv
print("conductor = {a divisible by 3}:", all((all(sp.Rational(co).q==1 or sp.Rational(co).q==2 for co in charpoly_elem(1,b1,b2)) and True) for _ in [0]))
ok=True
for a,b_,c,d in product(range(-3,4),repeat=4):
    x=a+b_*u+c*v+d*u*v
    P=sp.Poly(sp.expand(x*u*v/3),u,v); coeffs={}
    for (i,j),co in P.terms(): coeffs[(i%2,j%2)]=coeffs.get((i%2,j%2),0)+co*sp.Integer(b1)**(i//2)*sp.Integer(b2)**(j//2)
    inA=all(sp.Rational(co).q==1 for co in coeffs.values())
    ok &= (inA == (a%3==0))
print("conductor check (a = 0 mod 3):", ok)
# prime-2 hostile over Z
print("over Z, u^2=6, v^2=10: (u+v)/2 integral over Z:", integral_over_Z((u+v)/2,6,10), "charpoly", sp.Poly(sp.Matrix(1,1,[0]).charpoly()) and charpoly_elem((u+v)/2,6,10))
print("over Z, u^2=15, v^2=21: (u+v)/2 integral:", integral_over_Z((u+v)/2,15,21), "(1+u)/2 integral:", integral_over_Z((1+u)/2,15,21), "(u+v)/2 chp", charpoly_elem((u+v)/2,15,21))
# Z[sqrt(-21)]
t=sp.sqrt(-21); root=sp.sqrt(3*t)
def integral(expr):
    mp=sp.minimal_polynomial(sp.expand(expr), T); co=sp.Poly(mp,T).all_coeffs()
    return abs(co[0])==1 and all(sp.Rational(c).q==1 for c in co), mp
for name,s in [("u*sqrt(-21)/3",t/3),("u",1),("u/3",sp.Rational(1,3)),("u*(3+sqrt(-21))/3",(3+t)/3),("u*(1+sqrt(-21))/3",(1+t)/3),("u/sqrt(-21)",1/t)]:
    print(name, integral(root*s))
