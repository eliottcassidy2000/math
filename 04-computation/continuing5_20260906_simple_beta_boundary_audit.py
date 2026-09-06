"""Independent rational-series and algebraic-number audit of simple beta roots.
No producer imports. The analytic argument supplies the all-model quantifiers.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import permutations
import sympy as S
import sys
sys.stdout.reconfigure(newline='\n')
GATES=0
def need(ok,label):
    global GATES
    GATES+=1
    if not ok:raise ArithmeticError(label)
v,u,x,y,z,r=S.symbols('v u x y z r')
B=v**5-13*v**4+55*v**3-x*v**2+y*v-z
C=v**4-12*v**3+45*v**2-S.Rational(2,3)*x*v+S.Rational(3,7)*y
D=v**4-11*v**3+36*v**2-S.Rational(5,12)*x*v+y/7
den=1-13*u+55*u**2-x*u**3+y*u**4-z*u**5
# A finite geometric inverse, rather than the forward moment recurrence.
inverse=S.Poly(sum((1-den)**k for k in range(6)),u)
inverse=sum(inverse.nth(k)*u**k for k in range(6))
def moments(A):
    rev=S.Poly(A,v)
    num=sum(rev.nth(4-k)*u**k for k in range(5))
    P=S.Poly(S.expand(num*inverse),u)
    return [S.expand(P.nth(k)) for k in range(6)]
mu,nu=moments(C),moments(D)
def determinant(a):
    n=len(a);out=0
    for p in permutations(range(n)):
        parity=sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
        term=(-1)**parity
        for i,j in enumerate(p):term*=a[i][j]
        out+=term
    return S.expand(out)
dC2=determinant([[mu[i+j+1] for j in range(2)] for i in range(2)])
dC3=determinant([[mu[i+j+1] for j in range(3)] for i in range(3)])
dD3=determinant([[nu[i+j] for j in range(3)] for i in range(3)])
need(S.expand(dC2-(x-75)/3)==0,'shifted C variance identity')
need(S.expand(dC3-(-x*(x-75)**2/27+S.Rational(15,7)*(x-75)*y-S.Rational(16,49)*y**2+(x-75)*z/3))==0,'full shifted C3 identity')
need(S.expand(dD3-(-S.Rational(49,144)*x**2+S.Rational(269,4)*x-S.Rational(18,7)*y-3132))==0,'full ordinary D3 identity')
need(dD3.subs({x:75,y:0,z:0})==-S.Rational(37,16),'double-zero boundary contradiction')

# Solve the two common-node equations independently as a linear system.
solution=S.solve([C.subs(v,r),D.subs(v,r)],(x,y))
need(len(solution)==2,'unique x,y at nonzero common node')
xr=S.Rational(24,7)*r**3-36*r**2+108*r
yr=3*r**4-28*r**3+63*r**2
need(S.expand(solution[x]-xr)==0 and S.expand(solution[y]-yr)==0,'common-node coefficient map')
derivative=S.factor(S.diff(B,v).subs(v,r).subs(solution))
need(S.expand(derivative-S.Rational(4,7)*r**2*(2*r**2-14*r+21))==0,'double-root equation with zero case retained')
rel=2*r**2-14*r+21
remainder=lambda expr:S.rem(S.Poly(expr,r),S.Poly(rel,r)).as_expr()
need(S.expand(remainder(xr)-(126-12*r))==0,'root-field x reduction')
need(S.expand(remainder(yr)-(S.Rational(735,4)-49*r))==0,'root-field y reduction')
need(S.expand(remainder(dD3.subs({x:xr,y:yr}))-(5*r-S.Rational(75,4)))==0,'native D residue norm on repeated branch')
lo=(7-S.sqrt(7))/2;hi=(7+S.sqrt(7))/2
need(S.simplify(yr.subs(r,hi))==-S.Rational(49,4)*(2*S.sqrt(7)-1),'upper branch y exact expression')
need(bool(S.simplify(yr.subs(r,hi))<0),'upper branch violates nonnegative elementary y')
need(S.simplify((5*r-S.Rational(75,4)).subs(r,lo))==-S.Rational(5,4)-S.Rational(5,2)*S.sqrt(7),'lower branch norm exact expression')
need(bool((-S.Rational(5,4)-S.Rational(5,2)*S.sqrt(7))<0),'lower branch violates PSD')

G=v**4-13*v**3+55*v**2-x*v+y
Ct=v**3-S.Rational(45,4)*v**2+S.Rational(75,2)*v-S.Rational(5,12)*x
Dt=v**3-S.Rational(32,3)*v**2+S.Rational(197,6)*v-S.Rational(23,72)*x
need(S.expand(C-S.Rational(3,7)*G-S.Rational(4,7)*v*Ct)==0,'C zero-atom peeling identity')
need(S.expand(D-S.Rational(1,7)*G-S.Rational(6,7)*v*Dt)==0,'D zero-atom peeling identity')
need(S.Rational(3,7)+S.Rational(4,7)==1 and S.Rational(1,7)+S.Rational(6,7)==1,'positive normalized masses')

# Native one-channel hostile, with its actual root lists and failed moment.
host={x:75,y:0,z:0}
need(S.expand(B.subs(host)-v**2*(v-3)*(v-5)**2)==0,'hostile B exact factorization')
need(S.expand(C.subs(host)-v*(v-2)*(v-5)**2)==0,'hostile C exact factorization')
bb=(0,0,3,5,5);cc=(0,2,5,5)
need(all(bb[i]<=cc[i]<=bb[i+1] for i in range(4)),'one-channel weak interlacing hostile')
need(dD3.subs(host)<0,'other channel rejects hostile')

# Nonvacuity is independently reconstructed from exact central sign brackets.
def variation(seq):
    signs=[S.sign(a) for a in seq if a!=0]
    return sum(signs[i]!=signs[i+1] for i in range(len(signs)-1))
def sturm_count(P,a,b):
    chain=S.sturm(P,v)
    return variation([q.subs(v,a) for q in chain])-variation([q.subs(v,b) for q in chain])
for zz in (0,1):
    center={x:84,y:35,z:zz}
    need(sturm_count(B.subs(center),0,8)==(4 if zz==0 else 5),'center positive B roots with zero boundary separated')
    need(S.gcd(B.subs(center),S.diff(B,v).subs(center))==1,'center B has no repeated root')
    for A in (C,D):
        need(sturm_count(A.subs(center),0,8)==4,'actual center interlacer real roots')
        # Reconstruct all nine center moments by the geometric-series method.
        num=sum(S.Poly(A,v).nth(4-k)*u**k for k in range(5))
        inv=sum((1-den.subs(center))**k for k in range(9))
        poly=S.Poly(S.expand(num.subs(center)*inv),u)
        vals=[poly.nth(k) for k in range(9)]
        for size in range(1,6):
            need(determinant([[vals[i+j] for j in range(size)] for i in range(size)])>0,'nonvacuous full positive moment reference')
print('All-model repeated-root alternatives: upper branch y<0; lower branch D/B H3<0; zero branch C-shifted then D contradiction.')
print('Exact zero-atom masses C3/7,D1/7; quotient geometry and original carriers kept distinct.')
print('One-channel repeated-root hostile and genuine full center checked.')
print('PASS',GATES,'always-active exact gates')
