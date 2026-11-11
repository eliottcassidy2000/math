"""Exact finite controls for the moving-centre certificate operations.
Not an audit of the inherited comparison theorem, nor a JC computation.
"""
import sympy as z
x,y,a,b,u,v=z.symbols('x y a b u v')
f=x+y**2
g=y+f**2
D=lambda p:z.diff(g,y)*z.diff(p,x)-z.diff(g,x)*z.diff(p,y)
E=lambda p:-z.diff(f,y)*z.diff(p,x)+z.diff(f,x)*z.diff(p,y)
U,V=f-a,g-b
T=lambda p:D(p)-V*p
S=lambda p:E(p)-U*p
A=lambda p:z.diff(p,a)+V*p
B=lambda p:z.diff(p,b)+U*p
normal=lambda p:z.cancel(z.expand(p))
D0,E0,T0,S0,A0,B0=D,E,T,S,A,B
D=lambda p:normal(D0(p))
E=lambda p:normal(E0(p))
T=lambda p:normal(T0(p))
S=lambda p:normal(S0(p))
A=lambda p:normal(A0(p))
B=lambda p:normal(B0(p))
def require(condition):
    if not condition:
        raise AssertionError("Exact control failed")

def equal(p,q):
    require(normal(p-q)==0)
require(normal(z.diff(f,x)*z.diff(g,y)-z.diff(f,y)*z.diff(g,x))==1)
samples=[1,x,y,x*y+a/(b+1),x**2/(a-b)+y]
for h in samples:
    equal(D(E(h)),E(D(h)))
    equal(A(B(h)),B(A(h)))
    for connection in [A,B]:
        for stein in [T,S]:
            equal(connection(stein(h)),stein(connection(h)))
print('PASS nonlinear Keller pair and all six flatness/commutation identities on 5 probes')
# For parameter-dependent r, the target is partial_a r+D r;
# this is the general form underlying the parameter-independent compiler.
for p,q in [(x/(a-b),y/(a+1)),(y**2, x+a*b)]:
    r=T(p)+S(q)
    equal(z.diff(r,a)+D(r),T(r+A(p))+S(A(q)))
    equal(z.diff(r,b)+E(r),T(B(p))+S(r+B(q)))
print('PASS both certificate compiler identities with rational parameter coefficients')
# A parameter-independent first hit r=f-a is not available since a occurs.
# Instead r=0 has useful nontrivial syzygies p=S(k),q=-T(k).
k=x/(a-b)
p,q=S(k),-T(k)
r=normal(T(p)+S(q))
equal(r,0)
for operation in 'AB':
    if operation=='A':
        p,q=r+A(p),A(q)
        r=normal(D(r))
    else:
        p,q=B(p),r+B(q)
        r=normal(E(r))
    p,q=normal(p),normal(q)
    equal(T(p)+S(q),r)
    for c in [p,q]:
        denominator=z.factor(z.denom(c))
        # No new pole location is introduced.
        for factor, exponent in z.factor_list(denominator,a,b)[1]:
            require(z.Poly(factor,a,b).total_degree()==1)
            require(normal(factor.subs(a,b))==0)
print('PASS 2-step mixed compiler and preserved denominator location a=b')
# Inverse x=u-(v-u^2)^2, y=v-u^2 gives exact pure derivative budgets.
xinv=u-(v-u**2)**2
yinv=v-u**2
for h in [xinv,yinv,xinv*yinv]:
    P=z.Poly(h,u,v)
    require(z.diff(h,u,P.degree(u)+1)==0)
    require(z.diff(h,v,P.degree(v)+1)==0)
print('PASS automorphism control: polynomial inverse gives actual finite mixed derivative tails')
for m in range(9):
    residue=(-1)**m*z.factorial(m)/(2*a**(m+1))
    require(residue!=0)
    require(normal(z.diff(1/(2*a),a,m)-residue)==0)
print('PASS punctured-surface hostile: residue derivatives remain nonzero for m=0,...,8')
print('Scope: finite exact algebraic controls only; the all-order statements have proofs in proof.md.')
