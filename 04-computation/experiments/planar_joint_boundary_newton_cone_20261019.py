#!/usr/bin/env python3
"""Exact finite controls for the joint-boundary Newton completion."""
import json
import math
import sympy as sp
x,t,s,L,z=sp.symbols('x t s L z')
q=1+x*x*t
checks=0

def check(ok,label):
    global checks
    if not bool(ok): raise RuntimeError(label)
    checks+=1

def clean(f):return sp.cancel(sp.expand(f))
def zero(f,label):check(clean(f)==0,label)
def weight(f):
    f=sp.expand(f)
    if f==0:return math.inf
    return min(2*j-i for (i,j),c in sp.Poly(f,x,t).terms() if c!=0)
def torder(f):
    if f==0:return math.inf
    return min(j for (i,j),c in sp.Poly(sp.expand(f),x,t).terms() if c!=0)
def ev(f):return sp.cancel(f.subs({x:-1/s,t:-s*s},simultaneous=True))
def delta(H,f):return sp.expand(sp.diff(f,x)*sp.diff(H,t)-sp.diff(f,t)*sp.diff(H,x))
def bd(v,f):return clean(v*sp.diff(f,s))
def partial(op,f,time=1,n=2):
    out=f;term=f
    for j in range(1,n+1):
        term=op(term);out+=sp.Rational(time)**j*term/sp.factorial(j)
    return clean(out)

Rs=[x**i*t**j for j in range(3) for i in range(2*j+3)]
Rs += [1+x+x*x,x*x*q,x**4*t+x*t*t]
records=[]
inputs=[sp.Integer(1),x,t,x*x,x*t,t*t]
for k,R in enumerate(Rs):
    H=t*t*q*R
    check(weight(R)>=-2 and weight(H)>=2,f'Newton cone {k}')
    v=clean(s**4*ev(R))
    for a in inputs:
        d=delta(H,a)
        check(weight(d)>=weight(a)+1,f'weight growth {k}:{a}')
        check(torder(d)>=torder(a)+1,f'source growth {k}:{a}')
        zero(ev(d)-bd(v,ev(a)),f'boundary transport {k}:{a}')
    zero(ev(delta(H,q)),f'boundary ideal {k}')
    records.append({'R':str(R),'weight_R':weight(R),'weight_H':weight(H),'boundary_velocity':str(v)})

# Different-invariant finite partial operators match their boundary operators.
for wi,word in enumerate([[(x,1),(x*x,1)],[(x*x,1),(x,-1)],[(1+x+x*x,2),(x*x*q,-1)]]):
    for f in [x,t,q]:
        source=f;boundary=ev(f)
        for R,time in word:
            H=t*t*q*R;v=clean(s**4*ev(R))
            source=partial(lambda z:delta(H,z),source,time)
            boundary=partial(lambda z:bd(v,z),boundary,time)
        zero(ev(source)-boundary,f'finite word transport {wi}:{f}')
        if f==q:zero(boundary,f'finite word boundary preservation {wi}')
        else:zero(source.subs(t,0)-(x if f==x else 0),f'finite word line {wi}:{f}')

# Sparse finite controls of the completion: stabilization below weight bounds.
def support_F(N):return {(2*n+1,2*n):1 for n in range(N+1)}
def support_G(N):return {(n,n):(-1)**n for n in range(N+1)}
def product(A,B):
    out={}
    for (i,j),a in A.items():
        for (k,l),b in B.items():
            ij=(i+k,j+l);out[ij]=out.get(ij,0)+a*b
    return {ij:c for ij,c in out.items() if c}
def low(A,M):return {ij:c for ij,c in A.items() if 2*ij[1]-ij[0]<=M}
def derivative(A,axis):
    out={}
    for (i,j),c in A.items():
        power=i if axis==0 else j
        if power:
            ij=(i-1,j) if axis==0 else (i,j-1)
            out[ij]=power*c
    return out
def boundary_dict(A):
    out={}
    for (i,j),c in A.items():
        w=2*j-i;out[w]=out.get(w,0)+(-1)**(i+j)*c
    return {w:c for w,c in out.items() if c}
def convolution(A,B):
    out={}
    for i,a in A.items():
        for j,b in B.items():out[i+j]=out.get(i+j,0)+a*b
    return {w:c for w,c in out.items() if c}
for M in range(-1,10):
    N=M+5
    F=support_F(N);G=support_G(N)
    P=product(F,G)
    check(low(P,M)==low(product(support_F(N+1),support_G(N+1)),M),f'weight truncation stability {M}')
    for axis in [0,1]:
        check(low(derivative(P,axis),M)==low(derivative(product(support_F(N+1),support_G(N+1)),axis),M),f'derivative continuity control {M}:{axis}')
    check(boundary_dict(P)==convolution(boundary_dict(F),boundary_dict(G)),f'finite evaluation multiplicativity {M}')
for i in range(16):
    for j in range(16):
        w=2*j-i
        check(j>=math.ceil(w/2),f'weight implies source order {i},{j}')

# The family strictly extends K[xt,t]: x and x^2 have i>j.
for R in [x,x*x]:
    check(weight(R)>=-2,'strict extension admissible')
    check(any(i>j for (i,j),c in sp.Poly(R,x,t).terms()),'strict extension outside old ring')
# Old family monomial A(s,tau)=s^a tau^b maps to x^a t^(a+b).
for a in range(5):
    for b in range(5):check(a<=2*(a+b)+2,'old cone included')

# Exact common-chart denominator for the new R=x^2 example.
common_H=clean((t*t*q*x*x).subs({x:s/z,t:z},simultaneous=True))
zero(common_H-s*s*(s*s+z)/z,'x squared common Hamiltonian denominator')
common_ds=clean(z*sp.diff(common_H,z))
common_dt=clean(-z*sp.diff(common_H,s))
zero(common_ds+s**4/z,'x squared common s component')
zero(common_dt+4*s**3+2*s*z,'x squared common tau component')
zero(common_ds.subs(z,-s*s)-s*s,'x squared common boundary velocity')

# Edge velocity s^2 has the stated Möbius scalar-time formal flow.
v=s*s;term=s;series=s
for n in range(1,7):
    term=bd(v,term);series+=L**n*term/sp.factorial(n)
zero(series-sum(L**n*s**(n+1) for n in range(7)),'Mobius boundary jet')
zero(s**4*ev(x*x)-s*s,'x squared boundary velocity')
zero(s**4*ev(x)+s**3,'x boundary velocity')

# R=x^3 marks failure of this filtration proof, not a counterexample.
H=t*t*q*x**3
check(weight(H)==1,'outside cone Hamiltonian order')
check(weight(delta(H,t))==weight(t),'outside cone no weight gain')
zero(s**4*ev(x**3)+s,'outside cone linear velocity')
f=s
for n in range(1,7):
    f=bd(-s,f);zero(f-(-1)**n*s,f'outside cone no boundary order gain {n}')
# Missing q allows a genuine polynomial shear.
zero(delta(t*t,x)-2*t,'missing boundary shear')
zero(delta(t*t,delta(t*t,x)),'missing boundary termination')
check(ev(delta(t*t,q))!=0,'missing boundary failure')

print(json.dumps({'status':'PASS','checks':checks,'R_controls':records,
 'input_polynomials_per_R':len(inputs),'finite_word_controls':3,
 'partial_exponential_order':2,'weight_stability_bounds':list(range(-1,10)),
 'scope':'Finite controls of the written universal theorem; no unrestricted Jacobian claim.'},indent=2,sort_keys=True))
