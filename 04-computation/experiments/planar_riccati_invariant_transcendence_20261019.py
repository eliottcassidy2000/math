#!/usr/bin/env python3
"""Exact finite controls for the formal Riccati/transcendence proof."""
import json
import math
import sympy as sp
x,t,s,e,z,L,Y=sp.symbols('x t s e z L Y')
checks=0

def check(ok,label):
    global checks
    if not bool(ok):raise RuntimeError(label)
    checks+=1

def clean(f):return sp.cancel(sp.expand(f))
def zero(f,label):check(clean(f)==0,label)
def trunc(f,N=10):
    f=sp.Poly(sp.expand(f),s,e)
    return sp.Add(*(c*s**i*e**j for (i,j),c in f.terms() if i+2*j<=N))
def delta(f,g):return sp.expand(sp.diff(f,e)*(s*s-e)*sp.diff(g,s))
def inverse_unit(D,N):
    d=trunc(D-1,N);term=sp.Integer(1);out=term
    for _ in range(N):term=trunc(-term*d,N);out+=term
    return trunc(out,N)
def CS(f,time,N):
    a=time*sp.diff(f,e)
    C=trunc(sum(a**(2*n)*e**n/sp.factorial(2*n) for n in range(N//2+1)),N)
    S=trunc(sum(a**(2*n+1)*e**n/sp.factorial(2*n+1) for n in range(N//2+1)),N)
    return C,S

q=1+x*x*t
E=x*x*t*t*q
source_s=x*t
xinverse=(e-s*s)/s**3
tinverse=s**4/(e-s*s)
zero((source_s).subs({x:xinverse,t:tinverse},simultaneous=True)-s,'inverse s')
zero(E.subs({x:xinverse,t:tinverse},simultaneous=True)-e,'inverse E')
zero(xinverse.subs({s:source_s,e:E},simultaneous=True)-x,'inverse x')
zero(tinverse.subs({s:source_s,e:E},simultaneous=True)-t,'inverse t')
source_ds=sp.diff(source_s,x)*sp.diff(E,t)-sp.diff(source_s,t)*sp.diff(E,x)
zero(source_ds-(source_s**2-E),'literal Riccati equation')
zero(xinverse*s**3+s*s-e,'source x algebraic identity')
zero(s**4+tinverse*s*s-tinverse*e,'source t algebraic identity')

fs=[e,e*e,e+e*e,2*e**3-e]
records=[]
N=10
for f in fs:
    # The literal source Hamiltonian belongs to the preceding Newton cone.
    R=clean((f.subs(e,E)-f.subs(e,0))/(t*t*q))
    check(sp.denom(R)==1,f'Newton polynomial {f}')
    check(all(2*j-i>=-2 for (i,j),c in sp.Poly(R,x,t).terms()),f'Newton support {f}')
    for time in [1,-1,2]:
        C,S=CS(f,time,N)
        zero(trunc(C*C-e*S*S-1,N),f'formal determinant {f}:{time}')
        phi=trunc((s*C-e*S)*inverse_unit(C-s*S,N),N)
        term=s;flow=s
        for n in range(1,N+1):
            term=trunc(delta(f,term),N)
            flow+=sp.Rational(time)**n*term/sp.factorial(n)
        zero(trunc(phi-flow,N),f'formal Mobius flow {f}:{time}')
        zero(trunc(phi.subs(s,0)*C+e*S,N),f'specialization formula {f}:{time}')
        g=sp.expand(-phi.subs({s:0,e:z*z},simultaneous=True)/z)
        U=time*z*sp.diff(f,e).subs(e,z*z)
        P=sp.diff(U,z)
        check(P!=0 and P.is_polynomial(z),f'nonzero polynomial logarithmic derivative {f}:{time}')
        residual=sp.Poly(sp.expand(sp.diff(g,z)-P*(1-g*g)),z)
        check(all(i[0]>N-2 or c==0 for i,c in residual.terms()),f'Riccati boundary jet {f}:{time}')
        zero(g.subs(z,0),f'boundary initial value {f}:{time}')
        records.append({'f':str(f),'time':time,'U':str(U),'P':str(P),'weight_cutoff':N})

# Formal-time differential identities, with time left symbolic.
f=e+e*e
C,S=CS(f,L,8)
zero(trunc(sp.diff(C,L)-sp.diff(f,e)*e*S,8),'C formal time derivative')
zero(trunc(sp.diff(S,L)-sp.diff(f,e)*C,8),'S formal time derivative')

# Every pair from the same invariant family commutes.
for f in fs:
    for g in fs:
        zero(delta(f,delta(g,s))-delta(g,delta(f,s)),f'commuting pair {f}:{g}')
for word in [[(e,1),(e,-1)],[(e*e+e,1),(e*e,-1),(e,-1)],[(e,1),(e*e,1)]]:
    F=sum(time*f for f,time in word)
    zero(sum(time*delta(f,s) for f,time in word)-delta(F,s),'effective polynomial derivation')
    check((sp.diff(F,e)==0)==(sp.degree(F,e)<=0 or F==0),'constant cancellation criterion')

# Rational logarithmic derivatives have no polynomial part at infinity.
for a in range(6):
    for b in range(6):
        r=(z**a+1)/(z**b-2)
        logder=sp.cancel(sp.diff(r,z)/r)
        num,den=sp.fraction(logder)
        check(num==0 or sp.degree(num,z)<sp.degree(den,z),f'rational norm infinity {a}:{b}')

# Symbolic quadratic norm trace identity, including derivative of the generator.
a,b,d,ap,bp,dp,u=sp.symbols('a b d ap bp dp u')
r=a+b*u;rc=a-b*u
rp=ap+bp*u+b*dp*u/(2*d)
rcp=ap-bp*u-b*dp*u/(2*d)
num,den=sp.fraction(sp.cancel(rp/r+rcp/rc))
trace=sp.expand(num).subs(u**2,d)/sp.expand(den).subs(u**2,d)
Nnorm=a*a-b*b*d
Nprime=2*a*ap-2*b*bp*d-b*b*dp
zero(trace-Nprime/Nnorm,'quadratic logarithmic norm trace')

# Removing a common s factor preserves a nonzero specialized relation.
for k in range(1,5):
    for Q in [Y-s-e,Y*Y-1-s-e,(1+e+s)*Y-e-s]:
        relation=s**k*Q
        check(sp.expand(relation/s**k).subs(s,0)!=0,f'specialization common factor {k}:{Q}')

print(json.dumps({'status':'PASS','checks':checks,'formal_Riccati_controls':records,
 'commuting_family_polynomials':[str(f) for f in fs],
 'scope':'Finite exact controls; transcendence follows from the written specialization and norm proofs, not finite truncations.'},indent=2,sort_keys=True))
