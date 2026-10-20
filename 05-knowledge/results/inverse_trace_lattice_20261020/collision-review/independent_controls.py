#!/usr/bin/env python3
"""Small exact controls for the collision-gap certificate and singular image."""
from pathlib import Path
import json
import sympy as S
x,y,a,s,t,u,v,z=S.symbols('x y a s t u v z')
labels=[]
def check(c,label):
    if not bool(c):
        raise RuntimeError(label)
    labels.append(label)
def dd(f):
    return S.cancel((f.subs({x:a,y:s})-f.subs({x:a,y:t}))/(s-t)).expand()
# Actual polynomial Keller automorphisms: verify the universal identity exactly.
for r,n in [(2,2),(2,3),(3,2)]:
    P=x+y**r
    Q=y+P**n
    J=S.expand(S.diff(P,x)*S.diff(Q,y)-S.diff(P,y)*S.diff(Q,x))
    U,V=dd(P),dd(Q)
    numerator=S.expand(S.diff(P,x).subs({x:a,y:s})*V-S.diff(Q,x).subs({x:a,y:s})*U-J)
    W,rem=S.div(numerator,s-t,s)
    check(J==1,f'constant Jacobian r={r},n={n}')
    check(S.expand(U.subs(t,s)-S.diff(P,y).subs({x:a,y:s}))==0,f'divided derivative P r={r},n={n}')
    check(S.expand(V.subs(t,s)-S.diff(Q,y).subs({x:a,y:s}))==0,f'divided derivative Q r={r},n={n}')
    check(rem==0,f'collision quotient W is polynomial r={r},n={n}')
    check(S.expand((s-t)*W-numerator)==0,f'unit certificate identity r={r},n={n}')
    check(list(S.groebner([U,V],s,t,a).polys)[0].as_expr()==1,f'automorphism collision scheme empty r={r},n={n}')
    graph=z-u+(v-u**n)**r
    check(S.expand(graph.subs({u:P,v:Q,z:x}))==0,f'actual Keller graph identity r={r},n={n}')
    check(S.diff(graph,z)==1,f'Keller graph has no singular or RZ locus r={r},n={n}')
    Ru=S.diff(graph,u).subs({u:P,v:Q,z:x}); Rv=S.diff(graph,v).subs({u:P,v:Q,z:x})
    check(S.expand(Ru*S.diff(P,x)+Rv*S.diff(Q,x)+1)==0 and S.expand(Ru*S.diff(P,y)+Rv*S.diff(Q,y))==0,f'graph chain-rule identities r={r},n={n}')
# Finite immersed nodal surface normalization: persistent all-a collisions without Keller.
P=y*y-1
Q=x+y*(y*y-1)
R=(v-z)**2-u*u*(u+1)
check(S.expand(R.subs({u:P,v:Q,z:x}))==0,'nodal family actual graph')
check(S.cancel((Q-x)/P)==y,'nodal family graph generically birational')
J=S.expand(S.diff(P,x)*S.diff(Q,y)-S.diff(P,y)*S.diff(Q,x))
check(J==-2*y,'nodal family is not Keller')
U,V=dd(P),dd(Q)
check(U==s+t and V==s*s+s*t+t*t-1,'nodal family collision equations')
I=S.groebner([U,V],s,t,a,domain=S.QQ)
check(I.reduce(S.expand((s-t)**2/4-1))[1]==0,'gap is a unit on nodal collisions')
check(S.expand((s-t)**2/4-V+3*U**2/4)==1,'diagonal-collision comaximality certificate')
check(I.reduce(s*s-1)[1]==0,'collision image has u=0')
check(I.reduce(s*(s*s-1))[1]==0,'collision image has v=a')
check(U.subs({s:1,t:-1})==V.subs({s:1,t:-1})==0,'one collision over every a')
singular=S.groebner([R,S.diff(R,u),S.diff(R,v),S.diff(R,z)],u,v,z)
check([p.as_expr() for p in singular.polys]==[u,v-z],'singular ideal is exactly the image line')
check(S.gcd(2*y,3*y*y-1)==1,'normalization remains immersive')
check(R.subs({u:-1,v:0,z:0})==S.diff(R,z).subs({u:-1,v:0,z:0})==0 and S.diff(R,u).subs({u:-1,v:0,z:0})!=0,'non-Keller RZ locus can contain smooth graph points')
# Dominance alone can omit a parameter value; affine finiteness is not projective properness.
check((a*t-1).subs(a,0)==-1,'dominant hyperbola projection misses zero')
check(S.cancel((a*t-1).subs(t,1/a))==0,'hyperbola projection contains every nonzero parameter')
# A unit test sees nilpotents too, rather than merely forgetting their structure.
N=S.groebner([u*u,v*v],u,v)
check(N.reduce(S.expand((1+u+v)*(1-u-v+2*u*v)-1))[1]==0,'unit lifts across nonreduced structure')
out={'status':'PASS','checks':len(labels),'labels':labels,'scope':'Exact small controls. No noninvertible Keller map or general classification is constructed.'}
Path(__file__).with_name('results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'checks':len(labels),'status':'PASS'},sort_keys=True))
