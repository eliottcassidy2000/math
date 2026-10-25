#!/usr/bin/env python3
"""Exact algebra supporting LOWER_BOUND.md; no assertions disappear under -O."""
from fractions import Fraction as F
import json
from pathlib import Path
import sympy as s
checks = []
def require(ok, label):
    if not bool(ok): raise RuntimeError(label)
    checks.append(label)
def zero(f, label): require(s.expand(f)==0, label)
def order(f, v, a):
    q=s.Poly(s.expand(f),v); n=0
    while q.eval(a)==0:
        q,r=s.div(q,s.Poly(v-a,v));require(r.is_zero, 'exact valuation division')
        n+=1
    return n
X,t,z,w,u=s.symbols('X t z w u');d=1-t*t
f=(X-d)*(X+d)*(X-d*(t+1))
zero(f.subs(X,-d*z)-d**3*(1-z*z)*(1+t+z),'birational identity')
require(s.expand(f.subs(X,d*z)-d**3*(1-z*z)*(1+t+z))!=0,'hostile wrong sign detected')
zero(f.subs(X,d*t)-d**4,'section membership')
require(s.expand(f.subs(X,d*t)-d**2)!=0,'hostile wrong section Y detected')
a2=-d*(t+1);a4=-d*d;a6=d**3*(t+1)
zero(f-(X**3+a2*X**2+a4*X+a6),'Weierstrass expansion')
c4=16*(a2*a2-3*a4);c6=-64*a2**3+288*a2*a4-864*a6
Delta=16*s.discriminant(f,X)
zero(c4-16*d*d*(t*t+2*t+4),'c4')
zero(c6+64*(t-2)*(t-1)**3*(t+1)**4*(t+4),'c6')
zero(Delta-64*d**6*t*t*(t+2)**2,'discriminant')
zero(c4**3-c6**2-1728*Delta,'c invariant identity independent discriminant route')
rows=[]
for a,expect,typ,rank in [(0,(0,0,2),'I2',1),(-2,(0,0,2),'I2',1),(1,(2,3,6),'I0*',4),(-1,(2,4,6),'I0*',4)]:
    vals=tuple(order(p,t,a) for p in [c4,c6,Delta])
    require(vals==expect,'valuation '+str(a)); require(vals[2]<12,'minimality '+str(a))
    rows.append({'place':a,'valuations':vals,'fiber':typ,'root_rank':rank})
c4i=s.cancel(u**8*c4.subs(t,1/u));c6i=s.cancel(u**12*c6.subs(t,1/u));Di=s.cancel(u**24*Delta.subs(t,1/u))
for a,weight,expected in [(a2,4,-u*(u*u-1)*(u+1)),(a4,8,-u**4*(u*u-1)**2),(a6,12,u**5*(u*u-1)**3*(u+1))]:
    zero(s.cancel(u**weight*a.subs(t,1/u))-expected,'infinity coefficient weight '+str(weight))
vals=tuple(order(p,u,0) for p in [c4i,c6i,Di]);require(vals==(2,3,8),'infinity valuations')
require(vals[2]<12,'infinity minimality');rows.append({'place':'infinity','valuations':vals,'fiber':'I2*','root_rank':6})
require(2+sum(r['root_rank'] for r in rows)==18,'trivial lattice rank')
require(sum(r['valuations'][2] for r in rows)==24,'Euler number')
require(Delta.subs(t,2)!=0,'smooth specialization')
zero(f.subs({t:2,X:u-3})-(u**3-36*u),'integral short specialization')

def on_curve(P):return P is None or P[1]**2==P[0]**3-36*P[0]
def add(P,Q):
    if P is None:return Q
    if Q is None:return P
    x,y=P;xx,yy=Q
    if x==xx and y==-yy:return None
    m=(3*x*x-36)/(2*y) if P==Q else (yy-y)/(xx-x)
    a=m*m-x-xx;b=m*(x-a)-y
    return a,b
Q=(F(-3),F(9));require(on_curve(Q),'specialized point')
Q2=add(Q,Q);require(Q2==(F(25,4),F(-35,8)),'exact tangent doubling')
def v2rat(a):
    def v2(n):
        n=abs(n);k=0
        if n==0:raise ValueError('infinite valuation')
        while n%2==0:n//=2;k+=1
        return k
    return v2(a.numerator)-v2(a.denominator)
zero(s.cancel(((3*u*u-36)**2/(4*(u**3-36*u))-2*u)-(u*u+36)**2/(4*u*(u*u-36))),'duplication identity')
R=Q;double_vals=[]
for n in range(1,7):
    R=add(R,R);require(on_curve(R),'double lies on curve '+str(n))
    require(v2rat(R[0])==-2*n,'valuation induction positive control '+str(n))
    double_vals.append({'exponent':n,'v2_x':v2rat(R[0])})
T=(F(0),F(0));Qpos=add(Q,T);require(Qpos==(F(12),F(36)),'Session57 two-torsion translation')
require(add(Qpos,Qpos)==Q2,'same doubled point under two-torsion translation')
a=(Q[0]*Q[0]-36)/Q[1];b=12*Q[0]/Q[1];h=(Q[0]*Q[0]+36)/Q[1]
require((a,b,h)==(F(-3),F(-4),F(5)),'signed triangle')
require(a*a+b*b==h*h and a*b/2==6,'right triangle and area')
require(4*Q2[0]==25 and [25-24,25,25+24]==[1,25,49],'Session57 fixed-gap center')
require([6*r-56 for r in [19,20]]==[58,64],'hourglass consequence')
result={'status':'PASS','checks':len(checks),'fibers':rows,'trivial_lattice_rank':18,'specialization':{'curve':'Y^2=u^3-36u','point':['-3','9'],'double':['25/4','-35/8'],'doubling_valuation_controls':double_vals},'proved_outer_picard_lower_bound':19,'lower_bound_snapshot_outer_possibilities':[19,20],'lower_bound_snapshot_hourglass_possibilities':[58,64],'scope':'Algebra certificate plus cited/inherited geometry; no exact-rank upper bound claimed by this script.'}
text=json.dumps(result,indent=2,sort_keys=True)+'\n'
Path(__file__).with_name('LOWER_BOUND_CHECKS.json').write_text(text)
print(text,end='')
