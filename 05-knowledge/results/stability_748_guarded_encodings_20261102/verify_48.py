#!/usr/bin/env python3
"""Exact two-regime certificate; all checks remain active under -O."""
from fractions import Fraction as F
from math import factorial
from itertools import product,combinations
import json

def need(ok,message):
    if not ok:raise RuntimeError(message)

def main():
    R=24;c=F(99,100);t=F(1,6);d=F(5,6);lam=F(1,30);a=F(3,4)
    need(sum(F(4**j,factorial(j)) for j in range(8))>50,'initial exponential')
    need(F(49,200)>F(4,R),'V>1 bootstrap')
    steps=[(F(49,200),F(7,20)),(F(7,20),F(13,20)),(F(13,20),F(47,50)),(F(47,50),c)]
    mass=lambda z,r:z*F(r,4)-1-F(1,2*r)
    boot=[]
    for prev,nxt in steps:
        m=mass(prev,R);T=sum(m**j/factorial(j) for j in range(13))
        slack=T*(1-nxt-F(1,4)**R)-1
        need(m>0 and 1-nxt-F(1,4)**R>0 and slack>0,'bootstrap endpoint')
        boot.append({'from':str(prev),'to':str(nxt),'mass':str(m),'taylor_degree':12,'slack':str(slack)})
    g=lambda r:(c-4*t)*r-4-F(2,r)
    need(1-4*t==F(1,3)>0 and c-4*t>0,'slot monotonicity')
    need(g(R)-11*(1-4*t)==F(1,100)>0,'twelve heavy edges')
    Q=lambda r:F(r)+F(33,2)+F(8,r)
    tails=[Q(R)*a**(R-1),Q(R)*a**12/(R-12)]
    need(all(x<c-d for x in tails) and c-d==F(47,300),'convex tail endpoints')
    need(t*R==4 and d==5*t and lam==t/5 and t-4*lam==F(1,30)>0,'final dual constants')
    G=lambda u,v:max(u*v,u*v/2+v*v/4)-v+1-F(4,5)*u
    local=0
    for i in range(101):
        u=F(i,100)
        vstar=2-u if u<=F(2,3) else 2*u
        minimum=u*(F(1,5)-u/4) if u<=F(2,3) else 2*(u-F(7,10))**2+F(1,50)
        need(G(u,vstar)==minimum>=0,'exact local minimizer')
        for j in range(251):
            need(G(u,F(j,100))>=minimum,'local lower envelope')
            local+=1
    u=F(17,24);v=2*u
    hostile=max(u*v,u*v/2+v*v/4)-v+1-F(5,6)*u
    need(hostile==F(-1,288),'ratio24/5 scalar hostile')
    # Direct six-edge graphs independently confirm the max invoice and
    # its summation, not a sum of the two incompatible local bounds.
    E=list(combinations(range(4),2));graph_controls=0
    for weights in product([F(0),F(1,8),F(1,4)],repeat=6):
        q=dict(zip(E,weights));remaining={e for e in E if q[e]};groups=[]
        while remaining:
            e=max(remaining,key=lambda e:(q[e],e));members={f for f in remaining if set(e)&set(f)}
            remaining-=members;w=q[e];S=sum(q[f] for f in members-{e})
            energy=sum(q[f]*q[z] for f,z in combinations(members,2) if set(f)&set(z))
            need(energy>=max(w*S,w*S/2+S*S/4),'literal max invoice')
            groups.append((w,S,energy))
        W=sum(q[e]*q[f] for e,f in combinations(E,2) if set(e)&set(f))
        light=[(w,S) for w,S,en in groups if w<=t]
        D=sum(S for w,S,en in groups);L=len(light);U=sum(w for w,S in light)
        need(W>=t*D-lam*(d*L-4*U),'whole disjoint dual')
        graph_controls+=1
    uniform=0
    for r in range(24,301):
        need(g(r)>11*(1-4*t),'finite slot control')
        need(Q(r)*a**(r-1)<c-d and Q(r)*a**12/(r-12)<c-d,'finite tail controls')
        if lam*r<1:need((1-lam*r)*F(4,r)<=t-4*lam,'finite final inequality')
        uniform+=1
    need(23*t<4 and g(23)<11*(1-4*t),'order47 rejected by this certificate')
    print(json.dumps({'status':'passed','threshold':48,'remaining_open_orders':[10,47],
        'bootstrap':boot,'slot_slack':'1/100','tail_endpoints':[str(x) for x in tails],
        'tail_allowance':'47/300','local_checks':local,'whole_graph_controls':graph_controls,
        'finite_R_controls':uniform,'scalar_ratio24over5_hostile':str(hostile),
        'negative_order':47,'scope':'The continuum minimization, monotone bootstrap and convex endpoints prove all n>=48. Finite controls are supplementary; no independent-agent audit claimed.'},indent=2,sort_keys=True))

if __name__=='__main__':main()
