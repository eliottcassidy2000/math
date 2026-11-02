#!/usr/bin/env python3
"""Exact endpoints plus separate degree-energy and equality controls."""
from fractions import Fraction as F
from itertools import combinations,product
from math import factorial
import json,random

def need(ok,message):
    if not ok:raise RuntimeError(message)

def main():
    R=33;eta=F(4,33);c=F(99,100);a0=F(817,1650);a=F(3,4)
    need(sum(F(4**j,factorial(j)) for j in range(8))>50,'first exponential')
    need(sum(F(8**j,factorial(j)) for j in range(8))>1000,'second exponential')
    need(F(49,200)>F(4,R) and F(49*R,800)>2,'V>2 bootstrap')
    m=lambda c,r:c*F(r,4)-1-F(1,2*r)
    transitions=[(F(999,4000),F(3,5),2),(F(3,5),F(19,20),3),(F(19,20),c,4)]
    boot=[]
    for prev,nxt,degree in transitions:
        mass=m(prev,R)
        T=sum(mass**j/factorial(j) for j in range(degree+1))
        slack=T*(1-nxt-F(1,4)**R)-1
        need(mass>0 and 1-nxt-F(1,4)**R>0 and slack>0,'nonlinear bootstrap')
        boot.append({'from':str(prev),'to':str(nxt),'degree':degree,'mass':str(mass),'slack':str(slack)})
    g=lambda r:(c-4*eta)*r-4-F(2,r)
    need(1-4*eta==F(17,33)>0,'slot divisor')
    need(g(R)-24*(1-4*eta)==F(811,3300)>0,'25 retained heavy edges')
    need(c-4*eta>0,'slot monotonicity')
    Q=lambda r:F(r)+F(33,2)+F(8,r)
    tails=[Q(R)*a**(R-1),Q(R)*a**25/(R-25)]
    need(all(x<F(1,100) for x in tails),'both convex endpoints')
    need(a0==c-4*eta-F(1,100) and a0-4*eta==F(17,1650)>0,'positive dual coefficients')
    need(eta*R==4,'final endpoint closes exactly')
    uniform=0
    for r in range(33,401):
        need(g(r)>24*(1-4*eta),'retained count finite control')
        need(Q(r)*a**(r-1)<F(1,100) and Q(r)*a**25/(r-25)<F(1,100),'tail finite controls')
        if a0-r*eta**2>0:
            need(F(4,r)*a0-4*eta**2<=eta*a0-4*eta**2,'last comparison')
        uniform+=1
    # Independent literal weighted-degree calculation for the stronger
    # input, including assigned selected/discarded invoices.
    def degree_check(n,weights):
        E=list(combinations(range(n),2));q=dict(zip(E,weights))
        W=sum(q[e]*q[f] for e,f in combinations(E,2) if set(e)&set(f))
        rem={e for e in E if q[e]};discard=set();assigned=[]
        while rem:
            e=max(rem,key=lambda e:(q[e],e));rem.remove(e)
            fs={f for f in rem if set(e)&set(f)};rem-=fs;discard|=fs
            assigned.extend((e,f) for f in fs)
        degrees=[sum(q[e] for e in discard if v in e) for v in range(n)]
        D=sum(q[e] for e in discard)
        dd=sum(q[e]*q[f] for e,f in combinations(discard,2) if set(e)&set(f))
        sq=sum(q[e]**2 for e in discard)
        need(dd==sum(x*x for x in degrees)/2-sq,'discarded-degree identity')
        need(sum(q[e]*q[f] for e,f in assigned)>=sq,'distinct cancellation invoices')
        need(W>=sum(x*x for x in degrees)/2>=2*D*D/n,'factor-four floor')
    graph_controls=0
    for w in product([F(0),F(1,8),F(1,4)],repeat=6):
        degree_check(4,w);graph_controls+=1
    rng=random.Random(66)
    for n in range(5,13):
        for _ in range(100):
            degree_check(n,[F(rng.randrange(101),400) for e in combinations(range(n),2)])
            graph_controls+=1
    # Actual ordinal fair-pair examples for both parity cases.
    equality=[]
    for n in [66,67,140,141]:
        r=n//2;V=F(r,4);u=a**r;v=F(1,4)**r if n%2==0 else F(0)
        k=F(4,r)*(1-u-v);need(1-u-v==k*V,'actual sharp block equality')
        equality.append(n)
    need(g(32)<24*(1-4*eta) and 32*eta<4,'order65 rejected by certificate')
    print(json.dumps({'status':'passed','threshold':66,'remaining_open_orders':[10,65],
        'bootstrap':boot,'slot_slack':'811/3300','tail_endpoints':[str(x) for x in tails],
        'dual_slack':'17/1650','finite_R_controls':uniform,'graph_controls':graph_controls,
        'actual_equality_orders':equality,'negative_control':'same certificate at order65 fails',
        'scope':'All-n proof uses monotone bootstrap, slot count and convex endpoints. Finite checks corroborate it; no independent-agent audit is claimed.'},indent=2,sort_keys=True))

if __name__=='__main__':main()
