#!/usr/bin/env python3
"""A second algebraic path checks the mass gain and literal wedge accounting."""
from fractions import Fraction as F
from itertools import combinations,product
import json,random

def need(ok,message):
    if not ok:raise RuntimeError(message)

def main():
    d=F(2491,10000);t=d/8;lam=d/64
    # Minimize the quadratic in S directly, rather than reading its square
    # decomposition. The residual after subtracting tw/4 is nonnegative.
    need(lam*d==t*t and 4*lam==t/2,'coefficient identities')
    points=0
    for w in [t*F(i,100) for i in range(101)]:
        Sstar=2*t-w
        need(Sstar>=0,'actual minimum domain')
        def residual(S):return w*S/2+S*S/4-t*S+lam*(d-4*w)-t*w/4
        need(residual(Sstar)==w*(t-w)/4>=0,'minimum residual')
        for S in [F(i,100) for i in range(51)]:
            need(residual(S)==residual(Sstar)+(S-Sstar)**2/4,'quadratic recovery')
            points+=1
    def check(n,weights):
        edges=list(combinations(range(n),2));q=dict(zip(edges,weights))
        W=sum(q[e]*q[f] for e,f in combinations(edges,2) if len(set(e+f))==3)
        rem={e for e in edges if q[e]};groups=[];spent=set()
        while rem:
            e=max(rem,key=lambda e:(q[e],e))
            group={f for f in rem if set(e)&set(f)};rem-=group
            w=q[e];S=sum(q[f] for f in group-{e})
            pairs={tuple(sorted((f,g))) for f,g in combinations(group,2) if set(f)&set(g)}
            need(not spent&pairs,'no repeated wedge');spent|=pairs
            energy=sum(q[f]*q[g] for f,g in pairs)
            need(energy>=w*S and energy>=w*S/2+S*S/4,'strong local invoices')
            groups.append((w,S))
        D=sum(S for w,S in groups);light=[(w,S) for w,S in groups if w<=t]
        U=sum(w for w,S in light);L=len(light)
        need(W>=t*D-lam*(d*L-4*U)+t*U/4,'summed U-retaining dual')
    exhaustive=0
    for weights in product([F(0),t/2,F(1,4)],repeat=6):
        check(4,weights);exhaustive+=1
    rng=random.Random(514);sampled=0
    for n in range(5,13):
        for _ in range(100):
            check(n,[F(rng.randrange(1001),4000) for e in combinations(range(n),2)])
            sampled+=1
    print(json.dumps({'status':'passed','quadratic_checks':points,
        'complete_K4_controls':exhaustive,'random_controls':sampled,'seed':514,
        'scope':'Separate local algebra and wedge accounting, not an independent-agent audit or a finite-sample global proof.'},indent=2,sort_keys=True))

if __name__=='__main__':main()
