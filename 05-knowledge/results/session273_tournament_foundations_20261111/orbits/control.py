#!/usr/bin/env python3
"""Exact small controls for the finite-horizon orbit law, not a new census."""
from itertools import permutations
from fractions import Fraction as F
from collections import defaultdict
import json

def pairs(n): return [(i,j) for i in range(n) for j in range(i+1,n)]
def restrict(code,n,k):
    ps=pairs(n); return sum(((code>>ps.index(p))&1)<<a for a,p in enumerate(pairs(k)))
def relabel(code,n,p):
    ps=pairs(n); ix={e:i for i,e in enumerate(ps)}; z=0
    for bit,(i,j) in enumerate(ps):
        u,v=p[i],p[j]; b=(code>>ix[tuple(sorted((u,v)))])&1
        z|=(b if u<v else 1-b)<<bit
    return z

def fixed_completion(n,k,t,p):
    """Independent signed pair-orbit union; count completions fixed by p."""
    ps=pairs(n); ix={e:i for i,e in enumerate(ps)}; seen=set(); free=0
    prescribed={ix[e]:(t>>i)&1 for i,e in enumerate(pairs(k))}
    for a in range(len(ps)):
        if a in seen: continue
        u=a; parity=0; forced=None
        while u not in seen:
            seen.add(u)
            if u in prescribed:
                val=prescribed[u]^parity
                if forced is not None and val!=forced:return 0
                forced=val
            i,j=ps[u]; i,j=p[i],p[j]
            parity^=i>j;u=ix[tuple(sorted((i,j)))]
        if parity:return 0
        free+=forced is None
    return 2**free

def run():
    report=[]
    for n in range(1,6):
        perms=list(permutations(range(n))); size=2**len(pairs(n)); g=[]; canon=[]
        for t in range(size):
            r=[relabel(t,n,p) for p in perms]
            g.append(r.count(t));canon.append(min(r))
        a=len(set(canon)); b=F(size,len(perms)); delta=F(a,b)-1
        assert sum(g)==a*len(perms)
        # Full-history terminal change of measure, conditional martingale.
        H=[]
        for k in range(n+1):
            bins=defaultdict(int)
            for t,aut in enumerate(g): bins[restrict(t,n,k)]+=aut
            denom=2**(len(pairs(n))-len(pairs(k)))
            h={s:F(v,denom) for s,v in bins.items()};H.append(h)
            for s,v in h.items():
                assert v>=1
                fixed=sum(fixed_completion(n,k,s,p) for p in perms)
                assert v==F(fixed,denom)
        for k in range(n):
            for t,h in H[k].items():
                children=[s for s in H[k+1] if restrict(s,k+1,k)==t]
                assert sum(H[k+1][s] for s in children)==h*2**k
                assert sum(F(1,2**k)*H[k+1][s]/h for s in children)==1
        tv=sum(abs(F(aut,sum(g))-F(1,size)) for aut in g)/2
        assert tv<=delta/(1+delta)
        row={'n':n,'a':a,'delta':str(delta),'TV_class_lift_vs_fair':str(tv)}
        if n==4:
            cyclic_codes={2,5} # pair order 01,02,12
            mass=sum(aut for t,aut in enumerate(g) if restrict(t,n,3) in cyclic_codes)
            assert F(mass,sum(g))==F(1,4)
            row['deleted_cyclic_probability']='1/4'
            row['uniform_three_class_cyclic_probability']='1/2'
        report.append(row)
    print(json.dumps({'status':'PASS','scope':'all labelled tournaments n<=5; conditional horizons through5','results':report},indent=2))
if __name__=='__main__':run()
