#!/usr/bin/env python3
"""Independent exact event-cell oracle and constructive conductor audit."""
from fractions import Fraction as F
from math import floor, lcm
from random import Random
import json

def require(p, msg):
    if not p: raise RuntimeError(msg)

def bad(t,v,p,r):
    x=v*t+p
    return abs(x-floor(x+F(1,2)))<r

def events(a,b,vs,ps,r):
    out={a,b}
    for v,p in zip(vs,ps):
        for n in range(floor(v*a+p)-2,floor(v*b+p)+3):
            for sign in [-1,1]:
                x=(n+sign*r-p)/v
                if a<=x<=b:out.add(x)
    return sorted(out)

def cell_measure(a,b,v,p,r):
    es=events(a,b,[v],[p],r)
    return sum((y-x for x,y in zip(es,es[1:]) if bad((x+y)/2,v,p,r)),F(0))

def primitive(x,h):
    n=floor(x)
    return h*n+min(x-n,h)

def measured(a,b,v,p,r):
    return (primitive(v*b+p+r,2*r)-primitive(v*a+p+r,2*r))/v

def component(t,v,p,r):
    require(bad(t,v,p,r),'not in component')
    n=floor(v*t+p+F(1,2))
    return (n-r-p)/v,(n+r-p)/v

counts={'inputs':0,'admissible':0,'rejected':0,'overlap_repairs':0,'strict_density':0,'direct':0,'halvings':0,'negative_times':0,'oracle_queries':0}

def run(a,b,vs,ps,r):
    counts['inputs']+=1
    if a<0:counts['negative_times']+=1
    ms=[measured(a,b,v,p,r) for v,p in zip(vs,ps)]
    for v,p,m in zip(vs,ps,ms):
        require(m==cell_measure(a,b,v,p,r),'measure mismatch')
        counts['oracle_queries']+=1
    W=b-a
    if any(m>=W for m in ms) or sum(ms)>W:
        counts['rejected']+=1
        return
    counts['admissible']+=1
    safe=lambda t: not any(bad(t,v,p,r) for v,p in zip(vs,ps))
    require(any(safe(t) for t in events(a,b,vs,ps,r)),'existence failure')
    if safe(a):counts['direct']+=1;return
    D=W-sum(ms); O=None
    if D==0:
        i=next(i for i,(v,p) in enumerate(zip(vs,ps)) if bad(a,v,p,r))
        c,e=component(a,vs[i],ps[i],r)
        require(a<e<b,'interior endpoint failure')
        if safe(e):counts['direct']+=1;return
        j=next(j for j,(v,p) in enumerate(zip(vs,ps)) if bad(e,v,p,r))
        require(j!=i,'same component bad at open endpoint')
        d,f=component(e,vs[j],ps[j],r)
        O=max(a,c,d),min(b,e,f)
        require(O[0]<O[1],'overlap empty')
        counts['overlap_repairs']+=1
    else:counts['strict_density']+=1
    def integral(x,y):
        ans=y-x-sum(measured(x,y,v,p,r) for v,p in zip(vs,ps))
        if O:ans+=max(F(0),min(y,O[1])-max(x,O[0]))
        # independent full-event geometric integration, including O endpoints
        es=set(events(x,y,vs,ps,r))
        if O:es.update(t for t in O if x<t<y)
        es=sorted(es);exact=F(0)
        for q,t in zip(es,es[1:]):
            mid=(q+t)/2
            density=1-sum(bad(mid,v,p,r) for v,p in zip(vs,ps))
            if O and O[0]<mid<O[1]:density+=1
            require(density<=int(safe(mid)),'density not dominated')
            exact+=(t-q)*density
        require(ans==exact,'integral mismatch')
        return ans
    Q=lcm(r.denominator,*(p.denominator for p in ps))
    P=max(v.numerator for v in vs)
    sep=F(1,Q*Q*P*P)
    es=events(a,b,vs,ps,r)
    actual=[t for t in es if t not in (a,b)]
    require(all(y-x>=sep for x,y in zip(actual,actual[1:])),'separation wrong')
    x,y=a,b
    while True:
        require(integral(x,y)>0,'lost positivity')
        if safe(x) or safe(y):return
        require(y-x>=sep,'terminal endpoints unsafe')
        mid=(x+y)/2
        if integral(x,mid)>0:y=mid
        else:x=mid
        counts['halvings']+=1

# equality with isolated safe points; equality with positive overlap;
# invalid all-bad input; negative times; reduced rational speeds.
for ps in ([F(0),F(1,2)],[F(0),-F(1,8)]):
    run(F(0),F(1),[F(1),F(1)],ps,F(1,4))
run(F(0),F(1,10),[F(1)],[F(0)],F(1,4))
rng=Random(12320261031)
for _ in range(1200):
    k=rng.randrange(1,7)
    a=F(rng.randrange(-12,13),rng.randrange(1,8))
    b=a+F(rng.randrange(1,16),rng.randrange(1,9))
    vs=[F(rng.randrange(1,16),rng.randrange(1,5)) for _ in range(k)]
    ps=[F(rng.randrange(-15,16),rng.randrange(1,10)) for _ in range(k)]
    r=F(1,rng.randrange(3,17))
    run(a,b,vs,ps,r)
require(counts['overlap_repairs']>0,'missing overlap branch')
require(counts['halvings']>0,'missing bisection branch')
require(counts['rejected']>0,'missing rejection control')
print(json.dumps({'status':'PASS','scope':'independent exact finite geometry; integer and rational speeds','counts':counts},sort_keys=True,indent=2))
