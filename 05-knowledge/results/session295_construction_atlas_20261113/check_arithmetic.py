"""Exact arithmetic controls for weighted finite differences and observer fan."""
from itertools import product
from math import factorial, gcd, prod
from functools import reduce
import json

def fval(f,x):
    return sum(a*x**i for i,a in enumerate(f))
def v(x,p):
    assert x>0
    a=0
    while x%p==0:
        x//=p
        a+=1
    return a
def brute(f,n):
    vals=[fval(f,j) for j in range(n+1)]
    dif=vals[:]
    raw=[]; weighted=[]
    for k in range(1,n+1):
        dif=[y-x for x,y in zip(dif,dif[1:])]
        raw.append(dif[0])
        weighted.append(factorial(k-1)*dif[0])
    rg=reduce(gcd,raw)
    wg=reduce(gcd,weighted)
    if rg!=reduce(gcd,vals[1:]) or wg!=reduce(gcd,(factorial(k-1)*vals[k] for k in range(1,n+1))):
        raise ValueError('integer basis')
    return rg,wg

def mul(f,g):
    h=[0]*(len(f)+len(g)-1)
    for i,a in enumerate(f):
        for j,b in enumerate(g):
            h[i+j]+=a*b
    return h

def observation_fan(fs,p):
    def visit(fixed):
        exponent=max((t+1 for i,t in fixed),default=0)
        j=next((j for j in range(1,p**exponent+1)
                if all(v(fval(fs[i],j),p)==t for i,t in fixed)),None)
        if j is None:
            return set()
        out={j}
        vals=[v(fval(f,j),p) for f in fs]
        for i,a in enumerate(vals):
            if any(i==i0 for i0,t in fixed):
                continue
            for t in range(a):
                out.update(visit(fixed+[(i,t)]))
        return out
    return sorted(visit([]))

f8=[0,0,246,0,232,0,34,0,1]
f3=[0,2,0,1]
fb=[0,0,4,0,1]
fs=[f8,f3,fb]
checks=0
results={}
for p in [2,3,5,7]:
    js=observation_fan(fs,p)
    results[p]=js
    for ms in product(range(4),repeat=len(fs)):
        if not sum(ms):
            continue
        f=[1]
        for g,m in zip(fs,ms):
            for _ in range(m):
                f=mul(f,g)
        raw,ctx=brute(f,len(f)-1)
        vals={j:sum(m*v(fval(g,j),p) for g,m in zip(fs,ms)) for j in js}
        expected=(min(vals.values()),min(a+v(factorial(j-1),p) for j,a in vals.items()))
        if expected!=(v(raw,p),v(ctx,p)):
            raise ValueError((p,ms,expected,raw,ctx))
        checks+=1
for n in range(2,9):
    for a in range(0,15,2):
        f=[0]*(n+1); f[-1]=1; f[n-2]=a
        if not f[0]:
            brute(f,n)
            checks+=1
print(json.dumps({'status':'pass','cases':checks,'common_fan_probe_indices':results},sort_keys=True))
