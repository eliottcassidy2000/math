#!/usr/bin/env python3
"""Test explicit 2-colourings of N against monochromatic FS({a,b,c}), a<b<c 3-smooth."""
import itertools, sys
from sympy import factorint
def smooth3(N):
    out=[]; a=1
    while a<=N:
        b=a
        while b<=N:
            out.append(b); b*=3
        a*=2
    return sorted(out)
def v(n,p):
    c=0
    while n%p==0: n//=p; c+=1
    return c
def sigma(n): return (v(n,2)+v(n,3))%2
def liou(n): return sum(factorint(n).values())%2
def popc(n): return bin(n).count('1')%2
rules={
 'sigma=(v2+v3) mod 2': sigma,
 'Liouville Omega mod 2': liou,
 'v2 mod 2': lambda n: v(n,2)%2,
 'sigma on S, popcount off S': None,
}
N=int(sys.argv[1]) if len(sys.argv)>1 else 2**30
S=smooth3(N); Sset=set(S)
def mk_hybrid(off):
    return lambda n: sigma(n) if n in Sset else off(n)
rules['sigma on S, popcount off S']=mk_hybrid(popc)
rules['sigma on S, 0 off S']=mk_hybrid(lambda n:0)
rules['sigma on S, (n mod 3==2) off S']=mk_hybrid(lambda n: int(n%3==2))
rules['sigma on S, v2(n) mod 2 off S']=mk_hybrid(lambda n: v(n,2)%2)
rules['sigma on S, (sigma(n)+[n odd&n%3!=0]) off S']=mk_hybrid(lambda n: (sigma(n)+ (1 if (n%2 and n%3) else 0))%2)
for name,f in rules.items():
    cache={}
    def col(n):
        if n not in cache: cache[n]=f(n)
        return cache[n]
    bad=[]; tot=0
    for a,b,c in itertools.combinations(S,3):
        if a+b+c>N: continue
        tot+=1
        fs={a,b,c,a+b,a+c,b+c,a+b+c}
        cs={col(x) for x in fs}
        if len(cs)==1:
            bad.append((a,b,c))
            if len(bad)>=5: break
    print(f"{name:45s}: first mono FS(a,b,c): {bad[:5]}  (checked {tot})", flush=True)
