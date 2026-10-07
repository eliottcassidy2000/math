#!/usr/bin/env python3
"""Folkman m=3 with 3-smooth generators, 2 colours:
 (i) SAT with colour fixed to sigma=Omega mod 2 on S (only non-S numbers free)
 (ii) explicit greedy rules."""
import itertools, sys, bisect
from pysat.solvers import Solver
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

def sat_fixed(N, fix_sigma=True):
    S=smooth3(N); Sset=set(S)
    s=Solver(name='cadical195'); nc=0
    if fix_sigma:
        for x in S: s.add_clause([x if sigma(x) else -x])
    for a,b,c in itertools.combinations(S,3):
        if a+b+c>N: continue
        fs=sorted({a,b,c,a+b,a+c,b+c,a+b+c})
        s.add_clause(fs); s.add_clause([-x for x in fs]); nc+=1
    ok=s.solve(); s.delete(); return ok,nc
for e in [8,10,12,14,16,18,20]:
    ok,nc=sat_fixed(2**e,True)
    print(f"(i) sigma fixed on S, N=2^{e}: {'SAT' if ok else 'UNSAT'} ({nc} configs)",flush=True)
    if not ok: break

# (ii) greedy rules
N=int(sys.argv[1]) if len(sys.argv)>1 else 2**36
S=smooth3(4*N); Sset=set(S)
def P(n): return S[bisect.bisect_right(S,n)-1]
def Qn(n): return S[bisect.bisect_left(S,n)]
def greedy_terms(n):
    k=0
    while n>0:
        n-=P(n); k+=1
    return k
rules={
 'R1 sigma|S, 1-sigma(P(n)) off S': lambda n: sigma(n) if n in Sset else 1-sigma(P(n)),
 'R2 sigma|S, 1-sigma(Q(n)) off S': lambda n: sigma(n) if n in Sset else 1-sigma(Qn(n)),
 'R3 sigma(P(n)) + (#greedy terms-1) mod 2': lambda n: (sigma(P(n))+greedy_terms(n)-1)%2,
 'R4 (#greedy terms) mod 2': lambda n: greedy_terms(n)%2,
}
SN=smooth3(N)
for name,f in rules.items():
    cache={}
    def col(n):
        r=cache.get(n)
        if r is None: r=cache[n]=f(n)
        return r
    bad=[]; tot=0
    for a,b,c in itertools.combinations(SN,3):
        if a+b+c>N: continue
        tot+=1
        ca=col(a)
        if col(b)!=ca or col(c)!=ca: continue
        if col(a+b)==ca and col(a+c)==ca and col(b+c)==ca and col(a+b+c)==ca:
            bad.append((a,b,c))
            if len(bad)>=4: break
    print(f"(ii) {name:42s}: mono FS(a,b,c) found: {bad[:4]} (checked {tot}, N={N})",flush=True)
