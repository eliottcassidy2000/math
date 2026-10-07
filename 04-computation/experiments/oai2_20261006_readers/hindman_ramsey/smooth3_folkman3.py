#!/usr/bin/env python3
"""Explicit 2-colourings built from representations as sums of 3-smooth numbers."""
import itertools, sys, bisect
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
N=int(sys.argv[1]); S=smooth3(N); Sset=set(S)
# top2[n] = largest b in S with n-b in S (n-b>=1), for non-S n that are 2-sums
top2={}; low2={}
for i,a in enumerate(S):
    for b in S[i:]:
        n=a+b
        if n>N: break
        if n in Sset: continue
        if n not in top2 or b>top2[n]: top2[n]=b
        if n not in low2 or a<low2[n]: low2[n]=a
def rule_factory(kind):
    def f(n):
        if n in Sset: return sigma(n)
        if n in top2:
            if kind=='top': return 1-sigma(top2[n])
            if kind=='low': return 1-sigma(low2[n])
            if kind=='topxlow': return (sigma(top2[n])+sigma(n-top2[n])+1)%2
        return 0
    return f
for kind in ['top','low','topxlow']:
    f=rule_factory(kind); cache={}
    def col(n):
        r=cache.get(n)
        if r is None: r=cache[n]=f(n)
        return r
    bad=[]; tot=0
    for a,b,c in itertools.combinations(S,3):
        if a+b+c>N: continue
        tot+=1
        ca=col(a)
        if col(b)!=ca or col(c)!=ca: continue
        if col(a+b)==ca and col(a+c)==ca and col(b+c)==ca and col(a+b+c)==ca:
            bad.append((a,b,c))
            if len(bad)>=6: break
    print(f"rule {kind:8s}: mono FS(a,b,c): {bad[:6]} (checked {tot})",flush=True)
