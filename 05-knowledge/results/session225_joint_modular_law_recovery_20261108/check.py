#!/usr/bin/env python3
"""Portable exact regeneration, Python standard library only.
Run: python3 check.py [q ...]; defaults 3 5 7 9 15 21.
No file or network access. Larger q use O(q^3) state memory.
"""
import sys
from fractions import Fraction as R
from math import gcd
from collections import defaultdict, Counter

def run(q):
    assert q >= 1 and q % 2 == 1
    S = [(D,E,U,d) for d in range(q) for E in sorted({0,d})
         for D in range(q) for U in range(q)]
    ix = {s:i for i,s in enumerate(S)}
    def O(s):
        D,E,U,d=s
        return D,0,2*U%q,d
    def Z(s):
        D,E,U,d=s
        return (D-E+U)%q,d,(2*U-d)%q,d
    def J(s):
        D,E,U,d=s
        return (D-E+U)%q,(D-E)%q,(D+U)%q,D
    edges = [[ix[O(s)],ix[Z(s)]] for s in S]
    rev = [[] for _ in S]
    for i,es in enumerate(edges):
        for j in es: rev[j].append(i)
    assert len(S)==q*q*(2*q-1)
    assert all(len(es)==2 for es in rev)
    seen=set(); order=[]
    for i in range(len(S)):
        if i in seen: continue
        seen.add(i); stack=[(i,0)]
        while stack:
            x,k=stack.pop()
            if k==2:
                order.append(x); continue
            stack.append((x,k+1)); y=edges[x][k]
            if y not in seen:
                seen.add(y); stack.append((y,0))
    cid={}; comps=[]
    for i in reversed(order):
        if i in cid: continue
        c=len(comps); comp=[]; stack=[i]; cid[i]=c
        while stack:
            x=stack.pop(); comp.append(x)
            for y in rev[x]:
                if y not in cid:
                    cid[y]=c; stack.append(y)
        comps.append(comp)
    assert all(cid[i]==cid[j] for i,es in enumerate(edges) for j in es)
    for comp in comps:
        D,E,U,d=S[comp[0]]
        g=gcd(gcd(d,U),q); c=D%g; Q=q//g
        fs=Counter((S[i][0]-S[i][1]+S[i][2])%q for i in comp)
        assert fs=={r:len(comp)//Q for r in range(q) if r%g==c}
        e0=[S[i] for i in comp if S[i][1]==0]
        vs=Counter(((D-c)//g-pow(2,-1,Q)*(U//g))%Q
                   for D,E,U,d in e0) if Q>1 else Counter({0:len(e0)})
        assert len(vs)==Q and set(vs.values())=={len(e0)//Q}
    def project(dist,doJ=False):
        mass=defaultdict(R)
        for s,w in dist.items():
            if doJ: s=J(s)
            for z in (O(s),Z(s)): mass[cid[ix[z]]]+=w/2
        out={}
        for c,w in mass.items():
            for i in comps[c]: out[S[i]]=w/len(comps[c])
        return out
    start=project({(1%q,0,2%q,1%q):R(1)})
    assert len(start)==(2*q*q if q>1 else 1)
    assert set(start.values())=={R(1,len(start))}
    first=project(start,True)
    for D,E,U,d in S:
        g=gcd(gcd(d,U),q)
        w=R(g,q**3)*(R(3,4)*(D%g==0)+R(1,4)*((D+1)%g==0))
        expected=w/(2 if d else 1)
        assert first.get((D,E,U,d),0)==expected
    second=project(first,True)
    actual=[R(0) for _ in range(q)]
    for (D,E,U,d),w in second.items(): actual[(D-E+U)%q]+=q*w
    numer=[0]*q
    # Integer accumulation of the four branches: denominator 16*q^3.
    for d in range(q):
        for U in range(q):
            g=gcd(gcd(d,U),q)
            for D in range(q):
                h=gcd(gcd(D,U),q)
                a=g*(3*(D%g==0)+((D+1)%g==0))
                for r in range(q):
                    numer[r]+=a*h*(3*(r%h==0)+((r+d)%h==0))
    theory=[R(n,16*q**3) for n in numer]
    assert theory==actual and sum(actual)==q
    print("q=%d states=%d SCCs=%d PASS first-density, second-law, every-SCC output and line; L(0)=%s"
          %(q,len(S),len(comps),actual[0]))
    if q==63:
        chosen=[c for c in comps if S[c[0]][3]==21
                and gcd(gcd(S[c[0]][2],21),q)==1]
        assert sorted(map(len,chosen))==[756]*6
        print("q=63 d=21 g=1: six SCCs of size 756 PASS")
if __name__=="__main__":
    for q in map(int,sys.argv[1:] or [3,5,7,9,15,21]): run(q)
