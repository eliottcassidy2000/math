#!/usr/bin/env python3
"""Exact independent interval/cell/candidate-clock audits; no third-party packages."""
from fractions import Fraction as F
from itertools import combinations_with_replacement
from math import gcd
import json

C=tuple(range(1,12))+(13,)
checks=0

def check(condition,label):
    global checks
    checks+=1
    if not condition: raise RuntimeError(str(label))

def distance(v,t):
    x=(v*t)%1
    return min(x,1-x)

def survivor(V,r):
    out=[(F(0),F(1))]
    for v in V:
        teeth=[((k+r)/v,(k+1-r)/v) for k in range(v)]
        new=[]; i=j=0
        while i<len(out) and j<len(teeth):
            a,b=out[i]; c,d=teeth[j]
            l,h=max(a,c),min(b,d)
            if l<=h:
                if new and l<=new[-1][1]:new[-1]=(new[-1][0],max(h,new[-1][1]))
                else:new.append((l,h))
            if b<d:i+=1
            elif d<b:j+=1
            else:i+=1;j+=1
        out=new
    return out

def cell_survivor(V,r):
    walls={F(0),F(1)}
    for v in V:
        for k in range(v):
            walls.update(((k+r)/v,(k+1-r)/v))
    w=sorted(walls)
    out=[]
    def add(a,b):
        if out and a<=out[-1][1]:out[-1]=(out[-1][0],max(b,out[-1][1]))
        else:out.append((a,b))
    for i,a in enumerate(w):
        if all(distance(v,a)>=r for v in V):add(a,a)
        if i+1<len(w):
            b=w[i+1]
            if all(distance(v,(a+b)/2)>=r for v in V):add(a,b)
    return out

def top_windows(r):
    I=[((1+r)/13,(1-r)/11),((2+r)/5,(3-r)/7)]
    return sorted(I+[(1-b,1-a) for a,b in I])

def maximum(V):
    # Every local maximum of a finite minimum of tent waves has active
    # opposite slopes v and -w, with v=w allowed for a tent's own apex.
    # Their intersection has t=a/(v+w). Score directly with exact integers.
    V=tuple(sorted(set(V))); best=F(0); times=set(); visits=0
    for d in sorted({v+w for v,w in combinations_with_replacement(V,2)}):
        for a in range(d):
            mn=d
            for v in V:
                res=v*a%d; mn=min(mn,res,d-res)
                if mn*best.denominator < best.numerator*d:break
            else:
                value=F(mn,d)
                if value>best:best=value;times={F(a,d)}
                elif value==best:times.add(F(a,d))
            visits+=1
    return best,sorted(times),visits

def germ_coefficient(V,r,S):
    B=F(0)
    for a,b in S:
        P=[v for v in V if (v*a)%1==r]
        N=[v for v in V if (v*b)%1==1-r]
        check(bool(P) and bool(N),('owners',V,r,a,b))
        B+=F(1,max(P))+F(1,max(N))
    return B

def formula(v):
    if v%12:
        M=F(1,12); times=[F(a,12) for a in (1,5,7,11)]
        cs={1:13,5:5,7:7,11:11}
        B=2*sum((F(1,max(c,v) if v%12==u else c) for u,c in cs.items()),F(0))
    elif v in (12,24):
        M=F(1,14);times=[F(a,14) for a in (1,3,5,9,11,13)];B=F(23324,6435)
    else:
        m=v//12;M=F(m,12*m+5);t=F(5*m+2,12*m+5)
        times=sorted([t,1-t]);B=2*(F(1,5)+F(1,v))
    return M,times,B

def main():
    global checks
    r0=F(1,14);rho=F(1,12)
    expected=sorted(top_windows(r0)+[(F(a,14),F(a,14)) for a in (1,3,5,9,11,13)])
    check(survivor(C,r0)==expected,'core full bottom certificate')
    check(cell_survivor(C,r0)==expected,'independent full wall/cell bottom certificate')
    # Every core speed has fixed, matching floor across each bottom window.
    for a,b in top_windows(r0):
        for v in C:
            check((v*a).numerator//(v*a).denominator==(v*b).numerator//(v*b).denominator,('fixed phase branch',v,a,b))
    for j in range(1,41):
        r=r0+(rho-r0)*F(j,40)
        check(survivor(C,r)==top_windows(r),('all core windows',r))
        check(cell_survivor(C,r)==top_windows(r),('independent core cells',r))
    check(maximum(C)[:2]==(rho,[F(a,12) for a in (1,5,7,11)]),'core global maximum')
    visits=0;rows=[]
    tails=list(range(1,481))+[997,1009,1012,1200]
    for v in tails:
        V=tuple(sorted(set(C+(v,))))
        M,times,B=formula(v)
        got,points,n=maximum(V);visits+=n
        check(got==M,('global maximum',v,got,M))
        check(points==times,('all maximizers',v,points,times))
        S=survivor(V,M)
        check(S==[(t,t) for t in times],('closed maximum survivor',v))
        check(germ_coefficient(V,M,S)==B,('terminal slope',v,B))
        H=max([2]+[(a+b)//gcd(a,b) for a,b in combinations_with_replacement(V,2)])
        eps=F(1,2*M.denominator*H)
        below=survivor(V,M-eps)
        check(sum((b-a for a,b in below),F(0))==eps*B,('actual relaxed mass',v))
        if v%12==0:rows.append({'v':v,'M':str(M),'maximizers':[str(t) for t in times],'B':str(B)})
    for m in (10**3,10**6,10**12,10**30+7):
        r=F(m,12*m+5);t=F(5*m+2,12*m+5)
        check(min(distance(v,t) for v in C+(12*m,))==r,('large m witness',m))
        check(t==F(5,12)-(rho-r)/5,('large m outer endpoint',m))
        check(12*m*(rho-r)/5==r,('large m upper bound equality',m))
    report={'status':'passed','checks':checks,'tail_count':len(tails),'exhaustive_tails':'1..480 plus 997,1009,1012,1200','candidate_clock_visits':visits,'core_algorithms':['exact interval intersection','all exact wall and midpoint cells'],'maximum_algorithm':'all a/(v+w), including v=w; integer modular scores','all_m_proof':'RESONANT_TAIL.md','resonant_rows':rows,'unrestricted_LRC14':'OPEN'}
    print(json.dumps(report,indent=2,sort_keys=True))

if __name__=='__main__':main()
