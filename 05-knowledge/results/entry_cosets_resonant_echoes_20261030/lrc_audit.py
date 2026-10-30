#!/usr/bin/env python3
"""Exact safe-window intersections, independent pair-clock maxima, echo controls."""
from fractions import Fraction as F
from pathlib import Path
from math import gcd,floor,ceil
from itertools import product
import importlib.util,json,random
ROOT=Path(__file__).resolve().parents[3]
pin=ROOT/'05-knowledge/results/session114_continuation_symmetry_20261030/lrc/check_resonant_tail.py'
spec=importlib.util.spec_from_file_location('pair_clocks',pin);pair=importlib.util.module_from_spec(spec);spec.loader.exec_module(pair)
def need(b,msg):
    if not b:raise RuntimeError(str(msg))
def safe(V,r):
    intervals=[(F(0),F(1))]
    for v in V:
        out=[]
        for L,R in intervals:
            for k in range(max(0,ceil(v*L-1+r)),min(v-1,floor(v*R-r))+1):
                l=max(L,F(k+r,v));u=min(R,F(k+1-r,v))
                if l<=u:out.append((l,u))
        intervals=out
    return sorted(set(intervals))
def params(q,m,J):
    s=next(i for i in range(2,q) if gcd(i,q)==1);a=m;b=0;rows=[]
    for j in range(J+1):
        D=q*a+s;c=a-s*b;rho=F(c,D);x=F(q*b+1,D)
        rows.append(dict(a=a,b=b,c=c,D=D,rho=rho,x=x))
        a,b=q*a+s,q*b+1
    return s,rows
def main():
    exact=[];structural=0;pair_checks=0
    for q in range(3,41):
        s=next(i for i in range(2,q) if gcd(i,q)==1)
        for m in (s//2+1,s//2+3):
            _,rows=params(q,m,8);V=list(range(1,q))+[q+1]
            limit=F(m*(q-1),q*(m*(q-1)+s));threshold=F(m,(q+1)*m+s)
            need(limit>threshold and limit>=F(1,q+2),'common entry chamber')
            for j,row in enumerate(rows):
                V.append(q*row['a']);rho=row['rho'];x=row['x']
                aa=pow(s,-1,q);t=(F(aa)-x)/q
                need(min(pair.distance(v,t) for v in V)==rho,'literal all-speed witness')
                need(rho>limit and rho>F(1,q+2),'strict top chamber')
                need(all(rho<=old['a']*x-old['b']<F(1,2) for old in rows[:j+1]),'unwrapped joint-tail phases')
                structural+=1
                if q in (3,4,5,6,8,11,12) and j<=2 and m==s//2+1:
                    points=sorted({t,1-t});S=safe(V,rho)
                    need(S==[(p,p) for p in points],'literal complete maximizing set')
                    eps=F(1,8*max(V)**2)
                    need(not safe(V,rho+eps),'no larger event maximum')
                    # Independent opposite-slope clock enumeration on a smaller universe.
                    if q<=6 or j<=1:
                        got,pts,n=pair.maximum(tuple(V));need(got==rho and pts==points,'independent global maximum');pair_checks+=1
                    exact.append(dict(q=q,m=m,j=j,speeds=list(V),M=str(rho),points=list(map(str,points))))
    changing=changing_global=0
    for q in (3,4,5,6,8,11):
        s=next(t for t in range(2,q) if gcd(q,t)==1);m=s//2+1
        for seq in product((1,2,7),repeat=3):
            aa=m;bb=0;P=1;H=F(1);V=list(range(1,q))+[q+1,q*aa];states=[(aa,bb)]
            for K in seq:
                aa,bb=K*(q*aa+s),K*(q*bb+1);P*=q*K;H+=F(1,P);V.append(q*aa);states.append((aa,bb))
                rho=F(m,q*m+s*H);x=H/(q*m+s*H)
                need(rho==F(aa-s*bb,q*aa+s),'changing coefficient identity')
                need(all(rho<=a0*x-b0<F(1,2) for a0,b0 in states),'changing old tails stay safe')
                t=(F(pow(s,-1,q))-x)/q
                need(min(pair.distance(v,t) for v in V)==rho,'changing actual witness');changing+=1
            if q<=4:
                need(safe(V,rho)==[(z,z) for z in sorted({t,1-t})],'changing complete maximum set')
                need(not safe(V,rho+F(1,8*max(V)**2)),'changing no larger maximum');changing_global+=1
    hostile=[]
    for V,want in [((1,6),F(3,7)),((1,6,7),F(4,13)),((1,2,3,12),F(3,13))]:
        got,pts,n=pair.maximum(V);need(got==want,'balance hostile exact maximum')
        need(safe(V,want)==[(x,x) for x in pts],'hostile independent windows')
        hostile.append(dict(speeds=V,M=str(got),points=list(map(str,pts))))
    rng=random.Random(20261030);arrival_checks=0
    for q in range(3,15):
        s=next(i for i in range(2,q) if gcd(i,q)==1)
        for trial in range(10):
            M=tuple(sorted(rng.sample(range(1,61),3)));V=tuple(range(1,q))+(q+1,)+tuple(q*m for m in M)
            for lam in (F(1,4),F(1,2),F(3,4),F(1)):
                rho=F(1,q+2)+lam*(F(1,q)-F(1,q+2));T=safe(M,rho)
                accepted=bool(T) and q*rho+s*T[0][0]<=1
                need(bool(safe(V,rho))==accepted,'joint first-arrival equivalence');arrival_checks+=1
    tau=safe((1,3),F(3,10))[0][0]
    need(tau==F(13,30)>F(3,10),'individual arrivals do not synchronize')
    out=dict(status='all-parameter proof in LRC.md; these are exact finite controls',
             structural_cases=structural,complete_window_cases=len(exact),independent_pair_clock_cases=pair_checks,
             changing_generator_checks=changing,changing_complete_windows=changing_global,
             joint_first_arrival_checks=arrival_checks,joint_arrival_hostile=str(tau),exact=exact,hostiles=hostile)
    Path(__file__).with_name('lrc_audit.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='exact'},indent=2))
if __name__=='__main__':main()
