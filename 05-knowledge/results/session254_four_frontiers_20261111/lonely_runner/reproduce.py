#!/usr/bin/env python3
"""Exact rational wall-union tests of Session254's new shifted-pack claims."""
from fractions import Fraction as F
from math import floor, ceil
import json, random
from pathlib import Path

def safe_measure(v, phases, delta, lo=F(0), hi=F(1)):
    intervals=[]
    for a,p in zip(v,phases):
        for k in range(floor(a*lo+p-delta)-1,ceil(a*hi+p+delta)+2):
            l=max(lo,(k-p-delta)/a); r=min(hi,(k-p+delta)/a)
            if l<r: intervals.append((l,r))
    intervals.sort(); total=F(0); end=lo
    for l,r in intervals:
        if r>end: total+=r-max(l,end); end=r
    return hi-lo-total

def independent_measure(v, phases, delta):
    walls={F(0),F(1)}
    for a,p in zip(v,phases):
        for k in range(floor(p)-1,ceil(a+p)+2):
            for s in (-1,1):
                x=(k-p+s*delta)/a
                if 0<x<1: walls.add(x)
    walls=sorted(walls); total=F(0)
    for l,r in zip(walls,walls[1:]):
        t=(l+r)/2
        if all(min((a*t+p)%1,1-(a*t+p)%1)>=delta for a,p in zip(v,phases)):
            total+=r-l
    return total

def check(condition,label):
    if not condition: raise RuntimeError(label)

def main():
    rng=random.Random(254); bound=F(311,18750); delta=F(2,25)
    carriers=[list(range(1,8)),[1,2,4,8,16,32,64],[2,3,5,7,11,13,17],list(range(100,107)),[1,2,2,2,2,2,2],[2,2,2,2,2,2,3]]
    rows=[]; count=0
    for a in carriers:
        best=None; bestphase=None
        phasesets=[[F(j,7) for j in range(7)],[F(0)]*7]
        phasesets += [[F(rng.randrange(91),91) for _ in a] for _ in range(24)]
        for phases in phasesets:
            p=safe_measure(a,phases,delta); count+=1
            check(p>bound,'uniform rational safe-mass lower bound')
            if best is None or p<best: best,bestphase=p,phases
        check(independent_measure(a,bestphase,delta)==best,'independent midpoint replay')
        rows.append({'slopes':a,'tested_phase_rows':len(phasesets),'minimum_measure':str(best),'phases':[str(x) for x in bestphase]})
    hostile=[F(j,7) for j in range(7)]
    check(safe_measure([1]*7,hostile,F(1,14))==0,'critical repeated-slope tiling hostile')
    check(safe_measure([1]*7,hostile,delta)==0,'robust repeated-slope hostile')
    # Six core speeds, their exact protected interval about t=1/7,
    # and a centered block of length 1/h within it.
    core=list(range(1,7)); M=6; t0=F(1,7)
    attachment=[]
    for h,a,b in [(42,list(range(1,8)),[0]*7),(175,list(range(1,8)),[-3,-2,-1,0,1,2,3]),(350,list(range(1,8)),[3,-3,2,-2,1,-1,0]),(175,[1,2,2,2,2,2,2],[0,-3,-2,-1,1,2,3])]:
        B=max(map(abs,b)); lo=t0-F(1,2*h); hi=t0+F(1,2*h)
        check(h>=7*M and 3*h>=175*B,'attachment hypotheses')
        speeds=core+[h*x+y for x,y in zip(a,b)]
        p=safe_measure(speeds,[F(0)]*13,F(1,14),lo,hi)
        check(p>bound/h,'physical attached block mass')
        attachment.append({'h':h,'offsets':b,'speeds':speeds,'block':[str(lo),str(hi)],'safe_measure':str(p),'certified_lower_bound':str(bound/h)})
    report={'status':'FINITE-EXACT new-claim controls only','phase_rows':count,'threshold':'2/25','proved_lower_bound':str(bound),'carrier_tests':rows,'attachment_tests':attachment,'repeated_slope_hostile':'zero safe measure at 1/14 and empty safe set at 2/25','independent_midpoint_replays':len(rows)}
    print(json.dumps(report,indent=2,sort_keys=True))
if __name__=='__main__': main()
