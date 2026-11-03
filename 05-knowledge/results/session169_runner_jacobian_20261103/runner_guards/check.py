#!/usr/bin/env python3
"""Exact controls for the two-guard runner lift; standard library only."""
from fractions import Fraction as F
from itertools import combinations, product
from math import gcd
from functools import reduce
import json

GRAY=((0,0),(1,0),(1,1),(0,1))
DELTA=F(1,8)

def require(test, message):
    if not test: raise RuntimeError(message)

def bits(x):
    x=x%1
    # Closed semicircles rooted at 1/8 and 3/8; guards win ties.
    return tuple(int((x-a)%1 <= F(1,2)) for a in (F(1,8),F(3,8)))

def distance(x):
    x=x%1
    return min(x,1-x)

def eventset(v):
    return {F(k,8*v) for k in range(1,8*v,2)}

def v2(v):
    a=0
    while v%2==0: a+=1;v//=2
    return a

def divisors(v):
    return {d for d in range(1,v+1) if v%d==0}

def phi(v):
    return sum(gcd(k,v)==1 for k in range(1,v+1))

def formula_count(vs):
    groups={}
    for v in vs:
        a=v2(v);groups.setdefault(a,set()).update(divisors(v//(2**a)))
    return sum(4*2**a*sum(phi(d) for d in ds) for a,ds in groups.items())

def direct_state(vs,t): return tuple(bits(v*t) for v in vs)
def safe(state): return all(a or b for a,b in state)

def sweep(vs):
    es=sorted(set().union(*(eventset(v) for v in vs)))
    require(len(es)==formula_count(vs),'event count formula')
    bounds=[F(0)]+es+[F(1)]
    rank=[0]*len(vs)
    steps=0
    isolated=[]
    safe_events=[]
    strict_intervals=[]
    for j,e in enumerate(es,1):
        left=(bounds[j-1]+e)/2;right=(e+bounds[j+1])/2
        before=direct_state(vs,left);after=direct_state(vs,right)
        exact=direct_state(vs,e)
        require(before==tuple(GRAY[x%4] for x in rank),'pre-event Gray state')
        active=[i for i,v in enumerate(vs) if e in eventset(v)]
        for i in active: rank[i]+=1
        require(after==tuple(GRAY[x%4] for x in rank),'post-event Gray state')
        require(exact==tuple(tuple(a|b for a,b in zip(x,y)) for x,y in zip(before,after)),'closed event OR')
        middle=[list(x) for x in before]
        increasing=[];decreasing=[]
        for i in active:
            diff=[b for b in range(2) if before[i][b]!=after[i][b]]
            require(len(diff)==1,'one guard edge per active runner')
            b=diff[0]
            (increasing if after[i][b] else decreasing).append((i,b))
        for i,b in increasing: middle[i][b]=1;steps+=1
        require(tuple(map(tuple,middle))==exact,'serialized boundary state')
        for i,b in decreasing: middle[i][b]=0;steps+=1
        require(tuple(map(tuple,middle))==after,'serialized final state')
        require(safe(exact)==all(distance(v*e)>=DELTA for v in vs),'threshold detector')
        if safe(exact):
            safe_events.append(str(e))
            if not safe(before) and not safe(after):
                isolated.append(str(e))
                entering=[i for i,v in enumerate(vs) if (v*e)%1==DELTA]
                exiting=[i for i,v in enumerate(vs) if (v*e)%1==1-DELTA]
                require(entering and exiting,'opposing boundaries at isolated witness')
                require(all((vs[i]+vs[j])//gcd(vs[i],vs[j])%8==0 for i in entering for j in exiting),'resonant pair')
    for a,b in zip(bounds,bounds[1:]):
        t=(a+b)/2
        if all(distance(v*t)>DELTA for v in vs):strict_intervals.append([str(a),str(b)])
    require(rank==[4*v for v in vs],'deck translation')
    require(steps==4*sum(vs),'physical arc flip count')
    nonresonant=not any((a+b)//gcd(a,b)%8==0 for a,b in combinations(vs,2))
    require(not(nonresonant and safe_events and not strict_intervals),'strictness upgrade')
    return dict(speeds=vs,events=len(es),physical_flips=steps,deck_translation=rank,safe_events=safe_events,isolated_witnesses=isolated,strict_intervals=strict_intervals)

def hamilton_count(n,mask):
    out=[0]*n;j=0
    for a in range(n):
        for b in range(a+1,n):
            x,y=(a,b) if (mask>>j)&1 else (b,a)
            out[x]|=1<<y;j+=1
    dp={}
    for i in range(n):dp[1<<i,i]=1
    for s in range(1,1<<n):
        for last in range(n):
            val=dp.get((s,last),0)
            if not val:continue
            avail=out[last]&~s
            while avail:
                b=avail&-avail;avail-=b
                key=(s|b,b.bit_length()-1)
                dp[key]=dp.get(key,0)+val
    return sum(dp.get(((1<<n)-1,j),0) for j in range(n))

def tournament_mask(state):
    # Labels 0=A, 1=B, 2,... are runners; A->B, transitive runners.
    n=len(state)+2;mask=0;j=0
    for a in range(n):
        for b in range(a+1,n):
            directed=(state[b-2][a] if a<2 and b>=2 else 1)
            mask|=int(directed)<<j;j+=1
    return mask

def main():
    grid_checks=0
    for den in range(1,101):
        for num in range(den):
            x=F(num,den)
            require(any(bits(x))==(distance(x)>=DELTA),'direct detector grid')
            grid_checks+=1
    event_intersections=0
    for v in range(1,65):
        for w in range(v,65):
            want=eventset(gcd(v,w)) if v2(v)==v2(w) else set()
            require(eventset(v)&eventset(w)==want,'event gcd intersection')
            event_intersections+=1
    examples=[(1,2,3,4,5,6,7),(7,13),(1,7,13,91),(1,2,4,8),(1,3,4,5,7,13,18)]
    reports=[sweep(vs) for vs in examples]
    sweep_checks=0
    for n in range(1,5):
        for vs in combinations(range(1,13),n):
            sweep(vs);sweep_checks+=1
    q_checks=0
    for n in range(1,5):
        states=list(product(GRAY,repeat=n))
        masks=[tournament_mask(s) for s in states]
        require(len(set(masks))==4**n,'marked injection')
        # Neighbor calculation is complete, without all-pairs quadratic work.
        for s,m in zip(states,masks):
            neighbors=[]
            for i in range(n):
                for b in range(2):
                    z=[list(x) for x in s];z[i][b]^=1
                    z=tuple(map(tuple,z));neighbors.append(tournament_mask(z))
            require(len(set(neighbors))==2*n,'cube degree')
            require(all((m^x).bit_count()==1 for x in neighbors),'literal single flip')
            q_checks+=1
    profile_blind=[]
    for n in range(1,8):
        zero=((0,0),)*n;one=((1,1),)*n
        require(not safe(zero) and safe(one),'opposed predicate')
        hz=hamilton_count(n+2,tournament_mask(zero));ho=hamilton_count(n+2,tournament_mask(one))
        require(hz==ho==1,'transitive profile-blind witnesses')
        profile_blind.append(dict(runners=n,order=n+2,H=hz))
    require(reports[0]['isolated_witnesses']==['1/8','3/8','5/8','7/8'],'tight family witnesses')
    require(not reports[0]['strict_intervals'],'tight family no strict witness')
    require(reports[1]['events']==76,'7x13 event count')
    require(eventset(7)&eventset(13)=={F(1,8),F(3,8),F(5,8),F(7,8)},'7x13 common events')
    print(json.dumps(dict(status='PASS',grid_points=grid_checks,event_pair_checks=event_intersections,sweep_tuples=sweep_checks,cube_vertices_checked=q_checks,profile_blind=profile_blind,examples=reports),indent=2,sort_keys=True))

if __name__=='__main__':main()
