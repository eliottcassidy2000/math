#!/usr/bin/env python3
"""Consecutive-N event scan, independent of the full B-chain construction."""
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
import json
import check as core

def run():
    p=Path(__file__).resolve().parent
    inherited=json.loads((p/'receipt.json').read_text())
    controls={x['N']:x for x in inherited['root']['scales']}
    V={int(l.split()[0]) for l in (p/'frontier.txt').read_text().splitlines() if l and not l.startswith('#')}
    memo={s:(s,0) for s in V}
    def route(s):
        trail=[];seen=set();x=s
        while x not in memo:
            core.check(x not in seen and len(trail)<10000,'unresolved source cannot disappear')
            seen.add(x);y,t=core.edge(x);trail.append((x,t));x=y
        u,t=memo[x]
        for x,dt in reversed(trail):t+=dt;memo[x]=(u,t)
        return memo[s]
    last=65536;K=(2*last-1).bit_length()
    bins=[defaultdict(lambda:[0,0]) for _ in range(K+1)];E=[0]*(K+1)
    charges=defaultdict(lambda:[0,0]);S=0
    upper=((0,1),(-1,-1),(1,0));change=((1,1),(-1,0),(0,-1))
    keys=('binary','fourary','eightary','adaptive_1to3')
    peaks={k:(0,1,0) for k in keys};violations={k:0 for k in keys}
    for N in range(1,last+1):
        for m,delta in ((4*N-3,upper[(4*N-3)%3]),(4*N-1,upper[(4*N-1)%3]),(2*N-1,change[(2*N-1)%3])):
            s=m
            while s%8==5:s=(s-1)//4
            u,t=route(s);a,b=delta;old=core.Q(charges[s]);charges[s][0]+=a;charges[s][1]+=b
            S+=core.Q(charges[s])-old;sign=1 if t%2==0 else -1
            for j in range(K+1):
                v=bins[j][u,((s-1)//2)%(1<<j)];old=core.Q(v)
                v[0]+=sign*a;v[1]+=sign*b;E[j]+=core.Q(v)-old
        core.check(E[-1]==S,'source norm at separation')
        # Stop at the source's own first separating scale; do not permit empty
        # padding levels to change the controller's declared length.
        stop=(2*N-1).bit_length();profile=E[:stop+1];B=core.budgets(profile)
        if N in controls:
            C=controls[N]
            core.check(S==C['S'] and profile==C['E'] and B==C['budgets'],'event engine versus independent B-chain control')
        if N>=4096:
            for k in keys:
                d=B[k];v,s,_=peaks[k]
                if d*s>v*S:peaks[k]=(d,S,N)
                if 2*d>S:violations[k]+=1
    return dict(status='PASS',first_N=4096,last_N=last,consecutive_scales=last-4096+1,
                violations_of_half_source=violations,
                peaks={k:dict(debt=d,source=S,N=N,ratio=str(Fraction(d,S))) for k,(d,S,N) in peaks.items()},
                scope='complete first arrivals only; selected finite depth and clock cuts are in receipt.json')

if __name__=='__main__':print(json.dumps(run(),indent=2,sort_keys=True))
