#!/usr/bin/env python3
"""Three phase lifts of the actual source, with identical first-arrival routes.

Exact integer coefficients and Gram matrices.  The other two phase rows are
hostile controls, not other Collatz trajectories. Every route has a cap;
failure to enter is reported, never silently discarded.
"""
import argparse
import json
from pathlib import Path
import numpy as np


def seed(x):
    while x % 8 == 5:
        x = (x - 1) // 4
    return x


def edge(x):
    a = 2 if x % 8 == 1 else 1
    y = (3*x + 1) >> a
    t = a
    while y % 8 == 5:
        y = (y - 1)//4
        t += 2
    return y, t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--max-power', type=int, default=19)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    here = Path(__file__).parent
    V = sorted(int(s.split()[0]) for s in (here/'first_arrival_20261018/V_closed_Y1000.txt').read_text().splitlines()
               if s.strip() and not s.startswith('#'))
    assert len(V) == 627 and all(edge(x)[0] in V for x in V)
    vi = {x:i for i,x in enumerate(V)}
    memo = {x:(vi[x], 0, 0) for x in V}
    def route(x):
        path=[]; seen=set()
        while x not in memo:
            if x in seen or len(path) >= 10000:
                raise RuntimeError(('unresolved capped route', x))
            seen.add(x); y,t=edge(x); path.append((x,t)); x=y
        u,t,d=memo[x]
        for x,dt in reversed(path):
            t+=dt; d+=1; memo[x]=(u,t,d)
        return memo[path[0][0]] if path else (u,t,d)
    nmax=2**args.max_power
    src=np.array([x for x in range(1,4*nmax,2) if x%8!=5],dtype=np.int64)
    routes=np.array([route(int(x)) for x in src],dtype=np.int64)
    maxclock=int(routes[:,1].max())+1
    # Independent literal verification of a distributed sample of seed paths.
    checked=0
    for i in np.linspace(0,len(src)-1,min(1000,len(src)),dtype=int):
        x=int(src[i]); t=0; d=0
        while x not in vi:
            x,dt=edge(x); t+=dt; d+=1
        assert tuple(routes[i])==(vi[x],t,d)
        checked+=1
    rows=[]
    for k in range(3,args.max_power+1):
        N=2**k
        sources=src[src<4*N]
        rr=routes[:len(sources)]
        c=np.zeros((len(sources),3),dtype=np.int64)
        # B^j(s) contributes its actual lower/upper-band weight, with the
        # residue rotated by r only in the two auxiliary rows.
        m=sources.copy(); active=np.ones(len(m),dtype=bool)
        while active.any():
            inds=np.flatnonzero(active)
            for r in range(3):
                res=(m[inds]+r)%3
                w=np.where(res==2,2,-1)
                low=np.array([0,-3,3],dtype=np.int64)[res]
                c[inds,r]+=np.where(m[inds]<2*N,low,w)
            m[inds]=4*m[inds]+1
            active=m<4*N
        assert np.all(c.sum(axis=1)==0)
        arrival=np.zeros((len(V),3),dtype=np.int64)
        signed=c*np.where(rr[:,1]%2==0,1,-1)[:,None]
        for r in range(3): np.add.at(arrival[:,r],rr[:,0],signed[:,r])
        byclock=np.zeros((len(V)*maxclock,3),dtype=np.int64)
        for r in range(3): np.add.at(byclock[:,r],rr[:,0]*maxclock+rr[:,1],signed[:,r])
        G=arrival.T@arrival; H=byclock.T@byclock; E=c.T@c
        assert np.all(G.sum(axis=0)==0) and np.all(H.sum(axis=0)==0)
        row=dict(N=N,source_gram=E.tolist(),arrival_gram=G.tolist(),clock_diagonal_gram=H.tolist(),
                 energies_per_N=(np.diag(G)/N).tolist(),clock_energies_per_N=(np.diag(H)/N).tolist(),
                 max_depth=int(rr[:,2].max()),max_clock=int(rr[:,1].max()))
        rows.append(row)
        print(N,row['energies_per_N'],row['clock_energies_per_N'],flush=True)
    result=dict(status='FINITE-EXACT',universe=f'all primitive seeds < 4*2^{args.max_power}',
                frontier=V,route_cap=10000,independent_route_checks=checked,rows=rows,
                scope='Auxiliary source weights rotate residues; routes, ports, and integer clocks are unchanged. No universal entry claim.')
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
