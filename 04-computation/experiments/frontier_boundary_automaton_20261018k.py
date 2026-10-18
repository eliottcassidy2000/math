#!/usr/bin/env python3
"""Complete fixed-frontier source sweep, with a three-scale stabilization audit.

The accompanying proof shows that all event orders above R repeat under
scaling by four, and all source vectors repeat after three such event blocks.
This turns the finite maximum below 64R into an all-N boundary bound. There
is no first-arrival or Collatz convergence calculation in this program.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from source_pair_tensor_20261018k import coefficients


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
    here=Path(__file__).parent
    V=np.array([int(s.split()[0]) for s in (here/'first_arrival_20261018/V_closed_Y1000.txt').read_text().splitlines()
                if s.strip() and not s.startswith('#')],dtype=np.int64)
    assert len(V)==627 and all(s%8!=5 for s in V)
    R=1 << int(V.max()).bit_length();end=256*R
    events={}
    for i,s in enumerate(V):
        m=int(s);j=0
        while m<4*end:
            for kind,n in (('E',m//4+1),('L',m//2+1)):
                if n<end:events.setdefault(n,[]).append((i,kind,j,m))
            m=4*m+1;j+=1
    sorted_events=sorted(events)
    c=np.zeros((len(V),3),dtype=np.int64)
    kap=np.array([0,-3,3]);lam=np.array([-1,-1,2])
    maximum=(-10**9,None);minimum=(10**9,None)
    profiles=[[] for _ in range(4)];states=[];witness=None;literal_checks=0
    sample=set(sorted_events[::max(1,len(sorted_events)//1024)])
    for n in sorted_events:
        if n in (R,4*R,16*R,64*R):raise AssertionError('unexpected block boundary event')
        for i,kind,j,m in events[n]:
            residues=(m+np.arange(3))%3
            c[i]+=lam[residues] if kind=='E' else kap[residues]-lam[residues]
        T=int((c*c).sum());E=int(c[:,0]@c[:,0]);D=3*E-T
        if n<64*R:
            if D>maximum[0]:maximum=(D,n);witness=c.tolist()
            if D<minimum[0]:minimum=(D,n)
        if n in sample:
            assert np.array_equal(c,coefficients(V,n));literal_checks+=len(V)
        if n>=R:
            k=0
            while n>=R*4**(k+1):k+=1
            assert len(events[n])==1
            i,kind,j,m=events[n][0]
            profiles[k].append((int(V[i]),kind,j,n,D))
            if k==0:states.append(c.copy())
            if k==3:
                assert np.array_equal(c,states[len(profiles[k])-1])
    assert all(len(p)==2*len(V) for p in profiles)
    for k in range(3):
        for old,new in zip(profiles[k],profiles[k+1]):
            s,kind,j,n,D=old
            assert new[:3]==(s,kind,j+1)
            assert new[3]==4*n-(2 if kind=='E' else 1)
    assert [x[-1] for x in profiles[0]]==[x[-1] for x in profiles[3]]
    for n in (R,4*R,16*R,64*R):
        assert np.array_equal(coefficients(V,64*n),coefficients(V,n))
    assert np.array_equal(np.array(witness),coefficients(V,maximum[1]))
    out=dict(status='FINITE-EXACT certificate for the proved event-stabilization reduction',
             frontier_size=len(V),largest_seed=int(V.max()),R=R,
             complete_transient_and_period_end=64*R,audit_end=end,
             event_count=len(sorted_events),events_per_stable_block=[len(p) for p in profiles],
             literal_source_checks=literal_checks,three_block_state_checks=len(states),
             maximum_D0=maximum[0],first_maximum_N=maximum[1],minimum_D0=minimum[0],first_minimum_N=minimum[1],
             block_D0_ranges=[[min(x[-1] for x in p),max(x[-1] for x in p)] for p in profiles],
             stable_profiles=profiles[:3],
             scope='Boundary source only. No bound on transported arrivals follows without the OPEN orientation inequality.')
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(out,indent=2)+'\n')
    print({k:v for k,v in out.items() if k!='stable_profiles'})


if __name__=='__main__':main()
