#!/usr/bin/env python3
"""Exact source/sibling tensors and source-preserving arrival controls.

All three rows use the same original seeds, first-arrival ports and integer
clocks. The rotated rows are controls, never replacement Collatz sources.
"""
import argparse
import json
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from frontier_phase_20261018j import edge


def coefficients(s, N):
    c = np.zeros((len(s), 3), dtype=np.int64)
    m = s.copy()
    while np.any(m < 4*N):
        ii = np.flatnonzero(m < 4*N)
        for r in range(3):
            rs = (m[ii]+r) % 3
            c[ii, r] += np.where(m[ii] < 2*N, np.array([0,-3,3])[rs], np.array([-1,-1,2])[rs])
        m[ii] = 4*m[ii]+1
    assert np.all(c.sum(axis=1) == 0)
    return c


def shell_tensor():
    """Independent rational summation of all shell residues, not fitted data."""
    U = lambda r: [(-1,-1,2)[(r+j) % 3] for j in range(3)]
    D = lambda r,L: [sum((0,-3,3)[(r+j+i) % 3] for i in range(L)) for j in range(3)]
    total = [[F(0) for _ in range(3)] for _ in range(3)]
    tables = []
    for L in (1,2,3):
        A = [[sum(U(r)[i]*D(2*r+1,L)[j] for r in range(3)) for j in range(3)] for i in range(3)]
        B = [[sum(D(r,L)[i]*U(2*r+1)[j] for r in range(3)) for j in range(3)] for i in range(3)]
        for i in range(3):
            for j in range(3):
                total[i][j] += F(4**(3-L),63)*(F(A[i][j],12)+F(B[i][j],6))
        tables.append(dict(L_mod_3=L % 3, upper_lower=A, lower_upper=B))
    assert total == [[F(t,7) for t in row] for row in ((5,-1,-4),(-1,-4,5),(-4,5,-1))]
    return dict(tables=tables, limit=[[str(t) for t in row] for row in total])


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--max-power', type=int, default=20)
    ap.add_argument('--small', type=int, default=1024)
    ap.add_argument('--out', type=Path, required=True)
    args=ap.parse_args()
    here=Path(__file__).parent
    V=[int(s.split()[0]) for s in (here/'first_arrival_20261018/V_closed_Y1000.txt').read_text().splitlines()
       if s.strip() and not s.startswith('#')]
    assert len(V)==627 and all(edge(s)[0] in V for s in V)
    vi={s:i for i,s in enumerate(V)}
    memo={s:(i,0,0) for s,i in vi.items()}
    def route(s):
        x=s; seen=set(); path=[]
        while x not in memo:
            if x in seen or len(path)>=10000:
                raise RuntimeError(('capped/unresolved; not discarded',s,x))
            seen.add(x); y,t=edge(x); path.append((x,t));x=y
        u,t,d=memo[x]
        for x,dt in reversed(path):
            t+=dt;d+=1;memo[x]=(u,t,d)
        return memo[s]
    # Independent scalar source calculation, including finite shell endpoints.
    literal_checks=0
    small_max_error=0
    for N in range(1,args.small+1):
        xs=np.arange(1,2*N,8,dtype=np.int64)
        X=coefficients(xs,N);Y=coefficients(2*xs+1,N)
        for arr,ss in ((X,xs),(Y,2*xs+1)):
            for ii,s in enumerate(ss):
                expected=[0,0,0];m=int(s)
                while m<4*N:
                    for r in range(3): expected[r]+=((0,-3,3) if m<2*N else (-1,-1,2))[(m+r)%3]
                    m=4*m+1
                assert expected==arr[ii].tolist();literal_checks+=1
        S=X.T@Y
        small_max_error=max(small_max_error,int(np.abs(7*S-N*np.array([[5,-1,-4],[-1,-4,5],[-4,5,-1]])).max()))
    Ns=sorted(set([1204,4096,5997]+[2**k for k in range(3,args.max_power+1)]))
    rows=[];route_checks=0;packet_checks=0
    for N in Ns:
        ss=np.array([s for s in range(1,4*N,2) if s%8!=5],dtype=np.int64)
        c=coefficients(ss,N); rr=np.array([route(int(s)) for s in ss],dtype=np.int64)
        E=c.T@c
        xs=np.arange(1,2*N,8,dtype=np.int64)
        X=coefficients(xs,N);Y=coefficients(2*xs+1,N)
        S=X.T@Y
        valid=np.array([int(x) not in vi and int(2*x+1) not in vi for x in xs])
        Sv=X[valid].T@Y[valid]
        packet=E-Sv-Sv.T
        # Literal packet construction, not only subtraction of covariance.
        index={int(s):i for i,s in enumerate(ss)}
        raw=c.copy(); pair_ids=[]
        for x in xs[valid]:
            i,j=index[int(x)],index[int(2*x+1)]
            assert tuple(rr[i,[0,2]])==tuple(rr[j,[0,2]]) and rr[j,1]==rr[i,1]+1
            raw[i]=c[i]-c[j];raw[j]=0;pair_ids.append((i,j));packet_checks+=1
        assert np.array_equal(raw.T@raw,packet)
        signed=c*np.where(rr[:,1]%2==0,1,-1)[:,None]
        hist=np.zeros((int(rr[:,2].max())+1,len(V),3),dtype=np.int64)
        for r in range(3):np.add.at(hist[:,:,r],(rr[:,2],rr[:,0]),signed[:,r])
        pref=np.cumsum(hist,axis=0)
        grams=np.einsum('dur,dus->drs',pref,pref)
        trace=np.trace(grams,axis1=1,axis2=2)
        oriented=3*grams[:,0,0]-trace
        # Equal native depth makes the whole polynomial packet survive every
        # first-arrival prefix. Keep t and t+1 as separate integer exponents.
        paired_signed=raw*np.where(rr[:,1]%2==0,1,-1)[:,None]
        hist2=np.zeros_like(hist)
        for r in range(3):np.add.at(hist2[:,:,r],(rr[:,2],rr[:,0]),paired_signed[:,r])
        assert np.array_equal(hist2,hist)
        for ix in np.linspace(0,len(ss)-1,min(128,len(ss)),dtype=int):
            x=int(ss[ix]);t=d=0
            while x not in vi:
                x,dt=edge(x);t+=dt;d+=1
            assert tuple(rr[ix])==(vi[x],t,d);route_checks+=1
        # This is a diagnostic, not a pre-asserted universal cone.
        row=dict(N=N,source_gram=E.tolist(),sibling_tensor=S.tolist(),
                 valid_sibling_tensor=Sv.tolist(),packet_source_gram=packet.tolist(),
                 arrival_gram=grams[-1].tolist(),prefix_grams=grams.tolist(),
                 max_oriented_defect=int(oriented.max()),max_oriented_depth=int(oriented.argmax()),
                 final_oriented_defect=int(oriented[-1]),
                 max_actual_prefix=int(grams[:,0,0].max()),
                 prefix_bound_packet_margin=int(packet[0,0]-grams[:,0,0].max()),
                 max_depth=int(rr[:,2].max()),max_clock=int(rr[:,1].max()),
                 sibling_pairs=int(valid.sum()),frontier_exceptions=int((~valid).sum()))
        rows.append(row)
        print(N,'max oriented defect',row['max_oriented_defect'],'at',row['max_oriented_depth'],
              'prefix packet margin',row['prefix_bound_packet_margin'],flush=True)
    out=dict(status='FINITE-EXACT',shell_proof=shell_tensor(),small_universe=f'all 1 <= N <= {args.small}',
             literal_source_checks=literal_checks,small_max_scaled_tensor_error=small_max_error,
             route_checks=route_checks,packet_checks=packet_checks,rows=rows,
             scope='Exact source and clock-preserving packet identities; no uniform concentration or entry inference.')
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
