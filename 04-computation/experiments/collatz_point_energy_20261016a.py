#!/usr/bin/env python3
"""Held-out point refinements of frozen first-reset maps, exact integer tests.

Default universe k=8,14,22. The k22 run tests all stages of the inherited
default-bin stopping policy, not a new controller chosen using the point data.
No nonROOT cycles or source deletions are allowed. C++ only sweeps the same
annular coefficients; small runs compare the prior full point-matrix method.
"""
import argparse
import ctypes
import hashlib
import json
import subprocess
import tempfile
from pathlib import Path
import numpy as np
from collatz_compact_resets_20261015f import Graph, need, row_norms

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'05-knowledge/results'
SOURCE=Path(__file__).with_suffix('.cpp')


def squares(a):
    """Exact sum without overflowing a NumPy accumulator."""
    m=int(np.max(np.abs(a)))
    need(m*m < 2**62, 'individual int64 square bound')
    chunk=min(1<<19, max(1, (2**62-1)//max(1,m*m)))
    return sum(int(np.sum(a[i:i+chunk]*a[i:i+chunk])) for i in range(0,len(a),chunk))


def frozen(k):
    if k>=19:
        p=DATA/f'collatz_compact_resets_20261015f_k{k}_r1.json'
        r=json.loads(p.read_text())
    else:
        p=DATA/'collatz_survivor_packets_20261015e.json'
        r=next(s for s in json.loads(p.read_text())['first_reset'] if s['k']==k)
    return r, dict(path=str(p.relative_to(ROOT)), sha256=hashlib.sha256(p.read_bytes()).hexdigest())


def run(k, kernel):
    X=1<<k; b=1<<(k//2); g=Graph(8*X)
    old, provenance=frozen(k)
    cols,inv=np.unique(g.endpoint[1:],return_inverse=True)
    N=np.arange(X,2*X,dtype=np.int64)
    u0=4*N-3*((2*N+3)//6)+(4*N+1)//6-(2*N+1)//6
    fine=[squares(u0)]; q0=u0.reshape(-1,b).sum(axis=1)
    records=[dict(stage=0, fine_energy=fine[0],
                  coarse_numerator=sum(int(v)**2 for v in q0),
                  variance_numerator=b*fine[0]-sum(int(v)**2 for v in q0),
                  mixture_numerator=0, denominator=b)]
    del N,u0,q0
    for j in range(1,len(old['stages'])+1):
        g.extend(cols);g.close(cols,inv)
        cols,inv=np.unique(g.endpoint[1:],return_inverse=True)
        v,_=row_norms(k,g,inv,len(cols))
        need(v==old['stages'][j-1]['row_norms'], 'same frozen default rows')
        labels=np.empty(8*X+1,dtype=np.int64);labels[0]=0;labels[1:]=inv
        values=np.full(8*X+1,-1,dtype=np.int64);values[2::3]=2
        values[(g.clock&1)!=0]*=-1;values[0]=0
        points=np.empty(X,dtype=np.int64)
        code=kernel(X,len(cols),labels,values,points)
        need(code==0, 'literal halving reduction and all four-Lipschitz increments')
        if k<=8:
            from paired_bin_refinement_20261015g import point_rows
            direct=point_rows(k,g,inv,len(cols))
            need(np.array_equal(np.abs(direct).sum(axis=1),points), 'independent matrix sweep')
        f=squares(points); fine.append(f)
        q=[int(x) for x in points.reshape(-1,b).sum(axis=1)]
        c=sum(x*x for x in v); qq=sum(x*x for x in q)
        variance=b*f-qq; mixture=qq-c
        need(variance>=0 and mixture>=0 and variance+mixture==b*f-c,
             'exact two-price refinement invoice')
        need(3*variance<=4*b*X*(b*b-1),'scalar variance bound')
        # Direct, non-radical sharp Banach path bound, per physical bin.
        # L=4; m_b=(b^2-1)/(3b), t_b=(7b^4-15b^2+8)/(60b^2).
        need(15*b*(b*f-c) <= 40*b*(b*b-1)*sum(v) +
             4*X*(7*b**4-15*b*b+8), 'uniform refinement budget')
        clocks=g.internal_clock_audit()
        g.replay(range(1,8*X+1,max(1,X//32)))
        rec=dict(stage=j, fine_energy=f, coarse_numerator=c, denominator=b,
                 variance_numerator=variance, mixture_numerator=mixture,
                 internal_clock_edges_checked=clocks,max_native_clock=int(g.clock.max()),
                 terminals=len(cols),point_norm_sha256=hashlib.sha256(points.tobytes()).hexdigest())
        records.append(rec)
        print(json.dumps(dict(k=k,stage=j,fine_over_X2=f/X**2,
                              coarse_over_X2=c/b/X**2,mixture_over_X2=mixture/b/X**2)),flush=True)
        del labels,values,points
    pairs=[]
    for j in range(len(fine)-1):
        t=min(j+2,len(fine)-1)
        pairs.append(dict(start=j,end=t,old=fine[j],new=fine[t],
                          passes=fine[t]<=49*X*X or 2*fine[t]<=fine[j]))
    return dict(k=k,bin_width=b,frozen_input=provenance,stages=records,paired_tests=pairs,
                failures=[p for p in pairs if not p['passes']])


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--k',nargs='+',type=int,default=[8,14,22])
    args=ap.parse_args()
    with tempfile.TemporaryDirectory(prefix='collatz-point-') as tmp:
        lib=Path(tmp)/'point.dylib'
        subprocess.run(['clang++','-std=c++17','-O2','-shared','-fPIC',str(SOURCE),'-o',str(lib)],check=True)
        dll=ctypes.CDLL(str(lib));kernel=dll.point_norms
        arr=np.ctypeslib.ndpointer(dtype=np.int64,ndim=1,flags='C_CONTIGUOUS')
        kernel.argtypes=[ctypes.c_int64,ctypes.c_int64,arr,arr,arr];kernel.restype=ctypes.c_int
        for k in args.k:
            r=run(k,kernel)
            out=DATA/f'collatz_point_energy_20261016a_k{k}.json'
            out.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n')
            print('RESULT',k,'failures',r['failures'],'sha256',hashlib.sha256(out.read_bytes()).hexdigest(),flush=True)
