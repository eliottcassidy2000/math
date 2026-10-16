#!/usr/bin/env python3
"""Same native maps, all dyadic physical bins, exact finite refinement test."""
import hashlib
import json
from pathlib import Path
import numpy as np
from collatz_compact_resets_20261015f import Graph, need
from collatz_reset_debt_analysis_20261015f import initial_norms

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'05-knowledge/results'


def point_rows(k,g,inv,ncols):
    X=1<<k
    labels=np.empty(8*X+1,dtype=np.int64);labels[0]=0;labels[1:]=inv
    vals=np.full(8*X+1,-1,dtype=np.int64);vals[2::3]=2
    vals[(g.clock&1)!=0]*=-1;vals[0]=0
    row=np.zeros(ncols,dtype=np.int64)
    np.add.at(row,labels[X+1:4*X+1],vals[X+1:4*X+1])
    result=np.empty((X,ncols),dtype=np.int64)
    for j,N in enumerate(range(X,2*X)):
        result[j]=row
        for n in range(4*N+1,4*N+5):row[labels[n]]+=vals[n]
        row[labels[N+1]]-=vals[N+1]
    return result


def run(k):
    X=1<<k;g=Graph(8*X);cols,inv=np.unique(g.endpoint[1:],return_inverse=True)
    old=json.loads((DATA/'collatz_survivor_packets_20261015e.json').read_text())
    rec=next(r for r in old['first_reset'] if r['k']==k)
    N=np.arange(X,2*X,dtype=np.int64)
    l0=4*N-3*((2*N+3)//6)+(4*N+1)//6-(2*N+1)//6
    histories={p:[[int(v) for v in l0.reshape(-1,1<<p).sum(axis=1)]] for p in range(k+1)}
    invoices=[]
    need(histories[k//2][0]==initial_norms(k),'initial same source weights')
    for j in range(1,len(rec['stages'])+1):
        g.extend(cols);g.close(cols,inv)
        cols,inv=np.unique(g.endpoint[1:],return_inverse=True)
        points=point_rows(k,g,inv,len(cols))
        point_norms=[int(v) for v in np.abs(points).sum(axis=1)]
        fine=sum(v*v for v in point_norms)
        for lo in range(0,X-1,128):
            hi=min(lo+128,X-1)
            need(bool(np.all(np.abs(points[lo+1:hi+1]-points[lo:hi]).sum(axis=1)<=4)),
                 'exact three-label halving reduction gives four-Lipschitz rows')
        for p in range(k+1):
            b=1<<p
            rows=points.reshape(X//b,b,len(cols)).sum(axis=1)
            v=[int(v) for v in np.abs(rows).sum(axis=1)]
            histories[p].append(v)
            q=[sum(point_norms[a:a+b]) for a in range(0,X,b)]
            variance=b*fine-sum(a*a for a in q)
            mixture=sum(a*a-c*c for a,c in zip(q,v))
            need(variance>=0 and mixture>=0,'two nonnegative refinement prices')
            need(variance+mixture==b*fine-sum(c*c for c in v),'exact refinement invoice')
            need(3*variance<=4*b*X*(b*b-1),'four-Lipschitz scalar variance bound')
            invoices.append(dict(stage=j,bin_exponent=p,denominator=b,
                                 variance_numerator=variance,mixture_numerator=mixture))
        need(histories[k//2][-1]==rec['stages'][j-1]['row_norms'],'default bins reproduce frozen native rows')
        g.internal_clock_audit()
        g.replay(range(1,8*X+1,max(1,X//32)))
        del points,rows
    result=[];failures=[]
    for p,h in histories.items():
        b=1<<p;e=[sum(v*v for v in row) for row in h]
        for j in range(len(e)-1):
            end=min(j+2,len(e)-1)
            ok=e[end]<=49*b*X*X or 2*e[end]<=e[j]
            row=dict(bin_exponent=p,start=j,end=end,old_numerator=e[j],new_numerator=e[end],denominator=b,passes=ok)
            result.append(row)
            if not ok:failures.append(row)
    # Every represented graph has the same scalar source interval at every p.
    for j in range(len(rec['stages'])+1):
        previous=None
        for p in range(k+1):
            b=1<<p;e=sum(v*v for v in histories[p][j])
            if previous is not None:need(e<=2*previous,'energy decreases under equal-bin coarsening')
            previous=e
    return dict(k=k,policy_stop_stages=len(rec['stages']),comparisons=result,failures=failures,
                refinement_invoices=invoices)


if __name__=='__main__':
    for L in range(1,21):
        fine=2*(L+1)**2;coarse=2*L*L
        need(fine-coarse==4*L+2,'large envelope gap with increment norm only two')
        for z in [-1,1]:
            need((L+z)**2+(L-z)**2-coarse==2,'one shared terminal label has only scalar variance cost')
    out=DATA/'paired_bin_refinement_20261015g.json'
    result=dict(status='FINITE-EXACT same maps, all dyadic bin widths, k8,10,12,14',scales=[run(k) for k in [8,10,12,14]])
    out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    for r in result['scales']:print('k',r['k'],'comparisons',len(r['comparisons']),'failures',r['failures'])
    print('sha256',hashlib.sha256(out.read_bytes()).hexdigest())
