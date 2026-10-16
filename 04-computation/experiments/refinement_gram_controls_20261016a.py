#!/usr/bin/env python3
"""Exact controls for sharp Banach refinement and native Gram remainder.

Universe: scalar bounded-increment paths of lengths 1..6; explicit sharp
l1-valued paths of lengths 1..24; native first-reset graphs k=8,10. No
spectral floating-point test or heuristic terminal labels are used.
"""
from fractions import Fraction as F
from itertools import product
import hashlib
import json
from pathlib import Path
import numpy as np
from collatz_compact_resets_20261015f import Graph, need
from paired_bin_refinement_20261015g import point_rows

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'05-knowledge/results'


def absnorm(a):
    return sum(abs(x) for x in a)


def budget(rows,L):
    b=len(rows);s=absnorm([sum(r[e] for r in rows) for e in range(len(rows[0]))])
    fine=sum(absnorm(r)**2 for r in rows)
    m=F(b*b-1,3*b);t=F(7*b**4-15*b*b+8,60*b*b)
    bound=F(s*s,b)+2*L*m*s+L*L*b*t
    return fine,bound


def algebra_controls():
    count=0
    for b in range(1,7):
        for steps in product(range(-2,3),repeat=b-1):
            for initial in [-2,0,3]:
                rows=[[initial]]
                for x in steps: rows.append([rows[-1][0]+x])
                fine,bound=budget(rows,2)
                need(fine<=bound,'scalar hostile census'); count+=1
    sharp=[]
    for b in range(1,25):
        for baseline in [0,1,13]:
            rows=[[F(baseline)]+[4*(F(i>j)-F(b-j,b)) for j in range(1,b)]
                  for i in range(1,b+1)]
            need(all(absnorm([y-x for x,y in zip(rows[i],rows[i+1])])==4
                     for i in range(b-1)), 'sharp path increment')
            fine,bound=budget(rows,4)
            need(fine==bound,'sharp l1 bound including nonzero mean')
            sharp.append(dict(b=b,baseline=baseline,energy=[fine.numerator,fine.denominator]))
    # Actual native odd columns 3 and 5 coalesce at 5 with clocks 1 and 0.
    old=np.array([[1,0],[0,1]],dtype=np.int64)
    new=(old@np.array([-1,1],dtype=np.int64)).reshape(2,1)
    def invoice(v):
        b=len(v);r=v.sum(axis=0)
        envelope=b*sum(int(x)**2 for x in np.abs(v).sum(axis=1))-int(np.abs(r).sum())**2
        J=b*(v.T@v)-np.outer(r,r)
        return envelope,int(np.abs(J).sum())
    need(invoice(old)==(0,4) and invoice(new)==(4,4),
         'envelope refinement can increase while Gram remainder does not')
    need((3*3+1)//2==5,'native hostile clock')
    return dict(scalar_paths=count,sharp_paths=sharp,
                hostile=dict(columns=[3,5],literal_clocks=[1,0],
                             old_envelope_gap=0,new_envelope_gap=2,
                             old_gram_remainder_norm=2,new_gram_remainder_norm=2))


def native(k):
    X=1<<k;b=1<<(k//2);g=Graph(8*X)
    frozen=json.loads((DATA/'collatz_survivor_packets_20261015e.json').read_text())
    stop=len(next(r for r in frozen['first_reset'] if r['k']==k)['stages'])
    cols,idx,inv=np.unique(g.endpoint[1:],return_index=True,return_inverse=True)
    records=[];previous=None
    for stage in range(1,stop+1):
        reps=idx+1;oldclock=g.clock[reps].copy()
        g.extend(cols);g.close(cols,inv)
        transport_end=g.endpoint[reps].copy();extra=g.clock[reps]-oldclock
        cols,idx,inv=np.unique(g.endpoint[1:],return_index=True,return_inverse=True)
        points=point_rows(k,g,inv,len(cols));rows=points.reshape(-1,b,len(cols)).sum(axis=1)
        fine=points.T@points;coarse=rows.T@rows;J=b*fine-coarse
        need(np.array_equal(J,J.T),'symmetric exact Gram remainder')
        # Independent Haar outer-product decomposition; all denominators divide b.
        haar=np.zeros_like(J);charge=0
        for length in [1<<p for p in range(1,k//2+1)]:
            small=points.reshape(-1,length,len(cols))
            d=small[:,:length//2].sum(axis=1)-small[:,length//2:].sum(axis=1)
            haar+=(b//length)*(d.T@d)
            charge+=(b//length)*sum(int(v)**2 for v in np.abs(d).sum(axis=1))
        need(np.array_equal(J,haar),'exact positive Haar decomposition')
        jnorm=int(np.abs(J).sum())
        need(jnorm<=charge and 3*charge<=4*b*X*(b*b-1),'uniform Gram refinement price')
        if previous is not None:
            # Congruence through the *same* literal source paths, independently
            # accumulated in ordered target-pair fibres.
            target=np.searchsorted(cols,transport_end)
            sign=1-2*(extra.astype(np.int64)&1)
            pushed=np.zeros_like(J)
            for i in range(len(target)):
                np.add.at(pushed[target[i]],target,previous[i]*sign[i]*sign)
            need(np.array_equal(pushed,J),'native/refinement commuting square')
            need(jnorm<=previous_norm,'entrywise Gram remainder cannot increase')
        previous=J;previous_norm=jnorm
        cf=int(np.abs(coarse).sum());ff=int(np.abs(fine).sum())
        need(abs(b*ff-cf)<=jnorm,'refined Gram budget transfer')
        records.append(dict(stage=stage,columns=len(cols),denominator=b,
                            gram_remainder_numerator=jnorm,haar_budget_numerator=charge,
                            coarse_gram_numerator=cf,fine_gram_norm=ff))
        g.internal_clock_audit();g.replay(range(1,8*X+1,max(1,X//32)))
    return dict(k=k,bin_width=b,stages=records)


if __name__=='__main__':
    result=dict(status='FINITE-EXACT controls, not a uniform reset theorem',
                algebra=algebra_controls(),native=[native(k) for k in [8,10]])
    out=DATA/'refinement_gram_controls_20261016a.json'
    out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print('scalar',result['algebra']['scalar_paths'],'sharp',len(result['algebra']['sharp_paths']))
    print('native',result['native'])
    print('sha256',hashlib.sha256(out.read_bytes()).hexdigest())
