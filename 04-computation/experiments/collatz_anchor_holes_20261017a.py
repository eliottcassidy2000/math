#!/usr/bin/env python3
"""Exact hole moments of source-authenticated independent anchor controllers.

Predeclared continuation: stage floor(k/2)+2+(anchor index mod 3), as in
the preceding repair, and thresholds beta=3,4,5,6. No ROOT oracle.
"""
import argparse
import hashlib
import json
from pathlib import Path
import collatz_interval_atlas_repair_20261016e as R

checks=0
def need(ok,label):
    global checks
    checks+=1
    if not ok:raise RuntimeError(label)

def square_sum(n):return n*(n+1)*(2*n+1)//6

def interval_moments(lo,hi,anchors):
    if not anchors:return None
    A=sorted(anchors)
    need(len(A)==len(set(A)) and lo<=A[0]<=A[-1]<hi,'literal anchor rows')
    left,right=A[0]-lo,hi-1-A[-1]
    gaps=[b-a for a,b in zip(A,A[1:])]
    R1=left*(left+1)//2+right*(right+1)//2+sum(l*l//4 for l in gaps)
    R2=square_sum(left)+square_sum(right)+sum((l**3+(2*l if l%2==0 else -l))//12 for l in gaps)
    return dict(count=len(A),left_tail=left,right_tail=right,max_interior_gap=max(gaps,default=0),
                gap_cube_sum=sum(l**3 for l in gaps),R1=R1,R2=R2)

def moments(X,anchors):return interval_moments(X,2*X,anchors)

def graded(r,parts):
    X=r['X'];need(X%parts==0,'equal physical slabs')
    chosen=[]
    for lo in range(X,2*X,X//parts):
        hi=lo+X//parts;bank=[(a,B) for a,B,s in r['anchors'] if lo<=a<hi]
        need(bool(bank),'each slab has an authenticated anchor')
        best=None
        for T in sorted({B for a,B in bank}):
            S=[a for a,B in bank if B<=T];m=interval_moments(lo,hi,S)
            score=(hi-lo)*T*T+8*T*m['R1']+16*m['R2']
            if best is None or score<best['score']:
                best=dict(lo=lo,hi=hi,T=T,score=score,anchors=S,**m)
        chosen.append(best)
    total=sum(x['score'] for x in chosen)
    return dict(parts=parts,score=total,passes_7X=total<=49*X*X,slabs=chosen)

def profile(r):
    X=r['X'];out=[]
    for beta in (3,4,5,6):
        good=[a for a,B,stage in r['anchors'] if B*B<=beta*beta*X]
        m=moments(X,good)
        if m:
            # beta X + 4 sqrt(R2) <= 7 X, tested without radicals.
            m.update(beta=beta,passes_7X=16*m['R2']<=(7-beta)**2*X*X,
                     good_anchors=good)
        out.append(m)
    adaptive=[]
    for b in (2**(r['k']//2-i) for i in range(4)):
        E=r['stages'][-1]['energy'];R2=X*(b-1)*(2*b-1)//6
        slack=49*X*X-E-16*R2
        adaptive.append(dict(bin_width=b,anchor_count=X//b,common_stage=r['stages'][-1]['stage'],
                             fine_energy=E,distance_R2_upper=R2,
                             passes_7X=slack>=0 and 64*E*R2<=slack*slack,
                             status='PROVED selection from authenticated finite row energy; selected rows not stored'))
    return dict(k=r['k'],X=X,base_stage=r['base_stage'],profiles=out,
                graded=[graded(r,p) for p in (1,2,4,8)],
                adaptive_selection=adaptive,
                native_cone_energy=r.get('cone_energy',r.get('native_cone_energy')),anchor_energy=r['anchor_energy'],
                anchors=r['anchors'],stages=r['stages'],native_extension_audits=r['native_extension_audits'])

def controls():
    from itertools import combinations
    for X in range(1,13):
        for size in range(1,min(X,3)+1):
            for A in combinations(range(X,2*X),size):
                m=moments(X,A);dist=[min(abs(n-a) for a in A) for n in range(X,2*X)]
                need(sum(dist)==m['R1'] and sum(d*d for d in dist)==m['R2'],'independent direct distance sums')
    # Same anchor count and same interior gaps, different unobserved boundary.
    X=1024;b=16
    centred=list(range(X+256,X+768,b));left=list(range(X,X+512,b))
    a,c=moments(X,centred),moments(X,left)
    need(a['count']==c['count'] and a['max_interior_gap']==c['max_interior_gap'],'matched misleading summaries')
    need(a['R2']!=c['R2'],'boundary tails change actual coverage')
    plateau=[]
    # Abstract Lipschitz hostile, not a Collatz trajectory: a plateau has
    # energy O(X^2), while every single global threshold costs Omega(X^9/4).
    for m in (1,2):
        Xp=2**(8*m);height=2**(5*m);length=2**(6*m);ramp=height//4
        vals=[4*t for t in range(ramp)]+[height]*length+[4*t for t in range(ramp-1,-1,-1)]
        energy=sum(v*v for v in vals)
        need(energy<=2*Xp*Xp,'plateau native scalar energy')
        # Threshold T=height pays X*height^2; T<height leaves a gap of
        # at least length+1 between selected rows and pays its cubic moment.
        gap=length+1;gap_R2=(gap**3+(2*gap if gap%2==0 else -gap))//12
        floor=min(Xp*height*height,16*gap_R2)
        need(floor>=2**(18*m),'all global thresholds have superquadratic score')
        plateau.append(dict(m=m,X=Xp,height=height,length=length,energy=energy,
                            all_threshold_score_lower_bound=floor))
    from itertools import product
    for values in product(range(3),repeat=8):
        selected=[lo+min(range(4),key=lambda t:values[lo+t]) for lo in (0,4)]
        Q2=sum(4*values[a]**2 for a in selected)
        distance2=sum((n-selected[n//4])**2 for n in range(8))
        need(Q2<=sum(v*v for v in values) and distance2<=8*3*7//6,
             'adaptive minimum selection versus independent full fine energy')
    return dict(X=X,centred=a,left_cluster=c,abstract_plateau=plateau)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--k',type=int,nargs='*',default=[])
    ap.add_argument('--receipt',type=Path);args=ap.parse_args()
    hostile=controls();runs=[];inherited=None
    if args.receipt:
        raw=args.receipt.read_bytes();old=json.loads(raw)
        inherited=dict(sha256=hashlib.sha256(raw).hexdigest(),native_checks=old.get('native_checks',old.get('checks')))
        runs.extend(profile(r) for r in old['runs'])
    for k in args.k:runs.append(profile(R.run(k)))
    print(json.dumps(dict(status='FINITE-EXACT; no uniform deficit conclusion',runs=runs,
          boundary_hostile=hostile,inherited_receipt=inherited,checks=checks,native_checks=R.CHECKS,atlas_checks=R.A.CHECKS),
          sort_keys=True,separators=(',',':')))
