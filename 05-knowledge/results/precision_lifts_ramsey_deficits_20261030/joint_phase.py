#!/usr/bin/env python3
"""Exact common-clock two-speed attachment via integer Euclidean floor sums."""
from fractions import Fraction as F
from math import ceil, floor, gcd
from pathlib import Path
from itertools import combinations
import importlib.util, json, random

PREV=Path(__file__).parent.parent/'terminal_fees_phase_covers_20261030/phase_cover.py'
spec=importlib.util.spec_from_file_location('phase_cover',PREV)
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
need=old.need

def floor_sum(n,m,a,b):
    """Sum floor((a*i+b)/m), 0<=i<n; exact for signed a,b."""
    need(n>=0 and m>0,'floor-sum domain')
    ans=0
    while True:
        qa,a=divmod(a,m);qb,b=divmod(b,m)
        ans+=qa*n*(n-1)//2+qb*n
        y=a*n+b
        if y<m:return ans
        n,b=divmod(y,m);m,a=a,m

def one_tooth(left,right,r,v):
    lo=ceil(v*left-1+r);hi=floor(v*right-r)
    count=max(0,hi-lo+1)
    return count,(max(left,(lo+r)/v) if count else None)

def pair_interval(left,right,r,u,v,witness=True):
    """Count all intersecting native tooth pairs on a closed interval."""
    r=F(r);left=F(left);right=F(right)
    need(0<r<=F(1,2) and left<=right and u>0 and v>0,'pair domain')
    if u>v:u,v=v,u
    low=ceil(u*left-1+r);high=floor(u*right-r)
    if low>high:return 0,None
    first=ceil(u*left-r);last=floor(u*right-1+r)
    total=0;t=None
    # Only the first and last intersecting u-teeth can be clipped.
    for k in sorted({low,high}):
        if first<=k<=last:continue
        a=max(left,(k+r)/u);b=min(right,(k+1-r)/u)
        if a>b:continue
        c,s=one_tooth(a,b,r,v);total+=c
        if t is None and c:t=s
    if first<=last:
        s,d=r.numerator,r.denominator
        M=d*u;A=d*v;B=v*(d-s)-s*u;C=v*s-u*(d-s)-1
        def count_prefix(n):
            return floor_sum(n,M,A,A*first+B)-floor_sum(n,M,A,A*first+C)
        n=last-first+1;c=count_prefix(n);need(c>=0,'negative pair count')
        total+=c
        if witness and t is None and c:
            lo,hi=1,n
            while lo<hi:
                mid=(lo+hi)//2
                if count_prefix(mid):hi=mid
                else:lo=mid+1
            k=first+lo-1
            _,t=one_tooth((k+r)/u,(k+1-r)/u,r,v)
    return total,t if witness else None

def attach_pair(intervals,r,u,v):
    for a,b in intervals:
        count,t=pair_interval(a,b,r,u,v)
        if count:return t
    return None

def pair_bank(core,r):
    """Finite proof obligation for all unordered positive attached-speed pairs."""
    need(0<r<=F(1,4),'two-tail conductor requires r<=1/4')
    intervals=old.safe_intervals(core,r)
    W=max((b-a for a,b in intervals),default=F(0));need(W>0,'positive core width')
    # For r<=1/6, a full u-period has <=2 safe pieces, one of length
    # >=(1-2r)/(2u)>=2r/u. For 1/6<r<=1/4, fit a full safe tooth.
    B=ceil((1 if r<=F(1,6) else 2-2*r)/W);rows=[];queries=0
    for u in range(1,B):
        joint=old.safe_intervals((*core,u),r)
        W2=max((b-a for a,b in joint),default=F(0))
        need(W2>0,('point-only or empty branch needs a separate residue certificate',core,u))
        cutoff=ceil(2*r/W2);bad=[]
        for v in range(u,cutoff):
            t=attach_pair(intervals,r,u,v)
            ref=old.attachment(joint,r,v)
            need((ref is not None)==(t is not None),('independent one-tail bank reference',core,u,v))
            if t is None:bad.append(v)
            else:need(all(old.dist(x*t)>=r for x in (*core,u,v)),'bank witness')
            queries+=1
        rows.append(dict(smaller_speed=u,larger_speed_tail=cutoff,rejected_larger_speeds=bad))
    return dict(core=core,radius=str(r),core_components=len(intervals),both_large_cutoff=B,
                finite_queries=queries,branches=rows,
                all_pairs_accepted=not any(x['rejected_larger_speeds'] for x in rows))

def run():
    floor_tests=0
    for n in range(12):
        for m in range(1,14):
            for a,b in ((-19,-23),(-3,8),(0,0),(7,19),(25,-6)):
                need(floor_sum(n,m,a,b)==sum((a*i+b)//m for i in range(n)),'floor sum')
                floor_tests+=1
    rng=random.Random(34357);queries=0
    for _ in range(2500):
        a,b=sorted((F(rng.randrange(41),40),F(rng.randrange(41),40)))
        r=rng.choice((F(1,14),F(1,4),F(1,3),F(1,2)))
        u,v=rng.randrange(1,60),rng.randrange(1,60)
        literal=[]
        for k in range(ceil(u*a-1+r),floor(u*b-r)+1):
            for l in range(ceil(v*a-1+r),floor(v*b-r)+1):
                lo=max(a,(k+r)/u,(l+r)/v);hi=min(b,(k+1-r)/u,(l+1-r)/v)
                if lo<=hi:literal.append((lo,hi))
        count,t=pair_interval(a,b,r,u,v)
        need(count==len(literal),('literal tooth count',a,b,r,u,v,count,len(literal)))
        need((t is not None)==bool(literal),'witness existence')
        if t is not None:need(a<=t<=b and min(old.dist(u*t),old.dist(v*t))>=r,'witness')
        queries+=1
    for core in combinations(range(1,7),2):
        for r in (F(1,4),F(1,3),F(1,2)):
            ints=old.safe_intervals(core,r)
            for u in range(1,17):
                for v in range(u,17):
                    t=attach_pair(ints,r,u,v);ref=old.safe_intervals((*core,u,v),r)
                    need(bool(ref)==(t is not None),('core comparison',core,r,u,v))
                    if t is not None:need(min(old.dist(s*t) for s in (*core,u,v))>=r,'core witness')
                    queries+=1
    hostile=old.safe_intervals((1,2),F(1,4))
    need(old.attachment(hostile,F(1,4),3) is not None and old.attachment(hostile,F(1,4),12) is not None,'individual success')
    need(attach_pair(hostile,F(1,4),3,12) is None,'joint hostile rejected')
    # Thirteen moving speeds: eleven-speed cores plus two arbitrary tails.
    cores=[tuple(range(1,12)),tuple(26*j for j in range(1,12)),
           (1,2,3,5,8,13,21,34,55,89,144),(2,3,5,7,11,13,17,19,23,29,31)]
    banks=[]
    for c in cores:
        banks.append(pair_bank(c,F(1,14)))
        print('pair bank',c,'queries',banks[-1]['finite_queries'],flush=True)
    need(all(b['all_pairs_accepted'] for b in banks),'an infinite pair bank fails')
    large=[]
    ints=old.safe_intervals(cores[0],F(1,14))
    for digits in (100,1000):
        u=10**digits+39;v=10**digits+123
        t=attach_pair(ints,F(1,14),u,v)
        need(t is not None and min(old.dist(s*t) for s in (*cores[0],u,v))>=F(1,14),'large pair')
        large.append(dict(input_digits=len(str(u)),witness=str(t),actual_minimum=str(min(old.dist(s*t) for s in (*cores[0],u,v)))))
    out=dict(status='FINITE-EXACT checks of proved two-tail compiler and infinite fixed-core banks',
             floor_sum_tests=floor_tests,literal_queries=queries,hostile='(1,2)+3+12 rejected at 1/4',
             all_positive_pair_banks=banks,large_inputs=large)
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(floor_sum_tests=floor_tests,literal_queries=queries,
                         infinite_pair_banks=len(banks),bank_queries=sum(b['finite_queries'] for b in banks),
                         large_input_digits=[x['input_digits'] for x in large]),indent=2))

if __name__=='__main__':run()
