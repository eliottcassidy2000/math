#!/usr/bin/env python3
"""Source-atom interval certificates for fixed-controller annular cancellation.

Every atom keeps its original integer source, copy number, endpoint and full
shortcut clock. Matching changes with the physical row N. No ROOT oracle.
"""
import argparse
import copy
from collections import defaultdict
import hashlib
import heapq
import json
from pathlib import Path
import numpy as np
from collatz_compact_resets_20261015f import Graph, T

CHECKS=0
def check(ok, label):
    global CHECKS
    CHECKS+=1
    if not ok:raise RuntimeError(label)
def w(n):return 2 if n%3==2 else -1
def sign(n,clocks):return (1 if w(n)>0 else -1)*(-1)**int(clocks[n])
def tokens(n):return [(n,i) for i in range(abs(w(n)))]
def bounds(n,X):return max(X,(n+3)//4),min(2*X,n)

def matching(X,ends,clocks):
    free=defaultdict(lambda:[set(),set()]); queues=defaultdict(lambda:[[],[]])
    partner={}; start={}; records=[]
    active=set();changes=set()
    def key(a):return int(ends[a[0]]),int(sign(a[0],clocks)<0)
    def put(a):
        e,s=key(a);free[e][s].add(a);heapq.heappush(queues[e][s],a);changes.add(e)
    def add_source(n):
        for a in tokens(n):active.add(a);put(a)
    def remove_source(n,N):
        for a in tokens(n):
            check(a in active,'removed atom is active')
            active.remove(a)
            if a in partner:
                b=partner.pop(a);partner.pop(b)
                pair=tuple(sorted((a,b)));lo=start.pop(pair)
                if lo<N:records.append([list(pair[0]),list(pair[1]),lo,N])
                put(b)
            else:
                e,s=key(a);free[e][s].remove(a);changes.add(e)
    def close(N):
        for e in sorted(changes):
            pos,neg=free[e]
            while pos and neg:
                pq,nq=queues[e]
                while pq[0] not in pos:heapq.heappop(pq)
                while nq[0] not in neg:heapq.heappop(nq)
                a=heapq.heappop(pq);b=heapq.heappop(nq);pos.remove(a);neg.remove(b)
                partner[a]=b;partner[b]=a;start[tuple(sorted((a,b)))]=N
        changes.clear()
    for n in range(X+1,4*X+1):add_source(n)
    close(X);norms=[4*X-len(partner)]
    for N in range(X,2*X-1):
        remove_source(N+1,N+1)
        for n in range(4*N+1,4*N+5):add_source(n)
        close(N+1)
        norms.append(4*(N+1)-len(partner))
    for pair,lo in sorted(start.items()):
        records.append([list(pair[0]),list(pair[1]),lo,2*X])
    records.sort()
    budget=4*X-2+sum(abs(w(n)) for n in range(X+1,2*X))
    check(len(records)<=budget and 3*budget<16*X,'sharp universal interval-piece budget')
    return records,norms

def audit(X,ends,clocks,records,claimed_norms,direct=True):
    occupancy=defaultdict(list);coverage=[0]*(X+1)
    for aa,bb,lo,hi in records:
        a,b=tuple(aa),tuple(bb);n,m=a[0],b[0]
        check(a in tokens(n) and b in tokens(m) and n!=m,'literal distinct source atom types')
        check(int(ends[n])==int(ends[m]) and sign(n,clocks)==-sign(m,clocks),
              'same integer terminal and opposite weighted clock sign')
        l1,h1=bounds(n,X);l2,h2=bounds(m,X)
        check(max(l1,l2)<=lo<hi<=min(h1,h2),'pair interval stays in both source annuli')
        occupancy[a].append((lo,hi));occupancy[b].append((lo,hi))
        coverage[lo-X]+=1;coverage[hi-X]-=1
    for intervals in occupancy.values():
        intervals.sort()
        check(all(a[1]<=b[0] for a,b in zip(intervals,intervals[1:])),
              'each original unit of source mass used at most once in each physical row')
    cur=0;norms=[]
    for k in range(X):
        cur+=coverage[k];norms.append(4*(X+k)-2*cur)
    check(norms==claimed_norms,'coverage independently recovers every claimed norm')
    if direct:
        # Independent terminal event sweep; no matching data used.
        row=defaultdict(int);current=0
        def push(n,c):
            nonlocal current
            e=int(ends[n]);old=row[e];row[e]+=c*w(n)*(-1)**int(clocks[n])
            current+=abs(row[e])-abs(old)
        for n in range(X+1,4*X+1):push(n,1)
        for k,N in enumerate(range(X,2*X)):
            check(current==norms[k],'optimality against literal terminal row')
            if X<=256:check(current==sum(abs(v) for v in row.values()),'small independent full-vector norm')
            if N+1<2*X:
                push(N+1,-1)
                for n in range(4*N+1,4*N+5):push(n,1)
    return sum(n*n for n in norms)

def fixed_matching_hostile():
    paths={}
    for n in [3,4,9]:
        x=n;path=[x]
        while x!=1:
            check(len(path)<30,'declared finite hostile completed')
            x=T(x);path.append(x)
        paths[n]=path
    check([len(paths[n])-1 for n in [3,4,9]]==[5,2,13],'literal hostile clocks')
    check([w(n)*(-1)**(len(paths[n])-1) for n in [3,4,9]]==[1,-1,1],
          'original hostile weights')
    # N=2 needs all capacity 3--4; N=3 needs all capacity 4--9.
    check([n for n in paths if 2<n<=8]==[3,4], 'first restricted annular row')
    check([n for n in paths if 3<n<=12]==[4,9], 'second restricted annular row')
    return dict(sources=[3,4,9],paths=paths,rows=[2,3],
                message='restricted three-source subledger; a fixed source-capacity matching cannot attain both zero rows')

def rejects(fn,label):
    try:fn()
    except (RuntimeError,AssertionError):
        check(True,label);return
    raise RuntimeError('hostile was accepted: '+label)

def mutation_controls(X,g,records,norms):
    bad=copy.deepcopy(records);bad[0][0][1]=9
    rejects(lambda:audit(X,g.endpoint,g.clock,bad,norms,False),'invalid source-copy rejection')
    bad=copy.deepcopy(records);bad.append(bad[0])
    rejects(lambda:audit(X,g.endpoint,g.clock,bad,norms,False),'duplicated mass rejection')
    bad=copy.deepcopy(records);bad[0][2]=X-1
    rejects(lambda:audit(X,g.endpoint,g.clock,bad,norms,False),'wrong physical interval rejection')
    source=records[0][0][0];old=int(g.clock[source]);g.clock[source]=old+2
    # Parity-only tests cannot detect this. Literal native replay must do so.
    rejects(lambda:g.replay([source]),'same-parity false integer clock rejection')
    g.clock[source]=old

def sharp_format_controls():
    out=[]
    for X in range(1,65):
        ends=[1]*(8*X+1);clocks=[0]*(8*X+1)
        records,norms=matching(X,ends,clocks)
        check(audit(X,ends,clocks,records,norms)==0,'unrestricted periodic-weight null has zero row norms')
        budget=4*X-2+sum(abs(w(n)) for n in range(X+1,2*X))
        check(len(records)==budget,'sharpness of the record bound for unrestricted terminal labels')
        out.append([X,budget])
    return dict(scope='arbitrary signed terminal data, NOT an authenticated native controller',cases=out)

def run(k,stages):
    X=1<<k;g=Graph(8*X);out=[];last=None
    for stage in range(stages+1):
        if stage:
            cols,inv=np.unique(g.endpoint[1:],return_inverse=True)
            g.extend(cols);g.close(cols,inv)
        records,norms=matching(X,g.endpoint,g.clock)
        energy=audit(X,g.endpoint,g.clock,records,norms)
        if k==6 and stage==1:mutation_controls(X,g,records,norms)
        if last is not None:check(all(a>=b for a,b in zip(last,norms)),'every row contracts under native extension')
        last=norms
        g.replay(range(1,8*X+1))
        encoded=json.dumps(records,separators=(',',':')).encode()
        rec=dict(stage=stage,vertices=int(len(set(g.endpoint[1:]))),pieces=len(records),
                 energy=energy,energy_over_X2=energy/X**2,
                 norm_sha256=hashlib.sha256(json.dumps(norms).encode()).hexdigest(),
                 interval_sha256=hashlib.sha256(encoded).hexdigest(),
                 max_clock=int(g.clock.max()))
        if k<=8:
            rec.update(intervals=records,norms=norms,
                       endpoints=[int(x) for x in g.endpoint],clocks=[int(x) for x in g.clock])
        out.append(rec)
    pairs=[]
    for j in range(len(out)-2):
        pairs.append(dict(start=j,end=j+2,passes=out[j+2]['energy']<=49*X*X or
                          2*out[j+2]['energy']<=out[j]['energy']))
    return dict(k=k,X=X,piece_budget=4*X-2+sum(abs(w(n)) for n in range(X+1,2*X)),
                stages=out,paired_tests=pairs)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--k',nargs='+',type=int,default=[6,8,10])
    ap.add_argument('--stages',type=int,default=8);a=ap.parse_args()
    r=dict(status='FINITE-EXACT; no uniform arithmetic bound is inferred',
           hostile=fixed_matching_hostile(),sharp_format=sharp_format_controls(),
           runs=[run(k,a.stages) for k in a.k],checks=CHECKS)
    print(json.dumps(r,indent=2,sort_keys=True))
