#!/usr/bin/env python3
"""Exact checks of PB1--PB6. Python standard library only; no numerical floats."""
from itertools import combinations, permutations, product
from functools import lru_cache
from math import factorial
from pathlib import Path
import json


def paths(a):
    n = len(a)
    return [p for p in permutations(range(n)) if all(a[u][v] for u, v in zip(p,p[1:]))]


def count(a):
    n = len(a)
    if not n: return 1
    out = [tuple(j for j in range(n) if a[i][j]) for i in range(n)]
    full = (1 << n)-1
    @lru_cache(None)
    def suffix(mask,v):
        if mask == full: return 1
        return sum(suffix(mask | 1<<w,w) for w in out[v] if not mask>>w&1)
    return sum(suffix(1<<v,v) for v in range(n))


def tau(a):
    return [[int(i != j and not a[j][i]) for j in range(len(a))] for i in range(len(a))]


def strong(a):
    n=len(a)
    for s in range(n):
        seen={s}; todo=[s]
        while todo:
            i=todo.pop()
            for j in range(n):
                if a[i][j] and j not in seen:seen.add(j);todo.append(j)
        if len(seen)!=n:return False
    return True


def cycle(n):
    a=[[0]*n for _ in range(n)]
    for i in range(n):a[i][(i+1)%n]=1
    return a


def family(h,m):
    # h=2: two cycles, with first of length3 and second of length m.
    if h==2:
        a=[[0]*(m+2) for _ in range(m+2)]
        for p in ([0,1,2],[0]+list(range(3,m+2))):
            for u,v in zip(p,p[1:]+p[:1]):a[u][v]=1
        return a
    # h-cycle 1..h, hub0; I={1}, O={2,...,h}; m-cycle at hub.
    n=h+m;a=[[0]*n for _ in range(n)]
    p=list(range(1,h+1))
    for u,v in zip(p,p[1:]+p[:1]):a[u][v]=1
    a[1][0]=1
    for i in range(2,h+1):a[0][i]=1
    p=[0]+list(range(h+1,n))
    for u,v in zip(p,p[1:]+p[:1]):a[u][v]=1
    return a


def main():
    results={'status':'FINITE-EXACT','parity':[],'oriented_bridges':[], 'strong_fixed_count':[]}
    for n in range(1,5):
        es=[(i,j) for i in range(n) for j in range(n) if i!=j]
        # Every possible linear forest, represented by its directed edge mask.
        edge_index={e:k for k,e in enumerate(es)}
        forests=set()
        for p in permutations(range(n)):
            pe=[edge_index[e] for e in zip(p,p[1:])]
            for subset in range(1<<len(pe)):
                forests.add(sum(1<<e for k,e in enumerate(pe) if subset>>k&1))
        for bits in product((0,1),repeat=len(es)):
            a=[[0]*n for _ in range(n)]
            for (i,j),b in zip(es,bits):a[i][j]=b
            ha=len(paths(a));at=tau(a)
            assert ha==count(a)
            assert (ha-count(at))%2==0
            amask=sum(b<<k for k,b in enumerate(bits))
            reciprocity=sum((-1)**f.bit_count()*factorial(n-f.bit_count())
                            for f in forests if amask&f==f)
            assert reciprocity==count(at)
            assert tau(at)==a
        results['parity'].append({'order':n,'digraphs':2**len(es)})
    for n in range(1,5):
        pairs=list(combinations(range(n),2))
        for states in product((0,1,2),repeat=len(pairs)):
            o=[[0]*n for _ in range(n)];missing=[]
            for (i,j),s in zip(pairs,states):
                if s==0:missing.append((i,j))
                elif s==1:o[i][j]=1
                else:o[j][i]=1
            l=tau(o); m=len(missing);coeff=[0]*(1<<m)
            for p in paths(l):
                used={frozenset(e) for e in zip(p,p[1:])}
                mask=sum(1<<k for k,e in enumerate(missing) if frozenset(e) in used)
                coeff[mask]+=1
            totals=[];allsets=[]
            for bits in product((0,1),repeat=m):
                t=[r[:] for r in o]
                for (i,j),b in zip(missing,bits):t[i][j]=b;t[j][i]=1-b
                ps=set(paths(t));totals.append(len(ps));allsets.append(ps)
            assert coeff[0]==count(o)
            assert sum(coeff)==count(l)
            assert sum(totals)==sum(c*2**(m-s.bit_count()) for s,c in enumerate(coeff))
            ns=[sum(c for u,c in enumerate(coeff) if u&s==s) for s in range(1<<m)]
            assert sum((-1)**s.bit_count()*c for s,c in enumerate(ns))==count(o)
            assert set.intersection(*allsets)==set(paths(o))
            assert set.union(*allsets)==set(paths(l))
            assert count(o)<=min(totals)<=max(totals)<=count(l)
        results['oriented_bridges'].append({'order':n,'graphs':3**len(pairs)})
    for h in (2,3,4,7,21):
        for m in (3,4,8):
            a=family(h,m)
            assert strong(a)
            assert all(not(a[i][j] and a[j][i]) for i,j in combinations(range(len(a)),2))
            assert count(a)==h
            if len(a)<=8:assert len(paths(a))==h
            results['strong_fixed_count'].append({'h':h,'cycle_length':m,'order':len(a),'actual':count(a)})
    # Controls: raw H can change parity under case2 when 00 is admitted.
    a=[[0,1,0],[1,0,1],[0,0,0]]
    b=[[a[i][j]*(1-a[j][i]) for j in range(3)] for i in range(3)]
    assert count(a)==1 and count(b)==0
    results['case2_parity_counterexample']={'raw':a,'raw_h':1,'decoded_h':0}
    output=Path(__file__).with_name('checks.json')
    output.write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps(results,indent=2))

if __name__=='__main__':main()
