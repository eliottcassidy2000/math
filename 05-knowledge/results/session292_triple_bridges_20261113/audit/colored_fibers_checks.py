#!/usr/bin/env python3
"""Reviewer extension: n5 three colors, and exact new-port refinement gates."""
from itertools import combinations,product
from collections import Counter
from pathlib import Path
import importlib.util,json
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('colors',ROOT/'triples/colored_fibers.py')
f=importlib.util.module_from_spec(spec);spec.loader.exec_module(f)
def check(ok,msg):
    if not ok:raise RuntimeError(msg)
def cycle_signature(a,n):
    ans=[]
    for triple in combinations(range(n),3):
        ans.append(all(sum((a[v]>>w)&1 for w in triple if w!=v)==1 for v in triple))
    return tuple(ans)
def criterion(t,cols,bits):
    kind,mask,kids=t
    if not all(criterion(k,cols,bits) for k in kids):return False
    const=lambda m:len({bits[i] for i in range(len(bits)) if m>>i&1})<=1
    if kind=='prime' and f.color(mask,cols) is not None and not const(mask):return False
    if kind=='linear':
        old=None;run=0
        for k in kids:
            cur=f.color(k[1],cols)
            if cur is not None and cur==old:run|=k[1]
            else:
                if not const(run):return False
                old=cur;run=k[1] if cur is not None else 0
        if not const(run):return False
    return True
n=5; pairs=list(combinations(range(n),2));data=[]
for code in range(1<<len(pairs)):
    a=f.rows(n,code)
    data.append((code,cycle_signature(a,n),f.tree(a,(1<<n)-1)))
tested=0
for cols in product(range(3),repeat=n):
    mask=sum(1<<j for j,(u,v) in enumerate(pairs) if cols[u]!=cols[v])
    counts=Counter((c3,code&mask) for code,c3,t in data)
    for code,c3,t in data:
        check(f.predicted(t,cols)==counts[c3,code&mask],str((code,cols)))
        tested+=1
refinements=0
for n in range(1,5):
    for code in range(1<<(n*(n-1)//2)):
        t=f.tree(f.rows(n,code),(1<<n)-1)
        for cols in product(range(2),repeat=n):
            old=f.predicted(t,cols)
            for bits in product(range(2),repeat=n):
                new=f.predicted(t,list(zip(cols,bits)))
                check(new<=old,'refinement monotonicity')
                check((new==old)==criterion(t,cols,bits),'exact refinement criterion')
                refinements+=1
print(json.dumps({'status':'PASS','all_three_color_cases_at_n5':tested,'refinement_cases_n1_to_n4':refinements,
 'scope':'Direct C3/cross-color grouping; producer strong-module tree routine reused, independent cycle test and refinement criterion.'},indent=2))
