#!/usr/bin/env python3
"""Independent recursive port assignment and direct relative-fibre enumeration."""
from itertools import combinations,product
from collections import Counter
from pathlib import Path
import importlib.util,json
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('colors',ROOT/'triples/colored_fibers.py')
f=importlib.util.module_from_spec(spec);spec.loader.exec_module(f)
def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def color_tree(t,bits):
    kind,mask,kids=t
    if kind=='leaf':bits[(mask&-mask).bit_length()-1]=0;return
    for child in kids:color_tree(child,bits)
    if kind=='linear':
        bit=0
        for child in kids:
            if child[0]=='leaf':bits[(child[1]&-child[1]).bit_length()-1]=bit;bit^=1
            else:bit=0
    elif all(child[0]=='leaf' for child in kids):
        bits[(kids[0][1]&-kids[0][1]).bit_length()-1]=1
    need(len({bits[i] for i in range(len(bits)) if mask>>i&1})==2,'every nontrivial node mixed')
result=[]
for n in range(1,6):
    pairs=list(combinations(range(n),2));triples=list(combinations(range(n),3));data=[]
    for code in range(1<<len(pairs)):
        a=f.rows(n,code)
        c3=tuple(all(sum((a[v]>>w)&1 for w in tri if w!=v)==1 for v in tri) for tri in triples)
        tree=f.tree(a,(1<<n)-1);bits=[0]*n;color_tree(tree,bits)
        data.append((code,c3,tuple(bits)))
    lookup={}
    for bits in product(range(2),repeat=n):
        mask=sum(1<<i for i,(u,v) in enumerate(pairs) if bits[u]!=bits[v])
        lookup[bits]=(mask,Counter((c3,code&mask) for code,c3,_ in data))
    for code,c3,bits in data:
        mask,counts=lookup[bits];need(counts[c3,code&mask]==1,'engineered port singles out orientation')
    result.append(dict(n=n,all_tournaments=len(data),singleton_fibres=len(data)))
print(json.dumps({'status':'PASS','counts':result,'scope':'All labelled tournaments n<=5, direct C3/cross-color grouping with independently implemented port coloring; no isomorphism or H preservation claim.'},indent=2))
