#!/usr/bin/env python3
"""Independent order-ten stable-cell classification audit. No producer imports."""
import argparse, itertools as it, json
from collections import Counter
from pathlib import Path

def need(ok,msg):
    if not ok: raise RuntimeError(msg)
def cyc(n): return tuple(sum(1<<((i+k)%n) for k in range(1,(n+1)//2)) for i in range(n))
def relabel(rows,p):
    out=[0]*len(rows)
    for i,row in enumerate(rows): out[p[i]]=sum(1<<p[j] for j in range(len(rows)) if row>>j&1)
    return tuple(out)
def stable(rows):
    c=(0,)*len(rows)
    while True:
        sig=tuple((c[i],)+tuple(sum(rows[i]>>j&1 for j in range(len(rows)) if c[j]==k) for k in sorted(set(c))) for i in range(len(rows)))
        names={s:k for k,s in enumerate(sorted(set(sig)))}; d=tuple(names[s] for s in sig)
        if len(set(c))==len(set(d)): return tuple(tuple(i for i in range(len(rows)) if d[i]==k) for k in sorted(set(d)))
        c=d
def strong(rows):
    n=len(rows)
    for start in range(n):
        seen=todo=1<<start
        while todo:
            bit=todo&-todo; todo^=bit; more=rows[bit.bit_length()-1]&~seen; todo|=more; seen|=more
        if seen!=(1<<n)-1: return False
    return True
def autos(rows):
    n=len(rows); cells=stable(rows); targets={v:c for c in cells for v in c}; order=sorted(range(n),key=lambda v:(len(targets[v]),v)); p=[-1]*n; out=[]
    def visit(k,used):
        if k==n: out.append(tuple(p)); return
        v=order[k]
        for w in targets[v]:
            if used>>w&1: continue
            if any((rows[v]>>u&1)!=(rows[w]>>p[u]&1) for u in order[:k]): continue
            p[v]=w; visit(k+1,used|1<<w); p[v]=-1
    visit(0,0); return tuple(out)
def error(rows,group):
    sizes={v:len(c) for c in stable(rows) for v in c}; seen=set(); total=0; orbits=[]
    for v in range(len(rows)):
        if v in seen: continue
        orbit={p[v] for p in group}; seen|=orbit; orbits.append(sorted(orbit)); total+=sizes[v]-len(orbit)
    return total,orbits
def graph(n,bits):
    rows=[0]*n
    for k,(i,j) in enumerate(it.combinations(range(n),2)):
        a,b=(i,j) if bits>>k&1 else (j,i); rows[a]|=1<<b
    return tuple(rows)
def partitions(n,lo=1):
    if not n: yield ()
    for a in range(lo,n+1,2):
        for tail in partitions(n-a,a): yield (a,)+tail

