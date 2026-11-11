#!/usr/bin/env python3
"""Exact chamber fibers of small labelled graphs, no external dependencies."""
from itertools import permutations
from collections import Counter, defaultdict
from math import isqrt, factorial
from pathlib import Path
import json

def need(b, msg):
    if not b: raise RuntimeError(msg)
def orient(p,E):
    pos={v:i for i,v in enumerate(p)}
    return sum((pos[u]<pos[v])<<i for i,(u,v) in enumerate(E))
def check(n,E):
    E=tuple(sorted(tuple(sorted(e)) for e in E)); es=set(E)
    fibers=defaultdict(list)
    for p in permutations(range(1,n+1)): fibers[orient(p,E)].append(p)
    single=[]
    for key,ps in fibers.items():
        pool=set(ps); seen={ps[0]}; stack=[ps[0]]
        while stack:
            p=stack.pop()
            for i in range(n-1):
                e=tuple(sorted(p[i:i+2]))
                q=p[:i]+(p[i+1],p[i])+p[i+2:]
                same=orient(q,E)==key
                need(same==(e not in es),'swap iff nonedge')
                if same and q not in seen: seen.add(q);stack.append(q)
        need(seen==pool,'fiber nonedge-swap connectivity')
        for p in ps:
            hp=all(tuple(sorted(p[i:i+2])) in es for i in range(n-1))
            need(hp==(len(ps)==1),'singleton iff Hamiltonian')
        if len(ps)==1:single.append(ps[0])
    hist=Counter(map(len,fibers.values()))
    need(sum(k*v for k,v in hist.items())==factorial(n),'partition mass')
    return {'n':n,'edges':E,'orientations':len(fibers),'fiber_size_histogram':dict(sorted(hist.items())), 'oriented_H':len(single),'singleton_witnesses':single}

def main():
    cases={
      'path4':(4,[(1,2),(2,3),(3,4)]),
      'star4':(4,[(1,2),(1,3),(1,4)]),
      'cycle4':(4,[(1,2),(2,3),(3,4),(1,4)]),
      'complete5':(5,[(i,j) for i in range(1,6) for j in range(i+1,6)]),
      'empty5':(5,[]),
    }
    for n in range(2,9):
        cases['Q'+str(n)]=(n,[(i,j) for i in range(1,n+1) for j in range(i+1,n+1) if isqrt(i+j)**2==i+j])
    out={name:check(*c) for name,c in cases.items()}
    need(out['path4']['oriented_H']==2 and out['star4']['oriented_H']==0,'tree hostile')
    receipt={'status':'PASS','universe':'All permutations of each listed graph; square-sum Q2..Q8 and named positive/hostile graphs','claims_checked':['Each fiber is connected by adjacent nonedge swaps','Singleton fibers equal oriented Hamiltonian paths','Exact fiber mass equals n!'],'graphs':out}
    Path(__file__).with_name('results.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:{x:v[x] for x in ['orientations','fiber_size_histogram','oriented_H']} for k,v in out.items()},indent=2))
if __name__=='__main__':main()
