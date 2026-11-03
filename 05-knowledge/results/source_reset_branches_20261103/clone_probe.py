"""Exact one-clone-ray test of the three adaptive H31 sources and continuations.

Uses the inherited Session155 clone-slot identity, not numerical curve fitting.
"""
from collections import defaultdict
from itertools import permutations
from math import comb
import json,sys
from pathlib import Path

def coefficients(T,r):
    c=defaultdict(int);old=[x for x in range(len(T)) if x!=r]
    for pi in permutations(old):
        t=((T[r]>>pi[0])&1)+((T[pi[-1]]>>r)&1);q=0
        for x,y in zip(pi,pi[1:]):
            good=((T[x]>>r)&1) and ((T[r]>>y)&1)
            t+=good
            if not (T[x]>>y)&1:
                if not good:break
                q+=1
        else:
            for j in range(q+1):c[t-j]+=(-1)**j*comb(q,j)
    return {j:v for j,v in sorted(c.items()) if v}

def main():
    if len(sys.argv)!=2:raise SystemExit('usage: clone_probe.py selected_sources.json')
    data=json.loads(Path(sys.argv[1]).read_text(encoding='utf-8'))
    records=[]
    for i,case in enumerate(data):
        for kind,graph in [('source',case['source'])]+[('child',c) for c in case['children']]:
            for v in range(len(case['source']['rows'])):
                coef=coefficients(graph['rows'],v)
                if sum(c*j for j,c in coef.items())!=graph['H']:raise RuntimeError('k1 control')
                records.append(dict(source=i,kind=kind,H=graph['H'],vertex=v,coefficients=coef))
    sources=[r for r in records if r['kind']=='source'];matches=[];branches=[]
    for i,a in enumerate(sources):
        for b in sources[:i]:
            if a['source']!=b['source'] and {j:c for j,c in a['coefficients'].items() if j}=={j:c for j,c in b['coefficients'].items() if j}:matches.append([a,b])
        for b in records:
            if b['kind']!='child' or a['source']!=b['source'] or a['vertex']!=b['vertex']:continue
            diff={j:b['coefficients'].get(j,0)-3*a['coefficients'].get(j,0) for j in set(a['coefficients'])|set(b['coefficients']) if j}
            if all(c==0 or j==1 for j,c in diff.items()):branches.append(dict(source=a,child=b,offset=diff.get(1,0)))
    result=dict(records=records,matches=matches,branches=branches)
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__=='__main__':main()
