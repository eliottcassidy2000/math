#!/usr/bin/env python3
"""Independent permutation/cut audit. No lane counting engine is imported."""
from itertools import permutations
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent

def need(q,msg):
    if not q:raise RuntimeError(msg)

def inflate(rows,O,lengths):
    labels=[(i,k) for i in range(len(rows)) for k in range(lengths.get(i,1))]
    new=[]
    for i,(u,k) in enumerate(labels):
        new.append(sum(1<<j for j,(v,l) in enumerate(labels) if i!=j and ((u==v and k<l) or (u!=v and rows[u]>>v&1))))
    selected={i for i,(v,k) in enumerate(labels) if v in O}
    need(not any(lengths.get(v,1)!=1 for v in O),'selected clone undefined')
    return new,selected

def count(rows,O):
    n=len(rows);h=w=u=0
    for perm in permutations(range(n)):
        bad=[]
        for j in range(n-1):
            if not rows[perm[j]]>>perm[j+1]&1:
                bad.append(j)
                if len(bad)==2:break
        if len(bad)>1:continue
        if not bad:
            h+=1
            u+=perm[0] in O
            cuts=range(n-1)
        else:cuts=bad
        for j in cuts:
            if perm[j+1] in O:
                u+=1
                w+=perm[j] in O
    return {'H':h,'W':w,'U':u,'delta':h-w}

def main():
    out={'method':'all labelled permutations, with at most one failed adjacency; cut that adjacency or any legal one to count ordered covers','bridge':[],'coupled':[]}
    for rows,O,v in [((90,28,9,16,36,79,30),{3,6},2),((54,28,24,17,32,78,31),{2,6},0)]:
        for m in (1,2,3):
            child,sel=inflate(rows,O,{v:m});r=count(child,sel)
            need(r['delta']==1-4*(3**m-3*2**m+3),'bridge literal mismatch')
            out['bridge'].append({'core':rows,'O':sorted(O),'port':v,'m':m,**r})
    rows=(6,4,8,19,39,15);O={0,4}
    for a,b in ((1,1),(2,1),(1,2),(2,2),(3,2),(2,3)):
        child,sel=inflate(rows,O,{1:a,5:b});r=count(child,sel)
        H=1+2*2**a+2**b+2*2**(a+b)
        U=-2*2**a+2*3**a-2*2**b+10*2**(a+b)+8*3**a*2**b
        need(r['H']==H and r['U']==U and r['delta']==1+2**a-2**b,'coupled literal mismatch')
        out['coupled'].append({'a':a,'b':b,**r})
    print(json.dumps(out,indent=2,sort_keys=True))
if __name__=='__main__':main()
