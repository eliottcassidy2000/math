#!/usr/bin/env python3
from itertools import permutations,product
from collections import Counter,defaultdict
import argparse, json
from pathlib import Path

def census(n):
    pairs=[(i,j) for i in range(n) for j in range(i+2,n)]
    ix={p:k for k,p in enumerate(pairs)}
    def rev(x):return sum(((x>>k)&1)<<ix[n-1-j,n-1-i] for k,(i,j) in enumerate(pairs))
    def canon(a):
        score=defaultdict(list)
        for i,row in enumerate(a):score[sum(row)].append(i)
        best=None
        for chunks in product(*(permutations(score[s]) for s in sorted(score))):
            p=sum(chunks,())
            val=tuple(a[p[i]][p[j]] for i in range(n) for j in range(i+1,n))
            if best is None or val<best:best=val
        return best
    masses=Counter();blue=Counter();classes=set();word_nodes=[]
    for x in range(1<<len(pairs)):
        a=[[0]*n for _ in range(n)]
        for i in range(n-1):a[i][i+1]=1
        for k,(i,j) in enumerate(pairs):
            v=(x>>k)&1;a[i][j]=v;a[j][i]=1-v
        c=canon(a);o=canon(list(map(list,zip(*a))));key=min(c,o)
        classes.add(c);masses[key]+=1;word_nodes.append(key)
        if rev(x)==x:blue[key]+=1
    nodes=[{'canonical':''.join(map(str,key)),'degree':d,'blue_degree':blue[key]} for key,d in sorted(masses.items())]
    full=(1<<len(pairs))-1
    blue_loops=sum(x<(x^full) and word_nodes[x]==word_nodes[x^full] and rev(x)==x for x in range(full+1))
    black_loops=sum(x<(x^full) and word_nodes[x]==word_nodes[x^full] and rev(x)!=x for x in range(full+1))
    return {'blue_loops':blue_loops,'black_loops':black_loops,'n':n,'A':len(classes),'V':len(nodes),'S':sum(z['blue_degree']>0 for z in nodes),'N':sum(masses.values()),'B':sum(blue.values()),'M2':sum(d*d for d in masses.values()),'nodes':nodes}
if __name__=='__main__':
    out=[census(n) for n in range(3,7)]
    parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true');args=parser.parse_args()
    path=Path(__file__).with_name('native_reference.json');encoded=json.dumps(out,indent=2,sort_keys=True)+'\n'
    if args.check:
        if path.read_text()!=encoded:raise RuntimeError('native reference replay differs')
    else:path.write_text(encoded)
    print(json.dumps([{k:v for k,v in z.items() if k!='nodes'} for z in out],sort_keys=True))
