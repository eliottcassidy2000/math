#!/usr/bin/env python3
"""rangeclass_parallel_kps_S128c33.py -- the staircase range classes as cut-open parallel classes:
does wiggly self-loop rate decompose as f(range class) + boundary-of-cut correction?"""
import sys
from itertools import permutations
sys.stdout.reconfigure(line_buffering=True)
for n in (5,6):
    tiles=[(x,y) for y in range(1,n-1) for x in range(n,y+1,-1) if x-y>=2]
    m=len(tiles)
    def build(t):
        B=[[False]*n for _ in range(n)]
        for k in range(2,n+1): B[k-1][k-2]=True
        for i,(x,y) in enumerate(tiles):
            if (t>>i)&1: B[x-1][y-1]=True
            else: B[y-1][x-1]=True
        return B
    def canon(B):
        best=None
        for p in permutations(range(n)):
            key=tuple(B[p[i]][p[j]] for i in range(n) for j in range(n))
            if best is None or key<best: best=key
        return best
    # cache canon per tiling
    can={t:canon(build(t)) for t in range(1<<m)}
    sl={i:0 for i in range(m)}
    for t in range(1<<m):
        for i in range(m):
            if can[t]==can[t^(1<<i)]: sl[i]+=1
    print("n=%d per-tile SL%%:"%n)
    byclass={}
    for i,(x,y) in enumerate(tiles):
        r=x-y
        pct=100.0*sl[i]/(1<<m)
        # distance from the cut: the base path is the cut-open cycle; wrap pair = {1,n}
        dcut=min(y-1, n-x)   # 0 = touches an end vertex of the path
        byclass.setdefault(r,[]).append((pct,dcut,(x,y)))
        print("   tile (%d,%d) range %d dcut %d : SL %.2f%%"%(x,y,r,dcut,pct))
    print("  class summary (range: mean%% spread%% | by dcut):")
    for r in sorted(byclass):
        v=[p for p,_,_ in byclass[r]]
        bd={}
        for p,d,_ in byclass[r]: bd.setdefault(d,[]).append(p)
        bys=" ".join("d%d:%s"%(d,["%.2f"%x for x in bd[d]]) for d in sorted(bd))
        print("   r=%d: mean %.2f spread %.2f | %s"%(r,sum(v)/len(v),max(v)-min(v),bys))
print("DONE")
