# Exact M4 ceiling at the KEY (iso-class) level, fast: dedupe by the realized arc-pattern
# BEFORE canon. For each residue multiset, forced arcs + tie pairs; enumerate the 2^t tie
# orientations that are REALIZABLE by a total order. Collect distinct arc-patterns, canon once.
from math import gcd
from itertools import combinations, permutations, combinations_with_replacement, product
MOD=14; UNITS=[a for a in range(1,MOD) if gcd(a,MOD)==1]
def depth(r): r%=MOD; return min(r,MOD-r)
WIN=[[0]*MOD for _ in range(MOD)]
for x in range(MOD):
    for y in range(MOD):
        wi=wj=0
        for a in UNITS:
            di=depth(x*a); dj=depth(y*a)
            if di>dj: wi+=1
            elif dj>di: wj+=1
        WIN[x][y]=1 if wi>wj else (-1 if wj>wi else 0)
m=5
PAIRS=[(i,j) for i in range(m) for j in range(i+1,m)]
PERMS=list(permutations(range(m)))
def canon_arc(arc):  # arc dict (i<j)->1 if i->j
    full=[[0]*m for _ in range(m)]
    for (i,j) in PAIRS:
        full[i][j]=arc[(i,j)]; full[j][i]=1-arc[(i,j)]
    best=None
    for p in PERMS:
        b=tuple(full[p[i]][p[j]] for i in range(m) for j in range(m) if i!=j)
        if best is None or b<best: best=b
    return best
def sigkey(key):
    full=[[0]*m for _ in range(m)]; idx=0
    for i in range(m):
        for j in range(m):
            if i!=j: full[i][j]=key[idx]; idx+=1
    H=0
    for p in PERMS:
        if all(full[p[k]][p[k+1]] for k in range(m-1)): H+=1
    c3=0
    for a,b,cc in combinations(range(m),3):
        if full[a][b] and full[b][cc] and full[cc][a]: c3+=1
        if full[a][cc] and full[cc][b] and full[b][a]: c3+=1
    sc=tuple(sorted(sum(full[i][j] for j in range(m) if j!=i) for i in range(m)))
    return (H,c3,sc)
arc_patterns=set()
for res in combinations_with_replacement(range(MOD),m):
    forced={}; tied=[]
    for (i,j) in PAIRS:
        w=WIN[res[i]][res[j]]
        if w==1: forced[(i,j)]=1
        elif w==-1: forced[(i,j)]=0
        else: tied.append((i,j))
    if not tied:
        arc_patterns.add(tuple(forced[p] for p in PAIRS)); continue
    # realizable orientations = those induced by SOME total order on the 5 slots
    for order in PERMS:
        arc=dict(forced)
        for (i,j) in tied:
            arc[(i,j)]=1 if order[i]<order[j] else 0
        arc_patterns.add(tuple(arc[p] for p in PAIRS))
# canon each distinct pattern
keys=set()
for pat in arc_patterns:
    arc={p:pat[idx] for idx,p in enumerate(PAIRS)}
    keys.add(canon_arc(arc))
print(f"distinct arc-patterns: {len(arc_patterns)}")
print(f"EXACT M4 ceiling: {len(keys)} of 12 iso classes reachable")
# 12 free classes
free={}
for bits in product([0,1],repeat=len(PAIRS)):
    arc={p:bits[idx] for idx,p in enumerate(PAIRS)}
    k=canon_arc(arc)
    free[k]=sigkey(k)
reached=[(free[k],k in keys) for k in free]
print("\nclass-by-class (signature, reached?):")
for s,r in sorted(reached): print(f"   {s}: {'REACHED' if r else 'NOT reached'}")
FORB=[(9,3,(1,1,2,3,3)),(13,4,(1,2,2,2,3)),(15,4,(1,2,2,2,3)),(15,5,(2,2,2,2,2))]
print("\nForbidden-signature reachability (by iso class):")
for fs in FORB:
    for k in [k for k,s in free.items() if s==fs]:
        print(f"   sig {fs} key {k}: {'REACHED' if k in keys else 'NOT reached'}")
