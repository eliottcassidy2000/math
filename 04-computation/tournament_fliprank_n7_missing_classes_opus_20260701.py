"""Fast vectorized n=7. Find the 2 classes missed by spans{1,3}, and the minimal config that recovers all 456."""
import numpy as np
from itertools import combinations, permutations
n=7; al=list(combinations(range(n),2)); m=len(al); idx={a:e for e,a in enumerate(al)}
pw=(1<<np.arange(m)).astype(np.int64)
perms=list(permutations(range(n)))
srcpos=np.empty((len(perms),m),np.int64); oflip=np.empty((len(perms),m),np.int64)
for pi,p in enumerate(perms):
    dest=np.empty(m,np.int64); flip=np.empty(m,np.int64)
    for e,(i,j) in enumerate(al):
        a,b=p[i],p[j]
        if a<b: dest[e]=idx[(a,b)]; flip[e]=0
        else: dest[e]=idx[(b,a)]; flip[e]=1
    sp=np.empty(m,np.int64); sp[dest]=np.arange(m); srcpos[pi]=sp; oflip[pi]=flip[sp]
def canon_batch(X):  # X: N x m bits -> N canonical ints
    N=X.shape[0]; best=np.full(N,1<<62,np.int64)
    for pi in range(len(perms)):
        OUT=X[:,srcpos[pi]]^oflip[pi]
        val=OUT@pw
        np.minimum(best,val,out=best)
    return best
def fiber_bits(free):  # transitive base, flip subsets of 'free' arc-indices
    k=len(free); N=2**k; X=np.zeros((N,m),np.int8)
    for s in range(N):
        for b in range(k):
            if (s>>b)&1: X[s,free[b]]=1
    return X
def scores(xbits):  # out-degrees sorted, for labeling
    A=np.zeros((n,n),int)
    for e,(i,j) in enumerate(al): 
        if xbits[e]==0: A[i,j]=1
        else: A[j,i]=1
    return tuple(sorted(A.sum(1)))
# FULL set via standard tiling (base Ham path span1 fixed, free all 15 tiles spans>=2)
tiles=[e for e,(i,j) in enumerate(al) if j-i>=2]
full=set(canon_batch(fiber_bits(tiles)).tolist())
print(f"standard tiling (free {len(tiles)} tiles) covers {len(full)} classes (should be 456)")
# spans {1,3} fixed
free13=[e for e,(i,j) in enumerate(al) if (j-i) not in {1,3}]
c13=canon_batch(fiber_bits(free13))
set13=set(c13.tolist())
print(f"spans{{1,3}} fixed, free {len(free13)} (k=11): covers {len(set13)}; MISSING {len(full-set13)} classes")
# identify missing by score sequence: map canonical int -> a representative's score seq
missing=full-set13
# build map from canonical int to score seq using standard tiling fiber
Xs=fiber_bits(tiles); cs=canon_batch(Xs)
rep={}
for r in range(Xs.shape[0]):
    if cs[r] in missing and cs[r] not in rep: rep[cs[r]]=scores(Xs[r])
print(f"  missing classes' score sequences: {list(rep.values())}")
# minimal fix: try freeing ONE more fixed arc (k=12) -> does it reach 456?
fixed13=[e for e,(i,j) in enumerate(al) if (j-i) in {1,3}]
best=None
for e in fixed13:
    fr=free13+[e]
    cc=set(canon_batch(fiber_bits(fr)).tolist())
    if len(cc)==456: 
        best=(e,al[e],(al[e][1]-al[e][0])); break
print(f"  freeing one more arc {best[1]} (span {best[2]}) -> covers all 456 at k=12" if best else "  no single extra free arc reaches 456 (need k>=13)")
