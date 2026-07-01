"""RESOLVE kappa(7)=11 vs 12. Best 11-free-set F=spans{2,4,5,6} (my 454-config, transitive base misses the
regular tournament). Question: does ANY base orientation of the 10 fixed span{1,3} arcs host all 456 classes?
Strategy: only bases hosting the 2 missing classes can work -> compute those candidate bases -> test full coverage."""
import numpy as np
from itertools import combinations, permutations
n=7; al=list(combinations(range(n),2)); m=len(al); idx={a:e for e,a in enumerate(al)}
pw=(1<<np.arange(m)).astype(np.int64); perms=list(permutations(range(n)))
srcpos=np.empty((len(perms),m),np.int64); oflip=np.empty((len(perms),m),np.int64)
permmat=np.array(perms)
for pi,p in enumerate(perms):
    dest=np.empty(m,np.int64); flip=np.empty(m,np.int64)
    for e,(i,j) in enumerate(al):
        a,b=p[i],p[j]
        if a<b: dest[e]=idx[(a,b)]; flip[e]=0
        else: dest[e]=idx[(b,a)]; flip[e]=1
    sp=np.empty(m,np.int64); sp[dest]=np.arange(m); srcpos[pi]=sp; oflip[pi]=flip[sp]
def canon_batch(X):
    best=np.full(X.shape[0],1<<62,np.int64)
    for pi in range(len(perms)):
        np.minimum(best,(X[:,srcpos[pi]]^oflip[pi])@pw,out=best)
    return best
fixed=[e for e,(i,j) in enumerate(al) if (j-i) in {1,3}]      # 10 arcs
free =[e for e,(i,j) in enumerate(al) if (j-i) in {2,4,5,6}]  # 11 arcs
tiles=[e for e,(i,j) in enumerate(al) if (j-i)>=2]            # standard tiling (all 456)
def fiber_bits(free_arcs, base=0):
    k=len(free_arcs); X=np.zeros((2**k,m),np.int8)
    base_bits=[(base>>e)&1 for e in range(m)]
    for e in range(m):
        if base_bits[e]: X[:,e]=1
    for s in range(2**k):
        for b in range(k):
            if (s>>b)&1: X[s,free_arcs[b]]^=1
    return X
# full 456 and the transitive-base spans{1,3} fiber
full=set(canon_batch(fiber_bits(tiles)).tolist())
Xs=fiber_bits(tiles); cs=canon_batch(Xs)
X13=fiber_bits(free,0); c13=canon_batch(X13); set13=set(c13.tolist())
missing=sorted(full-set13)
print(f"full={len(full)}, spans{{1,3}} transitive-base fiber={len(set13)}, missing={len(missing)}")
# representative bit-vector for each missing class + |Aut|
reps={}
for r in range(Xs.shape[0]):
    if cs[r] in missing and cs[r] not in reps: reps[cs[r]]=int((Xs[r].astype(np.int64)*pw).sum())
def aut_and_score(xint):
    A=np.zeros((n,n),np.int8)
    for e,(i,j) in enumerate(al):
        if (xint>>e)&1: A[j,i]=1
        else: A[i,j]=1
    aut=sum(1 for p in perms if (A[np.ix_(list(p),list(p))]==A).all())
    return aut, tuple(sorted(A.sum(1)))
for c in missing:
    a,sc=aut_and_score(reps[c]); print(f"  missing class score={sc} |Aut|={a}  (labeled reps={5040//a})")
# candidate bases = restriction-to-fixed of the ORBIT of each missing rep; a base must host BOTH missing
def compatible_bases(xint):
    A=np.zeros((n,n),np.int8)
    for e,(i,j) in enumerate(al):
        if (xint>>e)&1: A[j,i]=1
        else: A[i,j]=1
    bases=set()
    for p in perms:
        B=A[np.ix_(list(p),list(p))]
        xb=0
        for e,(i,j) in enumerate(al):
            b=1 if B[j,i]==1 else 0
            if b: xb|=1<<e
        # restrict to fixed arcs only
        base=0
        for e in fixed:
            if (xb>>e)&1: base|=1<<e
        bases.add(base)
    return bases
Bs=[compatible_bases(reps[c]) for c in missing]
cand=set.intersection(*Bs)
print(f"candidate bases hosting BOTH missing classes: {len(cand)}")
# test each candidate for full 456 coverage
found=None
for bi,base in enumerate(sorted(cand)):
    cov=len(set(canon_batch(fiber_bits(free,base)).tolist()))
    if cov==456: found=base; print(f"  base #{bi}: FULL 456 COVERAGE -> kappa(7)=11 ACHIEVED with free=spans{{2,4,5,6}}!"); break
if found is None:
    print(f"  tested all {len(cand)} candidate bases: NONE reaches 456 (best free-set spans{{2,4,5,6}} cannot do 11 under any base) => strong evidence kappa(7)>=12")
