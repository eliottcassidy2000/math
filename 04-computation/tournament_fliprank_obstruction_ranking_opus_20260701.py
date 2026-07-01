"""EXTEND: is the Paley heptagon the ROBUST obstruction? Over many 11-free-sets (transitive base) count how
often each class is missed; rank by |Aut|. Plus beta(Paley)=MFAS, and the max|Aut| obstruction sequence."""
import numpy as np, random
from itertools import combinations, permutations
n=7; al=list(combinations(range(n),2)); m=len(al); idx={a:e for e,a in enumerate(al)}
pw=(1<<np.arange(m)).astype(np.int64); perms=list(permutations(range(n))); P=len(perms)
srcpos=np.empty((P,m),np.int64); oflip=np.empty((P,m),np.int64); permL=[list(p) for p in perms]
for pi,p in enumerate(perms):
    dest=np.empty(m,np.int64); flip=np.empty(m,np.int64)
    for e,(i,j) in enumerate(al):
        a,b=p[i],p[j]
        if a<b: dest[e]=idx[(a,b)]; flip[e]=0
        else: dest[e]=idx[(b,a)]; flip[e]=1
    sp=np.empty(m,np.int64); sp[dest]=np.arange(m); srcpos[pi]=sp; oflip[pi]=flip[sp]
def canon_batch(X):
    best=np.full(X.shape[0],1<<62,np.int64)
    for pi in range(P): np.minimum(best,(X[:,srcpos[pi]]^oflip[pi])@pw,out=best)
    return best
def fiber(free_arcs, base=0):
    k=len(free_arcs); X=np.zeros((2**k,m),np.int8)
    for e in range(m):
        if (base>>e)&1: X[:,e]=1
    for s in range(2**k):
        for b in range(k):
            if (s>>b)&1: X[s,free_arcs[b]]^=1
    return X
def aut(xint):
    A=np.zeros((n,n),np.int8)
    for e,(i,j) in enumerate(al):
        if (xint>>e)&1: A[j,i]=1
        else: A[i,j]=1
    return sum(1 for p in permL if (A[np.ix_(p,p)]==A).all())
def mfas(xint):  # min backward arcs over all orderings = min feedback arc set
    A=np.zeros((n,n),np.int8)
    for e,(i,j) in enumerate(al):
        if (xint>>e)&1: A[j,i]=1
        else: A[i,j]=1
    best=99
    for p in permL:
        back=sum(1 for a in range(n) for b in range(a+1,n) if A[p[b],p[a]]==1)
        best=min(best,back)
    return best
# full class set + rep per class (from standard tiling)
tiles=[e for e,(i,j) in enumerate(al) if (j-i)>=2]
Xs=fiber(tiles); cs=canon_batch(Xs)
rep={}
for r in range(Xs.shape[0]):
    if cs[r] not in rep: rep[cs[r]]=int((Xs[r].astype(np.int64)*pw).sum())
full=set(rep); print(f"n=7: {len(full)} classes")
# miss-frequency over random 11-free-sets (transitive base)
random.seed(3); NF=25; miss=dict.fromkeys(full,0)
for t in range(NF):
    F=sorted(random.sample(range(m),11))
    cov=set(canon_batch(fiber(F,0)).tolist())
    for c in full-cov: miss[c]+=1
hard=sorted(full,key=lambda c:-miss[c])[:6]
print(f"Top obstruction classes over {NF} random 11-free-sets (transitive base):")
for c in hard:
    a=aut(rep[c]); b=mfas(rep[c]); print(f"   missed {miss[c]:2d}/{NF}  |Aut|={a:2d}  MFAS(beta)={b}")
# the Paley heptagon rep: rotational QR i->j iff (j-i)mod7 in {1,2,4}
Apal=np.zeros((n,n),np.int8)
for i in range(n):
    for d in (1,2,4): Apal[i,(i+d)%n]=1
xpal=0
for e,(i,j) in enumerate(al):
    if Apal[j,i]==1: xpal|=1<<e
print(f"\nPaley heptagon (QR Z/7): |Aut|={aut(xpal)}  MFAS(beta)={mfas(xpal)}  (min tile-flips from transitive)")
print(f"max|Aut| n=7 = 21 (Paley). klein-S72: max|Aut|=3,3,5,9 (n=3..6) => sequence 1,3,3,5,9,21 (n=2..7).")
print(f"k(n)-floor = 0,0,0,1,3 (n=3..7); the |Aut|=21 leap at n=7 is the extra bit (11->12).")
