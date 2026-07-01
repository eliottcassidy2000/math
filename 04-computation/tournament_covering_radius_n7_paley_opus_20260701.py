"""n=7 covering radius R(7)=max MFAS (optimized). MFAS(C)=min iso-dist to transitive orbit. Check if the
extremizer is the Paley heptagon (|Aut|=21) => depth-max and symmetry-max COINCIDE at the Paley prime."""
import numpy as np
from itertools import combinations, permutations
n=7; al=list(combinations(range(n),2)); m=len(al); idx={a:e for e,a in enumerate(al)}
pw=(1<<np.arange(m)).astype(np.int64); perms=list(permutations(range(n))); permL=[list(p) for p in perms]
srcpos=np.empty((len(perms),m),np.int64); oflip=np.empty((len(perms),m),np.int64)
for pi,p in enumerate(perms):
    dest=np.empty(m,np.int64); flip=np.empty(m,np.int64)
    for e,(i,j) in enumerate(al):
        a,b=p[i],p[j]
        if a<b: dest[e]=idx[(a,b)]; flip[e]=0
        else: dest[e]=idx[(b,a)]; flip[e]=1
    sp=np.empty(m,np.int64); sp[dest]=np.arange(m); srcpos[pi]=sp; oflip[pi]=flip[sp]
def canon_batch(X):
    best=np.full(X.shape[0],1<<62,np.int64)
    for pi in range(len(perms)): np.minimum(best,(X[:,srcpos[pi]]^oflip[pi])@pw,out=best)
    return best
tiles=[e for e,(i,j) in enumerate(al) if (j-i)>=2]
X=np.zeros((2**len(tiles),m),np.int8)
for s in range(2**len(tiles)):
    for b in range(len(tiles)):
        if (s>>b)&1: X[s,tiles[b]]=1
reps=sorted(set(canon_batch(X).tolist()))       # 456 canonical bit-ints
print(f"#classes={len(reps)}")
# transitive orbit = all linear orders as bit-ints (orbit of x=0)
def orbit_int(x):
    O=set()
    A=np.zeros((n,n),np.int8)
    for e,(i,j) in enumerate(al):
        if (x>>e)&1: A[j,i]=1
        else: A[i,j]=1
    for p in permL:
        B=A[np.ix_(p,p)]; y=0
        for e,(i,j) in enumerate(al):
            if B[j,i]==1: y|=1<<e
        O.add(y)
    return O
trans=list(orbit_int(0)); print(f"transitive orbit size={len(trans)}")
trans=np.array(trans,np.int64)
pct=np.array([bin(i).count("1") for i in range(1<<16)],np.int16)
def popcount(a): return pct[a&0xFFFF]+pct[(a>>16)&0xFFFF]
def aut(x):
    A=np.zeros((n,n),np.int8)
    for e,(i,j) in enumerate(al):
        if (x>>e)&1: A[j,i]=1
        else: A[i,j]=1
    return sum(1 for p in permL if (A[np.ix_(p,p)]==A).all())
mfas={}
for r in reps:
    mfas[r]=int(popcount(np.int64(r)^trans).min())
R=max(mfas.values())
ext=[r for r in reps if mfas[r]==R]
print(f"R(7)=max MFAS={R}; #extremizers={len(ext)}; their |Aut|={sorted({aut(r) for r in ext},reverse=True)}")
# Paley heptagon rep + its MFAS
Ap=np.zeros((n,n),np.int8)
for i in range(n):
    for d in (1,2,4): Ap[i,(i+d)%n]=1
xp=0
for e,(i,j) in enumerate(al):
    if Ap[j,i]==1: xp|=1<<e
# canonicalize Paley to compare
cp=int(canon_batch(np.array([[(xp>>e)&1 for e in range(m)]],np.int8))[0])
print(f"Paley heptagon: |Aut|={aut(cp)}, MFAS={mfas[cp]}  => Paley is a covering-radius extremizer? {mfas[cp]==R}")
print(f"max|Aut|(7)={max(aut(r) for r in reps)}  => at n=7 (Paley PRIME) depth-max & symmetry-max COINCIDE at Paley={mfas[cp]==R and aut(cp)==max(aut(r) for r in reps)}")
