"""Rigorous k(7) >= D(7) (iso-diameter). Compute eccentricities from transitive (=R) and Paley, and push to
the diametral pair. dist(class C, transitive)=MFAS(C); dist(C, Paley)=min over Paley's 240-orbit."""
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
def orbit_ints(xint):
    A=np.zeros((n,n),np.int8)
    for e,(i,j) in enumerate(al):
        if (xint>>e)&1: A[j,i]=1
        else: A[i,j]=1
    O=set()
    for p in permL:
        B=A[np.ix_(p,p)]; y=0
        for e,(i,j) in enumerate(al):
            if B[j,i]==1: y|=1<<e
        O.add(y)
    return np.array(sorted(O),np.int64)
tiles=[e for e,(i,j) in enumerate(al) if (j-i)>=2]
X=np.zeros((2**len(tiles),m),np.int8)
for s in range(2**len(tiles)):
    for b in range(len(tiles)):
        if (s>>b)&1: X[s,tiles[b]]=1
reps=np.array(sorted(set(canon_batch(X).tolist())),np.int64)
pct=np.array([bin(i).count("1") for i in range(1<<16)],np.int16)
def pcount(a): return pct[a&0xFFFF]+pct[(a>>16)&0xFFFF]
def dist_to_orbit(reps, orb):  # for each rep, min popcount(rep XOR o) over orbit
    out=np.empty(len(reps),np.int64)
    for i,r in enumerate(reps):
        out[i]=int(pcount(r^orb).min())
    return out
# Paley rep, orbit
Ap=np.zeros((n,n),np.int8)
for i in range(n):
    for d in (1,2,4): Ap[i,(i+d)%n]=1
xp=0
for e,(i,j) in enumerate(al):
    if Ap[j,i]==1: xp|=1<<e
pal_orb=orbit_ints(xp); trans_orb=orbit_ints(0)
dP=dist_to_orbit(reps,pal_orb); dT=dist_to_orbit(reps,trans_orb)
eccP=int(dP.max()); eccT=int(dT.max())
print(f"ecc(transitive)=R(7)={eccT}   ecc(Paley)={eccP}")
# farthest class from Paley -> compute ITS eccentricity (push diameter)
far=reps[int(np.argmax(dP))]; far_orb=orbit_ints(int(far))
dFar=dist_to_orbit(reps,far_orb); eccFar=int(dFar.max())
print(f"farthest-from-Paley class: dist={eccP}; ITS eccentricity={eccFar}")
D_lb=max(eccT,eccP,eccFar)
print(f"=> D(7) >= {D_lb}  (rigorous lower bound on the iso-diameter)")
print(f"   k(7) >= D(7) >= {D_lb}; info-floor=9, R(7)=7, k(7)=12(evidenced).  D beats floor? {D_lb>9}")
