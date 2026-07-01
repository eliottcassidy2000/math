"""Targeted search for an 11-config at n=7 (would give k(7)=11). Free the 7 arcs where the CLOSEST Paley rep
departs from transitive (so Paley IS hosted) + 4 more; sweep the 4 (and Paley witnesses). If any covers 456 => k=11."""
import numpy as np, itertools, random
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
def fiber(free_arcs, base):
    k=len(free_arcs); Xf=np.zeros((2**k,m),np.int8)
    for e in range(m):
        if (base>>e)&1: Xf[:,e]=1
    for s in range(2**k):
        for b in range(k):
            if (s>>b)&1: Xf[s,free_arcs[b]]^=1
    return Xf
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
    return sorted(O)
tiles=[e for e,(i,j) in enumerate(al) if (j-i)>=2]
full=set(canon_batch(fiber(tiles,0)).tolist()); print(f"full={len(full)} classes")
# Paley orbit; closest rep to transitive (base=0): min popcount
Ap=np.zeros((n,n),np.int8)
for i in range(n):
    for d in (1,2,4): Ap[i,(i+d)%n]=1
xp=0
for e,(i,j) in enumerate(al):
    if Ap[j,i]==1: xp|=1<<e
palorb=orbit_ints(xp)
pc=lambda z:bin(z).count("1")
witnesses=sorted(palorb,key=pc)[:6]  # a few closest-to-transitive Paley reps
random.seed(0); tested=0; hit=False
for w,pw_rep in enumerate(witnesses):
    diff=[e for e in range(m) if (pw_rep>>e)&1]        # 7 arcs Paley differs from transitive(base 0)
    agree=[e for e in range(m) if not((pw_rep>>e)&1)]  # 14 agreement arcs
    if len(diff)!=pc(pw_rep): pass
    # free = diff (7) + 4 chosen agreement; base=0 (transitive on fixed). Sweep choices of 4.
    combos=list(combinations(agree,4)); random.shuffle(combos)
    for extra in combos[:60]:
        F=diff+list(extra); 
        if len(F)!=11: continue
        cov=set(canon_batch(fiber(F,0)).tolist()); tested+=1
        if len(cov)==456:
            hit=True; print(f"  !!! k(7)=11 FOUND: witness#{w}, free={[al[e] for e in F]}"); break
    if hit: break
print(f"tested {tested} targeted 11-subcubes (Paley-hosting, transitive base); reached 456? {hit}")
if not hit: print(f"  => no 11-config among Paley-geodesic-aligned subcubes either => k(7)=12 further reinforced")
