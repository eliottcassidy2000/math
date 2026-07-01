"""TRUE k(n): allow ANY base orientation (not just transitive). For each free-set F, check if SOME fiber
(fixed orientation of the complement) covers all iso classes. Resolves whether the log2 bound is achievable."""
import numpy as np, math
from itertools import combinations, permutations
def setup(n):
    al=list(combinations(range(n),2)); m=len(al); idx={a:e for e,a in enumerate(al)}
    perms=list(permutations(range(n))); ptrans=[]
    for p in perms:
        dest=[0]*m; flip=[0]*m
        for e,(i,j) in enumerate(al):
            a,b=p[i],p[j]
            if a<b: dest[e]=idx[(a,b)]
            else: dest[e]=idx[(b,a)]; flip[e]=1
        ptrans.append((dest,flip))
    classof=np.empty(2**m,np.int32); cls={}
    for x in range(2**m):
        best=None
        for dest,flip in ptrans:
            y=0
            for e in range(m):
                b=((x>>e)&1)^flip[e]; y|=b<<dest[e]
            if best is None or y<best: best=y
        if best not in cls: cls[best]=len(cls)
        classof[x]=cls[best]
    return al,m,classof,len(cls)
def true_k(n,kmax=None):
    al,m,classof,ncls=setup(n); kmin=math.ceil(math.log2(ncls))
    if kmax is None: kmax=m
    for k in range(kmin,kmax+1):
        found=None
        for F in combinations(range(m),k):
            Fmask=0
            for e in F: Fmask|=1<<e
            subs=[]; s=Fmask
            while True:
                subs.append(s)
                if s==0: break
                s=(s-1)&Fmask
            subs=np.array(subs,np.int64)  # 2^k submasks of F
            nonF=[e for e in range(m) if e not in F]
            # bases = all values on nonF bits (F-bits 0): 2^(m-k)
            bb=[0]
            for e in nonF:
                bb=[b|(bit<<e) for b in bb for bit in (0,1)]
            bases=np.array(bb,np.int64)
            grid=bases[:,None]|subs[None,:]       # (#bases) x (2^k)
            cl=classof[grid]
            cls_sorted=np.sort(cl,axis=1)
            distinct=1+(np.diff(cls_sorted,axis=1)>0).sum(axis=1)
            if (distinct==ncls).any():
                bi=int(np.where(distinct==ncls)[0][0]); found=(F,int(bases[bi])); break
        if found is not None:
            return al,m,ncls,kmin,k,found
    return al,m,ncls,kmin,None,None
for n in [3,4,5,6]:
    al,m,ncls,kmin,k,found=true_k(n)
    F,base=found
    # decode base tournament orientation on nonF arcs -> is it acyclic/transitive?
    Farcs=[al[e] for e in F]
    print(f"n={n}: classes={ncls}, kmin={kmin}, TRUE k(n)={k}  (info bound? {k==kmin})   free-set F={Farcs}")
