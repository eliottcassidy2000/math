#!/usr/bin/env python3
"""
n8_invisible_census_klein_S316.py — klein-2026-07-16-S316

THE COMPLETE n=8 INVISIBLE-PAIR CENSUS (THM-918 follow-up, owner directive).

Method: every n=8 tournament contains an n=7 subtournament, so extending the 456
n=7 iso-class reps by all 2^7 arc patterns of a new vertex covers all 6880 classes
(A000568(8) = 6880, validated). For each of the 58,368 extension tournaments compute
the EXTENDED PANEL exactly:
    (cpA, cpL, H, sorted tau_in, sorted tau_out)
(cpK omitted — THM-924 proves it is a function of cpA). Bucket by panel; within each
bucket compute EXACT canonical forms (refined-partition perm search, not a hash);
distinct canonical forms in one bucket = panel-INVISIBLE class pairs.

Outputs: total classes (must be 6880); invisible groups on the extended panel
(prediction: the double-blind four cones, plus whatever non-cone pairs exist);
invisible count on the OLD one-eyed panel (cpA,cpL,H,tau_in) for comparison
(prediction: >= 27 manufactured + 4 = 31-ish); cpA-tie stats; distinct-cpK counts
n=4..8 (the symmetrization-collapse sequence 1,?,?,11,?).
"""
import itertools, sys, time
import numpy as np
from fractions import Fraction as Fr

t0=time.time()
def log(msg): print(f"[{time.time()-t0:7.1f}s] {msg}", flush=True)

# ---------- n=7 class reps (orbit sweep) ----------
def census(n):
    m=n*(n-1)//2
    pairs=list(itertools.combinations(range(n),2))
    pidx={pr:i for i,pr in enumerate(pairs)}
    remaps=[[(pidx[(min(g[u],g[v]),max(g[u],g[v]))],0 if g[u]<g[v] else 1)
             for (u,v) in pairs] for g in itertools.permutations(range(n))]
    seen=bytearray(1<<m); reps=[]
    for bits in range(1<<m):
        if seen[bits]: continue
        orb=set()
        for tab in remaps:
            out=0
            for i in range(m):
                b=(bits>>i)&1; j,fl=tab[i]; out|=((b^fl)<<j)
            orb.add(out)
        for t in orb: seen[t]=1
        reps.append(bits)
    return reps,pairs

log("n=7 orbit-sweep census...")
reps7,pairs7=census(7)
log(f"n=7 classes: {len(reps7)} (expect 456)")
assert len(reps7)==456

# ---------- build all 58368 extensions as adjacency stack ----------
N7=len(reps7)
A_all=np.zeros((N7*128,8,8),dtype=np.int64)
for ci,bits in enumerate(reps7):
    A7=np.zeros((7,7),dtype=np.int64)
    for i,(u,v) in enumerate(pairs7):
        if (bits>>i)&1: A7[u,v]=1
        else: A7[v,u]=1
    base=np.zeros((8,8),dtype=np.int64); base[:7,:7]=A7
    for p in range(128):
        M=base.copy()
        for u in range(7):
            if (p>>u)&1: M[u,7]=1
            else: M[7,u]=1
        A_all[ci*128+p]=M
NT=A_all.shape[0]
log(f"extension tournaments: {NT}")

# ---------- batched integer Faddeev–LeVerrier charpoly ----------
def charpoly_batch(M):
    """M: (N,8,8) int64 -> coeffs (N,9) int64 for x^8 + c1 x^7 + ... + c8"""
    N=M.shape[0]; n=8
    I=np.eye(n,dtype=np.int64)[None,:,:]
    coeffs=np.zeros((N,n+1),dtype=np.int64); coeffs[:,0]=1
    Nk=np.repeat(I,N,axis=0)
    for k in range(1,n+1):
        MN=np.matmul(M,Nk)
        tr=np.trace(MN,axis1=1,axis2=2)
        assert (tr%k==0).all(), "FL divisibility violated"
        ck=-(tr//k)
        coeffs[:,k]=ck
        Nk=MN+ck[:,None,None]*I
    return coeffs

log("cpA batch...")
cpA=charpoly_batch(A_all)
s_all=A_all.sum(axis=2)                       # out-degrees (N,8)
L_all=-A_all.copy()
L_all[:,np.arange(8),np.arange(8)]=s_all
log("cpL batch...")
cpL=charpoly_batch(L_all)

# ---------- batched H (Hamiltonian path counts) ----------
log("H batch DP...")
H_all=np.zeros(NT,dtype=np.int64)
BATCH=8192
order=sorted(range(1,256),key=lambda m:bin(m).count("1"))
for b0 in range(0,NT,BATCH):
    b1=min(NT,b0+BATCH); B=b1-b0
    A=A_all[b0:b1]
    dp=np.zeros((B,256,8),dtype=np.int64)
    for v in range(8): dp[:,1<<v,v]=1
    for mask in order:
        vs=[v for v in range(8) if (mask>>v)&1]
        for v in vs:
            col=dp[:,mask,v]
            if not col.any(): continue
            for u in range(8):
                if not (mask>>u)&1:
                    dp[:,mask|(1<<u),u]+=col*A[:,v,u]
    H_all[b0:b1]=dp[:,255,:].sum(axis=1)

# ---------- batched tau vectors ----------
log("tau batches...")
def tau_batch(L):
    """sorted in-arborescence vector via 8 batched 7x7 dets (float64 -> rint)"""
    outs=[]
    idx=np.arange(8)
    for r in range(8):
        keep=idx[idx!=r]
        minor=L[:,keep][:,:,keep].astype(np.float64)
        outs.append(np.rint(np.linalg.det(minor)).astype(np.int64))
    T=np.stack(outs,axis=1)
    T.sort(axis=1)
    return T
tin=tau_batch(L_all)
Ar=np.transpose(A_all,(0,2,1)).copy()
sr=Ar.sum(axis=2)
Lr=-Ar; Lr[:,np.arange(8),np.arange(8)]=sr
tout=tau_batch(Lr)

# exact-sample validation of float dets
def det_int(M):
    M=[row[:] for row in M]; n=len(M); sign,prev=1,1
    for k in range(n-1):
        if M[k][k]==0:
            piv=next((i for i in range(k+1,n) if M[i][k]!=0),None)
            if piv is None: return 0
            M[k],M[piv]=M[piv],M[k]; sign=-sign
        for i in range(k+1,n):
            for j in range(k+1,n):
                M[i][j]=(M[i][j]*M[k][k]-M[i][k]*M[k][j])//prev
        prev=M[k][k]
    return sign*M[-1][-1]
rng=np.random.default_rng(5)
for t in rng.integers(0,NT,60):
    L=L_all[t]
    exact=sorted(det_int([[int(L[u,v]) for v in range(8) if v!=r]
                          for u in range(8) if u!=r]) for r in range(8))
    assert list(tin[t])==exact, "float tau mismatch"
log("tau float-vs-Bareiss sample check OK")

# ---------- bucket by panels ----------
log("bucketing...")
ext_buckets={}; old_buckets={}
for t in range(NT):
    ka=(tuple(cpA[t]),tuple(cpL[t]),int(H_all[t]),tuple(tin[t]))
    ke=ka+(tuple(tout[t]),)
    ext_buckets.setdefault(ke,[]).append(t)
    old_buckets.setdefault(ka,[]).append(t)
log(f"extended-panel buckets: {len(ext_buckets)}; old-panel buckets: {len(old_buckets)}")

# ---------- exact canonical form ----------
pairs8=list(itertools.combinations(range(8),2))
pidx8={pr:i for i,pr in enumerate(pairs8)}
def encode(A):
    b=0
    for i,(u,v) in enumerate(pairs8):
        if A[u,v]: b|=(1<<i)
    return b
def canonical(t):
    A=A_all[t]
    lab=[int(x) for x in s_all[t]]
    for _ in range(3):
        newlab=[]
        for v in range(8):
            outs=tuple(sorted(lab[u] for u in range(8) if A[v,u]))
            ins=tuple(sorted(lab[u] for u in range(8) if A[u,v]))
            newlab.append((lab[v],outs,ins))
        comp={x:i for i,x in enumerate(sorted(set(newlab)))}
        lab2=[comp[x] for x in newlab]
        if lab2==lab: break
        lab=lab2
    cells={}
    for v in range(8): cells.setdefault(lab[v],[]).append(v)
    keys=sorted(cells)
    best=None
    from itertools import permutations,product
    pools=[list(permutations(cells[k])) for k in keys]
    npos=1
    for pl in pools: npos*=len(pl)
    for combo in product(*pools):
        gmap={}
        pos=0
        for cell_perm in combo:
            for v in cell_perm:
                gmap[v]=pos; pos+=1
        b=0
        for i,(u,v) in enumerate(pairs8):
            gu,gv=gmap[u],gmap[v]
            a=A[u,v]
            if gu<gv:
                if a: b|=(1<<pidx8[(gu,gv)])
            else:
                if not a: b|=(1<<pidx8[(gv,gu)])
        if best is None or b<best: best=b
    return best

log("canonicalizing within buckets...")
canon_cache={}
def cform(t):
    if t not in canon_cache: canon_cache[t]=canonical(t)
    return canon_cache[t]

total_classes=set()
ext_invisible=[]
done=0
for key,members in ext_buckets.items():
    forms={}
    for t in members:
        forms.setdefault(cform(t),t)
    total_classes.update(forms)
    if len(forms)>1:
        ext_invisible.append((key,forms))
    done+=1
    if done%2000==0: log(f"  {done}/{len(ext_buckets)} buckets")
log(f"TOTAL n=8 classes found: {len(total_classes)} (expect 6880)")
assert len(total_classes)==6880, f"class count {len(total_classes)} != 6880"

old_invisible=[]
for key,members in old_buckets.items():
    forms={}
    for t in members:
        forms.setdefault(cform(t),t)
    if len(forms)>1:
        old_invisible.append((key,forms))

def npairs(groups): return sum(len(f)*(len(f)-1)//2 for _,f in groups)
print()
print("="*80)
print(f"OLD one-eyed panel (cpA,cpL,H,tau_in): {len(old_invisible)} invisible groups, "
      f"{npairs(old_invisible)} pairs")
print(f"EXTENDED two-eyed panel (+tau_out):    {len(ext_invisible)} invisible groups, "
      f"{npairs(ext_invisible)} pairs")
print("="*80)
print()
print("EXTENDED-panel invisible groups (the census truth at n=8):")
for key,forms in sorted(ext_invisible,key=lambda kf:kf[0][2]):
    H=key[2]; tin_k=key[3]; tout_k=key[4]
    descr=[]
    for cf,t in forms.items():
        sc=tuple(sorted(int(x) for x in s_all[t]))
        sink="SINK" if 0 in sc else ("SRC" if 7 in sc else "core")
        descr.append(f"scores={sc} {sink}")
    print(f"  H={H:4d} tau_in={tin_k} tau_out={tout_k}")
    for d in descr: print(f"       {d}")
n_cone=sum(1 for key,forms in ext_invisible
           if all(0 in [int(x) for x in s_all[t]] for t in forms.values()))
n_noncone=len(ext_invisible)-n_cone
print()
print(f"cone groups (all members have a sink): {n_cone}; NON-cone groups: {n_noncone}")

# cpA-tie stats + distinct cpK counts
cpa_groups={}
for cf,t in {cf:t for key,forms in ext_buckets.items() for cf,t in
             [(cform(m),m) for m in forms.values()]}.items() if False else []:
    pass
# simpler: one representative t per class
class_rep={}
for key,members in ext_buckets.items():
    for t in members:
        class_rep.setdefault(cform(t),t)
cpa_tie={}
for cf,t in class_rep.items():
    cpa_tie.setdefault(tuple(cpA[t]),[]).append(cf)
gA=sum(1 for g in cpa_tie.values() if len(g)>1)
cA=sum(len(g) for g in cpa_tie.values() if len(g)>1)
print(f"n=8 cpA-tie groups: {gA} covering {cA} classes (of 6880); distinct cpA: {len(cpa_tie)}")
# distinct cpK via THM-924 fingerprint: cpK(y) = 2^7 [cpA((y-1)/2) + cpA((-y-1)/2)]
def cpk_fp(c):
    cp=[Fr(int(x)) for x in c]
    def ev(x):
        v=Fr(0)
        for ci in cp: v=v*x+ci
        return v
    return tuple(int(Fr(128)*(ev(Fr(y-1,2))+ev(Fr(-y-1,2)))) for y in range(9))
kset={}
for ca in cpa_tie: kset.setdefault(cpk_fp(ca),0)
print(f"n=8 distinct cpK values: {len(kset)}  (symmetrization collapse; n=7 had 11)")
log("DONE")
