# For the 21-word code in the 62-dim minimal literal transfer: BF units/pairs/triples, and for every cheap cyclic word
# (growth < rho), the bottleneck length L = least number of code units of a prefix of c^infinity with rank <= 1 (mod p).
import numpy as np, itertools, sys, pickle
import os as _o; tok=open(_o.environ.get('SPAN3_INTREP','../../../04-computation/experiments/ribbon_spectrum_20261104/span3_intrep.txt')).read().split(); pos=0
def nx():
    global pos; v=int(tok[pos]); pos+=1; return v
dim=nx(); nl=nx(); [nx() for _ in range(2*dim)]
Mi={}
for _ in range(nl):
    li=nx(); Mi[li]=np.array([[nx() for _ in range(dim)] for _ in range(dim)],dtype=np.int64)
p=33554393
code=[(i,a) for i in range(3) for a in range(1,8)]
CU=[]
for i,a in code:
    P=np.identity(dim,dtype=np.int64)
    for l in (0,)*i+(a,): P=(Mi[l]%p)@P%p
    CU.append(P)
def rank_mod(A):
    A=A.copy()%p; r=0; rows,cols=A.shape
    for c in range(cols):
        nz=np.nonzero(A[r:,c])[0]
        if len(nz)==0: continue
        pv=r+nz[0]
        if pv!=r: A[[r,pv]]=A[[pv,r]]
        inv=pow(int(A[r,c]),-1,p); A[r]=A[r]*inv%p
        col=A[:,c].copy(); col[r]=0; nzr=np.nonzero(col)[0]
        if len(nzr): A[nzr]=(A[nzr]-np.outer(col[nzr],A[r])%p)%p
        r+=1
        if r==rows: break
    return r
def prod(w):
    P=np.identity(dim,dtype=np.int64)
    for u in w: P=CU[u]@P%p
    return P
BF1=[u for u in range(21) if rank_mod(CU[u])>=2]
BF2=[(u,v) for u in BF1 for v in BF1 if rank_mod(CU[v]@CU[u]%p)>=2]
print('BF code units',len(BF1),'of 21; BF pairs',len(BF2))
BF2s=set(BF2); BF3=[]
for (u,v) in BF2:
    for w in BF1:
        if (v,w) in BF2s and rank_mod(CU[w]@CU[v]@CU[u]%p)>=2: BF3.append((u,v,w))
print('BF triples',len(BF3))
Mf={l:Mi[l].astype(float) for l in Mi}
def canon(c): return min(c[k:]+c[:k] for k in range(len(c)))
rho=float(sys.argv[1]) if len(sys.argv)>1 else 1.765
seen=set(); need={}
for L in range(1,5):
    for c in itertools.product(range(21),repeat=L):
        cc=canon(c)
        if cc in seen: continue
        seen.add(cc)
        P=np.identity(dim)
        for u in cc:
            i,a=code[u]
            for l in (0,)*i+(a,): P=Mf[l]@P
        ev=np.abs(np.linalg.eigvals(P)).max(); nv=sum(code[u][0]+1 for u in cc)
        if ev<=0 or ev**(1/nv)>=rho: continue
        # bottleneck length along c^infinity from each rotation start: max over rotations of the first rank<=1 prefix
        worst=0
        for r0 in range(len(cc)):
            seq=list(cc[r0:]+cc[:r0])*12; Q=np.identity(dim,dtype=np.int64); Lb=None
            for k,u in enumerate(seq):
                Q=CU[u]@Q%p
                if rank_mod(Q)<=1: Lb=k+1; break
            worst=max(worst,Lb if Lb else 999)
        need[cc]=worst
from collections import Counter
print('cheap cyclic words (< %.3f) of <= 4 units: %d; bottleneck length distribution (code units):'%(rho,len(need)),sorted(Counter(need.values()).items()))
pickle.dump({'BF1':BF1,'BF2':BF2,'BF3':BF3,'need':need},open('memneed_%s.pkl'%rho,'wb'))
