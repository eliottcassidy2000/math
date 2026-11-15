# Every cyclic word over the 21-unit code with <= 4 units and spectral growth < 1.77 per vertex:
# rank of P^64 modulo two primes (rank mod p >= 2 would certify a genuine bottleneck-free cheap cycle).
import numpy as np, itertools
import sys as _s
tok=open(_s.argv[1] if len(_s.argv)>1 else '../../../04-computation/experiments/ribbon_spectrum_20261104/span3_intrep.txt').read().split(); pos=0
def nx():
    global pos; v=int(tok[pos]); pos+=1; return v
dim=nx(); nl=nx(); [nx() for _ in range(2*dim)]
Mi={}
for _ in range(nl):
    li=nx(); Mi[li]=np.array([[nx() for _ in range(dim)] for _ in range(dim)],dtype=np.int64)
PR=[33554393,33554383]          # primes < 2^25: 62*p^2 < 2^63
code=[(i,a) for i in range(3) for a in range(1,8)]
def prod_mod(word,p):
    P=np.identity(dim,dtype=np.int64)
    for u in word:
        i,a=code[u]
        for l in (0,)*i+(a,): P=(Mi[l]%p)@P%p
    return P
def matpow(P,e,p):
    R=np.identity(dim,dtype=np.int64)
    while e:
        if e&1: R=R@P%p
        P=P@P%p; e>>=1
    return R
def rank_mod(A,p):
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
Mf={l:Mi[l].astype(float) for l in Mi}
def canon(c): return min(c[k:]+c[:k] for k in range(len(c)))
seen=set(); cheap=[]; total=0
for L in range(1,5):
    for c in itertools.product(range(len(code)),repeat=L):
        cc=canon(c)
        if cc in seen: continue
        seen.add(cc); total+=1
        P=np.identity(dim)
        for u in cc:
            i,a=code[u]
            for l in (0,)*i+(a,): P=Mf[l]@P
        ev=np.abs(np.linalg.eigvals(P)).max(); nv=sum(code[u][0]+1 for u in cc)
        if ev>0 and ev**(1.0/nv)<1.77: cheap.append((ev**(1.0/nv),cc))
print('cyclic code words of <= 4 units:',total,'; with spectral growth < 1.77 per vertex:',len(cheap))
genuine=0
for g,cc in cheap:
    rk=[rank_mod(matpow(prod_mod(cc,p),64,p),p) for p in PR]
    if max(rk)>=2:
        genuine+=1; print('  GENUINE bottleneck-free cheap cycle:',g,['0'*code[u][0]+str(code[u][1]) for u in cc],rk)
print('cheap cycles whose 64th power has rank >= 2 modulo some prime:',genuine)
