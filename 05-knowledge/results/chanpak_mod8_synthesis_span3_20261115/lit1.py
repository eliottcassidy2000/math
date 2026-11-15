import numpy as np, itertools, sys
import sys as _s
tok=open(_s.argv[1] if len(_s.argv)>1 else '../../../04-computation/experiments/ribbon_spectrum_20261104/span3_intrep.txt').read().split(); pos=0
def nx():
    global pos; v=int(tok[pos]); pos+=1; return v
dim=nx(); nl=nx(); init=np.array([nx() for _ in range(dim)],dtype=object); fin=np.array([nx() for _ in range(dim)],dtype=object)
M={}
for _ in range(nl):
    li=nx(); M[li]=np.array([[nx() for _ in range(dim)] for _ in range(dim)],dtype=object)
print('dim',dim,'letters',sorted(M))
p=2147483629
def rk(A):
    A=np.array(A,dtype=object)%p; A=A.astype(np.int64); r=0; rows,cols=A.shape
    for c in range(cols):
        nz=np.nonzero(A[r:,c])[0]
        if len(nz)==0: continue
        pv=r+nz[0]
        if pv!=r: A[[r,pv]]=A[[pv,r]]
        inv=pow(int(A[r,c]),-1,p); A[r]=(A[r].astype(object)*inv%p).astype(np.int64)
        col=A[:,c].copy(); col[r]=0; nzr=np.nonzero(col)[0]
        if len(nzr): A[nzr]=((A[nzr].astype(object)-np.outer(col[nzr].astype(object),A[r].astype(object)))%p).astype(np.int64)
        r+=1
        if r==rows: break
    return r
Mm={a:(M[a]%p).astype(object) for a in M}
for a in sorted(M):
    P=np.identity(dim,dtype=object); rs=[]
    for k in range(1,7):
        P=Mm[a].dot(P)%p; rs.append(rk(P))
    print('letter %s (back arcs %s): power ranks %s'%(a,[d for d in (1,2,3) if (a>>(d-1))&1],rs))
# H of a literal word: fin . M_w ... init  (first vertex letter 0 already in init?)
def H(word):
    v=init.copy()
    for a in word: v=M[a].dot(v)
    return int(fin.dot(v))
print('H checks: 3-cycle-ish words', [ (w,H(w)) for w in [(0,),(0,0),(2,),(0,2),(1,2),(0,4)]])
