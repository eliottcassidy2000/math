import numpy as np, itertools, sys
R='/Users/e/Documents/GitHub/math/04-computation/experiments/ribbon_spectrum_20261104/span3_intrep.txt'
tok=open(R).read().split(); pos=0
def nxt():
    global pos; v=tok[pos]; pos+=1; return int(v)
dim=nxt(); nl=nxt()
init=[nxt() for _ in range(dim)]; fin=[nxt() for _ in range(dim)]
MI={}
for _ in range(nl):
    li=nxt(); MI[li]=np.array([[nxt() for _ in range(dim)] for _ in range(dim)],dtype=np.int64)
masks=sorted(MI.keys())
PR=[33554393,33554383]
def rank_mod(A,p):
    A=A.copy()%p; r=0; rows,cols=A.shape
    for c in range(cols):
        piv=None
        for i in range(r,rows):
            if A[i,c]: piv=i;break
        if piv is None: continue
        A[[r,piv]]=A[[piv,r]]; inv=pow(int(A[r,c]),p-2,p); A[r]=(A[r]*inv)%p
        nz=np.nonzero(A[:,c])[0]
        for i in nz:
            if i!=r: A[i]=(A[i]-A[i,c]*A[r])%p
        r+=1
        if r==rows: break
    return r
I=np.identity(dim,dtype=np.int64)
found=[]; ev1=0; tot=0
for L in range(1,int(sys.argv[1])+1):
    for w in itertools.product(masks,repeat=L):
        tot+=1
        p=PR[0]
        M=I.copy()
        for li in w: M=(MI[li]@M)%p
        A=(M-I)%p
        r1=rank_mod(A,p)
        if r1==dim: continue
        ev1+=1
        r2=rank_mod((A@A)%p,p)
        if r2<r1:
            # confirm with second prime
            q=PR[1]; M2=I.copy()
            for li in w: M2=(MI[li]@M2)%q
            A2=(M2-I)%q; s1=rank_mod(A2,q); s2=rank_mod((A2@A2)%q,q)
            if s2<s1: found.append((w,dim-r1,r1-r2))
print('blocks tested',tot,' with eigenvalue 1:',ev1,' with a Jordan block at 1 (unipotent part):',len(found))
for f in found[:12]: print('  block masks',f[0],' geometric mult',f[1],' rank drop of (M-I)^2:',f[2])
