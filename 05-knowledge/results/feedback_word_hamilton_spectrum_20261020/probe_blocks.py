#!/usr/bin/env python3
from itertools import product,combinations
import numpy as np
from probe_affine import matrices,words

ms,u,v=matrices();ws=words(range(4),2)
blocks=[w for w in ws if w]
I=np.eye(len(u),dtype=np.int64)
B=[]
for w in blocks:
    mat=I.copy()
    for a in w:mat=mat@ms[a]
    B.append(mat)
R=[];C=[]
for w in ws:
    r=u.copy();c=v.copy()
    for a in w:r=r@ms[a]
    for a in reversed(w):c=ms[a]@c
    R.append(r);C.append(c)
R=np.array(R);C=np.array(C).T
base=R@C
hits=[]
for lam in [2,3,4,5,9]:
    D=[x-lam*I for x in B]
    V=[x@C for x in D]
    V0=[R@x for x in V]
    for a,b in combinations(range(len(blocks)),2):
        ok=(V0[a]!=V0[b])
        for w in words((a,b),2):
            rw=R
            for x in w:rw=rw@B[x]
            ok &= ((rw@V[a])==V0[a]) & ((rw@V[b])==V0[b])
            if not np.any(ok):break
        for i,j in zip(*np.nonzero(ok)):
            hits.append((ws[i],ws[j],blocks[a],blocks[b],lam,int(base[i,j]),int(V0[a][i,j]),int(V0[b][i,j])))
print('prefix/suffix <=2, blocks length1..2, slopes2,3,4,5,9, middle block depth<=2')
print('candidates',len(hits))
print(hits[:50])
