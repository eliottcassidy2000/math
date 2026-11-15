import os
HERE=os.path.dirname(os.path.abspath(__file__))
RIB=os.path.join(HERE,'..','..','..','..','04-computation','experiments','ribbon_spectrum_20261104')
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'ribbon_intrep_span2.py')).read().split("# test")[0])
import random
from math import gcd
from functools import reduce
idx={'0':0,'1':1,'2':2,'12':3}
def Mw(b):
    M=np.identity(dim,dtype=object)
    for t in b.split('.'): M=MATS[idx[t]].dot(M)
    return M
u='1.2.1.2.1.2.1.2'
M=Mw(u)
# rank-one: find c primitive integer column direction and r
col=[M[:,j] for j in range(dim) if any(M[:,j])][0]
g=reduce(gcd,[abs(int(v)) for v in col if v]); c=np.array([int(v)//g for v in col],dtype=object)
# r_j = M[i,j]/c[i] for an i with c[i]!=0
i0=[i for i in range(dim) if c[i]][0]
r=np.array([M[i0,j]//c[i0] for j in range(dim)],dtype=object)
assert all(M[i,j]==c[i]*r[j] for i in range(dim) for j in range(dim)), 'not rank one'
print('M_(1.2)^4 = c r^T with c=',list(c),' r=',list(r))
rng=random.Random(7); ok=0; ex=[]
for t in range(300):
    lp=rng.randint(1,9); ls=rng.randint(0,9)
    p=[0]+[rng.randrange(4) for _ in range(lp-1)]
    s=[rng.randrange(4) for _ in range(ls)]
    word=p+[idx[t2] for t2 in u.split('.')]+s
    # respect spans at start
    y=init.copy()
    for q,li in enumerate(word):
        a=L[li]; a=frozenset(z for z in a if z<=q); y=MATS[L.index(a)].dot(y)
    H=int(fin.dot(y))
    yp=init.copy()
    for q,li in enumerate(p):
        a=L[li]; a=frozenset(z for z in a if z<=q); yp=MATS[L.index(a)].dot(yp)
    Lp=int(r.dot(yp))
    ys=c.copy()
    for li in s: ys=MATS[li].dot(ys)
    Rs=int(fin.dot(ys))
    if H==Lp*Rs: ok+=1
    if len(ex)<6 and H>1: ex.append((H,Lp,Rs))
print('H(p u s) = L(p) R(s) on',ok,'of 300 random words; examples (H,L,R):',ex)
