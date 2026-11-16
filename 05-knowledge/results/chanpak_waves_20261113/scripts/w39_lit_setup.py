# Wave 39 (mac-mini-2026-11-15-frame). Run from a working directory; w39_lit_setup.py writes lit_base.pkl used by the others.
# Wave 39 prototype: literal span-3 transfer, minimal projection, floor-unit rank-one powers, dictionary sizes.
import sys, pickle, itertools, numpy as np
import os
RIB=os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', '04-computation', 'experiments', 'ribbon_spectrum_20261104')
sys.path.insert(0,RIB)
from ribbon import transfer_system
order,states,mats,final=transfer_system(3)
N=len(order)
def lt(b): return frozenset(g for g in (1,2,3) if (b>>(g-1))&1)
A=[]
for b in range(8):
    M=np.zeros((N,N),dtype=np.int64)
    for (j,i),m in mats[lt(b)].items(): M[j,i]+=m
    A.append(M)
e0=np.zeros(N,dtype=np.int64); e0[0]=1
v0=A[0]@e0
h=np.zeros(N,dtype=np.int64)
for j in final: h[j]=1
P1,P2=2147483629,2147483587
def rank_mod(M,p):
    M=[ [int(x)%p for x in row] for row in M]
    r=0; rows=len(M); cols=len(M[0]) if rows else 0
    for c in range(cols):
        piv=None
        for i in range(r,rows):
            if M[i][c]: piv=i;break
        if piv is None: continue
        M[r],M[piv]=M[piv],M[r]
        inv=pow(M[r][c],p-2,p)
        M[r]=[x*inv%p for x in M[r]]
        for i in range(rows):
            if i!=r and M[i][c]:
                f=M[i][c]; M[i]=[(a-f*b)%p for a,b in zip(M[i],M[r])]
        r+=1
    return r
# reachable span basis by BFS over words (raw vectors), greedy independence mod P1
basis=[];vecs=[v0.copy()];frontier=[v0.copy()]
def indep(cands,newv,p=P1):
    return rank_mod(cands+[list(newv)],p)>len(cands)
B=[]
for v in [v0]:
    B.append(list(v))
frontier=[v0]
while frontier:
    nf=[]
    for v in frontier:
        for b in range(8):
            w=A[b]@v
            if not w.any(): continue
            if indep(B,w): B.append(list(w)); nf.append(w)
    frontier=nf
r=len(B); print('reachable span dim',r)
Rb=np.array(B,dtype=np.int64).T   # N x r
# observable rows: h A_w
O=[list(h)];frontier=[h]
while frontier:
    nf=[]
    for u in frontier:
        for b in range(8):
            w=u@A[b]
            if indep([list(np.array(o)@Rb) for o in O], w@Rb):   # independence on reachable span
                O.append(list(w)); nf.append(w)
    frontier=nf
O=np.array(O,dtype=np.int64); print('observable rank on reachable span',O.shape[0])
units=['0.2','1.2','0.123','2.13','0.1.3','0.3.1','0.13.2','0.1.123','0.3.12','0.12.3','1.2.3','1.12.3','0.1.23']
def tob(s): return sum(1<<(int(c)-1) for c in s if c in '123')
U=[[tob(t) for t in u.split('.')] for u in units]
def prod(word):
    M=np.identity(N,dtype=np.int64)
    for b in word: M=A[b]@M
    return M
def mrank(word):   # rank in minimal rep
    return rank_mod((O@prod(word)@Rb).tolist(),P1)
for name,w in zip(units,U):
    ks=[mrank(w*k) for k in range(1,6)]
    print('unit',name,'letters',w,'ranks of powers 1..5',ks)
print('rank of empty letter powers',[mrank([0]*k) for k in range(1,5)])
pickle.dump({'A':A,'v0':v0,'h':h,'O':O,'Rb':Rb,'units':units,'U':U,'N':N},open('lit_base.pkl','wb'))
