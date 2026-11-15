import sys, pickle, numpy as np
import os
RIB=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','..','..','..','04-computation','experiments','ribbon_spectrum_20261104')
sys.path.insert(0,RIB)
from ribbon import transfer_system, count_with_matrices
from sympy import Matrix
W=3; order,states,mats,final=transfer_system(W)
print('profile states total',len(order),'letters',len(mats))
D=pickle.load(open('w30_hankel_13units.pkl','rb')); U=D['units']
def lets(s): return [frozenset(int(c) for c in t if c in '123') for t in s.split('.')]
# letters are frozensets; restrict back arcs by position handled by word() in w30: letter at position q keeps x<=q
N=len(order)
def letmat(a):
    A=np.zeros((N,N),dtype=object)
    for (j,i),m in mats[a].items(): A[j,i]+=m
    return A
# unit matrices assume position q>=3 (all back arcs allowed) -- check positions: first vertex letter 0 is empty set
LM={a:letmat(a) for a in mats}
def unitmat(u):
    A=np.identity(N,dtype=object)
    for a in lets(u): A=LM[a].dot(A)
    return A
A={u:unitmat(u) for u in U}
v0=np.zeros(N,dtype=object); v0[0]=1
v0=LM[frozenset()].dot(v0)   # first vertex
fin=np.zeros(N,dtype=object)
for j in final: fin[j]=1
# reachable states at unit boundaries
reach=set(np.nonzero(v0)[0]); frontier=list(reach)
while frontier:
    nf=[]
    for j in frontier:
        e=np.zeros(N,dtype=object); e[j]=1
        for u in U:
            for k in np.nonzero(A[u].dot(e))[0]:
                if k not in reach: reach.add(k); nf.append(k)
    frontier=nf
reach=sorted(reach)
print('reachable profile states at unit boundaries',len(reach))
# check count against Hankel rep on random words, using profile with positions>=3 (units start at position 1; letter positions small!)
import random
rng=random.Random(1)
Mh=D['red']; z0=D['z0']; finh=D['fin']
bad=0
for t in range(200):
    w=[rng.randrange(13) for _ in range(rng.randint(1,5))]
    v=v0.copy()
    for i in w: v=A[U[i]].dot(v)
    h1=fin.dot(v)
    z=z0
    for i in w: z=Mh[i]*z
    if h1!=(finh*z)[0]: bad+=1
print('mismatch vs Hankel (positions<3 ignored)',bad)
pickle.dump({'reach':reach,'A':{u:A[u][np.ix_(reach,reach)] for u in U},'v0':v0[reach],'fin':fin[reach]},open('profile_units.pkl','wb'))
# reachable span dimension
R=np.array([[int(x) for x in A[U[0]][np.ix_(reach,reach)].flatten()]])
