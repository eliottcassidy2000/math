import os
HERE=os.path.dirname(os.path.abspath(__file__))
RIB=os.path.join(HERE,'..','..','..','..','04-computation','experiments','ribbon_spectrum_20261104')
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'w29_lp_positive_cone.py')).read().split("def feasible(rho):")[0])
import math
# reachable directions: propagate actual vectors
z0f=np.array(z0.tolist(),dtype=float).flatten()
def coords(s,v):
    B=np.array(Bs[s].tolist(),dtype=float); x,*_=np.linalg.lstsq(B,v,rcond=None); return x
Mf=[np.array(M.tolist(),dtype=float) for M in red]
dirs={s:[] for s in rec}
# BFS over paths from start (initial unit), depth up to 14
from collections import deque
q=deque()
for j in range(len(units)):
    v=Mf[j]@z0f
    key,rk=canon(red[j]); s=states[key]
    q.append((s,v,1))
while q:
    s,v,dep=q.popleft()
    if s in rec: 
        x=coords(s,v); dirs[s].append(x/np.linalg.norm(x))
    if dep>=8: continue
    for j,t in adj[s]:
        w=Mf[j]@v; q.append((t,w/np.abs(w).sum(),dep+1))
# sector per state: angles
def ang(x): return math.atan2(x[1],x[0])
K={}
for s in rec:
    A=sorted(ang(x) for x in dirs[s])
    # sectors should be narrow; take min,max (assume no wrap)
    K[s]=(A[0],A[-1])
    print('state',s,'n',len(A),'sector',round(A[0],6),round(A[-1],6))
# invariance check (numeric) of extreme rays
def ray(a): return np.array([math.cos(a),math.sin(a)])
viol=0
for (s,j,t),X in T.items():
    Xf=np.array(X.tolist(),dtype=float)
    for a in K[s]:
        b=ang(Xf@ray(a))
        if b<K[t][0]-1e-9 or b>K[t][1]+1e-9: viol+=1
print('invariance violations (numeric):',viol)
def feasible(rho):
    idx={s:i for i,s in enumerate(rec)}; nv=2*len(rec); A=[]; b=[]
    for (s,j,t),X in T.items():
        Xf=np.array(X.tolist(),dtype=float)
        for a in K[s]:
            e=ray(a); row=np.zeros(nv)
            row[2*idx[t]:2*idx[t]+2]-= Xf@e; row[2*idx[s]:2*idx[s]+2]+= rho**Lw[j]*e
            A.append(row); b.append(0.0)
    for s in rec:
        for a in K[s]:
            e=ray(a); row=np.zeros(nv); row[2*idx[s]:2*idx[s]+2]=-e; A.append(row); b.append(-1.0)
    res=linprog(np.zeros(nv),A_ub=np.array(A),b_ub=np.array(b),bounds=[(None,None)]*nv,method='highs')
    return res.status==0
lo,hi=1.5,1.8
for _ in range(45):
    mid=(lo+hi)/2
    if feasible(mid): lo=mid
    else: hi=mid
print('LP certificate on reachable sectors: max rho =',lo)
