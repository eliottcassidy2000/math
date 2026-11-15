import os
HERE=os.path.dirname(os.path.abspath(__file__))
RIB=os.path.join(HERE,'..','..','..','..','04-computation','experiments','ribbon_spectrum_20261104')
import pickle, sys, numpy as np, itertools
from fractions import Fraction as Fr
from scipy.optimize import linprog
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'w29_plane_automaton.py')).read().split("print('plane states:'")[0].replace("print('units',units,'dim',r)","").replace("for u,M in zip(units,red): print(u, M.tolist())","").replace("print('z0',z0.T.tolist(),'fin',fin.tolist())",""))
from sympy import Matrix
N=len(order); Lw=[len(u.split('.')) for u in units]
rec=[s for s in range(N) if order[s][2]==2]
# recurrent part: states reachable from themselves
adj={s:[(j,trans[(s,j)]) for j in range(len(units)) if trans.get((s,j)) is not None] for s in range(N)}
def reach(s):
    seen=set(); st=[s]
    while st:
        x=st.pop()
        for j,t in adj[x]:
            if t not in seen: seen.add(t); st.append(t)
    return seen
rec=[s for s in range(N) if s in reach(s)]
print('recurrent states',rec)
Bs={}; 
for s in rec:
    C=order[s][1].columnspace(); Bs[s]=Matrix.hstack(*C[:2])
# 2x2 transfer per edge: M_u B_s = B_t T
T={}
for s in rec:
    for j,t in adj[s]:
        if t not in rec: continue
        X=(Bs[t].T*Bs[t]).inv()*Bs[t].T*(red[j]*Bs[s])
        assert Bs[t]*X==red[j]*Bs[s]
        T[(s,j,t)]=X
# positive cone of plane: x with B_s x >= 0 ; extreme rays in R^2
def rays(B):
    Bf=np.array(B.tolist(),dtype=float); cand=[]
    for i in range(Bf.shape[0]):
        a=Bf[i]
        if np.allclose(a,0): continue
        for v in ([-a[1],a[0]],[a[1],-a[0]]):
            v=np.array(v)
            if np.all(Bf@v>=-1e-12): cand.append(v/np.linalg.norm(v))
    # dedupe
    out=[]
    for v in cand:
        if not any(np.allclose(v,o) for o in out): out.append(v)
    return out
R={s:rays(Bs[s]) for s in rec}
print({s:len(R[s]) for s in rec})
def feasible(rho):
    idx={s:i for i,s in enumerate(rec)}; nv=2*len(rec)
    A=[]; b=[]
    for (s,j,t),X in T.items():
        Xf=np.array(X.tolist(),dtype=float)
        for e in R[s]:
            row=np.zeros(nv)
            # -(phi_t X e - rho^L phi_s e) <= 0
            row[2*idx[t]:2*idx[t]+2]-= Xf@e
            row[2*idx[s]:2*idx[s]+2]+= rho**Lw[j]*e
            A.append(row); b.append(0.0)
    for s in rec:
        for e in R[s]:
            row=np.zeros(nv); row[2*idx[s]:2*idx[s]+2]=-e; A.append(row); b.append(-1.0)
    res=linprog(np.zeros(nv),A_ub=np.array(A),b_ub=np.array(b),bounds=[(None,None)]*nv,method='highs')
    return res.status==0, (res.x if res.status==0 else None)
lo,hi=1.5,1.8
for _ in range(40):
    mid=(lo+hi)/2
    ok,_=feasible(mid)
    if ok: lo=mid
    else: hi=mid
print('LP certificate: max rho =',lo)
