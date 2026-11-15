import os
HERE=os.path.dirname(os.path.abspath(__file__))
RIB=os.path.join(HERE,'..','..','..','..','04-computation','experiments','ribbon_spectrum_20261104')
import pickle, itertools, sys, numpy as np
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'w29_plane_automaton.py')).read().split("print('plane states:'")[0].replace("print('units',units,'dim',r)","").replace("for u,M in zip(units,red): print(u, M.tolist())","").replace("print('z0',z0.T.tolist(),'fin',fin.tolist())",""))
N=len(order); L=[len(u.split('.')) for u in units]
Mf=[np.array(M.tolist(),dtype=float) for M in red]
# enumerate simple cycles (Johnson-lite via DFS, small graph)
adj={s:[(j,trans[(s,j)]) for j in range(len(units)) if trans.get((s,j)) is not None] for s in range(N)}
cycles=[]
def dfs(start,s,path,seen):
    for j,t in adj[s]:
        if t==start: cycles.append(path+[j])
        elif t not in seen and t>start: dfs(start,t,path+[j],seen|{t})
for s in range(N): dfs(s,s,[],{s})
res=[]
for c in cycles:
    M=np.identity(r)
    for j in c: M=Mf[j]@M
    rho=max(abs(np.linalg.eigvals(M))); ln=sum(L[j] for j in c)
    res.append((rho**(1/ln), ln, '|'.join(units[j] for j in c)))
res.sort()
print('simple cycles:',len(res))
for x in res[:12]: print(x)
