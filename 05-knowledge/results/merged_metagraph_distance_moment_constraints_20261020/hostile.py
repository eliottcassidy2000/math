from pathlib import Path
ROOT = Path(__file__).resolve().parent
import json,itertools
import numpy as np
R=json.load(open(ROOT/'audit.json'))[-1]
W=[np.array(x,dtype='int64') for x in R['W']]
Q=[np.array(x,dtype='int64') for x in R['Q']]
k=len(W[0]);mass=np.array(R['mass']);m=R['m']
for a,b,c,d in itertools.permutations(range(k),4):
    if a!=min(a,b,c,d):continue
    delta=np.zeros((k,k),dtype='int64')
    for i,j,v in [(a,b,2),(c,d,2),(a,c,-2),(b,d,-2)]:delta[i,j]=delta[j,i]=v
    if np.any(W[1]+delta<0) or np.any(W[2]-delta<0):continue
    changed=[v.copy() for v in W]
    changed[1]+=delta;changed[2]-=delta
    assert all(np.array_equal(a.sum(axis=1),b.sum(axis=1)) for a,b in zip(W,changed))
    assert all(np.array_equal(a.diagonal(),b.diagonal()) for a,b in zip(W,changed))
    assert all(np.array_equal(a%2,b%2) for a,b in zip(W,changed))
    assert np.array_equal(sum(W),sum(changed))
    assert np.array_equal(changed[0],W[0]) and np.array_equal(changed[m],W[m])
    Qm=Q[m]-2*delta
    bad=[]
    for u,v in itertools.combinations(range(k),2):
        det=int(Qm[u,u]*Qm[v,v]-Qm[u,v]**2)
        if det<0:bad.append([u,v,det,int(Qm[u,u]),int(Qm[v,v]),int(Qm[u,v])])
    if not bad:continue
    out={'n':R['n'],'codes':R['codes'],'mass':R['mass'],'delta_edges':[[a,b,2],[c,d,2],[a,c,-2],[b,d,-2]],'changes':[{ 'layer':s,'u':i,'v':j,'before':int(W[s][i,j]),'after':int(W[s][i,j]+sgn*delta[i,j])} for s,sgn in [(1,1),(2,-1)] for i,j in itertools.combinations(range(k),2) if delta[i,j]],'bad_principal_minors':bad,'top_walsh_before':Q[m].tolist(),'top_walsh_after':Qm.tolist()}
    json.dump(out,open(ROOT/'hostile.json','w'),indent=2)
    print(json.dumps(out,indent=2));break
else:print('NONE')
