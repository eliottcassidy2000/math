# memory-4 cones on triple states: K_s = cone( O A_u3 A_u2 A_u1 A_u0 e_j : u0 u1 u2 u3 BF, j reachable )
import pickle, numpy as np, sys, time
from scipy.optimize import linprog
from scipy.sparse import coo_matrix
P=pickle.load(open('profile_units.pkl','rb')); C=pickle.load(open('cone_data.pkl','rb'))
U=C['U']; n=len(U); Lw=C['Lw']; M=C['M']; states=C['states']; edges=C['edges']
A=[np.array(P['A'][u],dtype=np.int64) for u in U]; v0=np.array(P['v0'],dtype=np.int64); fin=np.array(P['fin'],dtype=np.int64)
J=[None,0,1,2,3,5,6,7]
O=np.array([fin if j is None else fin@A[j] for j in J],dtype=np.int64)
al4=set(pickle.load(open('alive_k4.pkl','rb'))['alive'])
sidx={s:i for i,s in enumerate(states)}; S=len(states)
def prim(v):
    g=0
    for a in v: g=np.gcd(g,abs(int(a)))
    return tuple(int(a)//g for a in v)
gens=[]
for s in states:
    u1,u2,u3=s
    B=O@A[u3]@A[u2]@A[u1]
    cols=set()
    for u0 in range(n):
        if (u0,u1,u2,u3) in al4:
            G=B@A[u0]
            for j in range(G.shape[1]):
                if np.any(G[:,j]): cols.add(prim(G[:,j]))
    gens.append(sorted(cols))
print('deep generators per state: mean %.1f max %d'%(np.mean([len(g) for g in gens]),max(len(g) for g in gens)))
# initial: BF 4-words
init4={}
for w in al4:
    x=A[w[3]]@A[w[2]]@A[w[1]]@A[w[0]]@v0
    init4.setdefault(sidx[w[1:]],[]).append(O@x)
def solve(rho,eta=1e-7):
    rows=[];cols=[];vals=[];b=[];r=0
    for (s,u,t) in edges:
        f=rho**Lw[u]
        for g in gens[s]:
            mg=[float(x) for x in M[u].dot(np.array(g,dtype=object))]; nrm=sum(abs(x) for x in mg)
            for k in range(8):
                if mg[k]!=0: rows.append(r);cols.append(8*t+k);vals.append(-mg[k])
                if g[k]!=0: rows.append(r);cols.append(8*s+k);vals.append(f*g[k])
            b.append(-eta*nrm); r+=1
    for s in range(S):
        for g in gens[s]:
            for k in range(8):
                if g[k]!=0: rows.append(r);cols.append(8*s+k);vals.append(float(g[k]))
            b.append(float(g[0])); r+=1
        for v in init4.get(s,[]):
            for k in range(8):
                if v[k]!=0: rows.append(r);cols.append(8*s+k);vals.append(-float(v[k]))
            b.append(-1.0); r+=1
    Am=coo_matrix((vals,(rows,cols)),shape=(r,8*S)).tocsr()
    res=linprog(np.zeros(8*S),A_ub=Am,b_ub=np.array(b),bounds=[(None,None)]*(8*S),method='highs')
    return res
for rho in [float(x) for x in sys.argv[1:]]:
    t0=time.time(); res=solve(rho)
    print('rho',rho,'status',res.status,'%.0fs'%(time.time()-t0)); sys.stdout.flush()
    if res.status==0:
        Q=2**34
        pickle.dump({'rho_lp':rho,'Q':Q,'Phi':[[int(round(res.x[8*s+k]*Q)) for k in range(8)] for s in range(S)],'states':states,'deep':True},open('cone_cert_deep_%s.pkl'%rho,'wb'))
