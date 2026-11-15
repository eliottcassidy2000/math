# Build exact integer generators of memory cones K_s = O A_{u3} A_{u2} A_{u1} R_+^108 for BF triples s
import pickle, numpy as np, random
P=pickle.load(open('profile_units.pkl','rb')); D=pickle.load(open('w30_hankel_13units.pkl','rb'))
U=D['units']; n=len(U)
A=[np.array(P['A'][u],dtype=np.int64) for u in U]
v0=np.array(P['v0'],dtype=np.int64); fin=np.array(P['fin'],dtype=np.int64)
J=[None,0,1,2,3,5,6,7]
O=np.array([fin if j is None else fin@A[j] for j in J],dtype=np.int64)   # 8 x 108
# consistency with the Hankel rep
Mh=[np.array(m.tolist(),dtype=object) for m in D['red']]; z0=np.array(D['z0'].tolist(),dtype=object).flatten()
rng=random.Random(5); bad=0
for t in range(300):
    w=[rng.randrange(n) for _ in range(rng.randint(0,6))]
    x=v0.copy(); z=z0.copy()
    for i in w: x=A[i]@x; z=Mh[i].dot(z)
    if any(int(a)!=b for a,b in zip(O@x,z)): bad+=1
print('O A_W v0 == M_W z0 mismatches:',bad)
al3=pickle.load(open('alive_k3.pkl','rb'))['alive']; al4=pickle.load(open('alive_k4.pkl','rb'))['alive']
states=sorted(al3); sidx={s:i for i,s in enumerate(states)}
edges=[(sidx[w[:3]],w[3],sidx[w[1:]]) for w in sorted(al4)]
print('states',len(states),'edges',len(edges))
def prim(v):
    g=0
    for a in v: g=np.gcd(g,abs(int(a)))
    return tuple(int(a)//g for a in v)
gens=[]
OA=[O@A[i] for i in range(n)]
for s in states:
    u1,u2,u3=s
    G=O@A[u3]@A[u2]@A[u1]
    cols=set()
    for j in range(G.shape[1]):
        c=G[:,j]
        if np.any(c): cols.add(prim(c))
    gens.append(sorted(cols))
print('generators per state: min %d max %d mean %.1f'%(min(map(len,gens)),max(map(len,gens)),np.mean(list(map(len,gens)))))
init=[O@A[s[2]]@A[s[1]]@A[s[0]]@v0 for s in states]
pickle.dump({'states':states,'edges':edges,'gens':gens,'init':[tuple(int(a) for a in v) for v in init],'U':U,
             'M':[np.array(m.tolist(),dtype=object) for m in D['red']],'Lw':[len(u.split('.')) for u in U]},open('cone_data.pkl','wb'))
