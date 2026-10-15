import json, numpy as np, sympy, itertools
W={int(k):v for k,v in json.load(open('k6p6_W_allprimes.json')).items()}
P=list(sympy.primerange(3,1000)); Pq=[p for p in P if p>=5]
M=np.array([(W[q]-2*q+6) if q%4==3 else (-W[q]-2*q+20) for q in Pq],dtype=float)
qs=np.array(Pq,dtype=float); idx=[P.index(q) for q in Pq]
forms=[]
for line in open('mf_basis2.out'):
    head,vals=line.split(':'); k,N,D,dim=map(int,head.split())
    if dim!=1: continue
    forms.append(((k,N,D),np.array(eval(vals.strip()),dtype=float)[idx]))
chis=[np.ones_like(qs),np.array([1 if q%4==1 else -1 for q in Pq],float),np.array([1 if q%8 in (1,7) else -1 for q in Pq],float),np.array([1 if q%8 in (1,3) else -1 for q in Pq],float)]
poly=[]
for c in chis: poly+=[qs**2*c,qs*c,c]
cands=[]
by={k:[(m,v) for m,v in forms if m[0]==k] for k in range(2,7)}
for (m,v) in by[3]: cands.append((f'{m}^2',[v*v]))
for (m1,v1),(m2,v2) in itertools.product(by[2],by[4]): cands.append((f'{m1}x{m2}',[v1*v2]))
for (m1,v1),(m2,v2) in itertools.combinations(by[3],2): cands.append((f'{m1}x{m2}',[v1*v2]))
for (m1,v1),(m2,v2) in itertools.combinations_with_replacement(by[2],2): cands.append((f'q {m1}x{m2}',[qs*v1*v2]))
res=[]
for name,cols in cands:
    A=np.column_stack(cols+poly); c,_,rk,_=np.linalg.lstsq(A,M,rcond=None); res.append((np.abs(M-A@c).max(),name,np.round(c[:1],4)))
res.sort(key=lambda t:t[0])
print(len(cands),'product candidates; best:')
for r in res[:8]: print(f'  {r[0]:.4g}  {r[1]}  lead coeff {r[2]}')
print('rational 2-power forms by weight:', {k:len(by[k]) for k in by})
