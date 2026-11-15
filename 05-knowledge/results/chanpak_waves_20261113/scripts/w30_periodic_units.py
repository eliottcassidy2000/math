import os
HERE=os.path.dirname(os.path.abspath(__file__))
RIB=os.path.join(HERE,'..','..','..','..','04-computation','experiments','ribbon_spectrum_20261104')
import pickle, itertools, numpy as np
D=pickle.load(open(os.path.join(HERE,'w30_hankel_13units.pkl'),'rb')); units=D['units']
M=[np.array(m.tolist(),dtype=float) for m in D['red']]
Lw=[len(u.split('.')) for u in units]
best=[]; seen=set()
def canon(c): return min(c[i:]+c[:i] for i in range(len(c)))
for k in range(1,6):
    for c in itertools.product(range(len(units)),repeat=k):
        cc=canon(c)
        if cc in seen: continue
        seen.add(cc)
        P=np.identity(8)
        for i in cc: P=M[i]@P
        ev=sorted(abs(np.linalg.eigvals(P)),reverse=True)
        if ev[0]<1e-9: continue
        # bottleneck test: rank of P^8 (normalized) via singular values
        Q=np.linalg.matrix_power(P/ev[0],8); sv=np.linalg.svd(Q,compute_uv=False)
        r=int((sv>1e-6*sv[0]).sum())
        if r<2 and ev[1]/ev[0]<1e-6: continue
        L=sum(Lw[i] for i in cc); g=ev[0]**(1/L)
        best.append((g,L,'|'.join(units[i] for i in cc),round(ev[1]/ev[0],4)))
best.sort()
print('bottleneck-free periodic unit words (<= 5 units), lowest growth per vertex:')
for b in best[:8]: print(b)
print('count',len(best))
