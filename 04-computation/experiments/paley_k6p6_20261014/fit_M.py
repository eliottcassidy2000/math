import json, numpy as np, sympy
from fractions import Fraction as F
W={int(k):v for k,v in json.load(open('k6p6_W_allprimes.json')).items()}
P=list(sympy.primerange(3,1000))
Pq=[p for p in P if p>=5]
# M_q from the exact identities: q=3 mod 4: W = M + 2q - 6 ; q=1 mod 4: W = -M - 2q + 20
M=np.array([ (W[q]-2*q+6) if q%4==3 else (-W[q]-2*q+20) for q in Pq ],dtype=float)
qs=np.array(Pq,dtype=float)
idx=[P.index(q) for q in Pq]
forms=[]
for line in open('mf_basis2.out'):
    head,vals=line.split(':'); k,N,D,dim=map(int,head.split())
    v=np.array(eval(vals.strip()),dtype=float)[idx]
    forms.append(((k,N,D,dim),v))
chis={'1':np.ones_like(qs),'m4':np.array([1 if q%4==1 else -1 for q in Pq],float),
      'p8':np.array([1 if q%8 in (1,7) else -1 for q in Pq],float),'m8':np.array([1 if q%8 in (1,3) else -1 for q in Pq],float)}
poly=[]
for c in chis.values(): poly += [qs**2*c, qs*c, c]
w5=[v for m,v in forms if m[0]==5]; w3q=[qs*v for m,v in forms if m[0]==3]
w4=[v for m,v in forms if m[0]==4]; w2q=[qs*v for m,v in forms if m[0]==2]
def fit(cols,label,mask=None):
    A=np.column_stack(cols); y=M
    if mask is not None: A=A[mask]; y=M[mask]
    c,_,rk,_=np.linalg.lstsq(A,y,rcond=None); err=np.abs(y-A@c).max()
    print(f'{label:34s} cols {A.shape[1]:3d} rank {rk:3d} rows {A.shape[0]} max|res| {err:.4g}'); return c,err
fit(poly,'twisted polys')
fit(w5+poly,'w5 + twisted polys')
fit(w5+w3q+poly,'w5 + q*w3 + twisted polys')
c,err=fit(w5+w3q+w4+w2q+poly,'w5+q w3+w4+q w2+polys')
