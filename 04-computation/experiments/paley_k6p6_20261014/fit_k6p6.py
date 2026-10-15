import json, re, itertools
import numpy as np
from fractions import Fraction as F
import sympy
d={int(k):v for k,v in json.load(open('k6p6_values.json')).items()}
P=[p for p in sympy.primerange(3,1501) if p%4==3]
W=np.array([d[q]//(q*(q-1)) for q in P],dtype=float)
qs=np.array(P,dtype=float)
forms=[]
for line in open('mf_basis.out'):
    head,vals=line.split(':')
    k,N,D,dim=map(int,head.split())
    v=eval(vals.strip())
    forms.append(((k,N,D,dim),np.array(v,dtype=float)))
print('forms:',len(forms),' primes:',len(P))
def fit(cols):
    A=np.column_stack(cols)
    c,res,rk,sv=np.linalg.lstsq(A,W,rcond=None)
    r=W-A@c
    return c,np.abs(r).max(),rk
poly=[qs**2,qs,np.ones_like(qs)]
best=[]
for (meta,v) in forms:
    c,err,rk=fit([v]+poly)
    best.append((err,meta,c))
best.sort(key=lambda t:t[0])
print('best single-form fits (max abs residual):')
for err,meta,c in best[:8]: print(f'  {meta}: max|res|={err:.3g}, coeffs={np.round(c,4)}')
# mixed: all weight-5 traces + q*weight-3 traces + poly
w5=[v for (m,v) in forms if m[0]==5]
w3q=[qs*v for (m,v) in forms if m[0]==3]
w4=[v for (m,v) in forms if m[0]==4]
for name,cols in [('all w5 + poly',w5+poly),('all w5 + q*w3 + poly',w5+w3q+poly),('w5+q*w3+w4+poly',w5+w3q+w4+poly)]:
    c,err,rk=fit(cols)
    print(f'{name}: {len(cols)} cols, rank {rk}, max|res| = {err:.3g}')
