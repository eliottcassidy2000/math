import json
import numpy as np
import sympy
d={int(k):v for k,v in json.load(open('k6p6_values.json')).items()}
P=[p for p in sympy.primerange(3,1501) if p%4==3]
W=np.array([d[q]//(q*(q-1)) for q in P],dtype=float)
qs=np.array(P,dtype=float)
forms=[]
for line in open('mf_basis.out'):
    head,vals=line.split(':')
    k,N,D,dim=map(int,head.split())
    forms.append(((k,N,D,dim),np.array(eval(vals.strip()),dtype=float)))
def nz(v): return np.abs(v).max()>0
f2=[(m,v) for m,v in forms if m[0]==2 and m[3]==1 and nz(v)]
f3=[(m,v) for m,v in forms if m[0]==3 and nz(v)]
f4=[(m,v) for m,v in forms if m[0]==4 and m[3]==1 and nz(v)]
print(len(f2),'rational weight-2,',len(f3),'weight-3 orbits,',len(f4),'rational weight-4')
poly=[qs**2,qs,np.ones_like(qs)]
res=[]
def tryfit(name,cols):
    A=np.column_stack(cols+poly)
    c,_,rk,_=np.linalg.lstsq(A,W,rcond=None)
    err=np.abs(W-A@c).max()
    res.append((err,name,np.round(c,6)))
for m2,v2 in f2:
    for m4,v4 in f4:
        tryfit(f'w2{m2[:3]} x w4{m4[:3]}',[v2*v4])
for (i,(m,v)) in enumerate(f3):
    tryfit(f'w3{m[:3]}^2',[v*v])
    for (m_,v_) in f3[i+1:]:
        tryfit(f'w3{m[:3]} x w3{m_[:3]}',[v*v_])
res.sort(key=lambda t:t[0])
for err,name,c in res[:10]: print(f'{err:12.4g}  {name}  {c}')
