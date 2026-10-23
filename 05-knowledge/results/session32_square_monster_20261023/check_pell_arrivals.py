from math import isqrt
from pathlib import Path
import json
D=lambda n:isqrt(2*n-1)-isqrt(n)
limit=1000000
actual=[]
for k in range(2,isqrt(limit)+1):
 n=k*k;silent=D(n)==D(n-1)
 a=isqrt(2*n-1);eq=a*a in (2*n-1,2*n-2)
 if silent!=eq:raise RuntimeError((k,'classification'))
 if silent:actual.append(n)
pred=set();families=[]
m,k=1,1
while k*k<=limit:
 if k>=2:pred.add(k*k);families.append({'n':k*k,'m':m,'k':k,'norm':-1})
 m,k=3*m+4*k,2*m+3*k
k,a=1,0
while k*k<=limit:
 if k>=2:pred.add(k*k);families.append({'n':k*k,'m':2*a,'k':k,'norm':-2})
 k,a=3*k+4*a,2*k+3*a
if actual!=sorted(pred):raise RuntimeError('Pell families')
for n in range(2,limit+1):
 change=D(n)-D(n-1)
 if change not in (-1,0,1):raise RuntimeError('jump')
 if isqrt(n)**2!=n and change<0:raise RuntimeError('non-square drop')
out={'status':'PASS exact integer checks','n_max':limit,'silent_square_births':actual,'Pell_witnesses':sorted(families,key=lambda r:r['n']),'scope':'Finite controls supplement the all-order Pell proof; no Hamiltonian criterion.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
