from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json
P=Path(__file__).resolve().parent
rows=json.loads((P/'independent-normal.json').read_text())['physical_profiles']
data={(t,j,m):(h,e,r)for t,j,m,h,e,r in rows}
def solve(A,b):
 a=[[F(x)for x in row]+[F(v)]for row,v in zip(A,b)];n=len(a)
 for c in range(n):
  p=next(i for i in range(c,n)if a[i][c]);a[c],a[p]=a[p],a[c]
  z=a[c][c];a[c]=[x/z for x in a[c]]
  for i in range(n):
   if i!=c:
    z=a[i][c];a[i]=[x-z*y for x,y in zip(a[i],a[c])]
 return[a[i][-1]for i in range(n)]
base=list(product((1,2),(1,2),(1,2,3)));grid=list(product((1,2),(1,2),(1,2,3)))
H=solve([[a**t*b**j*c**m for a,b,c in base]for t,j,m in grid],[data[t,j,m][0]for t,j,m in grid])
if H!=[F(1)if b==(1,1,1)else F(30)if b==(2,2,3)else F(0)for b in base]:raise RuntimeError('H tensor coefficients')
basis=list(product((1,2),(1,2)));grid=list(product((1,2),(1,2)))
E=solve([[a**t*b**j for a,b in basis]for t,j in grid],[data[t,j,1][1]for t,j in grid])
if E!=list(map(F,[-3,30,0,0])):raise RuntimeError('E tensor coefficients')
basis=list(product((1,2,3),(1,2)));grid=list(product((1,2,3),(1,2)))
Z=solve([[a**t*b**j for a,b in basis]for t,j in grid],[data[t,j,1][0]-data[t,j,1][1]+data[t,j,1][2]for t,j in grid])
if Z!=list(map(F,[5,-57,0,114,6,90])):raise RuntimeError('augmented tensor coefficients')
print(json.dumps(dict(status='PASS',H_coefficients=[int(x)for x in H],E_coefficients=[int(x)for x in E],augmented_coefficients=[int(x)for x in Z],method='Exact rational elimination on independently counted physical graphs'),indent=2,sort_keys=True))
