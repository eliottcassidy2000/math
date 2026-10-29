"""Direct finite-field y-value frequency count, independent of Legendre symbols."""
import json
from pathlib import Path
import sympy as S
from named_family import family
P=Path(__file__).resolve().parent;t=S.symbols('T');stored=json.loads((P/'k3_reduction_counts.json').read_text());out=[]
for entry in stored:
 p=entry['p'];freq=[]
 for b in range(p):
  f=[0]*p
  for y in range(p):f[(y*y+b*y)%p]+=1
  freq.append(f)
 row={'p':p,'direct_counts':{}}
 for i in [273,302,4239]:
  coeffs=[S.Poly(a,t)for a in family(i,t)[0]];total=0
  for z in range(p+1):
   a=[int(f.eval(z))%p for f in coeffs]if z<p else[int(f.nth(w))%p for f,w in zip(coeffs,[2,4,6,8,12])]
   a1,a2,a3,a4,a6=a;total+=1
   for x in range(p):total+=freq[(a1*x+a3)%p][(x**3+a2*x*x+a4*x+a6)%p]
  if total!=entry['counts'][str(i)]:raise ValueError('count mismatch')
  row['direct_counts'][str(i)]=total
 out.append(row)
print(json.dumps({'verified_primes':len(out),'counts':out},indent=2))
