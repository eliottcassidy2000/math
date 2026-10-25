import json,itertools,math
from fractions import Fraction as F
from pathlib import Path

def require(ok,msg):
 if not ok:raise RuntimeError(msg)
rows=[]
for r in range(8):
 integral=sum(F(2*((-2)**k)*math.comb(r,k),(k+2)*(k+3)) for k in range(r+1))
 expected=F(0) if r%2 else F(1,(r+1)*(r+3))
 require(integral==expected,'symbolic finite expansion')
 row={'r':r,'n':r+3,'integral':str(integral)}
 if r<=5:
  E=[(0,1)]+[(a,c) for a in (0,1) for c in range(2,r+2)]+[(0,r+2)]
  S=sum((-1)**sum(p[a]>p[b] for a,b in E) for p in itertools.permutations(range(r+3)))
  require(F(S,math.factorial(r+3))==integral,'direct signed order enumeration')
  row['order_sum']=S
 rows.append(row)
print(json.dumps(rows,indent=2))
