import sympy as s
from itertools import combinations
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent
n=7;B=s.Matrix([[0 if i==j else (1 if (j-i)%n<=n//2 else -1) for j in range(n)] for i in range(n)])
for tri in combinations(range(n),3):
 a,b,c=tri
 if B[a,b]==B[b,c]==B[c,a]:
  C=B.copy()
  for i,j in [(a,b),(b,c),(c,a)]:C[i,j]=-C[i,j];C[j,i]=-C[j,i]
  dd=list((C**6).diagonal())
  if len(set(dd))>1:B=C;break
M=B/7;I=s.eye(7);Q=M*M;G=-Q*(I+Q).inv();H=Q*(I-Q).inv();d=s.Matrix((H-G).diagonal()).T
v=14*G*d
print('B=',B);print('d=',list(d));print('limit nL=',list(v));print('decimal=',[float(x)for x in v]);print('sum=',sum(v))
json.dump({'B':list(map(list,B.tolist())),'limits':[str(x)for x in v],'decimal':[float(x)for x in v]},open(HERE/'lollipop_obstruction.json','w'),indent=2,default=int)
