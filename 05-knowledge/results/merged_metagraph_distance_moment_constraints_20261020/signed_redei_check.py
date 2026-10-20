from pathlib import Path
ROOT = Path(__file__).resolve().parent
from itertools import permutations, combinations
import json, math

def verify(rows):
 out=[]
 for data in rows:
  n=data['n'];m=data['m'];free=[(i,j) for i in range(n) for j in range(i+2,n)];perms=list(permutations(range(n)))
  for u,fib in enumerate(data['fibres']):
   t=fib[0];A=[[False]*n for _ in range(n)]
   for i in range(n-1):A[i][i+1]=True
   for k,(i,j) in enumerate(free):A[i][j]=bool((t>>k)&1);A[j][i]=not A[i][j]
   for i in range(n-1):A[i+1][i]=False
   H=R=a=0
   for p in perms:
    if all(A[p[i]][p[i+1]] for i in range(n-1)):
     H+=1;R+=(-1)**sum(p[i]>p[j] for i,j in combinations(range(n),2))
    a+=all(A[i][j]==A[p[i]][p[j]] for i,j in combinations(range(n),2))
   g=math.prod(1 if A[i][j] else -1 for i,j in combinations(range(n),2))
   eps=1 if len(fib)%2 else 2
   z=sum((-1)**x.bit_count() for x in fib)
   assert len(fib)==eps*H//a
   assert z*a==(-1)**m*eps*g*R
   out.append({'n':n,'node_code':data['codes'][u],'representative_mask':t,'H':H,'signed_R':R,'edge_sign_product':g,'aut':a,'sheets':eps,'top_z':z})
 return out
if __name__=='__main__':
 rows=json.load(open(ROOT/'audit.json'))
 out=verify(rows);json.dump(out,open(ROOT/'signed_redei_check.json','w'),indent=2)
 print('PASS:',len(out),'nodes, n=3..5')
