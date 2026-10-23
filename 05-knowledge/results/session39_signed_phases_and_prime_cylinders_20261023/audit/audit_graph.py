"""Independent cyclic-cell reconstruction, no sibling imports."""
from itertools import combinations, product
from collections import Counter
import json

def require(x):
 if not x:raise RuntimeError('independent graph audit failed')
def build(m,d):
 n=3*m;A=[[0]*n for _ in range(n)]
 for i in range(m):
  for x in range(3):A[3*i+x][3*i+(x+1)%3]=1
 for i,j in combinations(range(m),2):
  marked=(j==i+1 or (i==0 and j==m-1))
  phase=d if (i,j)==(0,m-1) else 0
  for x,y in product(range(3),repeat=2):
   z=int(y==(x+phase)%3) if marked else 1
   A[3*i+x][3*j+y]=z;A[3*j+y][3*i+x]=1-z
 return A

def data(A):
 n=len(A);c=[0]*n;e=[[0]*n for _ in A]
 for u,v,w in combinations(range(n),3):
  if A[u][v]==A[v][w]==A[w][u]:
   for x in (u,v,w):c[x]+=1
   for x,y in combinations((u,v,w),2):e[x][y]+=1;e[y][x]+=1
 E=[sorted((c[v],e[u][v]) for v in range(n) if A[u][v]) for u in range(n)]
 # Endpoint recurrence, backwards from last vertex rather than forwarding paths.
 D=[[0]*n for _ in range(1<<n)]
 for v in range(n):D[1<<v][v]=1
 for mask in range(1,1<<n):
  if mask&(mask-1)==0:continue
  for v in range(n):
   if mask>>v&1:
    prev=mask^(1<<v)
    D[mask][v]=sum(D[prev][u] for u in range(n) if prev>>u&1 and A[u][v])
 return {'E':E,'H':sum(D[-1]),'c':c,'scores':[sum(row) for row in A]}
rows=[data(build(5,d)) for d in range(3)]
require(rows[0]['E']==rows[1]['E']==rows[2]['E'])
for m,d in product((3,4,5),range(3)):
 A=build(m,d);R=lambda u:3*(m-1-u//3)+(-u%3)
 require(all(A[u][v]==A[R(v)][R(u)] for u,v in product(range(3*m),repeat=2)))
print(json.dumps({'vertexwise_E_equal':True,'H':[r['H'] for r in rows],'cell_scores':rows[0]['scores'][::3],'cell_c':rows[0]['c'][::3],'self_converse_cases':9},indent=2))
