"""Independent small rooted lift controls; no producer imports."""
from itertools import combinations,permutations
from collections import Counter
from functools import lru_cache
import json

def graph(g,d):
 n=3*g;A=[[0]*n for _ in range(n)]
 for i in range(g):
  for x in range(3):A[3*i+x][3*i+(x+1)%3]=1
 for i,j in combinations(range(g),2):
  for x in range(3):
   for y in range(3):
    z=int(y==(x+(d if (i,j)==(0,g-1) else 0))%3) if j==i+1 or (i,j)==(0,g-1) else 1
    A[3*i+x][3*j+y]=z;A[3*j+y][3*i+x]=1-z
 return A

def lift(A):
 n=len(A);B=[r+[1,0] for r in A]
 B += [[0]*n+[0,1],[1]*n+[0,0]]
 return B

@lru_cache(None)
def canon(bits,k):
 M=[[0]*k for _ in range(k)]
 for z,(i,j) in zip(bits,combinations(range(k),2)):M[i][j]=z;M[j][i]=1-z
 return min(tuple(M[p[i]][p[j]] for i,j in combinations(range(k),2)) for tail in permutations(range(1,k)) for p in [(0,)+tail])
def profiles(A,k):
 n=len(A);out=[]
 for u in range(n):
  row=[]
  for size in range(1,k+1):
   c=Counter()
   for S in combinations([v for v in range(n) if v!=u],size-1):
    order=(u,)+S;bits=tuple(A[order[i]][order[j]] for i,j in combinations(range(size),2))
    c[canon(bits,size)]+=1
   row.append(tuple(sorted(c.items())))
  out.append(tuple(row))
 return out
records=[]
for g in (3,4):
 pair=[graph(g,d) for d in (1,2)]
 anti=[3*(g-1-u//3)+(-u%3) for u in range(3*g)]
 for r in range(3):
  n=len(pair[0]);ps=[profiles(A,g) for A in pair]
  if ps[0]!=ps[1]:raise ValueError('rooted lift profile mismatch')
  for A in pair:
   if not all(A[u][v]==A[anti[v]][anti[u]] for u in range(n) for v in range(n)):raise ValueError('anti-map')
   if r:
    ds=[sum(row) for row in A]
    if ds.count(1)!=1 or ds.count(n-2)!=1 or ds[-2:]!=[1,n-2]:raise ValueError('recoverable extremes')
  records.append({'g':g,'lifts':r,'vertices':n,'cutoff':g,'actual_roots_checked':n,'status':'PASS'})
  pair=list(map(lift,pair));anti+= [n+1,n]
def incidence(A):
 n=len(A);c=[0]*n;e=[[0]*n for _ in range(n)]
 for u,v,w in combinations(range(n),3):
  if A[u][v]==A[v][w]==A[w][u]:
   for x in (u,v,w):c[x]+=1
   for x,y in ((u,v),(u,w),(v,w)):e[x][y]+=1;e[y][x]+=1
 return sorted(tuple(sorted((c[v],e[u][v]) for v in range(n) if A[u][v])) for u in range(n))
triple=[graph(5,d) for d in (0,1,2)];Ecases=[]
for r in range(3):
 for reverse in (False,True):
  rows=[incidence(list(map(list,zip(*A))) if reverse else A) for A in triple]
  if not rows[0]==rows[1]==rows[2]:raise ValueError('lifted E mismatch')
 Ecases.append({'vertices':len(triple[0]),'lifts':r,'all_three_exact_E_equal_both_orientations':True})
 triple=list(map(lift,triple))
print(json.dumps({'status':'PASS','rooted_records':records,'exact_E_controls':Ecases},indent=2))
