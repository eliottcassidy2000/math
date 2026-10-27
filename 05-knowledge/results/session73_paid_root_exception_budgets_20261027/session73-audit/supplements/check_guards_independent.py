from fractions import Fraction as F
from itertools import permutations,combinations
import json

def check(x,msg):
 if not x:raise RuntimeError(msg)
def rank(rows):
 A=[list(map(F,r)) for r in rows];i=0
 for j in range(len(A[0])):
  k=next((k for k in range(i,len(A)) if A[k][j]),None)
  if k is None:continue
  A[i],A[k]=A[k],A[i];z=A[i][j];A[i]=[v/z for v in A[i]]
  for k in range(i+1,len(A)):
   z=A[k][j]
   if z:A[k]=[a-z*b for a,b in zip(A[k],A[i])]
  i+=1
  if i==len(A):break
 return i

functional={i:1 for i in (8,9,11,12,13,40,41,43,44,56,57)}
rowsout=[]
for n in [3,4,5]:
 E=[(i,j) for i,j in combinations(range(n),2) if j-i>1];d=len(E);N=1<<d;paths=list(permutations(range(n)));guard={p:[] for p in paths}
 for state in range(N):
  T=[[False]*n for _ in range(n)]
  for i,j in combinations(range(n),2):
   T[i][j]=True if j==i+1 else bool(state>>E.index((i,j))&1);T[j][i]=not T[i][j]
  for p in paths:guard[p].append(int(all(T[a][b] for a,b in zip(p,p[1:]))))
 legal={p:g for p,g in guard.items() if any(g)}
 full={sum(1<<i for i,e in enumerate(E) if any(set(e)==set((a,b)) for a,b in zip(p,p[1:]))) for p in legal}
 coefficients=[]
 for p,g in legal.items():
  coef=[]
  for S in range(N):
   value=0;T=S
   while True:
    value+=(-1)**((S^T).bit_count())*g[T]
    if T==0:break
    T=(T-1)&S
   coef.append(value)
  coefficients.append(coef)
  if n==5:check(sum(coef[i]*v for i,v in functional.items())==0,'11-term coefficient relation')
 support={i for row in coefficients for i,x in enumerate(row) if x}
 r=rank(list(legal.values()));cr=rank(coefficients)
 check(r==cr=={3:2,4:8,5:44}[n],'independent rational ranks')
 if n==4:
  check(all((state&1)==1-legal[(2,0,1,3)][state]-legal[(3,1,2,0)][state] for state in range(N)),'n4 monomial witness')
  check(1 not in full,'support a absent')
 if n==5:
  check(8 in support and functional[8]==1,'n5 monomial separation')
 products={tuple(a*b for a,b in zip(g,h)) for g in legal.values() for h in legal.values()}
 pr=rank(list(products));check(pr==N,'two-guard full rank')
 rowsout.append({'n':n,'full_supports':len(full),'actual_coefficient_supports':len(support),'guard_rank_Q':r,'two_guard_rank_Q':pr,'nonzero_possible_paths':len(legal)})
print(json.dumps({'status':'PASS','scope':'independent literal adjacency guard tables, subset finite differences, and rational elimination n3..5','cases':rowsout,'n5_11_term_annihilator_on_all_guards':True,'sharp_depth_n7':'NOT COMPUTED'},indent=2,sort_keys=True))
