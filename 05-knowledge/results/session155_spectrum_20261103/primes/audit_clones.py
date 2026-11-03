def require(x):
 if not x: raise AssertionError("clone audit failed")
from itertools import permutations,product,combinations
from math import comb

def graph(n,c):
 A=[[False]*n for _ in range(n)]
 for i,(u,v) in enumerate(combinations(range(n),2)):
  A[u][v]=bool(c>>i&1);A[v][u]=not A[u][v]
 return A

def expand(A,r,k):
 old=[i for i in range(len(A)) if i!=r];labels=old+['r']*k;N=len(labels)
 return [[False if i==j else (i<j if labels[i]==labels[j]=='r' else A[r if labels[i]=='r' else labels[i]][r if labels[j]=='r' else labels[j]]) for j in range(N)] for i in range(N)]

def direct(A):
 n=len(A);return sum(all(A[p[i]][p[i+1]] for i in range(n-1)) for p in permutations(range(n)))

def slot_by_assignments(A,r,k):
 total=0
 for pi in permutations(i for i in range(len(A)) if i!=r):
  if not pi:total+=1;continue
  valid=[];forced=[]
  if A[r][pi[0]]:valid.append(0)
  for i,(u,v) in enumerate(zip(pi,pi[1:]),1):
   if A[u][r] and A[r][v]:valid.append(i)
   if not A[u][v]:forced.append(i)
  if A[pi[-1]][r]:valid.append(len(pi))
  for assignment in product(valid,repeat=k):
   if set(forced)<=set(assignment):total+=1
 return total
N=0
for n in range(1,5):
 for c in range(1<<(n*(n-1)//2)):
  A=graph(n,c)
  for r in range(n):
   for k in range(0,5):
    ex=expand(A,r,k);actual=direct(ex);computed=slot_by_assignments(A,r,k)
    if actual!=computed: raise AssertionError((n,c,r,k,actual,computed))
    N+=1
print('independent explicit slot assignments vs permutations, incl k=0:',N,'PASS')
for rows,H in [([62,60,120,16,96,8,43],35),([62,60,56,80,96,8,39],49)]:
 A=[[bool(rows[i]>>j&1) for j in range(7)] for i in range(7)]
 require(all(not A[i][i] for i in range(7)) and all(A[i][j]!=A[j][i] for i,j in combinations(range(7),2)))
 adj=lambda rem:[p for p in combinations([i for i in range(7) if i not in rem],3) if (A[p[0]][p[1]] and A[p[1]][p[2]] and A[p[2]][p[0]]) or (A[p[0]][p[2]] and A[p[2]][p[1]] and A[p[1]][p[0]])]
 require(direct(A)==H and not adj({5,6}) and all(adj({i}) for i in range(7)))
 R=[set([i])|{j for j in range(7) if A[i][j]} for i in range(7)]
 for _ in range(7):R=[set().union(*(R[j] for j in R[i])) for i in range(7)]
 require(all(len(x)==7 for x in R))
 print('H',H,'strong, deletion pair (5,6) transitive and no single deletion works')
