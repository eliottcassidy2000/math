import json
from verify_pathcovers import need,adj_triangle
from verify_seeds import circulate

def hp(A):
 n=len(A);dp=[{} for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for mask in range(1,1<<n):
  for v,c in dp[mask].items():
   for w in range(n):
    if not(mask>>w&1) and A[v][w]:
     d=dp[mask|1<<w];d[w]=d.get(w,0)+c
 return sum(dp[-1].values())

def fixhp(A,tau):
 n=len(A);p=[(i,tau[i]) for i in range(n) if i<tau[i]];z=next(i for i in range(n) if tau[i]==i);k=len(p)
 if not k:return 1
 dp=[{} for _ in range(1<<k)]
 for j,pair in enumerate(p):
  for v in pair:dp[1<<j][v]=1
 for mask in range(1,1<<k):
  for v,c in dp[mask].items():
   for j,pair in enumerate(p):
    if not mask>>j&1:
     for w in pair:
      if A[v][w]:
       d=dp[mask|1<<j];d[w]=d.get(w,0)+c
 return sum(c for v,c in dp[-1].items() if A[v][z])

def blow(A,B):
 q=len(A);m=len(B)
 return [[A[i//m][j//m] if i//m!=j//m else B[i%m][j%m] for j in range(q*m)] for i in range(q*m)]
results=[]
for q,S,H,D in [(5,{1,2},197030475,7857),(7,{1,2,4},123156134945235,5757723)]:
 A=blow(circulate(q,S),circulate(3,{1}))
 tau=[((-i)%q)*3+(-j)%3 for i in range(q) for j in range(3)]
 got=fixhp(A,tau);need(got==D,'independent pair-mask fixed Hamilton')
 if q==5:need(hp(A)==H,'independent15 Held Karp')
 results.append(dict(n=q*3,d1=got,H_checked=q==5))
A=[[0]];tau=[0]
for k in range(3):
 m=len(A);A=adj_triangle(A);tau=[((-i)%3)*m+tau[j] for i in range(3) for j in range(m)]
got=fixhp(A,tau);need(got==9782085285,'independent27 fixed Hamilton');results.append(dict(n=27,d1=got))
print(json.dumps(results,indent=2))
