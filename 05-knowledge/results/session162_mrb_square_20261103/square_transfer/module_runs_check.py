from itertools import permutations,combinations

def check(ok,msg):
 if not ok:raise AssertionError(msg)

def H(adj):
 n=len(adj); N=1<<n; dp=[[0]*n for _ in range(N)]
 for i in range(n):dp[1<<i][i]=1
 for mask in range(1,N):
  for v in range(n):
   if dp[mask][v]:
    for w in range(n):
     if not (mask>>w&1) and adj[v][w]:dp[mask|1<<w][w]+=dp[mask][v]
 return sum(dp[-1])

def covers(adj):
 n=len(adj); c=[0]*(n+1)
 for p in permutations(range(n)):
  for k in range(1,n+1):
   for cuts in combinations(range(1,n),k-1):
    if all(i in cuts or adj[p[i-1]][p[i]] for i in range(1,n)):c[k]+=1
 return c

def outer_coeff(adj):
 n=len(adj); other=list(range(1,n)); a=[0]*(n+1)
 for p in permutations(other):
  for k in range(1,n+1):
   for gs in combinations(range(n),k):
    G=set(gs)
    valid=all(i in G or adj[p[i-1]][p[i]] for i in range(1,n-1))
    valid=valid and all((g==0 or adj[p[g-1]][0]) and (g==n-1 or adj[0][p[g]]) for g in gs)
    if valid:a[k]+=1
 return a

def substitute(adj,B):
 m=len(adj);s=len(B);n=s+m-1; T=[[False]*n for _ in range(n)]
 for i in range(n):
  for j in range(n):
   if i==j:continue
   if i<s and j<s:T[i][j]=B[i][j]
   elif i<s:T[i][j]=adj[0][j-s+1]
   elif j<s:T[i][j]=adj[i-s+1][0]
   else:T[i][j]=adj[i-s+1][j-s+1]
 return T

B=[[False,True,False],[False,False,True],[True,False,False]]
cs=covers(B);check(cs==[0,3,6,6],('covers',cs))
seen21=[]
for m in range(1,6):
 pairs=list(combinations(range(m),2))
 for z in range(1<<len(pairs)):
  A=[[False]*m for _ in range(m)]
  for bit,(i,j) in enumerate(pairs):
   if z>>bit&1:A[j][i]=True
   else:A[i][j]=True
  ac=outer_coeff(A);direct=H(substitute(A,B));pred=sum(ac[k]*cs[k] for k in range(1,min(m,3)+1))
  check(direct==pred and direct%3==0,('outer',m,z,ac,pred,direct))
  if direct==21:seen21.append((m,z))
# B labels 7,8,9, nonuniform 10→7, 10→8, 9→10.
X=[[False]*4 for _ in range(4)]
for i,j in ((0,1),(1,2),(2,0),(3,0),(3,1),(2,3)):X[i][j]=True
c=covers(X);f2=2*c[1]+c[2];check((H(X),f2)==(5,32),('nonuniform X',c,f2))
# ordinal X plus 2 sink vertices, then outer singleton vertices a,b: X+T2→a→b→X+T2
s=6;T=[[False]*8 for _ in range(8)]
for i in range(4):
 for j in range(4):T[i][j]=X[i][j]
for i in range(6):
 for j in range(4,6):
  if i<j:T[i][j]=True
for i in range(s):T[i][6]=True;T[7][i]=True
T[6][7]=True
check(H(T)==133,('133 witness',H(T)))
print('PASS all outer tournaments through order 5; square-edge nonuniform X; direct 133; 21 cases',seen21)
